#include "nusift/engine/decay_engine.hpp"

#include <Eigen/SparseCore>
#include <algorithm>
#include <cmath>
#include <deque>
#include <exception>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include "cram/chain.hpp"
#include "cram/cram_solver.hpp"
#include "nusift/core/error.hpp"
#include "nusift/core/nuclide_name.hpp"
#include "nusift/nucdata/nuclear_data.hpp"
#include "nusift/nucdata/nuclear_data_internal.hpp"

// This and nucdata/nuclear_data.cpp are the only translation units that see cram or Eigen.
// Everything crossing out of here is plain std:: types, which is what keeps them off the
// public interface and out of every consumer's compile.

namespace nusift {
namespace {

constexpr const char* kModule = "decay engine";

[[noreturn]] void fail(const std::string& what) {
  throw NusiftError(tagged(kModule, what));
}

cram::CramOrder toCramOrder(CramOrder order) {
  return order == CramOrder::Order16 ? cram::CramOrder::CRAM16 : cram::CramOrder::CRAM48;
}

// Map the inventory onto chain indices, raising a specific error for a nuclide the store does
// not carry. Silently dropping it would understate every metric with no indication why.
std::vector<double> seedVector(const NuclearData& data, const Inventory& inventory) {
  std::vector<double> seed(static_cast<std::size_t>(data.size()), 0.0);
  for (const InventoryEntry& entry : inventory.entries()) {
    const int index = data.indexOfKey(entry.zaiKey);
    if (index < 0) {
      throw InputError(tagged(kModule, "the data store has no nuclide " +
                                           formatNuclideName(Zai::fromKey(entry.zaiKey)) +
                                           "; it cannot be decayed"));
    }
    seed[static_cast<std::size_t>(index)] += entry.atoms;
  }
  return seed;
}

// Indices forward-reachable from the seed, ascending.
//
// Reachability is read off the decay matrix's own sparsity rather than re-walking decay modes:
// A(j, i) != 0 means i produces j, which captures decay daughters, branching, and spontaneous
// fission products uniformly, with no chance of disagreeing with the matrix actually solved.
// The resulting set is closed under production, so restricting to it is exact.
std::vector<int> forwardClosure(const Eigen::SparseMatrix<double>& a,
                                const std::vector<double>& seed) {
  const int n = static_cast<int>(seed.size());
  std::vector<char> reached(static_cast<std::size_t>(n), 0);
  std::deque<int> queue;
  for (int i = 0; i < n; ++i) {
    if (seed[static_cast<std::size_t>(i)] != 0.0) {
      reached[static_cast<std::size_t>(i)] = 1;
      queue.push_back(i);
    }
  }
  while (!queue.empty()) {
    const int i = queue.front();
    queue.pop_front();
    // Column i holds everything i produces. Eigen's default storage is column-major, so this
    // iterates exactly that column.
    for (Eigen::SparseMatrix<double>::InnerIterator it(a, i); it; ++it) {
      const int j = static_cast<int>(it.row());
      if (it.value() != 0.0 && reached[static_cast<std::size_t>(j)] == 0) {
        reached[static_cast<std::size_t>(j)] = 1;
        queue.push_back(j);
      }
    }
  }
  std::vector<int> keep;
  keep.reserve(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    if (reached[static_cast<std::size_t>(i)] != 0) {
      keep.push_back(i);
    }
  }
  return keep;
}

// Restrict A to `keep`. Exact, because `keep` is closed under production: any entry dropped
// has a column outside the set, and those nuclides are identically zero for all time.
Eigen::SparseMatrix<double> restrict(const Eigen::SparseMatrix<double>& a,
                                     const std::vector<int>& keep, int fullSize) {
  std::vector<int> local(static_cast<std::size_t>(fullSize), -1);
  for (std::size_t k = 0; k < keep.size(); ++k) {
    local[static_cast<std::size_t>(keep[k])] = static_cast<int>(k);
  }
  const int m = static_cast<int>(keep.size());
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(static_cast<std::size_t>(a.nonZeros()));
  for (const int col : keep) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(a, col); it; ++it) {
      const int row = local[static_cast<std::size_t>(it.row())];
      if (row >= 0) {
        triplets.emplace_back(row, local[static_cast<std::size_t>(col)], it.value());
      }
    }
  }
  Eigen::SparseMatrix<double> reduced(m, m);
  reduced.setFromTriplets(triplets.begin(), triplets.end());
  return reduced;
}

// [[A, 0], [I, 0]] -- the augmented generator whose exponential carries the inventory in its
// top block and the exact cumulative integral in its bottom.
Eigen::SparseMatrix<double> augment(const Eigen::SparseMatrix<double>& a) {
  const int n = static_cast<int>(a.rows());
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(static_cast<std::size_t>(a.nonZeros()) + static_cast<std::size_t>(n));
  for (int col = 0; col < a.outerSize(); ++col) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(a, col); it; ++it) {
      triplets.emplace_back(static_cast<int>(it.row()), static_cast<int>(it.col()), it.value());
    }
  }
  for (int i = 0; i < n; ++i) {
    triplets.emplace_back(n + i, i, 1.0);  // bottom-left identity: y' = x
  }
  Eigen::SparseMatrix<double> m(2 * n, 2 * n);
  m.setFromTriplets(triplets.begin(), triplets.end());
  return m;
}

// How many workers to actually run. Never more than there are times to solve, and never more
// than one when there is a single time -- spawning a thread to do one solve costs more than it
// saves. A request of 0 means "ask the hardware", which can itself report 0.
int resolveThreadCount(int requested, int timeCount) {
  int workers = requested;
  if (workers <= 0) {
    workers = static_cast<int>(std::thread::hardware_concurrency());
  }
  if (workers <= 0) {
    workers = 1;
  }
  return std::min(workers, std::max(1, timeCount));
}

void requireValidTimes(std::span<const double> times) {
  if (times.empty()) {
    throw InputError(tagged(kModule, "no times requested"));
  }
  for (std::size_t k = 0; k < times.size(); ++k) {
    if (!std::isfinite(times[k]) || times[k] < 0.0) {  // isfinite also rejects NaN
      throw InputError(tagged(
          kModule, "times must be non-negative and finite; got " + std::to_string(times[k])));
    }
    if (k > 0 && times[k] <= times[k - 1]) {
      throw InputError(tagged(kModule, "times must be strictly ascending; " +
                                           std::to_string(times[k]) + " follows " +
                                           std::to_string(times[k - 1])));
    }
  }
}

// Shared setup: seed, prune, and build the augmented matrix over the kept index space.
struct Prepared {
  std::vector<int> keep;           // chain indices retained, ascending
  std::vector<std::int64_t> keys;  // their ZAI keys, same order
  Eigen::SparseMatrix<double> augmented;
  Eigen::VectorXd seed;  // length keep.size()
};

Prepared prepare(const NuclearData& data, const Inventory& inventory, const DecayOptions& options) {
  const cram::DepletionChain& chain = chainOf(data);
  const int n = chain.size();
  if (n == 0) {
    fail("the data store contains no nuclides");
  }

  const std::vector<double> seedFull = seedVector(data, inventory);
  const Eigen::SparseMatrix<double> full = chain.decayMatrix();
  if (full.rows() != n || full.cols() != n) {
    fail("decay matrix size disagrees with the chain size");
  }

  Prepared out;
  if (options.prune) {
    out.keep = forwardClosure(full, seedFull);
  } else {
    out.keep.resize(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
      out.keep[static_cast<std::size_t>(i)] = i;
    }
  }
  if (out.keep.empty()) {
    throw InputError(tagged(kModule, "the inventory is empty, so there is nothing to decay"));
  }

  const int m = static_cast<int>(out.keep.size());
  out.keys.resize(static_cast<std::size_t>(m));
  out.seed = Eigen::VectorXd::Zero(m);
  for (int k = 0; k < m; ++k) {
    const int global = out.keep[static_cast<std::size_t>(k)];
    out.keys[static_cast<std::size_t>(k)] = data.zaiAt(global).key();
    out.seed(k) = seedFull[static_cast<std::size_t>(global)];
  }

  const Eigen::SparseMatrix<double> reduced = options.prune ? restrict(full, out.keep, n) : full;
  out.augmented = augment(reduced);
  return out;
}

}  // namespace

DecayResult decay(const NuclearData& data, const Inventory& inventory,
                  std::span<const double> times, const DecayOptions& options) {
  requireValidTimes(times);

  const Prepared prepared = prepare(data, inventory, options);
  const int m = static_cast<int>(prepared.keys.size());
  const int nT = static_cast<int>(times.size());

  DecayResult result;
  result.times.assign(times.begin(), times.end());
  result.nuclideKeys = prepared.keys;
  result.atoms.assign(static_cast<std::size_t>(nT) * static_cast<std::size_t>(m), 0.0);
  result.integratedAtoms.assign(static_cast<std::size_t>(nT) * static_cast<std::size_t>(m), 0.0);
  result.seedProvenance = inventory.provenance();

  Eigen::VectorXd z0 = Eigen::VectorXd::Zero(2 * m);
  z0.head(m) = prepared.seed;

  // One factorization per time, applied once. Nothing is shared between times, so the loop
  // parallelises directly -- each worker owns its own solver because cram documents
  // CramSolver as not thread-safe.
  const auto solveOne = [&](cram::CramSolver& solver, int k) {
    const double t = times[static_cast<std::size_t>(k)];
    const std::size_t base = static_cast<std::size_t>(k) * static_cast<std::size_t>(m);
    if (t <= 0.0) {
      // t = 0 is the seed itself with a zero integral. Solving would be harmless but wasteful,
      // and CRAM at dt = 0 is a needless round trip through 48 factorizations.
      for (int i = 0; i < m; ++i) {
        result.atoms[base + static_cast<std::size_t>(i)] = prepared.seed(i);
      }
      return;
    }
    solver.prepare(prepared.augmented, t);
    const Eigen::VectorXd z = solver.apply(z0);
    for (int i = 0; i < m; ++i) {
      result.atoms[base + static_cast<std::size_t>(i)] = z(i);
      result.integratedAtoms[base + static_cast<std::size_t>(i)] = z(m + i);
    }
  };

  const int workers = resolveThreadCount(options.threads, nT);
  if (workers <= 1) {
    cram::CramSolver solver(toCramOrder(options.order));
    for (int k = 0; k < nT; ++k) {
      solveOne(solver, k);
    }
    return result;
  }

  // Times are dealt round-robin rather than in contiguous blocks. Cost rises with t -- a
  // later time needs more of the matrix -- so contiguous blocks would leave the worker holding
  // the early times idle while the last one finishes.
  std::mutex failureMutex;
  std::exception_ptr failure;
  std::vector<std::thread> pool;
  pool.reserve(static_cast<std::size_t>(workers));
  for (int w = 0; w < workers; ++w) {
    pool.emplace_back([&, w]() {
      try {
        cram::CramSolver solver(toCramOrder(options.order));
        for (int k = w; k < nT; k += workers) {
          solveOne(solver, k);
        }
      } catch (...) {
        // Rethrown on the calling thread once every worker has stopped. Letting it escape here
        // would call std::terminate and lose the message entirely.
        const std::lock_guard<std::mutex> lock(failureMutex);
        if (!failure) {
          failure = std::current_exception();
        }
      }
    });
  }
  for (std::thread& worker : pool) {
    worker.join();
  }
  if (failure) {
    std::rethrow_exception(failure);
  }

  return result;
}

std::vector<double> intervalIntegral(const NuclearData& data, const Inventory& inventory, double t1,
                                     double t2, std::vector<std::int64_t>* keysOut,
                                     const DecayOptions& options) {
  if (!std::isfinite(t1) || !std::isfinite(t2) || t1 < 0.0 || t2 < 0.0) {
    throw InputError(tagged(kModule, "interval endpoints must be non-negative and finite"));
  }
  if (t2 <= t1) {
    throw InputError(tagged(kModule, "interval end must be after its start; got [" +
                                         std::to_string(t1) + ", " + std::to_string(t2) + "]"));
  }

  const Prepared prepared = prepare(data, inventory, options);
  const int m = static_cast<int>(prepared.keys.size());
  if (keysOut != nullptr) {
    *keysOut = prepared.keys;
  }

  Eigen::VectorXd z0 = Eigen::VectorXd::Zero(2 * m);
  z0.head(m) = prepared.seed;

  cram::CramSolver solver(toCramOrder(options.order));

  Eigen::VectorXd lower = Eigen::VectorXd::Zero(2 * m);
  if (t1 > 0.0) {
    solver.prepare(prepared.augmented, t1);
    lower = solver.apply(z0);
  } else {
    lower.head(m) = prepared.seed;
  }

  solver.prepare(prepared.augmented, t2);
  const Eigen::VectorXd upper = solver.apply(z0);

  // How much of G(t2) survives the subtraction, over the WORST nuclide. When the interval is
  // narrow relative to the time already elapsed, G(t1) and G(t2) agree to within rounding and
  // the difference is noise -- so measure that directly rather than guessing from the interval
  // width, which would depend on the decay constants involved.
  //
  // Taking the minimum over every nuclide, not over the ones that matter, means a single
  // saturated trace species is enough to send the whole solve down the re-solve path. On a
  // full chain that is not the rare case: it fires for nearly every interval with t1 > 0. The
  // cost is one extra factorization and the result is the more accurate of the two, so the
  // trade is deliberately in that direction -- but it is the common path, not the exception.
  double worstRetained = 1.0;
  for (int i = 0; i < m; ++i) {
    const double g2 = upper(m + i);
    const double g1 = lower(m + i);
    if (g2 > 0.0) {
      worstRetained = std::min(worstRetained, (g2 - g1) / g2);
    }
  }

  std::vector<double> integral(static_cast<std::size_t>(m), 0.0);

  constexpr double kCancellationFloor = 1e-8;
  if (worstRetained >= kCancellationFloor) {
    for (int i = 0; i < m; ++i) {
      integral[static_cast<std::size_t>(i)] = upper(m + i) - lower(m + i);
    }
    return integral;
  }

  // Re-solve the interval directly: restart the augmented system at t1 with a zeroed
  // integrator block, so the bottom block accumulates only over [t1, t2]. Exact, and it
  // costs one additional factorization in the only regime where it is needed.
  Eigen::VectorXd restart = Eigen::VectorXd::Zero(2 * m);
  restart.head(m) = lower.head(m);
  solver.prepare(prepared.augmented, t2 - t1);
  const Eigen::VectorXd direct = solver.apply(restart);
  for (int i = 0; i < m; ++i) {
    integral[static_cast<std::size_t>(i)] = direct(m + i);
  }
  return integral;
}

}  // namespace nusift
