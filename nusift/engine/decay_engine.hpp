#pragma once
/**
 * @file
 * @brief Decays an inventory to a set of times, producing atoms and their exact time integrals.
 * @ingroup engine
 */
//
// THE AUGMENTED SYSTEM. Under pure decay the matrix A is constant, so augmenting it
//
//     M = [ A  0 ]        z0 = [ n0 ]        exp(M t) z0 = [ exp(A t) n0            ]
//         [ I  0 ]             [ 0  ]                      [ int_0^t exp(A tau) n0  ]
//
// yields, from ONE solve at t, both the inventory and the EXACT cumulative time integral --
// no quadrature, no error that depends on how finely the time grid was chosen. The bottom
// block is what every time-integrated metric is built from: total decays over an interval,
// integrated exposure, dose accrued between two cooling times.
//
// Because the integral comes for free with the inventory, there is no separate "plain A"
// path to maintain and no configuration in which NuSIFT computes one without the other.
//
#include <span>
#include <string>

#include "nusift/engine/decay_result.hpp"
#include "nusift/engine/inventory.hpp"

namespace nusift {

class NuclearData;

// The CRAM approximation order. Declared here rather than reusing cram::CramOrder because
// cram/cram.hpp includes Eigen, and no public NuSIFT header may.
enum class CramOrder {
  Order16 = 16,  // cheaper; adequate for screening
  Order48 = 48,  // the default, and what any reported number should use
};

struct DecayOptions {
  CramOrder order = CramOrder::Order48;

  // Restrict the solve to nuclides forward-reachable from the seed. This is exact, not an
  // approximation: the reachable set is closed under production, so nuclides outside it are
  // identically zero for all time and contribute nothing to any metric. It is also the single
  // largest performance lever, because sparse LU cost grows superlinearly -- a single-nuclide
  // seed touches tens of nuclides rather than thousands -- and it keeps the output legible by
  // omitting thousands of permanently-zero rows.
  bool prune = true;

  // Worker threads for the per-time solves. 0 means "as many as the hardware reports"; 1
  // forces the serial path.
  //
  // The solves are embarrassingly parallel: each time gets its own factorization of the same
  // matrix and shares nothing with the others. cram documents CramSolver as not thread-safe,
  // so each worker holds its own -- the pole tables it reads are immutable constants, and the
  // per-pole factorizations are what must not be shared.
  //
  // Factorizing dominates the cost, so the speedup is close to linear until memory bandwidth
  // bites. A single-time run stays serial regardless, since spawning a thread to do one solve
  // costs more than it saves.
  int threads = 0;
};

// Decay `inventory` to every time in `times` (seconds, ascending, non-negative).
//
// Throws InputError if the times are unsorted or negative, or if the inventory names a
// nuclide the data store does not carry.
DecayResult decay(const NuclearData& data, const Inventory& inventory,
                  std::span<const double> times, const DecayOptions& options = {});

// The exact per-nuclide integral \int_{t1}^{t2} n(tau) dtau, in atom-seconds, over the same
// pruned index space `decay()` would produce for this seed.
//
// Computed as G(t2) - G(t1) from two augmented solves, EXCEPT when that difference would lose
// its significant digits to cancellation -- both integrals large and nearly equal. There the
// interval is re-solved directly from t1, which is exact and costs one more solve. Without
// that guard the naive difference silently returns noise, and it does so precisely in the
// regime a user asking for a narrow late-time window cares about.
//
// The guard is judged over the worst-affected nuclide in the chain, so on a full chain the
// re-solve is the ordinary path for any t1 > 0 rather than an exception -- expect three
// factorizations, not two. Costing one solve to keep the accurate path is the intended trade;
// it is worth knowing it is the one usually taken.
std::vector<double> intervalIntegral(const NuclearData& data, const Inventory& inventory, double t1,
                                     double t2, std::vector<std::int64_t>* keysOut,
                                     const DecayOptions& options = {});

}  // namespace nusift
