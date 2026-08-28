//
// A measurement harness, not a feature and not a test.
//
// Two questions gate whether adjoint sensitivity is a subcommand, a batch job, or a bad idea,
// and neither can be answered from cram's headers:
//
//   COST         rateSensitivities() is quadrature over the forward and adjoint trajectories,
//                not a matrix exponential: roughly (subIntervals + 2*endRefinements) pieces x
//                gaussPoints nodes x two directions of solve, per schedule interval. That is
//                ~160 solves where a forecast point costs one. Three orders of magnitude is
//                the estimate; this measures it.
//
//   CONDITIONING cram tuned the SensitivityOptions defaults for "a thermal pin with a few
//                intervals of days to months". A decay chain spans microsecond isomers to
//                primordials, so both trajectories carry boundary layers far thinner than
//                that tuning anticipated.
//
// WHAT IT FOUND, since a harness whose comment still states the refuted hypothesis is worse
// than no comment. The conditioning problem is real -- a single [0,30 d] interval at the
// default endRefinements gets dR/dlambda wrong by 37% for the Cs-137/Ba-137m pair -- but
// splitting the schedule is NOT the fix. Refinement is: endRefinements = 12 on one interval
// lands at 1e-7, while an 8-point log grid at the default 6 only reaches 8e-3 for four times
// the work. cram's own rule predicts this exactly (h / 2^k below the shortest removal time:
// 30 d / 4 / 2^12 = 158 s, against Ba-137m's 153 s half-life), so the lever is the one cram
// documents, applied with NuSIFT's dynamic range in mind rather than a pin's.
//
// Correctness is checked against central finite differences on lambda, which is the only
// reference that does not assume the thing being measured.
//
// The third question, added once cost stopped being the interesting one: WHICH PARAMETER
// CLASS carries the uncertainty. The answer is that the two are evenly matched on physics --
// RSS elasticity 0.227 for half-lives against 0.236 for independent yields -- so the whole
// difference between them lies in the evaluated uncertainties, not in the sensitivity. They
// also implicate disjoint nuclides: yields point at the short-lived precursors that were
// seeded, half-lives at the nuclides actually present at the response time.
//
// Usage: nusift_sensitivity_spike [logIntervals]
//
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <string>
#include <vector>

#include <Eigen/SparseCore>

#include "cram/adjoint.hpp"
#include "cram/chain.hpp"
#include "cram/cram.hpp"
#include "cram/cram_solver.hpp"
#include "nusift/core/nuclide.hpp"
#include "nusift/core/nuclide_name.hpp"
#include "nusift/engine/inventory.hpp"
#include "nusift/nucdata/nuclear_data.hpp"
#include "nusift/nucdata/nuclear_data_internal.hpp"
#include "nusift/seed/seed_fission.hpp"

namespace {

using nusift::NuclearData;
using nusift::Zai;

constexpr double kDay = 86400.0;

double seconds(std::chrono::steady_clock::time_point a, std::chrono::steady_clock::time_point b) {
  return std::chrono::duration<double>(b - a).count();
}

// The pruned problem: the forward closure of the seed under production, as decay_engine's
// prepare() computes it. Reproduced here rather than exported from the library, because the
// point of a spike is to leave the library alone until the measurement says what to build.
struct Problem {
  std::vector<int> keep;          // chain indices, ascending
  Eigen::SparseMatrix<double> a;  // decay matrix restricted to keep
  Eigen::VectorXd n0;             // seed over keep
  Eigen::VectorXd w;              // response weight over keep: lambda, i.e. activity in Bq
};

std::vector<int> forwardClosure(const Eigen::SparseMatrix<double>& full,
                                const std::vector<double>& seed) {
  const int n = static_cast<int>(full.rows());
  std::vector<char> seen(static_cast<std::size_t>(n), 0);
  std::vector<int> stack;
  for (int i = 0; i < n; ++i) {
    if (seed[static_cast<std::size_t>(i)] != 0.0) {
      seen[static_cast<std::size_t>(i)] = 1;
      stack.push_back(i);
    }
  }
  while (!stack.empty()) {
    const int col = stack.back();
    stack.pop_back();
    for (Eigen::SparseMatrix<double>::InnerIterator it(full, col); it; ++it) {
      const int row = static_cast<int>(it.row());
      if (row != col && it.value() != 0.0 && seen[static_cast<std::size_t>(row)] == 0) {
        seen[static_cast<std::size_t>(row)] = 1;
        stack.push_back(row);
      }
    }
  }
  std::vector<int> keep;
  for (int i = 0; i < n; ++i) {
    if (seen[static_cast<std::size_t>(i)] != 0) {
      keep.push_back(i);
    }
  }
  return keep;
}

Eigen::SparseMatrix<double> restrictTo(const Eigen::SparseMatrix<double>& full,
                                       const std::vector<int>& keep, int n) {
  std::vector<int> pos(static_cast<std::size_t>(n), -1);
  for (std::size_t k = 0; k < keep.size(); ++k) {
    pos[static_cast<std::size_t>(keep[k])] = static_cast<int>(k);
  }
  const int m = static_cast<int>(keep.size());
  std::vector<Eigen::Triplet<double>> t;
  for (int k = 0; k < m; ++k) {
    const int col = keep[static_cast<std::size_t>(k)];
    for (Eigen::SparseMatrix<double>::InnerIterator it(full, col); it; ++it) {
      const int row = pos[static_cast<std::size_t>(it.row())];
      if (row >= 0) {
        t.emplace_back(row, k, it.value());
      }
    }
  }
  Eigen::SparseMatrix<double> out(m, m);
  out.setFromTriplets(t.begin(), t.end());
  return out;
}

Problem build(const NuclearData& data, const cram::DepletionChain& chain,
              const nusift::Inventory& inv) {
  const int n = chain.size();
  std::vector<double> seedFull(static_cast<std::size_t>(n), 0.0);
  for (const nusift::InventoryEntry& e : inv.entries()) {
    const int i = data.indexOfKey(e.zaiKey);
    if (i >= 0) {
      seedFull[static_cast<std::size_t>(i)] = e.atoms;
    }
  }

  int dropped = 0;
  const Eigen::SparseMatrix<double> full = chain.decayMatrix(&dropped);
  Problem p;
  p.keep = forwardClosure(full, seedFull);
  const int m = static_cast<int>(p.keep.size());
  p.a = restrictTo(full, p.keep, n);
  p.n0 = Eigen::VectorXd::Zero(m);
  p.w = Eigen::VectorXd::Zero(m);
  for (int k = 0; k < m; ++k) {
    const int g = p.keep[static_cast<std::size_t>(k)];
    p.n0(k) = seedFull[static_cast<std::size_t>(g)];
    p.w(k) = data.decayConstant(g);  // R = <lambda, n(T)>, activity in Bq
  }
  return p;
}

// R = <w, n(T)> by one forward solve: the reference the adjoint is checked against, and what
// finite differences perturb.
double responseAt(const Eigen::SparseMatrix<double>& a, const Eigen::VectorXd& n0,
                  const Eigen::VectorXd& w, double t) {
  cram::CramSolver solver(cram::CramOrder::CRAM48);
  solver.prepare(a, t);
  return w.dot(solver.apply(n0));
}

// A log grid over [t0, T] as interval widths -- the schedule NuSIFT already builds for a
// forecast, reused verbatim as the sensitivity schedule.
std::vector<double> logSchedule(double t0, double t, int points) {
  std::vector<double> dts;
  dts.reserve(static_cast<std::size_t>(points));
  const double l0 = std::log(t0);
  const double l1 = std::log(t);
  double prev = 0.0;
  for (int k = 0; k < points; ++k) {
    const double f = static_cast<double>(k + 1) / static_cast<double>(points);
    const double next = std::exp(l0 + (l1 - l0) * f);
    dts.push_back(next - prev);
    prev = next;
  }
  return dts;
}

struct SensitivityRun {
  Eigen::VectorXd dRdLambda;  // indexed like the chain: the IMPLICIT term only
  Eigen::VectorXd nT;         // indexed like the chain: n(T), for the explicit term
  Eigen::VectorXd dRdn0;      // indexed like the chain: adj.nStar.front(), the seed sensitivity
  double marchSeconds = 0.0;
  double quadratureSeconds = 0.0;
  int solves = 0;
};

// dR/dlambda over the whole chain from one forward march and one adjoint march.
//
// S comes back in the pruned space and is scattered into chain space, because
// decayConstantSensitivities() contracts it against a dA/dlambda derived from the full chain.
// Entries outside the closure contract to zero, which is the right answer rather than an
// accident: a nuclide the seed cannot reach has n == 0 for all time, so every entry its lambda
// touches multiplies zero.
SensitivityRun run(const cram::DepletionChain& chain, const Problem& p,
                   const std::vector<double>& dts, const cram::SensitivityOptions& opts) {
  const std::vector<Eigen::SparseMatrix<double>> as(dts.size(), p.a);

  const auto t0 = std::chrono::steady_clock::now();
  const cram::DepletionResult fwd = cram::depleteLinear(as, dts, p.n0);
  const cram::AdjointResult adj = cram::adjointDeplete(as, dts, p.w);
  const auto t1 = std::chrono::steady_clock::now();

  const std::vector<Eigen::SparseMatrix<double>> s =
      cram::rateSensitivities(fwd, adj, as, dts, opts);
  const auto t2 = std::chrono::steady_clock::now();

  const int n = chain.size();
  std::vector<Eigen::SparseMatrix<double>> big;
  big.reserve(s.size());
  for (const Eigen::SparseMatrix<double>& sk : s) {
    std::vector<Eigen::Triplet<double>> t;
    for (int col = 0; col < sk.outerSize(); ++col) {
      for (Eigen::SparseMatrix<double>::InnerIterator it(sk, col); it; ++it) {
        t.emplace_back(p.keep[static_cast<std::size_t>(it.row())],
                       p.keep[static_cast<std::size_t>(it.col())], it.value());
      }
    }
    Eigen::SparseMatrix<double> m(n, n);
    m.setFromTriplets(t.begin(), t.end());
    big.push_back(std::move(m));
  }

  SensitivityRun out;
  out.dRdLambda = cram::decayConstantSensitivities(chain, big);
  out.nT = Eigen::VectorXd::Zero(n);
  out.dRdn0 = Eigen::VectorXd::Zero(n);
  for (int k = 0; k < static_cast<int>(p.keep.size()); ++k) {
    const int g = p.keep[static_cast<std::size_t>(k)];
    out.nT(g) = fwd.n.back()(k);
    // dR/dn0 falls straight out of the adjoint march -- no quadrature, no convergence
    // parameter, nothing beyond the 20 ms the march already cost.
    out.dRdn0(g) = adj.nStar.front()(k);
  }
  out.marchSeconds = seconds(t0, t1);
  out.quadratureSeconds = seconds(t1, t2);
  const int pieces = opts.subIntervals + 2 * opts.endRefinements;
  out.solves = static_cast<int>(dts.size()) * pieces * opts.gaussPoints * 2;
  return out;
}

// Central difference on lambda_i with the chain rebuilt around the perturbation. setDecay()
// derives lambda from the half-life and ignores any decayConstant passed in, so the half-life
// is what moves: lambda(1+e) means T(1/(1+e)).
//
// The keep set is held fixed. A relative step of 1e-6 cannot change the topology, and
// re-deriving the closure would compare two different index spaces.
double finiteDifference(const NuclearData& data, const cram::DepletionChain& base, const Problem& p,
                        const Zai& zai, double t, double rel, bool moveWeight) {
  const cram::Zai cz = nusift::toCram(zai);
  const cram::DecayData* d = base.decay(cz);
  if (d == nullptr || d->decayConstant <= 0.0 || d->halfLife <= 0.0) {
    return 0.0;
  }
  const double lambda = d->decayConstant;
  const double h = lambda * rel;
  const int n = base.size();
  const int global = data.indexOf(zai);
  const auto at = std::lower_bound(p.keep.begin(), p.keep.end(), global);
  const bool inKeep = at != p.keep.end() && *at == global;
  const int local = static_cast<int>(at - p.keep.begin());

  double side[2] = {0.0, 0.0};
  for (int s = 0; s < 2; ++s) {
    const double lambdaPrime = lambda + (s == 0 ? -h : h);
    cram::DepletionChain c = base;
    cram::DecayData nd = *d;
    nd.halfLife = cram::kLn2 / lambdaPrime;
    c.setDecay(cz, nd);

    int dropped = 0;
    const Eigen::SparseMatrix<double> a = restrictTo(c.decayMatrix(&dropped), p.keep, n);
    // R = <w, n(T)> is cram's definition with w HELD FIXED, so its adjoint returns only the
    // implicit term dR/dlambda through A. For an activity response w IS lambda, so the total
    // derivative carries a second, explicit term. moveWeight selects which of the two
    // derivatives this measures.
    //
    // p.n0 is reused on both sides, which makes this the derivative AT FIXED ATOM INVENTORY.
    // That is the right one here: the seed is fissions * Y, a count of atoms, so it does not
    // move with lambda. It would NOT be right for an inventory specified in activity, where
    // n0 = A0/lambda and a third term dR/dn0 * (-A0/lambda^2) appears.
    Eigen::VectorXd w = p.w;
    if (inKeep && moveWeight) {
      w(local) = lambdaPrime;
    }
    side[s] = responseAt(a, p.n0, w, t);
  }
  return (side[1] - side[0]) / (2.0 * h);
}

}  // namespace

int main(int argc, char** argv) {
  const int logIntervals = argc > 1 ? std::atoi(argv[1]) : 8;

  const NuclearData data = NuclearData::open(NUSIFT_COMMITTED_STORE);
  const cram::DepletionChain& chain = nusift::chainOf(data);

  nusift::seed::FissionSeed fs;
  fs.fissile = Zai{92, 235, 0};
  fs.fissions = 2.9e24;  // about 20 kt
  const nusift::Inventory inv = nusift::seed::seedFromFission(data, fs);

  const Problem p = build(data, chain, inv);
  const double t = 30.0 * kDay;

  std::printf("chain %d nuclides, pruned to %d; response = total activity at 30 d\n", chain.size(),
              static_cast<int>(p.keep.size()));

  const auto b0 = std::chrono::steady_clock::now();
  const double r0 = responseAt(p.a, p.n0, p.w, t);
  const auto b1 = std::chrono::steady_clock::now();
  std::printf("R = %.6e Bq   (one forward solve: %.4f s)\n\n", r0, seconds(b0, b1));

  std::printf("COST AND SCHEDULE\n");
  const cram::SensitivityOptions defaults;
  struct Case {
    const char* label;
    std::vector<double> dts;
    cram::SensitivityOptions opts;
  };
  cram::SensitivityOptions refined = defaults;
  refined.endRefinements = 12;

  std::vector<Case> cases;
  cases.push_back({"single interval [0,T]", {t}, defaults});
  cases.push_back({"single, endRefinements 12", {t}, refined});
  cases.push_back({"log grid", logSchedule(kDay, t, logIntervals), defaults});

  std::vector<SensitivityRun> runs;
  for (const Case& c : cases) {
    const SensitivityRun r = run(chain, p, c.dts, c.opts);
    std::printf("  %-28s %3d interval(s)  march %7.3f s  quadrature %8.3f s  ~%6d solves\n",
                c.label, static_cast<int>(c.dts.size()), r.marchSeconds, r.quadratureSeconds,
                r.solves);
    runs.push_back(r);
  }
  std::printf("\n");

  // Two derivatives, not one. cram's R = <w, n(T)> holds w fixed, so its adjoint returns the
  // IMPLICIT term: how lambda moves the trajectory through A. But for an activity response the
  // weight IS lambda, so the physical derivative carries a second, EXPLICIT term
  // (dw_i/dlambda_i) * n_i(T) -- which for w = lambda is just n_i(T), and in general is
  // (w_i/lambda_i) * n_i(T), so an exposure weight works the same way.
  //
  // Checking both separately is the point: the first says whether cram is right, the second
  // says what NuSIFT would have to add on top of it.
  const Zai probes[] = {Zai{55, 137, 0}, Zai{56, 137, 1}, Zai{57, 140, 0},
                        Zai{40, 95, 0},  Zai{53, 131, 0}, Zai{38, 90, 0}};

  std::printf("IMPLICIT TERM ONLY -- cram's definition, w held fixed\n");
  std::printf("  %-9s %15s %15s", "nuclide", "FD (w fixed)", "adjoint");
  for (const Case& c : cases) {
    std::printf(" %11s", c.label[0] == 's' ? (c.opts.endRefinements == 12 ? "single/12" : "single")
                                           : "log grid");
  }
  std::printf("\n");
  for (const Zai& z : probes) {
    const int gi = data.indexOf(z);
    if (gi < 0) {
      continue;
    }
    const double fd = finiteDifference(data, chain, p, z, t, 1e-6, /*moveWeight=*/false);
    std::printf("  %-9s %15.6e %15.6e", nusift::formatNuclideName(z).c_str(), fd,
                runs.front().dRdLambda(gi));
    for (const SensitivityRun& r : runs) {
      const double v = r.dRdLambda(gi);
      std::printf(" %11.2e", fd == 0.0 ? 0.0 : std::abs(v - fd) / std::abs(fd));
    }
    std::printf("\n");
  }
  std::printf("\n");

  // The accurate run, not the log grid: see the table above. Both terms are printed because
  // for a nuclide held in secular equilibrium they cancel almost exactly -- its activity is
  // set by its parent feed rate, so moving its own lambda barely moves R -- and a relative
  // error against that near-zero difference is meaningless. Elasticity against R is the
  // quantity that stays interpretable.
  std::printf(
      "TOTAL DERIVATIVE at fixed atom inventory -- implicit + (w_i/lambda_i) * n_i(T), "
      "from single/12\n");
  std::printf("  %-9s %14s %14s %14s %14s %10s\n", "nuclide", "implicit", "explicit", "total",
              "FD (w moves)", "elasticity");
  for (const Zai& z : probes) {
    const int gi = data.indexOf(z);
    if (gi < 0) {
      continue;
    }
    const double fd = finiteDifference(data, chain, p, z, t, 1e-6, /*moveWeight=*/true);
    const SensitivityRun& r = runs[1];
    const double implicitTerm = r.dRdLambda(gi);
    const double explicitTerm = r.nT(gi);  // (w_i/lambda_i) * n_i(T); w = lambda for activity
    const double total = implicitTerm + explicitTerm;
    std::printf("  %-9s %14.5e %14.5e %14.5e %14.5e %10.3f\n", nusift::formatNuclideName(z).c_str(),
                implicitTerm, explicitTerm, total, fd, data.decayConstant(gi) * total / r0);
  }
  std::printf("\n");

  // --- Which parameter class actually carries the uncertainty --------------
  //
  // Elasticities are dimensionless, so a half-life elasticity and a yield elasticity are
  // directly comparable. The yield one costs nothing extra: the seed is n0_i = fissions * Y_i,
  // so Y_i and n0_i are proportional and their elasticities are IDENTICAL -- and dR/dn0 is
  // adj.nStar.front(), which the forward+adjoint march already produced in 20 ms. The
  // quadrature that dR/dlambda needs is 30x that and buys the other class.
  //
  // Combining a class in quadrature under one representative relative uncertainty gives
  // sigma_R/R = sigma_rel * sqrt(sum e_i^2), so the RSS of the elasticities is the part that
  // depends on the physics and sigma_rel is the part that depends on the evaluation. The store
  // carries no uncertainties yet, which is exactly why this reports an RSS to be multiplied
  // rather than a sigma.
  //
  // And an RSS is where it stops even when they land. That expression is g^T Sigma g with
  // Sigma taken DIAGONAL and every relative uncertainty taken EQUAL, and neither holds: yields
  // are correlated through the mass and charge balance they were fitted under -- "independent
  // yield" is a yield type, pre-decay, not a statistical claim -- and branchings normalized to
  // sum to one are anti-correlated by construction. What follows compares sensitivity NORMS
  // between two parameter classes, which is a real and useful thing to measure. It is not an
  // error budget, and no multiplication by a sigma turns it into one.
  const SensitivityRun& acc = runs[1];  // single/12, the accurate one
  std::vector<double> n0Full(static_cast<std::size_t>(chain.size()), 0.0);
  for (int k = 0; k < static_cast<int>(p.keep.size()); ++k) {
    n0Full[static_cast<std::size_t>(p.keep[static_cast<std::size_t>(k)])] = p.n0(k);
  }

  std::vector<std::pair<double, int>> byLambda;
  std::vector<std::pair<double, int>> byYield;
  double rssLambda = 0.0;
  double rssYield = 0.0;
  for (int i = 0; i < chain.size(); ++i) {
    const double lambda = data.decayConstant(i);
    // total derivative: implicit + explicit, the explicit term being n_i(T) for w = lambda
    const double total = lambda > 0.0 ? acc.dRdLambda(i) + acc.nT(i) : 0.0;
    const double eL = lambda > 0.0 ? lambda * total / r0 : 0.0;
    const double eY = acc.dRdn0(i) * n0Full[static_cast<std::size_t>(i)] / r0;
    if (eL != 0.0) {
      byLambda.emplace_back(std::abs(eL), i);
      rssLambda += eL * eL;
    }
    if (eY != 0.0) {
      byYield.emplace_back(std::abs(eY), i);
      rssYield += eY * eY;
    }
  }
  // For a LINEAR system the seed elasticities are not merely comparable, they PARTITION the
  // response: R = <w, exp(AT) n0> = sum_i n0_i <w, exp(AT) e_i>, and dR/dn0_i is exactly that
  // inner product. So sum_i n0_i * dR/dn0_i == R and the elasticities sum to 1. If that holds
  // numerically, seed attribution is an exact decomposition with a total and a covered
  // fraction -- the same shape as every other ranking in the tool -- not a sensitivity that
  // needs uncertainty data to mean anything.
  double partition = 0.0;
  for (int i = 0; i < chain.size(); ++i) {
    partition += acc.dRdn0(i) * n0Full[static_cast<std::size_t>(i)] / r0;
  }
  std::printf("seed elasticities sum to %.12f (exact partition of R would be 1)\n\n", partition);

  std::sort(byLambda.begin(), byLambda.end(), std::greater<>());
  std::sort(byYield.begin(), byYield.end(), std::greater<>());
  rssLambda = std::sqrt(rssLambda);
  rssYield = std::sqrt(rssYield);

  std::printf("ELASTICITY BY PARAMETER CLASS (dimensionless, %% per %%)\n");
  std::printf("  %-24s %-24s\n", "half-life (lambda)", "independent yield (Y)");
  for (std::size_t k = 0; k < 8 && k < byLambda.size() && k < byYield.size(); ++k) {
    const int li = byLambda[k].second;
    const int yi = byYield[k].second;
    const double lambda = data.decayConstant(li);
    const double eL = lambda * (acc.dRdLambda(li) + acc.nT(li)) / r0;
    const double eY = acc.dRdn0(yi) * n0Full[static_cast<std::size_t>(yi)] / r0;
    std::printf("  %-9s %13.4f   %-9s %13.4f\n", nusift::formatNuclideName(data.zaiAt(li)).c_str(),
                eL, nusift::formatNuclideName(data.zaiAt(yi)).c_str(), eY);
  }
  std::printf("\n  RSS over %4d nuclides: lambda %8.4f   |   over %4d: yield %8.4f\n",
              static_cast<int>(byLambda.size()), rssLambda, static_cast<int>(byYield.size()),
              rssYield);
  std::printf("  ratio yield/lambda (equal relative uncertainty): %.2f\n\n",
              rssLambda == 0.0 ? 0.0 : rssYield / rssLambda);

  // Illustrative only -- the store carries no uncertainties, so these are stand-ins chosen to
  // show how the two classes combine, not values read from an evaluation.
  const double sigmaHalfLife = 0.005;  // 0.5%, typical for a well-studied fission product
  const double sigmaYield = 0.20;      // 20%, typical for an independent yield
  std::printf("  with illustrative sigma_T = %.1f%% and sigma_Y = %.0f%%:\n", sigmaHalfLife * 100.0,
              sigmaYield * 100.0);
  std::printf("    sigma_R/R from half-lives %.4f%%   from yields %.4f%%   yields dominate %.0fx\n",
              sigmaHalfLife * rssLambda * 100.0, sigmaYield * rssYield * 100.0,
              (sigmaHalfLife * rssLambda) == 0.0
                  ? 0.0
                  : (sigmaYield * rssYield) / (sigmaHalfLife * rssLambda));
  std::printf("\n");

  return 0;
}
