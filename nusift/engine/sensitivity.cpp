#include "nusift/engine/sensitivity.hpp"

#include <Eigen/SparseCore>
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "cram/adjoint.hpp"
#include "cram/chain.hpp"
#include "nusift/core/error.hpp"
#include "nusift/core/nuclide.hpp"
#include "nusift/core/nuclide_name.hpp"
#include "nusift/engine/decay_engine_internal.hpp"
#include "nusift/nucdata/nuclear_data.hpp"
#include "nusift/nucdata/nuclear_data_internal.hpp"

// The fourth translation unit that sees cram and Eigen, and the second to include
// cram/adjoint.hpp. It reaches further into cram than any other: rateSensitivities() and
// decayConstantSensitivities() have no counterpart on the forward path.

namespace nusift {
namespace {

constexpr const char* kModule = "sensitivity";

cram::CramOrder toCramOrder(CramOrder order) {
  return order == CramOrder::Order16 ? cram::CramOrder::CRAM16 : cram::CramOrder::CRAM48;
}

// A log schedule over (0, time], as interval widths. Log rather than linear because a decay
// response varies over decades and a uniform schedule spends its intervals where nothing is
// happening. The first edge is placed a decade below the last rather than at zero, which a log
// grid cannot reach.
std::vector<double> logSchedule(double time, int intervals) {
  if (intervals <= 1) {
    return {time};
  }
  std::vector<double> widths;
  widths.reserve(static_cast<std::size_t>(intervals));
  const double last = std::log(time);
  const double first = last - std::log(10.0) * intervals;
  double previous = 0.0;
  for (int k = 0; k < intervals; ++k) {
    const double fraction = static_cast<double>(k + 1) / static_cast<double>(intervals);
    const double edge = std::exp(first + (last - first) * fraction);
    widths.push_back(edge - previous);
    previous = edge;
  }
  return widths;
}

}  // namespace

int chooseEndRefinements(double intervalSeconds, int subIntervals, double shortestRemovalSeconds,
                         int cap) {
  if (!(intervalSeconds > 0.0) || !(shortestRemovalSeconds > 0.0) || subIntervals < 1) {
    return 0;
  }
  // cram's rule, stated in its own header: the smallest piece is h / 2^k with h the uniform
  // piece length, and k is raised until that falls below the shortest removal time in the
  // problem. Solved rather than searched, since it is a logarithm.
  const double h = intervalSeconds / subIntervals;
  if (h <= shortestRemovalSeconds) {
    return 0;
  }
  const double needed = std::log2(h / shortestRemovalSeconds);
  return std::min(cap, static_cast<int>(std::ceil(needed)));
}

DecaySensitivities decaySensitivities(const NuclearData& data, const Inventory& inventory,
                                      double time, const ResponseSpec& spec,
                                      const SensitivityOptions& sensitivity,
                                      const DecayOptions& options) {
  if (!std::isfinite(time) || !(time > 0.0)) {
    throw InputError(tagged(kModule,
                            "a decay-constant sensitivity needs a positive time; at t = 0 the "
                            "response is the seed and no half-life has had a chance to act on it"));
  }
  if (!unitSuitsDomain(spec.unit, Domain::Instant)) {
    throw InputError(tagged(kModule, std::string(unitName(spec.unit)) +
                                         " describes a total accrued over a window; this is the "
                                         "derivative of an INSTANTANEOUS response"));
  }
  if (spec.aggregate == Aggregate::GammaLine) {
    throw InputError(tagged(kModule, "a photon line has no decay constant of its own"));
  }
  if (sensitivity.scheduleIntervals < 1 || sensitivity.gaussPoints < 1 ||
      sensitivity.gaussPoints > 5 || sensitivity.subIntervals < 1) {
    throw InputError(tagged(kModule, "the quadrature options are outside what cram accepts"));
  }
  // cram's own hard limit, checked here so the message names the option the caller set rather
  // than surfacing from inside a quadrature they never asked about.
  if (sensitivity.maxEndRefinements < 0 || sensitivity.maxEndRefinements > 30) {
    throw InputError(tagged(kModule,
                            "the refinement ceiling is 0..30, which is cram's limit on end "
                            "refinements; got " +
                                std::to_string(sensitivity.maxEndRefinements)));
  }

  const engine_internal::Prepared prepared = engine_internal::prepare(data, inventory, options);
  const cram::DepletionChain& chain = chainOf(data);
  const int chainSize = chain.size();
  const int kept = static_cast<int>(prepared.keep.size());

  const std::vector<double> weights = responseWeights(data, spec);
  const std::vector<double> weightDerivatives = weightDecayDerivatives(data, spec);

  // The response weight over the pruned space, and the shortest removal time in it. The removal
  // time is 1/lambda -- the e-folding of the boundary layer the refinement has to resolve --
  // taken over the nuclides actually retained, because a microsecond isomer the seed cannot
  // reach imposes no boundary layer on this problem.
  Eigen::VectorXd w = Eigen::VectorXd::Zero(kept);
  double fastest = 0.0;
  for (int k = 0; k < kept; ++k) {
    const int index = prepared.keep[static_cast<std::size_t>(k)];
    w(k) = weights[static_cast<std::size_t>(index)];
    fastest = std::max(fastest, data.decayConstant(index));
  }

  DecaySensitivities result;
  result.time = time;
  result.shortestRemovalSeconds = fastest > 0.0 ? 1.0 / fastest : 0.0;

  const std::vector<double> widths = logSchedule(time, sensitivity.scheduleIntervals);
  const double shortestInterval = *std::min_element(widths.begin(), widths.end());
  const int refinements =
      chooseEndRefinements(shortestInterval, sensitivity.subIntervals,
                           result.shortestRemovalSeconds, sensitivity.maxEndRefinements);
  result.endRefinements = refinements;
  // The cap binding means the smallest piece is still above the shortest removal time, which is
  // exactly the condition cram warns produces a silently wrong answer.
  result.refinementCapped =
      refinements >= sensitivity.maxEndRefinements &&
      chooseEndRefinements(shortestInterval, sensitivity.subIntervals,
                           result.shortestRemovalSeconds,
                           sensitivity.maxEndRefinements + 1) > sensitivity.maxEndRefinements;

  cram::SensitivityOptions cramOptions;
  cramOptions.gaussPoints = sensitivity.gaussPoints;
  cramOptions.subIntervals = sensitivity.subIntervals;
  cramOptions.endRefinements = refinements;

  const std::vector<Eigen::SparseMatrix<double>> matrices(widths.size(), prepared.reduced);
  const cram::DepletionResult forward =
      cram::depleteLinear(matrices, widths, prepared.seed, toCramOrder(options.order));
  const cram::AdjointResult adjoint =
      cram::adjointDeplete(matrices, widths, w, toCramOrder(options.order));
  const std::vector<Eigen::SparseMatrix<double>> reduced = cram::rateSensitivities(
      forward, adjoint, matrices, widths, cramOptions, toCramOrder(options.order));

  // decayConstantSensitivities() contracts against a dA/dlambda derived from the WHOLE chain, so
  // the per-interval matrices are scattered back into chain space first. Entries outside the
  // closure contract to zero, which is the right answer rather than an accident: a nuclide the
  // seed cannot reach has n == 0 for all time.
  std::vector<Eigen::SparseMatrix<double>> scattered;
  scattered.reserve(reduced.size());
  for (const Eigen::SparseMatrix<double>& piece : reduced) {
    std::vector<Eigen::Triplet<double>> triplets;
    for (int col = 0; col < piece.outerSize(); ++col) {
      for (Eigen::SparseMatrix<double>::InnerIterator it(piece, col); it; ++it) {
        triplets.emplace_back(prepared.keep[static_cast<std::size_t>(it.row())],
                              prepared.keep[static_cast<std::size_t>(it.col())], it.value());
      }
    }
    Eigen::SparseMatrix<double> full(chainSize, chainSize);
    full.setFromTriplets(triplets.begin(), triplets.end());
    scattered.push_back(std::move(full));
  }
  const Eigen::VectorXd implicitTerm = cram::decayConstantSensitivities(chain, scattered);

  const Eigen::VectorXd& atomsAtTime = forward.n.back();
  const Eigen::VectorXd& seedImportance = adjoint.nStar.front();
  result.response = w.dot(atomsAtTime);

  const int pieces = sensitivity.subIntervals + 2 * refinements;
  result.solves = static_cast<int>(widths.size()) * pieces * sensitivity.gaussPoints * 2;

  double variance = 0.0;
  double covered = 0.0;
  double totalMagnitude = 0.0;

  for (int k = 0; k < kept; ++k) {
    const int index = prepared.keep[static_cast<std::size_t>(k)];
    const double lambda = data.decayConstant(index);
    if (!(lambda > 0.0)) {
      continue;
    }

    DecaySensitivity one;
    one.key = prepared.keys[static_cast<std::size_t>(k)];
    one.label = formatNuclideName(Zai::fromKey(one.key));
    one.decayConstant = lambda;
    one.halfLifeSeconds = data.halfLifeSeconds(index);

    one.implicit = implicitTerm(index);
    // The weight moving with lambda. Zero for a pack whose coefficient multiplies atoms or mass
    // rather than an activity, which weightDecayDerivatives() decides rather than assuming.
    one.explicitWeight = weightDerivatives[static_cast<std::size_t>(index)] * atomsAtTime(k);
    // The seed moving with lambda, for whatever part of this row was written as an activity:
    // n0 = A0/lambda, so dn0/dlambda = -n0/lambda.
    const double fromActivity = [&] {
      for (const InventoryEntry& entry : inventory.entries()) {
        if (entry.zaiKey == one.key) {
          return entry.atomsFromActivity;
        }
      }
      return 0.0;
    }();
    one.basis = seedImportance(k) * (-fromActivity / lambda);

    one.total = one.implicit + one.explicitWeight + one.basis;
    one.elasticity = result.response != 0.0 ? lambda * one.total / result.response : 0.0;

    const double sigma = data.halfLifeUncertainty(index);
    one.relativeUncertainty = one.halfLifeSeconds > 0.0 ? sigma / one.halfLifeSeconds : 0.0;
    one.sigmaContribution = std::abs(one.elasticity) * one.relativeUncertainty;

    totalMagnitude += std::abs(one.elasticity);
    if (one.relativeUncertainty > 0.0) {
      ++result.withUncertainty;
      covered += std::abs(one.elasticity);
      variance += one.sigmaContribution * one.sigmaContribution;
    } else {
      ++result.withoutUncertainty;
    }
    result.nuclides.push_back(std::move(one));
  }

  result.relativeNorm = std::sqrt(variance);
  result.coveredFraction = totalMagnitude > 0.0 ? covered / totalMagnitude : 0.0;

  const bool anyUncertainty = result.withUncertainty > 0;
  std::sort(result.nuclides.begin(), result.nuclides.end(),
            [anyUncertainty](const DecaySensitivity& a, const DecaySensitivity& b) {
              const double left = anyUncertainty ? a.sigmaContribution : std::abs(a.elasticity);
              const double right = anyUncertainty ? b.sigmaContribution : std::abs(b.elasticity);
              if (left != right) {
                return left > right;
              }
              return a.key < b.key;
            });
  return result;
}

}  // namespace nusift
