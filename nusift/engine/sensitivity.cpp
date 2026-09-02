#include "nusift/engine/sensitivity.hpp"

#include <Eigen/SparseCore>
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include <optional>
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

// The ENDF RTYP, spelled. Only the leading digit is named: a multi-step mode like 1.5 is a
// beta-minus followed by an alpha, and naming the first step is what identifies the branch a
// reader is looking at without inventing a vocabulary for every sequence.
std::string modeName(double rtyp) {
  switch (static_cast<int>(rtyp)) {
    case 0:
      return "gamma";
    case 1:
      return "beta-";
    case 2:
      return "beta+/EC";
    case 3:
      return "IT";
    case 4:
      return "alpha";
    case 5:
      return "neutron";
    case 6:
      return "SF";
    case 7:
      return "proton";
    default:
      break;
  }
  return "rtyp " + std::to_string(rtyp);
}

// The constrained covariance of one nuclide's branchings, contracted against its sensitivities.
//
//     C = D - (D u u^T D) / (u^T D u)
//
// so e^T C e = sum(e_i^2 d_i) - (sum(e_i d_i))^2 / sum(d_i) with d_i = sigma_i^2. Written in
// that reduced form rather than by forming C: the block is small, but the identity makes the
// two limits visible -- one mode gives exactly zero, and two modes give (e_1 - e_2)^2 times the
// harmonic-ish combination the constraint forces.
double constrainedVariance(const std::vector<double>& sensitivities,
                           const std::vector<double>& sigmas) {
  double weighted = 0.0;  // sum e_i^2 d_i
  double crossed = 0.0;   // sum e_i d_i
  double total = 0.0;     // sum d_i
  for (std::size_t i = 0; i < sensitivities.size(); ++i) {
    const double d = sigmas[i] * sigmas[i];
    weighted += sensitivities[i] * sensitivities[i] * d;
    crossed += sensitivities[i] * d;
    total += d;
  }
  if (!(total > 0.0)) {
    return 0.0;
  }
  // Clamped at zero: the expression is a variance and cannot be negative, but it is a difference
  // of two numbers much larger than the result whenever the sensitivities are close -- which the
  // constraint makes the common case, since equal sensitivities give exactly zero. Two or three
  // digits go there, which is immaterial to a variance quoted to two and is the same cancellation
  // the lambda terms show, met again in the algebra that makes the block cheap.
  return std::max(0.0, weighted - crossed * crossed / total);
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

  // --- branchings ------------------------------------------------------------
  //
  // cram contracts S against dA/dlambda for us and has no counterpart for branchings, but the
  // pattern is simpler: A(d, j) carries +lambda_j * b_m for mode m of nuclide j, and the
  // diagonal A(j, j) = -lambda_j does NOT move with b, because the modes sum to one and the
  // total removal rate is lambda whatever the split. So dA/db_m is a single entry and
  //
  //     dR/db_m = lambda_j * sum over intervals of S_k(daughter, j)
  //
  // The daughter comes from cram's own decayDaughter(), which is public precisely so a caller
  // contracting over the decay topology sees the same product the matrix was built with.
  Eigen::SparseMatrix<double> summed(chainSize, chainSize);
  for (const Eigen::SparseMatrix<double>& piece : scattered) {
    summed += piece;
  }

  double branchingVariance = 0.0;
  double branchingDiagonal = 0.0;
  for (int k = 0; k < kept; ++k) {
    const int index = prepared.keep[static_cast<std::size_t>(k)];
    const double lambda = data.decayConstant(index);
    if (!(lambda > 0.0)) {
      continue;
    }
    const Zai parent = Zai::fromKey(prepared.keys[static_cast<std::size_t>(k)]);
    const cram::DecayData* decay = chain.decay(toCram(parent));
    if (decay == nullptr || decay->modes.size() < 1) {
      continue;
    }

    BranchingBlock block;
    block.parentKey = parent.key();
    block.parent = formatNuclideName(parent);
    block.modes = static_cast<int>(decay->modes.size());

    std::vector<double> sensitivities;
    std::vector<double> sigmas;
    for (std::size_t m = 0; m < decay->modes.size(); ++m) {
      const cram::DecayMode& mode = decay->modes[m];
      BranchingSensitivity one;
      one.parentKey = parent.key();
      one.parent = block.parent;
      one.mode = modeName(mode.rtyp);
      one.branching = mode.branching;
      one.sigma = data.modeBranchingUncertainty(index, static_cast<int>(m));

      // Spontaneous fission produces from the yield table rather than a single daughter, so its
      // branching moves a whole column. Left out of this contraction rather than approximated by
      // one product: it belongs with the yield covariance, which is the next piece of work.
      const std::optional<cram::Zai> daughter =
          cram::DepletionChain::decayDaughter(toCram(parent), mode);
      if (!daughter.has_value()) {
        one.daughter = mode.isFission ? "fission products" : "(none tracked)";
        result.branchings.push_back(std::move(one));
        continue;
      }
      const Zai product = fromCram(*daughter);
      one.daughter = formatNuclideName(product);
      const int productIndex = data.indexOfKey(product.key());
      if (productIndex >= 0) {
        one.total = lambda * summed.coeff(productIndex, index);
        one.elasticity = result.response != 0.0 ? one.branching * one.total / result.response : 0.0;
      }

      if (one.sigma > 0.0) {
        sensitivities.push_back(one.total);
        sigmas.push_back(one.sigma);
      }
      result.branchings.push_back(std::move(one));
    }

    if (!sensitivities.empty()) {
      for (std::size_t i = 0; i < sensitivities.size(); ++i) {
        const double term = sensitivities[i] * sigmas[i];
        block.diagonalVariance += term * term;
      }
      block.constrainedVariance = constrainedVariance(sensitivities, sigmas);
      branchingVariance += block.constrainedVariance;
      branchingDiagonal += block.diagonalVariance;
      result.branchingBlocks.push_back(std::move(block));
    }
  }

  result.branchingNorm =
      result.response != 0.0 ? std::sqrt(branchingVariance) / std::abs(result.response) : 0.0;
  result.branchingNormDiagonal =
      result.response != 0.0 ? std::sqrt(branchingDiagonal) / std::abs(result.response) : 0.0;
  std::sort(result.branchingBlocks.begin(), result.branchingBlocks.end(),
            [](const BranchingBlock& a, const BranchingBlock& b) {
              if (a.constrainedVariance != b.constrainedVariance) {
                return a.constrainedVariance > b.constrainedVariance;
              }
              return a.parentKey < b.parentKey;
            });

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
