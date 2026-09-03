#include "nusift/triage/yield_uncertainty.hpp"

#include <algorithm>
#include <cmath>
#include <string>
#include <unordered_map>
#include <vector>

#include "nusift/core/error.hpp"
#include "nusift/core/nuclide.hpp"
#include "nusift/core/nuclide_name.hpp"
#include "nusift/engine/adjoint_engine.hpp"
#include "nusift/engine/inventory.hpp"
#include "nusift/nucdata/fission_yield.hpp"
#include "nusift/nucdata/nuclear_data.hpp"

namespace nusift {
namespace {

constexpr const char* kModule = "yield uncertainty";

[[noreturn]] void fail(const std::string& what) {
  throw InputError(tagged(kModule, what));
}

}  // namespace

YieldUncertainty yieldUncertainty(const NuclearData& data, const seed::FissionSeed& seed,
                                  const YieldCorrelation& correlation, double time,
                                  const ResponseSpec& spec, const DecayOptions& options) {
  if (!(time >= 0.0)) {
    fail("a response time cannot be negative");
  }
  // Both refusals are responseUncertainty's, for the same two reasons: an integrated response
  // needs the integrated adjoint, and a photon line has no seed to carry a yield.
  if (!unitSuitsDomain(spec.unit, Domain::Instant)) {
    fail("an error bar here is on an INSTANTANEOUS response; " + std::string(unitName(spec.unit)) +
         " describes a total accrued over a window, which needs the integrated adjoint");
  }
  if (spec.aggregate == Aggregate::GammaLine) {
    fail("a photon line has no seed, so it has no yield uncertainty either");
  }

  // seedFromFission does the lookup, the refusal, and the provenance; this needs the SET as
  // well, for the sigmas, and asks for it the same way so the two cannot disagree about which
  // energy was landed on.
  const FissionYieldSet* set = data.fissionYields().nearest(seed.fissile, seed.incidentEnergyEv);
  const Inventory inventory = seed::seedFromFission(data, seed);
  if (set == nullptr) {
    fail("no fission yields for " + formatNuclideName(seed.fissile));  // seedFromFission threw
  }

  std::unordered_map<std::int64_t, const FissionProduct*> byKey;
  byKey.reserve(set->products.size() * 2);
  for (const FissionProduct& product : set->products) {
    byKey.emplace(product.nuclide.key(), &product);
  }

  const std::vector<double> weight = responseWeights(data, spec);
  const SeedImportance importance = seedImportance(data, inventory, weight, time, options);

  YieldUncertainty result;
  result.metric = spec.metric;
  result.unit = spec.unit;
  result.time = time;
  result.fissions = seed.fissions;
  result.incidentEnergyEv = set->energyEv;
  result.provenance = correlation.provenance();
  result.storeLibrary = data.provenance().library;

  // d_i = |dR/dn0_i| * N_f * sigma_Y_i, the per-product standard deviation. Assembled over the
  // matched products in the matrix's own row order, so the contraction below is a dense
  // symmetric quadratic form rather than a hash lookup per pair.
  std::vector<int> matchedRows;
  std::vector<double> matchedTerms;
  matchedRows.reserve(importance.nuclideKeys.size());
  matchedTerms.reserve(importance.nuclideKeys.size());

  double diagonalVariance = 0.0;
  double matchedDiagonalVariance = 0.0;
  double covered = 0.0;
  double sigmaCovered = 0.0;
  bool anySigma = false;

  for (std::size_t i = 0; i < importance.nuclideKeys.size(); ++i) {
    const double atoms = importance.seedAtoms[i];
    if (!(atoms > 0.0)) {
      continue;
    }
    YieldContribution product;
    product.key = importance.nuclideKeys[i];
    product.label = formatNuclideName(Zai::fromKey(product.key));
    product.seedAtoms = atoms;
    product.importance = importance.importance[i];
    product.share = product.importance * atoms;

    const auto it = byKey.find(product.key);
    if (it != byKey.end()) {
      product.yield = it->second->yield;
      product.sigmaYield = std::max(it->second->yieldUncertainty, 0.0);
    }

    // The seed is n0 = N_f * Y with N_f taken exact, so the sigma on the seed is N_f * sigma_Y
    // and the whole propagation is the one already shipped with a different Sigma.
    const double term = product.importance * seed.fissions * product.sigmaYield;
    product.sigmaContribution = std::abs(term);
    diagonalVariance += term * term;

    result.response += product.share;
    if (product.sigmaYield > 0.0) {
      anySigma = true;
      ++result.productsWithSigma;
      sigmaCovered += product.share;
    } else {
      ++result.productsWithoutSigma;
    }

    const int row = correlation.indexOf(product.key);
    if (row >= 0) {
      product.correlated = true;
      ++result.productsMatched;
      covered += product.share;
      matchedRows.push_back(row);
      matchedTerms.push_back(term);
      matchedDiagonalVariance += term * term;
    } else {
      ++result.productsUnmatched;
    }
    result.products.push_back(std::move(product));
  }

  if (!anySigma) {
    fail("no staged yield for " + formatNuclideName(seed.fissile) +
         " carries an uncertainty, so there is nothing to propagate. A store staged before the "
         "uncertainty columns were carried has none at all -- check `nusift data info`");
  }
  if (matchedRows.empty()) {
    fail("\"" + correlation.provenance().path +
         "\" and this seed share no fission product. The "
         "matrix is per fissioning SYSTEM, so a " +
         formatNuclideName(seed.fissile) + " seed needs that system's own file");
  }

  // The contraction. Symmetric, so the off-diagonal is walked once and doubled -- which halves
  // the work and, more usefully, makes the diagonal and off-diagonal contributions separately
  // visible in the arithmetic.
  double correlated = matchedDiagonalVariance;
  const std::size_t n = matchedRows.size();
  for (std::size_t i = 0; i < n; ++i) {
    double offDiagonal = 0.0;
    for (std::size_t j = i + 1; j < n; ++j) {
      offDiagonal += correlation.at(static_cast<std::size_t>(matchedRows[i]),
                                    static_cast<std::size_t>(matchedRows[j])) *
                     matchedTerms[j];
    }
    correlated += 2.0 * matchedTerms[i] * offDiagonal;
  }

  result.varianceCorrelated = correlated;
  result.varianceNegative = correlated < 0.0;
  result.sigmaDiagonal = std::sqrt(diagonalVariance);
  result.sigmaDiagonalMatched = std::sqrt(matchedDiagonalVariance);
  result.sigmaCorrelated = correlated > 0.0 ? std::sqrt(correlated) : 0.0;
  result.varianceRatio = correlated > 0.0 && matchedDiagonalVariance > 0.0
                             ? correlated / matchedDiagonalVariance
                             : 0.0;
  result.smallestEigenvalue = correlation.smallestEigenvalue(matchedRows);
  result.productsUnusedInMatrix = static_cast<int>(correlation.size()) - result.productsMatched;

  if (result.response > 0.0) {
    result.relativeDiagonal = result.sigmaDiagonal / result.response;
    result.relativeDiagonalMatched = result.sigmaDiagonalMatched / result.response;
    result.relativeCorrelated = result.sigmaCorrelated / result.response;
    result.coveredFraction = covered / result.response;
    result.sigmaCoveredFraction = sigmaCovered / result.response;
  }

  for (YieldContribution& product : result.products) {
    product.varianceFraction =
        diagonalVariance > 0.0
            ? (product.sigmaContribution * product.sigmaContribution) / diagonalVariance
            : 0.0;
  }
  // Ranked by the DIAGONAL variance fraction, which is the only per-product share that exists:
  // a correlated variance does not decompose into per-product terms, since half of every
  // off-diagonal belongs to each of two products. The ordering says which yield to re-evaluate,
  // and that question is the diagonal's to answer.
  std::sort(result.products.begin(), result.products.end(),
            [](const YieldContribution& a, const YieldContribution& b) {
              if (a.varianceFraction != b.varianceFraction) {
                return a.varianceFraction > b.varianceFraction;
              }
              return a.key < b.key;
            });
  return result;
}

}  // namespace nusift
