#include "nusift/triage/uncertainty.hpp"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "nusift/core/error.hpp"
#include "nusift/core/nuclide.hpp"
#include "nusift/core/nuclide_name.hpp"
#include "nusift/engine/adjoint_engine.hpp"
#include "nusift/io/time_spec.hpp"
#include "nusift/nucdata/nuclear_data.hpp"

namespace nusift {
namespace {

constexpr const char* kModule = "uncertainty";

[[noreturn]] void fail(const std::string& what) {
  throw InputError(tagged(kModule, what));
}

// The sigma an assay stated for one nuclide, or zero. Entries are sorted by key, so this is a
// binary search rather than the scan an inline lambda would have been.
double sigmaOf(const Inventory& inventory, std::int64_t key) {
  const std::span<const InventoryEntry> entries = inventory.entries();
  const auto it = std::lower_bound(
      entries.begin(), entries.end(), key,
      [](const InventoryEntry& entry, std::int64_t k) { return entry.zaiKey < k; });
  return it != entries.end() && it->zaiKey == key ? it->sigmaAtoms : 0.0;
}

}  // namespace

ResponseUncertainty responseUncertainty(const NuclearData& data, std::span<const AssayGroup> assays,
                                        double epochSeconds, double time, const ResponseSpec& spec,
                                        const DecayOptions& options) {
  if (assays.empty()) {
    fail("no assays to propagate");
  }
  if (!(time >= 0.0)) {
    fail("a response time cannot be negative");
  }
  // The shares are instantaneous for the reason attributeToSeed's are: an integrated response
  // needs the integrated adjoint, which is a different solve and is not this one.
  if (!unitSuitsDomain(spec.unit, Domain::Instant)) {
    fail("an error bar here is on an INSTANTANEOUS response; " + std::string(unitName(spec.unit)) +
         " describes a total accrued over a window, which needs the integrated adjoint");
  }
  if (spec.aggregate == Aggregate::GammaLine) {
    fail("a photon line has no seed, so it has no seed uncertainty either");
  }

  const std::vector<double> weight = responseWeights(data, spec);

  // A row assayed as 0 +/- sigma has to reach the index space, or its stated uncertainty is
  // dropped without a word. See DecayOptions::sigmaSeedsClosure.
  DecayOptions closureOptions = options;
  closureOptions.sigmaSeedsClosure = true;

  ResponseUncertainty result;
  result.metric = spec.metric;
  result.unit = spec.unit;
  result.time = time;

  double variance = 0.0;
  double covered = 0.0;

  for (const AssayGroup& assay : assays) {
    if (!assay.dated && epochSeconds != 0.0) {
      fail("the assay \"" + assay.label + "\" states no date, so it cannot be carried to " +
           formatCalendarDate(epochSeconds) +
           ". Give its rows an `assayed` column, or ask about no epoch at all");
    }
    if (assay.dateSeconds > epochSeconds) {
      fail("the assay dated " + formatCalendarDate(assay.dateSeconds) +
           " is later than the epoch " + formatCalendarDate(epochSeconds) +
           "; carrying one backward is an inverse problem, not a solve");
    }
    const double carried = epochSeconds - assay.dateSeconds;

    // THE WHOLE TRICK, and it needs no new machinery. The importance of an assay carried
    // forward by `carried` is the ordinary adjoint run for T + carried, because the carry and
    // the response interval share one decay matrix and their exponentials commute. So this is
    // the same solve `attribute` runs, at a later time for an older sheet.
    const SeedImportance importance =
        seedImportance(data, assay.inventory, weight, time + carried, closureOptions);

    for (std::size_t i = 0; i < importance.nuclideKeys.size(); ++i) {
      const double atoms = importance.seedAtoms[i];
      // From the assay's OWN row, which is the only place a diagonal sigma is meaningful. A
      // nuclide the adjoint reaches but the sheet never listed carries none, and says so.
      const double sigma = sigmaOf(assay.inventory, importance.nuclideKeys[i]);
      // Two ways to be a seed here, and a row needs only one of them. Atoms make a share of R;
      // a stated sigma makes a term in its variance. A non-detect reported as 0 +/- MDA has the
      // second without the first, and is exactly the row whose uncertainty must not vanish --
      // dropping it would quote an error bar that silently ignores a measurement the sheet
      // made. Everything the closure merely REACHED -- a daughter the sheet never listed -- has
      // neither, and is not a seed at all.
      if (!(atoms > 0.0) && !(sigma > 0.0)) {
        continue;
      }
      SeedUncertainty seed;
      seed.key = importance.nuclideKeys[i];
      seed.label = formatNuclideName(Zai::fromKey(seed.key));
      seed.assay = assay.label;
      seed.carriedSeconds = carried;
      seed.seedAtoms = atoms;
      seed.importance = importance.importance[i];
      seed.share = seed.importance * atoms;

      seed.sigmaAtoms = sigma;
      seed.sigmaContribution = std::abs(seed.importance) * sigma;
      variance += seed.sigmaContribution * seed.sigmaContribution;

      result.response += seed.share;
      if (sigma > 0.0) {
        ++result.rowsWithSigma;
        covered += seed.share;
      } else {
        ++result.rowsWithoutSigma;
      }
      result.seeds.push_back(std::move(seed));
    }
  }

  result.sigma = std::sqrt(variance);
  result.relative = result.response > 0.0 ? result.sigma / result.response : 0.0;
  result.coveredFraction = result.response > 0.0 ? covered / result.response : 0.0;

  for (SeedUncertainty& seed : result.seeds) {
    seed.varianceFraction =
        variance > 0.0 ? (seed.sigmaContribution * seed.sigmaContribution) / variance : 0.0;
  }
  // Ranked by what they contribute to the VARIANCE, not to the answer: the question this
  // ordering exists to answer is which measurement to improve, and a row with a negligible
  // share of R can dominate its error bar. Ties break by key, as every ordering here does.
  std::sort(result.seeds.begin(), result.seeds.end(),
            [](const SeedUncertainty& a, const SeedUncertainty& b) {
              if (a.varianceFraction != b.varianceFraction) {
                return a.varianceFraction > b.varianceFraction;
              }
              return a.key < b.key;
            });
  return result;
}

}  // namespace nusift
