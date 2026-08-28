#include "nusift/triage/attribution.hpp"

#include <algorithm>
#include <cmath>
#include <string>
#include <unordered_set>
#include <vector>

#include "nusift/core/error.hpp"
#include "nusift/core/nuclide_name.hpp"
#include "nusift/engine/adjoint_engine.hpp"
#include "nusift/nucdata/nuclear_data.hpp"

namespace nusift {
namespace {

constexpr const char* kModule = "seed attribution";

[[noreturn]] void fail(const std::string& what) {
  throw InputError(tagged(kModule, what));
}

void requireAttributableSpec(const ResponseSpec& spec) {
  if (!unitSuitsDomain(spec.unit, Domain::Instant)) {
    fail(std::string("unit ") + unitName(spec.unit) +
         " is a time-integrated total, and seed shares are of an instantaneous response");
  }
  if (spec.aggregate == Aggregate::GammaLine) {
    fail(
        "seed shares are attributed to seeded nuclides, and a photon line has no seed -- rank "
        "by gamma line to see which lines carry the dose, or attribute by nuclide to see which "
        "seeds carry them");
  }
  if (spec.aggregate != Aggregate::Nuclide) {
    fail(std::string("seed shares are attributed by nuclide, not by ") +
         aggregateName(spec.aggregate) +
         " -- an inventory row names a nuclide, so that is what a share can name");
  }
}

}  // namespace

std::int64_t requireSeedPin(const NuclearData& data, const Inventory& inventory,
                            std::string_view text) {
  const Zai zai = requireNuclideName(text, tagged(kModule, "pin"));
  const std::int64_t key = zai.key();
  for (const InventoryEntry& entry : inventory.entries()) {
    if (entry.zaiKey == key && entry.atoms != 0.0) {
      return key;
    }
  }
  // Refused rather than pinned to a zero row, for the reason requirePin() refuses: a pin that
  // silently resolves to nothing reads as "this seed contributes nothing", when in fact it was
  // never seeded at all. Those are different answers and only one of them is true.
  const bool known = data.indexOf(zai) >= 0;
  fail("\"" + std::string(text) + "\" names " + formatNuclideName(zai) +
       (known ? ", which this inventory does not seed -- only seeded nuclides have a share"
              : ", which this data store does not carry"));
}

SeedAttribution attributeToSeed(const NuclearData& data, const Inventory& inventory, double time,
                                const ResponseSpec& spec, const RankRequest& request,
                                const DecayOptions& options) {
  requireAttributableSpec(spec);

  // The same weight vector the forward path would build for this metric, so the two
  // attributions of one number cannot disagree about what the number is.
  const std::vector<double> weight = responseWeights(data, spec);
  const SeedImportance imp = seedImportance(data, inventory, weight, time, options);

  SeedAttribution out;
  out.metric = spec.metric;
  out.unit = spec.unit;
  out.time = time;
  out.total = imp.response;
  out.seedProvenance = inventory.provenance();
  out.buildup = spec.geometry.buildup;

  // The forward solve the caveats need, paid only for exposure because activity carries
  // neither of them. The adjoint produces the response and its decomposition but never forms
  // n(T), and both warnings are properties of what is EMITTING at `time` -- the OTHER
  // attribution of this number. Reporting an attributed exposure without them would leave the
  // same figure warned about under `rank` and silent under `attribute`.
  if (spec.metric == Metric::Exposure) {
    const double at[] = {time};
    const DecayResult forward = decay(data, inventory, at, options);
    ExposureCaveats caveats = exposureCaveats(data, forward.nuclideKeys, forward.atomsAt(0), spec);
    out.unmodeledEnergyFraction = caveats.unmodeledEnergyFraction;
    out.meanOpticalDepth = caveats.meanOpticalDepth;
    out.unmodeledContinuum = std::move(caveats.unmodeledContinuum);
  }

  const std::unordered_set<std::int64_t> pinned(request.pinned.begin(), request.pinned.end());

  // Only seeded nuclides hold a share. Every other entry of the importance vector is a real
  // derivative -- it says what an atom placed there WOULD be worth -- but it multiplies a seed
  // of zero, so it is not part of this partition.
  //
  // A seed that WAS placed and is still worth nothing -- a stable nuclide, or a pure beta
  // emitter under an exposure metric -- is separated out rather than ranked. It has no place
  // in an ordering by value, and leaving it in `all` would spend a topN slot on it and then
  // count it as an omitted contributor when the cut fell above it; both mislead, because
  // showing the row would not move the total by an atom. It is still a real answer to a
  // question someone asked explicitly, so a pin brings it back below the ranking.
  std::vector<SeedShare> all;
  std::vector<SeedShare> inert;
  all.reserve(imp.nuclideKeys.size());
  for (std::size_t i = 0; i < imp.nuclideKeys.size(); ++i) {
    const double n0 = imp.seedAtoms[i];
    if (n0 == 0.0) {
      continue;
    }
    SeedShare share;
    share.key = imp.nuclideKeys[i];
    share.label = formatNuclideName(Zai::fromKey(share.key));
    share.seedAtoms = n0;
    share.importance = imp.importance[i];
    share.value = n0 * imp.importance[i];
    if (share.value == 0.0) {
      if (pinned.count(share.key) != 0) {
        share.pinned = true;
        inert.push_back(std::move(share));
      }
      continue;
    }
    all.push_back(std::move(share));
  }

  std::stable_sort(all.begin(), all.end(),
                   [](const SeedShare& a, const SeedShare& b) { return a.value > b.value; });

  double running = 0.0;
  for (std::size_t i = 0; i < all.size(); ++i) {
    all[i].fraction = out.total != 0.0 ? all[i].value / out.total : 0.0;
    running += all[i].fraction;
    all[i].cumulativeFraction = running;
    // Every row here carries value != 0, so every row holds a place. Rank 0 means "contributes
    // nothing and so holds no place in the ordering", matching Contributor::rank, and reaches
    // the output only on a pinned inert row.
    all[i].rank = static_cast<int>(i) + 1;
  }

  std::size_t cut = all.size();
  if (request.topN > 0) {
    cut = std::min(cut, static_cast<std::size_t>(request.topN));
  }
  if (request.coverage > 0.0) {
    for (std::size_t i = 0; i < all.size(); ++i) {
      if (all[i].cumulativeFraction >= request.coverage) {
        cut = std::min(cut, i + 1);
        break;
      }
    }
  }
  if (request.minFraction > 0.0) {
    for (std::size_t i = 0; i < all.size(); ++i) {
      if (all[i].fraction < request.minFraction) {
        cut = std::min(cut, i);
        break;
      }
    }
  }

  for (std::size_t i = 0; i < cut; ++i) {
    out.shares.push_back(all[i]);
  }
  std::size_t rankedShown = out.shares.size();
  for (std::size_t i = cut; i < all.size(); ++i) {
    if (pinned.count(all[i].key) != 0) {
      SeedShare row = all[i];
      row.pinned = true;
      out.shares.push_back(std::move(row));
      ++rankedShown;
    }
  }
  // Last, after the pinned tail, so the rankless rows sit at the bottom of the table where the
  // reader has already been told what a dash in the rank column means.
  for (SeedShare& row : inert) {
    out.shares.push_back(std::move(row));
  }

  // Summed over what was RETURNED rather than taken from the last cumulative, because a pin
  // below the cut makes those two different numbers and only the sum describes what the reader
  // can actually see. Ranking::coveredFraction is defined the same way for the same reason.
  for (const SeedShare& share : out.shares) {
    out.coveredFraction += share.fraction;
  }
  // Over the ranked set alone. `inert` was never in the ordering, so a seed that contributes
  // nothing is not something the reader is missing -- see the note on omittedCount.
  out.omittedCount = static_cast<int>(all.size() - rankedShown);
  return out;
}

}  // namespace nusift
