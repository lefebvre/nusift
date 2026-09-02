#include "nusift/engine/reconcile.hpp"

#include <algorithm>
#include <vector>

#include "nusift/core/error.hpp"
#include "nusift/core/nuclide.hpp"
#include "nusift/io/time_spec.hpp"

namespace nusift {
namespace {

constexpr const char* kModule = "reconcile";

}  // namespace

double latestAssayDate(std::span<const AssayGroup> groups) {
  if (groups.empty()) {
    throw InputError(tagged(kModule, "no assays to reconcile"));
  }
  double latest = groups.front().dateSeconds;
  for (const AssayGroup& group : groups) {
    latest = std::max(latest, group.dateSeconds);
  }
  return latest;
}

Reconciliation reconcile(const NuclearData& data, std::span<const AssayGroup> groups,
                         double epochSeconds, const DecayOptions& options) {
  if (groups.empty()) {
    throw InputError(tagged(kModule, "no assays to reconcile"));
  }

  Reconciliation result;
  result.epochSeconds = epochSeconds;

  std::vector<const AssayGroup*> ordered;
  ordered.reserve(groups.size());
  for (const AssayGroup& group : groups) {
    if (group.inventory.empty()) {
      throw InputError(tagged(kModule, "the assay dated " + formatCalendarDate(group.dateSeconds) +
                                           " has no rows, so there is nothing to carry forward"));
    }
    // The refusal the header argues for, named per assay rather than once for the set: which
    // assay is in the future of the epoch is the thing the caller has to fix.
    if (group.dateSeconds > epochSeconds) {
      throw InputError(tagged(
          kModule, "the assay dated " + formatCalendarDate(group.dateSeconds) +
                       " is later than the epoch " + formatCalendarDate(epochSeconds) +
                       ", and carrying an assay BACKWARD is not a decay solve but an inverse "
                       "problem: daughters measured at assay cannot be told from daughters grown "
                       "in since, so no unique earlier composition exists. Reconcile to " +
                       formatCalendarDate(latestAssayDate(groups)) + " or later"));
    }
    ordered.push_back(&group);
  }
  std::sort(ordered.begin(), ordered.end(), [](const AssayGroup* a, const AssayGroup* b) {
    return a->dateSeconds < b->dateSeconds;
  });

  result.spanSeconds = ordered.back()->dateSeconds - ordered.front()->dateSeconds;

  for (const AssayGroup* group : ordered) {
    AssayContribution contribution;
    contribution.label = group->label;
    contribution.dateSeconds = group->dateSeconds;
    contribution.carriedSeconds = epochSeconds - group->dateSeconds;
    contribution.nuclides = group->inventory.size();
    contribution.atomsAtAssay = group->inventory.totalAtoms();

    if (contribution.carriedSeconds == 0.0) {
      // Already at the epoch. Not solved for, because a solve over a zero interval is the
      // identity and paying for one would only add floating-point noise to numbers a user
      // supplied exactly.
      for (const InventoryEntry& entry : group->inventory.entries()) {
        result.inventory.addKey(entry.zaiKey, entry.atoms);
      }
      contribution.atomsAtEpoch = contribution.atomsAtAssay;
    } else {
      const DecayResult carried =
          decay(data, group->inventory, std::vector<double>{contribution.carriedSeconds}, options);
      const std::span<const double> atoms = carried.atomsAt(0);
      for (std::size_t i = 0; i < atoms.size(); ++i) {
        if (atoms[i] > 0.0) {
          result.inventory.addKey(carried.nuclideKeys[i], atoms[i]);
          contribution.atomsAtEpoch += atoms[i];
        }
      }
    }
    result.contributions.push_back(std::move(contribution));
  }

  return result;
}

}  // namespace nusift
