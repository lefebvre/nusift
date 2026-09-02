#include "nusift/triage/task_plan.hpp"

#include <algorithm>
#include <string>
#include <vector>

#include "nusift/core/error.hpp"
#include "nusift/triage/events.hpp"

namespace nusift {
namespace {

constexpr const char* kModule = "task plan";

[[noreturn]] void fail(const std::string& what) {
  throw InputError(tagged(kModule, what));
}

void validate(double startSeconds, std::span<const PlanLeg> legs) {
  if (!(startSeconds >= 0.0)) {
    fail("a plan cannot begin before the inventory exists");
  }
  if (legs.empty()) {
    fail("a plan needs at least one leg");
  }
  for (const PlanLeg& leg : legs) {
    if (leg.name.empty()) {
      fail(
          "every leg needs a name; \"leg 3 costs 60% of your dose\" is not an answer anybody can "
          "act on");
    }
    if (!(leg.durationSeconds > 0.0)) {
      fail("leg \"" + leg.name + "\" has no duration, so it is not a stretch of anything");
    }
    if (!(leg.occupancy >= 0.0) || leg.occupancy > 1.0) {
      fail("leg \"" + leg.name +
           "\" has an occupancy outside [0, 1]; it is the fraction of the "
           "leg actually spent in the field");
    }
  }
}

}  // namespace

TaskPlan runTaskPlan(const NuclearData& data, const Inventory& inventory, const ResponseSpec& spec,
                     double startSeconds, std::span<const PlanLeg> legs, double budget,
                     const DecayOptions& options) {
  validate(startSeconds, legs);

  TaskPlan plan;
  plan.unit = spec.unit;
  plan.startSeconds = startSeconds;
  plan.budget = budget > 0.0 ? budget : 0.0;

  double clock = startSeconds;
  for (const PlanLeg& leg : legs) {
    LegResult result;
    result.name = leg.name;
    result.startSeconds = clock;
    result.endSeconds = clock + leg.durationSeconds;
    result.occupancy = leg.occupancy;
    result.distanceM = leg.geometry.distanceM;
    result.isBreak = !(leg.occupancy > 0.0);

    if (!result.isBreak) {
      // The geometry is this leg's, and the atoms are not: a window's atom-seconds are the same
      // however far away anyone is standing, which is why a leg at 0.8 m and a leg at 3 m cost
      // one solve each rather than one model each.
      ResponseSpec legSpec = spec;
      legSpec.aggregate = Aggregate::Nuclide;
      legSpec.geometry = leg.geometry;

      std::vector<std::int64_t> keys;
      const std::vector<double> integral =
          intervalIntegral(data, inventory, result.startSeconds, result.endSeconds, &keys, options);
      const double accrued = buildIntervalResponse(data, keys, integral, result.startSeconds,
                                                   result.endSeconds, legSpec)
                                 .totals.front();
      result.accrued = leg.occupancy * accrued;
    }

    const double exposed = leg.durationSeconds * leg.occupancy;
    result.meanRate = exposed > 0.0 ? result.accrued / exposed : 0.0;

    plan.total += result.accrued;
    plan.elapsedSeconds += leg.durationSeconds;
    plan.exposedSeconds += exposed;
    result.cumulative = plan.total;
    plan.legs.push_back(std::move(result));
    clock = plan.legs.back().endSeconds;
  }
  plan.endSeconds = clock;

  for (LegResult& leg : plan.legs) {
    leg.fraction = plan.total > 0.0 ? leg.accrued / plan.total : 0.0;
    leg.cumulativeFraction = plan.total > 0.0 ? leg.cumulative / plan.total : 0.0;
  }

  // --- where the budget runs out ----------------------------------------------
  //
  // Not merely whether the plan fits. A plan that does not fit is only actionable if it says
  // where to stop, and that is the same root-find stayTime() runs -- on this leg's own window,
  // for the budget still unspent when the leg begins.
  if (plan.budget > 0.0 && plan.total > plan.budget) {
    plan.budgetSpent = true;
    double before = 0.0;
    for (std::size_t i = 0; i < plan.legs.size(); ++i) {
      const LegResult& leg = plan.legs[i];
      if (leg.cumulative <= plan.budget) {
        before = leg.cumulative;
        continue;
      }
      plan.spentInLeg = static_cast<int>(i);

      ResponseSpec legSpec = spec;
      legSpec.aggregate = Aggregate::Nuclide;
      legSpec.geometry = legs[i].geometry;
      // The budget left when this leg starts, divided by the occupancy, because what the
      // integral has to reach is the UNOCCUPIED accrual that scales to it.
      const double remaining = (plan.budget - before) / leg.occupancy;
      const StayTime stay = stayTime(data, inventory, legSpec, leg.startSeconds, remaining,
                                     leg.endSeconds - leg.startSeconds, options);
      // Bounded by construction: this leg's full accrual already exceeds what is left, which is
      // how it was chosen. The guard is for the pathological case rather than the expected one.
      plan.spentAtSeconds = stay.bounded ? leg.startSeconds + stay.durationSeconds : leg.endSeconds;
      break;
    }
  }

  return plan;
}

}  // namespace nusift
