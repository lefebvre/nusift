// A job as a sequence of legs.
//
// The arithmetic is the interval integral applied per leg, so it is checked against the closed
// form like every other integral here. What is worth testing beyond that is the composition: the
// legs abut, the sum is the total, a leg's own geometry is the one used for it and no other, and
// a break costs nothing at all.
//
#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "nusift/core/error.hpp"
#include "nusift/engine/decay_engine.hpp"
#include "nusift/nucdata/nuclear_data.hpp"
#include "nusift/triage/task_plan.hpp"
#include "synthetic_chain.hpp"

namespace nusift {
namespace {

constexpr double kLambda = 1.0e-5;

// One decaying nuclide with a stable terminator, so the total decay rate is the parent's alone
// and every accrual has a closed form.
NuclearData oneEmitter() {
  StoreArrays arrays = synth::linearChain({kLambda});
  synth::addLines(arrays, 0, {661657.0}, {0.9});
  return NuclearData::fromArrays(std::move(arrays));
}

Inventory seeded(double atoms = 1.0e20) {
  Inventory inv;
  inv.add(Zai{50, 100, 0}, atoms);
  return inv;
}

PlanLeg leg(const char* name, double seconds, double distanceM, double occupancy = 1.0) {
  PlanLeg one;
  one.name = name;
  one.durationSeconds = seconds;
  one.occupancy = occupancy;
  one.geometry.distanceM = distanceM;
  return one;
}

PlanLeg breakFor(const char* name, double seconds) {
  PlanLeg one;
  one.name = name;
  one.durationSeconds = seconds;
  one.occupancy = 0.0;
  return one;
}

TEST(TaskPlan, EachLegIsTheExactIntegralOverItsOwnWindow) {
  const NuclearData data = oneEmitter();
  const double atoms = 1.0e20;
  const double start = 1.0e4;

  ResponseSpec spec;
  spec.unit = Unit::Decays;

  const std::vector<PlanLeg> legs = {leg("in", 600.0, 3.0), leg("work", 1800.0, 1.0)};
  const TaskPlan plan = runTaskPlan(data, seeded(atoms), spec, start, legs);

  // Decays over [t1, t2] for a single nuclide is N0 (e^{-l t1} - e^{-l t2}), exactly.
  const auto decays = [&](double t1, double t2) {
    return atoms * (std::exp(-kLambda * t1) - std::exp(-kLambda * t2));
  };
  ASSERT_EQ(plan.legs.size(), 2u);
  EXPECT_NEAR(plan.legs[0].accrued, decays(start, start + 600.0), atoms * 1.0e-12);
  EXPECT_NEAR(plan.legs[1].accrued, decays(start + 600.0, start + 2400.0), atoms * 1.0e-12);
  // The legs abut: the second begins where the first ends, with no gap the plan did not name.
  EXPECT_DOUBLE_EQ(plan.legs[1].startSeconds, plan.legs[0].endSeconds);
  EXPECT_NEAR(plan.total, decays(start, start + 2400.0), atoms * 1.0e-12)
      << "one job split in two accrues what the whole job does";
  EXPECT_DOUBLE_EQ(plan.elapsedSeconds, 2400.0);
  EXPECT_DOUBLE_EQ(plan.exposedSeconds, 2400.0);
}

// The point of legs. Activity ignores the distance entirely, so this is the test that says the
// per-leg geometry is real: the same two legs cost the same in decays and different in roentgen.
TEST(TaskPlan, UsesEachLegsOwnDistanceAndNoOthers) {
  const NuclearData data = oneEmitter();

  ResponseSpec activity;
  activity.unit = Unit::Decays;
  ResponseSpec exposure;
  exposure.metric = Metric::Exposure;
  exposure.unit = Unit::Roentgen;

  // Two one-leg plans over the SAME window, differing only in distance. Comparing two legs of
  // one plan would not isolate the distance: the second leg starts where the first ended, and
  // over that gap the source decays by lambda*dt -- which the test above is about and this one
  // is not.
  PlanLeg near = leg("near", 60.0, 1.0);
  PlanLeg far = leg("far", 60.0, 2.0);
  near.geometry.airAttenuation = false;
  far.geometry.airAttenuation = false;

  const TaskPlan nearActivity = runTaskPlan(data, seeded(), activity, 0.0, {&near, 1});
  const TaskPlan farActivity = runTaskPlan(data, seeded(), activity, 0.0, {&far, 1});
  EXPECT_DOUBLE_EQ(nearActivity.total, farActivity.total)
      << "a count of decays does not care where anyone is standing";

  const TaskPlan nearExposure = runTaskPlan(data, seeded(), exposure, 0.0, {&near, 1});
  const TaskPlan farExposure = runTaskPlan(data, seeded(), exposure, 0.0, {&far, 1});
  // Inverse square with attenuation off is exact, so the far leg is a quarter of the near one
  // to floating point rather than approximately.
  EXPECT_NEAR(farExposure.total, nearExposure.total / 4.0, nearExposure.total * 1.0e-12)
      << "inverse square, per leg";
}

TEST(TaskPlan, ABreakAdvancesTheClockAndAccruesNothing) {
  const NuclearData data = oneEmitter();
  ResponseSpec spec;
  spec.unit = Unit::Decays;

  const std::vector<PlanLeg> withBreak = {leg("a", 600.0, 1.0), breakFor("lunch", 3600.0),
                                          leg("b", 600.0, 1.0)};
  const TaskPlan plan = runTaskPlan(data, seeded(), spec, 0.0, withBreak);

  EXPECT_TRUE(plan.legs[1].isBreak);
  EXPECT_DOUBLE_EQ(plan.legs[1].accrued, 0.0);
  EXPECT_DOUBLE_EQ(plan.legs[1].meanRate, 0.0);
  EXPECT_DOUBLE_EQ(plan.elapsedSeconds, 4800.0);
  EXPECT_DOUBLE_EQ(plan.exposedSeconds, 1200.0) << "a break earns no dose and costs no time in";
  // And it really did advance the clock: the leg after it starts an hour later, so it accrues
  // less than the identical leg before it.
  EXPECT_LT(plan.legs[2].accrued, plan.legs[0].accrued);
  EXPECT_DOUBLE_EQ(plan.legs[2].startSeconds, 4200.0);
}

TEST(TaskPlan, OccupancyScalesALegAndTheTimeThatEarnsIt) {
  const NuclearData data = oneEmitter();
  ResponseSpec spec;
  spec.unit = Unit::Decays;

  const TaskPlan full =
      runTaskPlan(data, seeded(), spec, 0.0, std::vector<PlanLeg>{leg("all", 600.0, 1.0, 1.0)});
  const TaskPlan half =
      runTaskPlan(data, seeded(), spec, 0.0, std::vector<PlanLeg>{leg("half", 600.0, 1.0, 0.5)});

  EXPECT_NEAR(half.total, 0.5 * full.total, full.total * 1.0e-12);
  EXPECT_DOUBLE_EQ(half.exposedSeconds, 300.0);
  EXPECT_DOUBLE_EQ(half.elapsedSeconds, 600.0) << "the clock runs whether or not anyone is there";
  // The mean rate is per second ACTUALLY spent there, so it is unchanged by occupancy -- which
  // is what makes the column comparable between legs of different occupancy.
  EXPECT_NEAR(half.legs[0].meanRate, full.legs[0].meanRate, full.legs[0].meanRate * 1.0e-12);
}

TEST(TaskPlan, SaysWhereABudgetRunsOutRatherThanOnlyThatItDoes) {
  const NuclearData data = oneEmitter();
  ResponseSpec spec;
  spec.unit = Unit::Decays;

  const std::vector<PlanLeg> legs = {leg("a", 600.0, 1.0), leg("b", 600.0, 1.0),
                                     leg("c", 600.0, 1.0)};
  const TaskPlan whole = runTaskPlan(data, seeded(), spec, 0.0, legs);

  // A budget that the first leg clears and the second does not.
  const double budget = whole.legs[0].accrued + 0.5 * whole.legs[1].accrued;
  const TaskPlan plan = runTaskPlan(data, seeded(), spec, 0.0, legs, budget);

  EXPECT_TRUE(plan.budgetSpent);
  EXPECT_EQ(plan.spentInLeg, 1);
  EXPECT_GT(plan.spentAtSeconds, plan.legs[1].startSeconds);
  EXPECT_LT(plan.spentAtSeconds, plan.legs[1].endSeconds);
  // Over ten minutes this source barely decays, so the budget is spent almost exactly halfway
  // through the leg -- which is the check that the located instant is the right one and not
  // merely inside the bracket.
  EXPECT_NEAR(plan.spentAtSeconds, plan.legs[1].startSeconds + 300.0, 5.0);

  // A budget the plan fits leaves the location unset rather than pointing at the last leg.
  const TaskPlan fits = runTaskPlan(data, seeded(), spec, 0.0, legs, 10.0 * whole.total);
  EXPECT_FALSE(fits.budgetSpent);
  EXPECT_EQ(fits.spentInLeg, -1);
}

// A budget exhausted exactly at a leg boundary is spent AT the boundary, not inside anything.
//
// The natural way to reach this is to read a leg's accrual off one run and ask the next what
// happens on precisely that budget, which round-trips exactly through the JSON report. The leg
// that spent it then clears the budget rather than exceeding it, so the search moves on to the
// next leg with nothing left to spend -- and a root-find for a duration on a budget of zero is
// not a question about a duration. Answering it as the end of the leg that did the spending is
// both the truthful instant and the one that needs no search.
TEST(TaskPlan, LocatesABudgetExhaustedExactlyAtALegBoundary) {
  const NuclearData data = oneEmitter();
  ResponseSpec spec;
  spec.unit = Unit::Decays;

  const std::vector<PlanLeg> legs = {leg("a", 600.0, 1.0), leg("b", 600.0, 1.0)};
  const TaskPlan whole = runTaskPlan(data, seeded(), spec, 0.0, legs);

  const TaskPlan plan = runTaskPlan(data, seeded(), spec, 0.0, legs, whole.legs[0].accrued);
  EXPECT_TRUE(plan.budgetSpent);
  EXPECT_EQ(plan.spentInLeg, 0) << "the leg that spent it, not the one that would have";
  EXPECT_DOUBLE_EQ(plan.spentAtSeconds, plan.legs[0].endSeconds);
}

// And the instant is the end of the leg that spent it even when the clock runs on before
// anything else accrues. A break between the two is the case that separates "the end of the last
// leg that cost something" from "the start of the next leg that does": only the first is when
// the budget actually ran out, and they differ by the whole length of the break.
TEST(TaskPlan, DoesNotChargeABoundaryExhaustionToTheFarSideOfABreak) {
  const NuclearData data = oneEmitter();
  ResponseSpec spec;
  spec.unit = Unit::Decays;

  const std::vector<PlanLeg> legs = {leg("a", 600.0, 1.0), breakFor("lunch", 1800.0),
                                     leg("b", 600.0, 1.0)};
  const TaskPlan whole = runTaskPlan(data, seeded(), spec, 0.0, legs);

  const TaskPlan plan = runTaskPlan(data, seeded(), spec, 0.0, legs, whole.legs[0].accrued);
  EXPECT_TRUE(plan.budgetSpent);
  EXPECT_EQ(plan.spentInLeg, 0);
  EXPECT_DOUBLE_EQ(plan.spentAtSeconds, plan.legs[0].endSeconds);
  EXPECT_LT(plan.spentAtSeconds, plan.legs[2].startSeconds) << "the break is not on the meter";
}

TEST(TaskPlan, FractionsPartitionTheTotal) {
  const NuclearData data = oneEmitter();
  ResponseSpec spec;
  spec.metric = Metric::Exposure;
  spec.unit = Unit::Roentgen;

  const std::vector<PlanLeg> legs = {leg("in", 120.0, 5.0), leg("work", 1200.0, 0.5, 0.8),
                                     breakFor("wait", 600.0), leg("out", 120.0, 5.0)};
  const TaskPlan plan = runTaskPlan(data, seeded(), spec, 1.0e5, legs);

  double sum = 0.0;
  double fractions = 0.0;
  for (const LegResult& one : plan.legs) {
    sum += one.accrued;
    fractions += one.fraction;
  }
  EXPECT_NEAR(sum, plan.total, plan.total * 1.0e-12);
  EXPECT_NEAR(fractions, 1.0, 1.0e-12);
  EXPECT_NEAR(plan.legs.back().cumulative, plan.total, plan.total * 1.0e-12);
  EXPECT_NEAR(plan.legs.back().cumulativeFraction, 1.0, 1.0e-12);
}

TEST(TaskPlan, RefusesALegThatIsNotOne) {
  const NuclearData data = oneEmitter();
  ResponseSpec spec;
  spec.unit = Unit::Decays;
  const Inventory inv = seeded();

  EXPECT_THROW(runTaskPlan(data, inv, spec, 0.0, {}), InputError);
  EXPECT_THROW(runTaskPlan(data, inv, spec, -1.0, std::vector<PlanLeg>{leg("a", 60.0, 1.0)}),
               InputError);
  EXPECT_THROW(runTaskPlan(data, inv, spec, 0.0, std::vector<PlanLeg>{leg("a", 0.0, 1.0)}),
               InputError);
  EXPECT_THROW(runTaskPlan(data, inv, spec, 0.0, std::vector<PlanLeg>{leg("", 60.0, 1.0)}),
               InputError);
  EXPECT_THROW(runTaskPlan(data, inv, spec, 0.0, std::vector<PlanLeg>{leg("a", 60.0, 1.0, 1.5)}),
               InputError);
  EXPECT_THROW(runTaskPlan(data, inv, spec, 0.0, std::vector<PlanLeg>{leg("a", 60.0, 1.0, -0.1)}),
               InputError);

  // A leg accrues a total, so a rate unit is the wrong dimension for it -- the same refusal the
  // task curve and the stay time make, from the same place.
  ResponseSpec rate;
  rate.unit = Unit::Becquerel;
  EXPECT_THROW(runTaskPlan(data, inv, rate, 0.0, std::vector<PlanLeg>{leg("a", 60.0, 1.0)}),
               InputError);
}

}  // namespace
}  // namespace nusift
