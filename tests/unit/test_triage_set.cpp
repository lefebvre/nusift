// The smallest set that holds a floor everywhere at once.
//
// The properties that matter here are not numerical. A covering answer is only worth anything if
// it actually covers -- at every constrained point, not on average -- and if it is reproducible,
// since a monitoring list that changed between machines would be worse than none. Both are
// checked directly. What is NOT claimed anywhere is minimality: greedy is not exact, and the
// tests assert what greedy guarantees rather than what an exact solver would.
//
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <set>
#include <string>
#include <vector>

#include "nusift/core/error.hpp"
#include "nusift/engine/decay_engine.hpp"
#include "nusift/nucdata/nuclear_data.hpp"
#include "nusift/triage/triage_set.hpp"
#include "synthetic_chain.hpp"

namespace nusift {
namespace {

// Four nuclides whose activities cross over: a short-lived one that dominates early and a
// long-lived one that dominates late, so no single time's top-N is the answer for the grid.
NuclearData fourNuclides() {
  StoreArrays arrays;
  arrays.provenance.version = 1;
  // Four independent chains at different mass numbers, so nothing feeds anything and the
  // activities are exactly the exponentials the test reasons about.
  const double lambdas[] = {1.0e-2, 1.0e-3, 1.0e-4, 1.0e-6};
  for (int i = 0; i < 4; ++i) {
    arrays.nuclideKey.push_back(Zai{50, 100 + i, 0}.key());
    arrays.halfLife.push_back(synth::halfLifeFor(lambdas[i]));
  }
  arrays.modeOffset.assign(5, 0);
  return NuclearData::fromArrays(std::move(arrays));
}

Inventory allFour(double atoms) {
  Inventory inv;
  for (int i = 0; i < 4; ++i) {
    inv.add(Zai{50, 100 + i, 0}, atoms);
  }
  return inv;
}

ResponseTable activityOver(const NuclearData& data, const Inventory& inv,
                           const std::vector<double>& times) {
  return buildResponse(data, decay(data, inv, times), ResponseSpec{});
}

CoverageRequirement require(const ResponseTable& table, double fraction, const char* label) {
  CoverageRequirement requirement;
  requirement.table = &table;
  requirement.fraction = fraction;
  requirement.label = label;
  return requirement;
}

// The property the whole thing exists for, checked against the tables rather than against the
// search's own bookkeeping: at EVERY constrained point, the chosen set holds the floor.
void expectCoversEverywhere(const TriageSet& set,
                            const std::vector<CoverageRequirement>& requirements) {
  std::set<std::int64_t> keys;
  for (const SetMember& member : set.members) {
    keys.insert(member.id.key);
  }
  for (const CoverageRequirement& requirement : requirements) {
    const ResponseTable& table = *requirement.table;
    for (int k = 0; k < table.timeCount(); ++k) {
      const double total = table.totals[static_cast<std::size_t>(k)];
      if (!(total > 0.0)) {
        continue;
      }
      double held = 0.0;
      for (int c = 0; c < table.contributorCount(); ++c) {
        if (keys.count(table.contributors[static_cast<std::size_t>(c)].key) != 0) {
          held += table.valuesAt(k)[static_cast<std::size_t>(c)];
        }
      }
      EXPECT_GE(held / total, requirement.fraction - 1.0e-12)
          << requirement.label << " at t=" << table.times[static_cast<std::size_t>(k)];
    }
  }
}

TEST(TriageSet, HoldsTheFloorAtEveryTimeAndNotOnlyOnAverage) {
  const NuclearData data = fourNuclides();
  const std::vector<double> times = {1.0, 1.0e2, 1.0e4, 1.0e6, 1.0e8};
  const ResponseTable table = activityOver(data, allFour(1.0e20), times);

  const std::vector<CoverageRequirement> requirements = {require(table, 0.9, "activity")};
  const TriageSet set = robustTriageSet(requirements);

  EXPECT_FALSE(set.members.empty());
  EXPECT_TRUE(set.shortfalls.empty());
  EXPECT_EQ(set.candidateCount, 4);
  EXPECT_EQ(set.constraintCount, static_cast<int>(times.size()));
  expectCoversEverywhere(set, requirements);

  // The binding point is a real point of the grid, and the margin there is the smallest one.
  EXPECT_GE(set.binding.achieved, set.binding.required - 1.0e-12);
  EXPECT_DOUBLE_EQ(set.binding.required, 0.9);
}

// A single time needs whatever that time needs; the whole grid needs the union of what every
// time needs, which is more. Growing the grid can only grow the set, never shrink it -- the
// constraints are added, never relaxed.
TEST(TriageSet, AskingAboutMoreTimesNeverAsksForFewerNuclides) {
  const NuclearData data = fourNuclides();
  const Inventory inv = allFour(1.0e20);

  const ResponseTable early = activityOver(data, inv, {1.0});
  const ResponseTable whole = activityOver(data, inv, {1.0, 1.0e2, 1.0e4, 1.0e6, 1.0e8});

  const std::vector<CoverageRequirement> one = {require(early, 0.9, "activity")};
  const std::vector<CoverageRequirement> many = {require(whole, 0.9, "activity")};

  EXPECT_GE(robustTriageSet(many).members.size(), robustTriageSet(one).members.size());
}

// Two requirements at once is the point of the type. A set chosen for one metric guarantees
// nothing about the other, and the joint set has to satisfy both -- which is checked against
// both tables independently.
TEST(TriageSet, SatisfiesEveryRequirementJointly) {
  const NuclearData data = fourNuclides();
  const std::vector<double> times = {1.0, 1.0e4, 1.0e8};
  // Two inventories of the same nuclides in different proportions stand in for two metrics:
  // what matters to the search is that the two tables disagree about who dominates.
  Inventory heavyOnFast;
  heavyOnFast.add(Zai{50, 100, 0}, 1.0e22);
  heavyOnFast.add(Zai{50, 103, 0}, 1.0e18);
  Inventory heavyOnSlow;
  heavyOnSlow.add(Zai{50, 100, 0}, 1.0e18);
  heavyOnSlow.add(Zai{50, 103, 0}, 1.0e22);

  const ResponseTable a = activityOver(data, heavyOnFast, times);
  const ResponseTable b = activityOver(data, heavyOnSlow, times);
  const std::vector<CoverageRequirement> both = {require(a, 0.95, "sheet A"),
                                                 require(b, 0.95, "sheet B")};

  const TriageSet joint = robustTriageSet(both);
  expectCoversEverywhere(joint, both);
  EXPECT_EQ(joint.constraintCount, 2 * static_cast<int>(times.size()));

  // And the joint set is at least as large as either alone, since it satisfies strictly more.
  const std::vector<CoverageRequirement> justA = {require(a, 0.95, "sheet A")};
  EXPECT_GE(joint.members.size(), robustTriageSet(justA).members.size());
}

// Re-expressing one requirement in a different unit must not change the answer.
//
// Nothing about the covering question depends on whether an activity is stated in Bq or in Ci,
// or on whether it sits beside an exposure whose numbers are twenty orders of magnitude smaller.
// The requirement is a FRACTION of its own total, so the constraint is already scale-free; what
// was not, before the gains were normalised, was the greedy's choice of who to take next. Adding
// raw shortfalls across requirements let whichever one happened to carry the larger numbers
// decide the whole list, and the same problem stated in another unit came back a different list.
//
// Scaling one table's inventory is exactly the change a unit conversion makes -- every value and
// its total move by one common factor -- and the factor here is Bq per Ci, so this is the
// conversion rather than a stand-in for it.
//
// The two slowest nuclides and a grid inside their decade: what is being compared is two runs of
// the same search, so both have to be reading real numbers. Ten thousand e-foldings of the fast
// nuclide would leave a total that is pure solver round-off, whose SIGN is arbitrary, and a
// negative total is a constraint with nothing to cover -- the runs would then differ over which
// constraints exist at all, which is a fact about the grid and not about the units.
TEST(TriageSet, DoesNotDependOnTheUnitsARequirementIsStatedIn) {
  const NuclearData data = fourNuclides();
  const std::vector<double> times = {1.0, 1.0e4, 1.0e5};
  const Zai faster{50, 102, 0};  // lambda 1e-4
  const Zai slower{50, 103, 0};  // lambda 1e-6

  Inventory heavyOnFaster;
  heavyOnFaster.add(faster, 1.0e22);
  heavyOnFaster.add(slower, 1.0e18);

  // The same second requirement twice, differing only by the factor its numbers are carried in.
  const double becquerelPerCurie = 3.7e10;
  Inventory heavyOnSlower;
  heavyOnSlower.add(faster, 1.0e18);
  heavyOnSlower.add(slower, 1.0e22);
  Inventory heavyOnSlowerRescaled;
  heavyOnSlowerRescaled.add(faster, 1.0e18 * becquerelPerCurie);
  heavyOnSlowerRescaled.add(slower, 1.0e22 * becquerelPerCurie);

  const ResponseTable a = activityOver(data, heavyOnFaster, times);
  const ResponseTable b = activityOver(data, heavyOnSlower, times);
  const ResponseTable bRescaled = activityOver(data, heavyOnSlowerRescaled, times);

  const std::vector<CoverageRequirement> asStated = {require(a, 0.95, "sheet A"),
                                                     require(b, 0.95, "sheet B")};
  const std::vector<CoverageRequirement> rescaled = {require(a, 0.95, "sheet A"),
                                                     require(bRescaled, 0.95, "sheet B")};

  const TriageSet first = robustTriageSet(asStated);
  const TriageSet second = robustTriageSet(rescaled);

  ASSERT_FALSE(first.members.empty());
  ASSERT_EQ(first.members.size(), second.members.size());
  for (std::size_t i = 0; i < first.members.size(); ++i) {
    EXPECT_EQ(first.members[i].id.key, second.members[i].id.key);
    EXPECT_EQ(first.members[i].order, second.members[i].order);
    // And the marginal value that earned each its place is a fraction of a required coverage,
    // not an amount in anyone's unit, so it too survives the conversion.
    EXPECT_NEAR(first.members[i].closedShortfall, second.members[i].closedShortfall,
                1.0e-9 * std::abs(first.members[i].closedShortfall));
  }
  expectCoversEverywhere(second, rescaled);
}

// A monitoring list that changed between runs would be worse than no list. Ties are broken by
// contributor key, so the same tables give the same list in the same order every time.
TEST(TriageSet, IsReproducibleDownToTheOrder) {
  const NuclearData data = fourNuclides();
  const std::vector<double> times = {1.0, 1.0e4, 1.0e8};
  const ResponseTable table = activityOver(data, allFour(1.0e20), times);
  const std::vector<CoverageRequirement> requirements = {require(table, 0.9, "activity")};

  const TriageSet first = robustTriageSet(requirements);
  const TriageSet second = robustTriageSet(requirements);
  ASSERT_EQ(first.members.size(), second.members.size());
  for (std::size_t i = 0; i < first.members.size(); ++i) {
    EXPECT_EQ(first.members[i].label, second.members[i].label);
    EXPECT_EQ(first.members[i].order, second.members[i].order);
  }
}

// A floor of 1.0 needs every contributor that holds anything, which is the boundary case where
// the covering answer and the whole table coincide.
TEST(TriageSet, TakesEverythingWhenTheFloorIsEverything) {
  const NuclearData data = fourNuclides();
  const ResponseTable table = activityOver(data, allFour(1.0e20), {1.0});
  const std::vector<CoverageRequirement> requirements = {require(table, 1.0, "activity")};

  const TriageSet set = robustTriageSet(requirements);
  EXPECT_EQ(static_cast<int>(set.members.size()), set.candidateCount);
  EXPECT_TRUE(set.shortfalls.empty());
}

// A time at which a metric is identically zero constrains no set. Treating it as unmeetable
// would make every exposure answer report a shortfall it does not have, since a fresh inventory
// of pure beta emitters has no photon output at all.
TEST(TriageSet, IsNotConstrainedByATimeWithNothingToCover) {
  // A pure beta emitter whose daughter carries the only photon line. At t = 0 the daughter does
  // not exist yet, so the exposure total is exactly zero there and non-zero afterwards.
  StoreArrays arrays = synth::linearChain({1.0e-3, 5.0e-4});
  synth::addLines(arrays, 1, {661657.0}, {0.9});
  const NuclearData data = NuclearData::fromArrays(std::move(arrays));

  Inventory inv;
  inv.add(Zai{50, 100, 0}, 1.0e20);

  ResponseSpec spec;
  spec.metric = Metric::Exposure;
  spec.unit = Unit::RoentgenPerHour;
  const ResponseTable table =
      buildResponse(data, decay(data, inv, std::vector<double>{0.0, 1.0e3}), spec);
  ASSERT_EQ(table.timeCount(), 2);
  ASSERT_DOUBLE_EQ(table.totals[0], 0.0);
  ASSERT_GT(table.totals[1], 0.0);

  const std::vector<CoverageRequirement> requirements = {require(table, 0.95, "exposure")};
  const TriageSet set = robustTriageSet(requirements);
  // The t = 0 constraint is satisfied by the empty set and reports no shortfall: nothing to
  // cover is covered. Only the later time asks for anything.
  EXPECT_EQ(set.members.size(), 1u);
  EXPECT_TRUE(set.shortfalls.empty()) << "nothing to cover is covered";
  expectCoversEverywhere(set, requirements);
}

TEST(TriageSet, RefusesAQuestionItCannotAnswer) {
  const NuclearData data = fourNuclides();
  const ResponseTable byNuclide = activityOver(data, allFour(1.0e20), {1.0});

  EXPECT_THROW(robustTriageSet({}), InputError);

  CoverageRequirement noTable;
  noTable.label = "nothing";
  EXPECT_THROW(robustTriageSet(std::vector<CoverageRequirement>{noTable}), InputError);

  // A shortfall reported against requirement 2 of 4 names nothing a reader can act on.
  EXPECT_THROW(robustTriageSet(std::vector<CoverageRequirement>{require(byNuclide, 0.9, "")}),
               InputError);
  EXPECT_THROW(robustTriageSet(std::vector<CoverageRequirement>{require(byNuclide, 0.0, "a")}),
               InputError);
  EXPECT_THROW(robustTriageSet(std::vector<CoverageRequirement>{require(byNuclide, 1.5, "a")}),
               InputError);

  // A nuclide column and a mass-chain column are not the same kind of thing, and a set mixing
  // them would answer neither question.
  ResponseSpec chainSpec;
  chainSpec.aggregate = Aggregate::MassChain;
  const ResponseTable byChain =
      buildResponse(data, decay(data, allFour(1.0e20), std::vector<double>{1.0}), chainSpec);
  EXPECT_THROW(robustTriageSet(std::vector<CoverageRequirement>{require(byNuclide, 0.9, "a"),
                                                                require(byChain, 0.9, "b")}),
               InputError);
}

}  // namespace
}  // namespace nusift
