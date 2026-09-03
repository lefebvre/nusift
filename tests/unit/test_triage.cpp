#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <numbers>
#include <sstream>
#include <string>
#include <vector>

#include "nusift/core/error.hpp"
#include "nusift/engine/decay_engine.hpp"
#include "nusift/exposure/air_coefficients.hpp"
#include "nusift/exposure/point_source.hpp"
#include "nusift/nucdata/coefficient_pack.hpp"
#include "nusift/nucdata/nuclear_data.hpp"
#include "nusift/triage/ranking.hpp"
#include "nusift/triage/response.hpp"
#include "nusift/units.hpp"
#include "synthetic_chain.hpp"

namespace nusift {
namespace {

// A response table built directly, so ranking can be checked on values chosen for the edge
// case rather than on whatever a decay happens to produce.
ResponseTable tableOf(const std::vector<double>& values, const std::vector<std::int64_t>& keys) {
  ResponseTable table;
  table.times = {0.0};
  double total = 0.0;
  for (std::size_t i = 0; i < values.size(); ++i) {
    table.contributors.push_back(ContributorId{keys[i], 0});
    table.labels.push_back("c" + std::to_string(keys[i]));
    table.flags.push_back(kFlagNone);
    total += values[i];
  }
  table.values = values;
  table.totals = {total};
  return table;
}

TEST(Ranking, OrdersByValueDescending) {
  const ResponseTable table = tableOf({3.0, 10.0, 5.0}, {1, 2, 3});
  const Ranking ranking = rank(table, 0, RankRequest{});

  ASSERT_EQ(ranking.contributors.size(), 3u);
  EXPECT_DOUBLE_EQ(ranking.contributors[0].value, 10.0);
  EXPECT_DOUBLE_EQ(ranking.contributors[1].value, 5.0);
  EXPECT_DOUBLE_EQ(ranking.contributors[2].value, 3.0);
  EXPECT_EQ(ranking.contributors[0].rank, 1);
  EXPECT_EQ(ranking.contributors[2].rank, 3);
}

TEST(Ranking, FractionsSumToOneAndCumulativeIsMonotone) {
  const ResponseTable table = tableOf({3.0, 10.0, 5.0, 2.0}, {1, 2, 3, 4});
  const Ranking ranking = rank(table, 0, RankRequest{});

  double sum = 0.0;
  double previous = 0.0;
  for (const Contributor& c : ranking.contributors) {
    sum += c.fraction;
    EXPECT_GE(c.cumulativeFraction, previous);
    previous = c.cumulativeFraction;
  }
  EXPECT_NEAR(sum, 1.0, 1e-12);
  EXPECT_NEAR(ranking.coveredFraction, 1.0, 1e-12);
  EXPECT_EQ(ranking.omittedCount, 0);
}

// Exactly-tied values are common with symmetric inputs. Without a specified secondary key
// their order depends on the sort implementation, and any golden test over them becomes
// flaky across platforms.
TEST(Ranking, TiesBreakDeterministicallyByAscendingKey) {
  const ResponseTable table = tableOf({5.0, 5.0, 5.0}, {30, 10, 20});
  const Ranking ranking = rank(table, 0, RankRequest{});

  ASSERT_EQ(ranking.contributors.size(), 3u);
  EXPECT_EQ(ranking.contributors[0].id.key, 10);
  EXPECT_EQ(ranking.contributors[1].id.key, 20);
  EXPECT_EQ(ranking.contributors[2].id.key, 30);
}

// The literal question "which contributors are 95% of the total" wants the SMALLEST prefix
// reaching 95%, not the largest one below it.
TEST(Ranking, CoverageReturnsTheMinimalPrefixReachingTheTarget) {
  // 50, 30, 15, 5 -> cumulative 0.50, 0.80, 0.95, 1.00
  const ResponseTable table = tableOf({50.0, 30.0, 15.0, 5.0}, {1, 2, 3, 4});

  RankRequest request;
  request.topN = 0;
  request.coverage = 0.95;
  const Ranking ranking = rank(table, 0, request);

  ASSERT_EQ(ranking.contributors.size(), 3u);
  EXPECT_NEAR(ranking.coveredFraction, 0.95, 1e-12);
  EXPECT_EQ(ranking.omittedCount, 1);

  // A target between two steps still stops at the first one that reaches it.
  request.coverage = 0.6;
  const Ranking tighter = rank(table, 0, request);
  EXPECT_EQ(tighter.contributors.size(), 2u);
  EXPECT_NEAR(tighter.coveredFraction, 0.80, 1e-12);
}

// A truncated ranking must never present itself as complete. coveredFraction and
// omittedCount are what make the difference between a top-10 worth 40% and one worth 99%
// visible at all.
TEST(Ranking, TruncatedRankingReportsWhatItLeftOut) {
  const ResponseTable table = tableOf({50.0, 30.0, 15.0, 5.0}, {1, 2, 3, 4});
  RankRequest request;
  request.topN = 2;
  const Ranking ranking = rank(table, 0, request);

  EXPECT_EQ(ranking.contributors.size(), 2u);
  EXPECT_NEAR(ranking.coveredFraction, 0.80, 1e-12);
  EXPECT_EQ(ranking.omittedCount, 2);
  // The total is over everything, not over the returned prefix.
  EXPECT_DOUBLE_EQ(ranking.total, 100.0);
}

TEST(Ranking, MinFractionTrimsTheTail) {
  const ResponseTable table = tableOf({50.0, 30.0, 15.0, 5.0}, {1, 2, 3, 4});
  RankRequest request;
  request.topN = 0;
  request.minFraction = 0.10;
  const Ranking ranking = rank(table, 0, request);

  ASSERT_EQ(ranking.contributors.size(), 3u);
  EXPECT_DOUBLE_EQ(ranking.contributors.back().value, 15.0);
}

TEST(Ranking, TopNLargerThanTheTableClamps) {
  const ResponseTable table = tableOf({1.0, 2.0}, {1, 2});
  RankRequest request;
  request.topN = 100;
  const Ranking ranking = rank(table, 0, request);
  EXPECT_EQ(ranking.contributors.size(), 2u);
  EXPECT_EQ(ranking.omittedCount, 0);
}

// An inventory of nothing but stable nuclides has no activity. That is a legitimate answer,
// not an error, and it must not divide by zero on the way to reporting it.
TEST(Ranking, ZeroTotalYieldsAnEmptyRankingRatherThanNaN) {
  const ResponseTable table = tableOf({0.0, 0.0}, {1, 2});
  const Ranking ranking = rank(table, 0, RankRequest{});
  EXPECT_TRUE(ranking.contributors.empty());
  EXPECT_DOUBLE_EQ(ranking.total, 0.0);
  EXPECT_DOUBLE_EQ(ranking.coveredFraction, 0.0);
}

TEST(Ranking, RejectsOutOfRangeTimeIndex) {
  const ResponseTable table = tableOf({1.0}, {1});
  EXPECT_THROW(rank(table, 5, RankRequest{}), NusiftError);
  EXPECT_THROW(rank(table, -1, RankRequest{}), NusiftError);
}

// --- pinning ----------------------------------------------------------------

// The point of a pin: a contributor below every cut still appears, and it says where it really
// stands rather than borrowing a place in the prefix it was excluded from.
TEST(Ranking, PinnedContributorSurvivesTheCutCarryingItsTrueRank) {
  const ResponseTable table = tableOf({50.0, 30.0, 15.0, 5.0}, {1, 2, 3, 4});
  RankRequest request;
  request.topN = 2;
  request.pinned = {4};
  const Ranking ranking = rank(table, 0, request);

  ASSERT_EQ(ranking.contributors.size(), 3u);
  EXPECT_FALSE(ranking.contributors[0].pinned);
  EXPECT_FALSE(ranking.contributors[1].pinned);

  const Contributor& pinned = ranking.contributors[2];
  EXPECT_TRUE(pinned.pinned);
  EXPECT_EQ(pinned.id.key, 4);
  EXPECT_EQ(pinned.rank, 4) << "the rank it holds in the full ordering, not the row it landed on";
  EXPECT_NEAR(pinned.fraction, 0.05, 1e-12);
  // Everything down to it, which is the whole table -- so the number says the rows above it
  // account for all but its own share.
  EXPECT_NEAR(pinned.cumulativeFraction, 1.0, 1e-12);
}

// Coverage is the sum of what the reader can see, so a pinned row counts toward it and out of
// the omitted count. Reporting the prefix cumulative instead would understate what was shown.
TEST(Ranking, PinnedRowCountsTowardCoverageAndOutOfWhatWasOmitted) {
  const ResponseTable table = tableOf({50.0, 30.0, 15.0, 5.0}, {1, 2, 3, 4});
  RankRequest request;
  request.topN = 2;
  request.pinned = {4};
  const Ranking ranking = rank(table, 0, request);

  EXPECT_NEAR(ranking.coveredFraction, 0.85, 1e-12);
  EXPECT_EQ(ranking.omittedCount, 1) << "only the unpinned, unranked contributor is missing";
}

// A pin that ranked on its own is not a second row. Anything else would double-count it in the
// coverage, and print a contributor twice under two different headings.
TEST(Ranking, PinnedContributorAlreadyInTheTopIsNotRepeated) {
  const ResponseTable table = tableOf({50.0, 30.0, 15.0, 5.0}, {1, 2, 3, 4});
  RankRequest request;
  request.topN = 2;
  request.pinned = {1, 1, 2};  // including one named twice
  const Ranking ranking = rank(table, 0, request);

  ASSERT_EQ(ranking.contributors.size(), 2u);
  EXPECT_FALSE(ranking.contributors[0].pinned);
  EXPECT_FALSE(ranking.contributors[1].pinned);
  EXPECT_NEAR(ranking.coveredFraction, 0.80, 1e-12);
}

// Pinning is not filtering: it reaches past the cut without moving it, so the rows that would
// have been shown are all still shown, in the order they were.
TEST(Ranking, PinnedRowsFollowTheRankedOnesInTheOrderTheyStand) {
  const ResponseTable table = tableOf({50.0, 30.0, 15.0, 4.0, 1.0}, {1, 2, 3, 4, 5});
  RankRequest request;
  request.topN = 1;
  request.pinned = {5, 3};  // deliberately given worst-first
  const Ranking ranking = rank(table, 0, request);

  ASSERT_EQ(ranking.contributors.size(), 3u);
  EXPECT_EQ(ranking.contributors[0].id.key, 1);
  EXPECT_EQ(ranking.contributors[1].id.key, 3) << "3rd comes before 5th, whatever order they were "
                                                  "pinned in";
  EXPECT_EQ(ranking.contributors[2].id.key, 5);
  EXPECT_EQ(ranking.contributors[1].rank, 3);
  EXPECT_EQ(ranking.contributors[2].rank, 5);
}

// The cuts that drop a tail have to be reached past as well, or `--min-fraction` and
// `--coverage` would each quietly defeat a pin the way `--top` does not.
TEST(Ranking, PinReachesPastCoverageAndMinFractionToo) {
  const ResponseTable table = tableOf({50.0, 30.0, 15.0, 5.0}, {1, 2, 3, 4});

  RankRequest byCoverage;
  byCoverage.topN = 0;
  byCoverage.coverage = 0.5;
  byCoverage.pinned = {4};
  const Ranking covered = rank(table, 0, byCoverage);
  ASSERT_EQ(covered.contributors.size(), 2u);
  EXPECT_EQ(covered.contributors[1].id.key, 4);

  RankRequest byFraction;
  byFraction.topN = 0;
  byFraction.minFraction = 0.20;
  byFraction.pinned = {4};
  const Ranking trimmed = rank(table, 0, byFraction);
  ASSERT_EQ(trimmed.contributors.size(), 3u);
  EXPECT_EQ(trimmed.contributors[2].id.key, 4);
}

// A pinned contributor that contributes nothing is a real answer, not a missing row: a pure
// beta emitter pinned in an exposure ranking is exactly this case. It holds no place in the
// ordering, and rank 0 says so rather than implying it came last.
TEST(Ranking, PinnedContributorThatContributesNothingHasNoRank) {
  const ResponseTable table = tableOf({50.0, 30.0, 0.0}, {1, 2, 3});
  RankRequest request;
  request.pinned = {3};
  const Ranking ranking = rank(table, 0, request);

  ASSERT_EQ(ranking.contributors.size(), 3u);
  const Contributor& pinned = ranking.contributors[2];
  EXPECT_TRUE(pinned.pinned);
  EXPECT_EQ(pinned.rank, 0);
  EXPECT_DOUBLE_EQ(pinned.value, 0.0);
  EXPECT_DOUBLE_EQ(pinned.fraction, 0.0);
  // It was never in the count of contributors, so it cannot come out of it.
  EXPECT_EQ(ranking.omittedCount, 0);
  EXPECT_NEAR(ranking.coveredFraction, 1.0, 1e-12);
}

// One emitter is several columns in a gamma-line table, and a pin names the emitter. Matching
// on the key alone would be wrong in the other direction too -- it must catch every line, not
// just the first.
TEST(Ranking, PinningAnEmitterInALineTablePinsEveryLineItEmits) {
  const std::int64_t emitter = Zai{56, 137, 1}.key();
  const std::int64_t other = Zai{55, 134, 0}.key();

  ResponseTable table;
  table.aggregate = Aggregate::GammaLine;
  table.times = {0.0};
  table.contributors = {ContributorId{other, other, 604700.0},
                        ContributorId{emitter, emitter, 661700.0},
                        ContributorId{emitter, emitter, 31800.0}};
  table.labels = {"Cs-134 604.7 keV", "Ba-137m 661.7 keV", "Ba-137m 31.8 keV"};
  table.flags.assign(3, kFlagNone);
  table.values = {60.0, 30.0, 10.0};
  table.totals = {100.0};

  RankRequest request;
  request.topN = 1;
  request.pinned = {emitter};
  const Ranking ranking = rank(table, 0, request);

  ASSERT_EQ(ranking.contributors.size(), 3u);
  EXPECT_EQ(ranking.contributors[0].label, "Cs-134 604.7 keV");
  EXPECT_EQ(ranking.contributors[1].label, "Ba-137m 661.7 keV");
  EXPECT_EQ(ranking.contributors[2].label, "Ba-137m 31.8 keV");
  EXPECT_TRUE(ranking.contributors[1].pinned);
  EXPECT_TRUE(ranking.contributors[2].pinned);
}

// --- naming a contributor to pin --------------------------------------------

ResponseTable tableRankedBy(Aggregate aggregate, const std::vector<std::int64_t>& keys) {
  ResponseTable table = tableOf(std::vector<double>(keys.size(), 1.0), keys);
  table.aggregate = aggregate;
  return table;
}

TEST(Pin, ResolvesANuclideHoweverItIsSpelled) {
  const std::int64_t cs137 = Zai{55, 137, 0}.key();
  const ResponseTable table = tableRankedBy(Aggregate::Nuclide, {cs137});

  EXPECT_EQ(requirePin(table, "Cs-137"), cs137);
  EXPECT_EQ(requirePin(table, "cs137"), cs137);
  EXPECT_EQ(requirePin(table, "  551370 "), cs137);
}

// A user who knows a nuclide name should not have to work out which isobar it sits in. The
// bare and qualified forms of the number itself both work, so a pin can say what it means.
TEST(Pin, ResolvesAMassChainFromANumberOrANuclide) {
  const ResponseTable table = tableRankedBy(Aggregate::MassChain, {137});

  EXPECT_EQ(requirePin(table, "137"), 137);
  EXPECT_EQ(requirePin(table, "A=137"), 137);
  EXPECT_EQ(requirePin(table, "a=137"), 137);
  EXPECT_EQ(requirePin(table, "Cs-137"), 137) << "the chain the nuclide belongs to";
  EXPECT_EQ(requirePin(table, "551370"), 137) << "and the same through the raw-key spelling";
}

TEST(Pin, ResolvesAnElementFromASymbolNumberOrNuclide) {
  const ResponseTable table = tableRankedBy(Aggregate::Element, {55});

  EXPECT_EQ(requirePin(table, "Cs"), 55);
  EXPECT_EQ(requirePin(table, "cs"), 55);
  EXPECT_EQ(requirePin(table, "55"), 55);
  EXPECT_EQ(requirePin(table, "Z=55"), 55);
  EXPECT_EQ(requirePin(table, "Cs-137"), 55);
}

// A pin that resolves to nothing is refused rather than answered with zeros. "Cs-137 is not in
// this chain" and "Cs-137 contributes nothing here" are different statements, and a row of
// zeros makes the first look like the second.
TEST(Pin, RefusesAContributorTheTableDoesNotCarry) {
  const ResponseTable table = tableRankedBy(Aggregate::Nuclide, {Zai{55, 137, 0}.key()});
  EXPECT_THROW(requirePin(table, "Co-60"), InputError);
  EXPECT_THROW(requirePin(table, "not-a-nuclide"), InputError);
  EXPECT_THROW(requirePin(table, ""), InputError);

  const ResponseTable chains = tableRankedBy(Aggregate::MassChain, {137});
  EXPECT_THROW(requirePin(chains, "140"), InputError);
  EXPECT_THROW(requirePin(chains, "A=9999"), InputError);
}

// --- response construction -------------------------------------------------

// Activity is lambda*N per nuclide. With a single nuclide the whole table reduces to a value
// anyone can check by hand.
TEST(Response, ActivityIsDecayConstantTimesAtoms) {
  const double lambda = 1.0e-3;
  const NuclearData data = NuclearData::fromArrays(synth::linearChain({lambda}));
  Inventory inv;
  inv.add(Zai{50, 100, 0}, 1.0e20);

  const DecayResult result = decay(data, inv, std::vector<double>{0.0});
  const ResponseTable table = buildResponse(data, result, ResponseSpec{});

  EXPECT_EQ(table.contributorCount(), 2);  // parent + stable daughter
  EXPECT_NEAR(table.totals[0], lambda * 1.0e20, lambda * 1.0e20 * 1e-12);
}

TEST(Response, CurieScalesTheWholeTable) {
  const double lambda = 1.0e-3;
  const NuclearData data = NuclearData::fromArrays(synth::linearChain({lambda}));
  Inventory inv;
  inv.add(Zai{50, 100, 0}, 1.0e20);
  const DecayResult result = decay(data, inv, std::vector<double>{0.0});

  ResponseSpec spec;
  spec.unit = Unit::Curie;
  const ResponseTable table = buildResponse(data, result, spec);
  EXPECT_NEAR(table.totals[0], lambda * 1.0e20 / 3.7e10, lambda * 1.0e20 / 3.7e10 * 1e-12);
}

// Aggregating by mass chain sums nuclide columns; it is not a different calculation, so the
// total must be identical to the per-nuclide one.
TEST(Response, MassChainAggregationPreservesTheTotal) {
  const NuclearData data = NuclearData::fromArrays(synth::linearChain({1.0e-3, 5.0e-4}));
  Inventory inv;
  inv.add(Zai{50, 100, 0}, 1.0e20);
  const DecayResult result = decay(data, inv, std::vector<double>{2000.0});

  const ResponseTable byNuclide = buildResponse(data, result, ResponseSpec{});

  ResponseSpec chainSpec;
  chainSpec.aggregate = Aggregate::MassChain;
  const ResponseTable byChain = buildResponse(data, result, chainSpec);

  // A linear chain at fixed A collapses to exactly one isobar.
  EXPECT_EQ(byChain.contributorCount(), 1);
  EXPECT_GT(byNuclide.contributorCount(), 1);
  EXPECT_NEAR(byChain.totals[0], byNuclide.totals[0], byNuclide.totals[0] * 1e-12);
  EXPECT_EQ(byChain.labels[0].rfind("A=100", 0), 0u) << byChain.labels[0];
}

// A mass chain named only by its number is not actionable. The label carries the dominant
// member so a reader knows what to look at.
TEST(Response, MassChainLabelNamesItsDominantMember) {
  const NuclearData data = NuclearData::fromArrays(synth::linearChain({1.0e-3, 5.0e-4}));
  Inventory inv;
  inv.add(Zai{50, 100, 0}, 1.0e20);
  const DecayResult result = decay(data, inv, std::vector<double>{0.0});

  ResponseSpec spec;
  spec.aggregate = Aggregate::MassChain;
  const ResponseTable table = buildResponse(data, result, spec);
  ASSERT_EQ(table.contributorCount(), 1);
  // At t=0 only the seeded nuclide has any activity, so it must be the one named.
  EXPECT_NE(table.labels[0].find("Sn-100"), std::string::npos) << table.labels[0];
}

TEST(Response, ElementAggregationGroupsByAtomicNumber) {
  const NuclearData data = NuclearData::fromArrays(synth::linearChain({1.0e-3, 5.0e-4}));
  Inventory inv;
  inv.add(Zai{50, 100, 0}, 1.0e20);
  const DecayResult result = decay(data, inv, std::vector<double>{2000.0});

  ResponseSpec spec;
  spec.aggregate = Aggregate::Element;
  const ResponseTable table = buildResponse(data, result, spec);

  // Sn, Sb, Te -- three distinct elements in this chain.
  EXPECT_EQ(table.contributorCount(), 3);
  EXPECT_EQ(table.labels[0].rfind("Sn", 0), 0u) << table.labels[0];
}

// Bq is a rate and cannot express an interval total; decays cannot express an instant.
// Catching that at the boundary is what stops a number being reported in the wrong
// dimension, which is invisible in a table of bare figures.
TEST(Response, RejectsAUnitThatDoesNotSuitTheDomain) {
  const NuclearData data = NuclearData::fromArrays(synth::linearChain({1.0e-3}));
  Inventory inv;
  inv.add(Zai{50, 100, 0}, 1.0e20);
  const DecayResult result = decay(data, inv, std::vector<double>{100.0});

  ResponseSpec spec;
  spec.unit = Unit::Decays;  // interval-only
  EXPECT_THROW(buildResponse(data, result, spec), InputError);

  std::vector<std::int64_t> keys;
  const std::vector<double> integral = intervalIntegral(data, inv, 0.0, 100.0, &keys);
  ResponseSpec rate;
  rate.unit = Unit::Becquerel;  // instant-only
  EXPECT_THROW(buildIntervalResponse(data, keys, integral, 0.0, 100.0, rate), InputError);
}

// Integrated activity is a count of decays. For a single nuclide over [0, t] that count is
// N0 (1 - e^{-lambda t}) -- every atom that decayed -- which is a check anyone can do on
// paper and which pins that lambda is applied to atom-seconds and not to atoms.
TEST(Response, IntegratedActivityCountsTheDecaysThatOccurred) {
  const double lambda = 1.0e-3;
  const double n0 = 1.0e20;
  const NuclearData data = NuclearData::fromArrays(synth::linearChain({lambda}));
  Inventory inv;
  inv.add(Zai{50, 100, 0}, n0);

  const double t = 2000.0;
  std::vector<std::int64_t> keys;
  const std::vector<double> integral = intervalIntegral(data, inv, 0.0, t, &keys);

  ResponseSpec spec;
  spec.unit = Unit::Decays;
  const ResponseTable table = buildIntervalResponse(data, keys, integral, 0.0, t, spec);

  const double expected = n0 * (1.0 - std::exp(-lambda * t));
  EXPECT_NEAR(table.totals[0], expected, expected * 1e-9);
  EXPECT_EQ(table.domain, Domain::Interval);
}

// --- exposure metric ---------------------------------------------------------

// A chain where one nuclide emits photons and the other does not. Exposure and activity must
// therefore rank it differently, which is the entire reason both metrics exist.
NuclearData chainWithOnePhotonEmitter() {
  StoreArrays arrays = synth::linearChain({1.0e-3, 5.0e-4});
  // Lines on the middle member only. The seeded parent is a pure beta emitter.
  synth::addLines(arrays, 1, {661657.0}, {0.9});
  return NuclearData::fromArrays(std::move(arrays));
}

TEST(ResponseExposure, RanksDifferentlyFromActivity) {
  const NuclearData data = chainWithOnePhotonEmitter();
  Inventory inv;
  inv.add(Zai{50, 100, 0}, 1.0e20);
  // Early, while the seeded parent still holds most of the activity. Left long enough and the
  // longer-lived daughter out-accumulates it, at which point both rankings would agree on the
  // leader and the test would prove nothing.
  const DecayResult result = decay(data, inv, std::vector<double>{100.0});

  ResponseSpec activitySpec;
  const ResponseTable byActivity = buildResponse(data, result, activitySpec);

  ResponseSpec exposureSpec;
  exposureSpec.metric = Metric::Exposure;
  exposureSpec.unit = Unit::RoentgenPerHour;
  const ResponseTable byExposure = buildResponse(data, result, exposureSpec);

  const Ranking activityRank = rank(byActivity, 0, RankRequest{});
  const Ranking exposureRank = rank(byExposure, 0, RankRequest{});

  // The seeded parent dominates activity but contributes no exposure at all, so the two
  // rankings do not merely reorder -- they have different memberships.
  ASSERT_FALSE(activityRank.contributors.empty());
  ASSERT_FALSE(exposureRank.contributors.empty());
  EXPECT_EQ(activityRank.contributors[0].label, "Sn-100");
  EXPECT_EQ(exposureRank.contributors[0].label, "Sb-100");
  EXPECT_EQ(exposureRank.contributors.size(), 1u)
      << "only the photon emitter should contribute exposure";
}

// Exposure is lambda * N * (exposure per becquerel), so it must equal the activity times the
// per-becquerel factor the physics layer computes independently.
TEST(ResponseExposure, EqualsActivityTimesTheperBecquerelFactor) {
  const NuclearData data = chainWithOnePhotonEmitter();
  Inventory inv;
  inv.add(Zai{50, 100, 0}, 1.0e20);
  const DecayResult result = decay(data, inv, std::vector<double>{2000.0});

  ResponseSpec exposureSpec;
  exposureSpec.metric = Metric::Exposure;
  exposureSpec.unit = Unit::RoentgenPerHour;
  const ResponseTable table = buildResponse(data, result, exposureSpec);

  const int emitter = data.indexOf(Zai{51, 100, 0});
  ASSERT_GE(emitter, 0);
  const int resultIdx = [&] {
    for (int i = 0; i < result.nuclideCount(); ++i) {
      if (result.nuclideKeys[static_cast<std::size_t>(i)] == Zai{51, 100, 0}.key()) {
        return i;
      }
    }
    return -1;
  }();
  ASSERT_GE(resultIdx, 0);

  const double activity = data.decayConstant(emitter) * result.atomsAt(0)[resultIdx];
  const double perBq =
      exposure::exposureRatePerBecquerel(data.lines(emitter), exposureSpec.geometry);
  EXPECT_NEAR(table.totals[0], activity * perBq, activity * perBq * 1e-10);
}

// Halving nothing but the distance must quarter every exposure value, since the geometry
// enters as 1/d^2 and the response layer applies the same coefficient the physics layer does.
TEST(ResponseExposure, ScalesWithGeometry) {
  const NuclearData data = chainWithOnePhotonEmitter();
  Inventory inv;
  inv.add(Zai{50, 100, 0}, 1.0e20);
  const DecayResult result = decay(data, inv, std::vector<double>{2000.0});

  ResponseSpec near;
  near.metric = Metric::Exposure;
  near.unit = Unit::RoentgenPerHour;
  near.geometry.airAttenuation = false;
  near.geometry.distanceM = 1.0;

  ResponseSpec far = near;
  far.geometry.distanceM = 2.0;

  const double atOne = buildResponse(data, result, near).totals[0];
  const double atTwo = buildResponse(data, result, far).totals[0];
  EXPECT_NEAR(atTwo, atOne / 4.0, atOne * 1e-12);
}

// The unit selects the QUANTITY on this metric, and the two quantities are not proportional.
// A sievert is ICRP 116 effective dose; a roentgen and a gray are air kerma. Whether they
// happen to agree is a fact about the photon's energy, not about the units.
//
// At 662 keV they nearly do -- ICRP's own Table A.2 puts effective dose at 1.02 times air
// kerma there, which is the coincidence that let the old air-kerma-as-sievert column look
// right for cesium and cobalt.
TEST(ResponseExposure, SievertIsEffectiveDoseAndNotAConvertedRoentgen) {
  const NuclearData data = chainWithOnePhotonEmitter();  // one 661.657 keV line
  Inventory inv;
  inv.add(Zai{50, 100, 0}, 1.0e20);
  const DecayResult result = decay(data, inv, std::vector<double>{2000.0});

  ResponseSpec roentgen;
  roentgen.metric = Metric::Exposure;
  roentgen.unit = Unit::RoentgenPerHour;
  roentgen.geometry.airAttenuation = false;
  ResponseSpec sievert = roentgen;
  sievert.unit = Unit::SievertPerHour;

  const double inR = buildResponse(data, result, roentgen).totals[0];
  const double inSv = buildResponse(data, result, sievert).totals[0];
  const double asIfConverted = inR * units::kGyPerR;
  EXPECT_NEAR(inSv / asIfConverted, 1.02, 0.01);
}

// And where the old column was wrong, it was wrong by a factor rather than by a rounding. At
// 20 keV air keeps absorbing strongly while a body's organs are shielded by everything in
// front of them, so effective dose is a seventh of the air kerma -- ICRP Table A.2 gives 0.130
// Sv/Gy in AP. The nuclide this describes is Am-241, whose spectrum is mostly this soft.
TEST(ResponseExposure, ASoftEmitterNoLongerReportsAirKermaAsASievert) {
  StoreArrays arrays = synth::linearChain({1.0e-3, 5.0e-4});
  synth::addLines(arrays, 1, {2.0e4}, {0.9});
  const NuclearData data = NuclearData::fromArrays(std::move(arrays));

  Inventory inv;
  inv.add(Zai{50, 100, 0}, 1.0e20);
  const DecayResult result = decay(data, inv, std::vector<double>{2000.0});

  ResponseSpec roentgen;
  roentgen.metric = Metric::Exposure;
  roentgen.unit = Unit::RoentgenPerHour;
  roentgen.geometry.airAttenuation = false;
  ResponseSpec sievert = roentgen;
  sievert.unit = Unit::SievertPerHour;

  const double inR = buildResponse(data, result, roentgen).totals[0];
  const double inSv = buildResponse(data, result, sievert).totals[0];
  EXPECT_NEAR(inSv / (inR * units::kGyPerR), 0.130, 0.130 * 0.02);
}

// The irradiation geometry is part of what a sievert means, so it has to reach the weight.
// Facing the source and facing away from it are different answers for the same field, and a
// table that ignored the field would report one of them under both names.
TEST(ResponseExposure, TheIrradiationGeometryReachesTheWeight) {
  const NuclearData data = chainWithOnePhotonEmitter();
  Inventory inv;
  inv.add(Zai{50, 100, 0}, 1.0e20);
  const DecayResult result = decay(data, inv, std::vector<double>{2000.0});

  ResponseSpec ap;
  ap.metric = Metric::Exposure;
  ap.unit = Unit::SievertPerHour;
  ResponseSpec pa = ap;
  pa.geometry.irradiation = exposure::Irradiation::PA;

  const double facing = buildResponse(data, result, ap).totals[0];
  const double away = buildResponse(data, result, pa).totals[0];
  EXPECT_GT(facing, away);
  // ICRP 116 Table A.1 at 0.662 MeV, AP over PA. Not to the last digit: the cesium line is
  // 661.657 keV, a few hundred eV below ICRP's grid point, so both columns are interpolated
  // and the ratio moves in the fifth figure.
  EXPECT_NEAR(facing / away, 3.17 / 2.62, 1.0e-3);
}

// --- a kernel pack as a metric ------------------------------------------------

CoefficientPack packFrom(const std::string& text) {
  std::istringstream in(text);
  return CoefficientPack::read(in, "test.csv");
}

// A kernel is evaluated against the nuclide's own lines and the spec's geometry, which is the
// same construction the built-in photon metrics use. With one line and a flat curve the whole
// weight is checkable by hand: lambda * intensity * kappa / (4 pi d^2).
TEST(ResponsePack, AKernelWeighsTheLinesWithTheGeometryInside) {
  StoreArrays arrays = synth::linearChain({1.0e-3, 5.0e-4});
  synth::addLines(arrays, 1, {6.0e5}, {0.5});
  const NuclearData data = NuclearData::fromArrays(std::move(arrays));
  Inventory inv;
  inv.add(Zai{50, 100, 0}, 1.0e20);

  const CoefficientPack pack = packFrom(
      "# pack: flat\n# version: 1\n# quantity: a dose quantity\n# unit: Sv/s\n"
      "# shape: kernel\n# applies: fluence\n# domain: instant\n# source: invented\n"
      "energy_ev,coefficient\n1.0e5,2.0e-16\n1.0e7,2.0e-16\n");
  const ResolvedPack resolved = resolvePack(pack, data, inv);

  ResponseSpec spec;
  spec.metric = Metric::Pack;
  spec.unit = Unit::PackDefined;
  spec.pack = &resolved;
  spec.geometry.distanceM = 1.0;
  spec.geometry.airAttenuation = false;

  const ResponseTable table =
      buildResponse(data, decay(data, inv, std::vector<double>{2000.0}), spec);
  const int emitter = data.indexOfKey(Zai{51, 100, 0}.key());
  ASSERT_GE(emitter, 0);

  const double atoms =
      decay(data, inv, std::vector<double>{2000.0}).atomsAt(0)[static_cast<std::size_t>(emitter)];
  const double expected =
      data.decayConstant(emitter) * atoms * 0.5 * 2.0e-16 / (4.0 * std::numbers::pi);
  EXPECT_NEAR(table.totals[0], expected, expected * 1.0e-12);
  EXPECT_EQ(table.unitLabel, "Sv/s");
  // Every line sits inside the table, so the curve spoke for all of the answer.
  EXPECT_DOUBLE_EQ(table.packCoverage[0], 1.0);
}

// What a kernel can fail to speak for is an ENERGY, not a nuclide, so coverage measures the
// share of the answer that came from lines inside the tabulated range. A spectrum sitting off
// the end of a published curve is not an answer that curve can give.
TEST(ResponsePack, AKernelReportsTheShareOfTheAnswerItHadToClamp) {
  StoreArrays arrays = synth::linearChain({1.0e-3, 5.0e-4});
  // One line inside the curve's range and one far below it, of equal intensity.
  synth::addLines(arrays, 1, {5.0e5, 1.0e3}, {0.5, 0.5});
  const NuclearData data = NuclearData::fromArrays(std::move(arrays));
  Inventory inv;
  inv.add(Zai{50, 100, 0}, 1.0e20);

  const CoefficientPack pack = packFrom(
      "# pack: flat\n# version: 1\n# quantity: a dose quantity\n# unit: Sv/s\n"
      "# shape: kernel\n# applies: fluence\n# domain: instant\n# source: invented\n"
      "energy_ev,coefficient\n1.0e5,1.0e-16\n1.0e7,1.0e-16\n");
  const ResolvedPack resolved = resolvePack(pack, data, inv);

  ResponseSpec spec;
  spec.metric = Metric::Pack;
  spec.unit = Unit::PackDefined;
  spec.pack = &resolved;
  spec.geometry.airAttenuation = false;

  const ResponseTable table =
      buildResponse(data, decay(data, inv, std::vector<double>{2000.0}), spec);
  // The curve is flat and the two lines have equal intensity, so the clamped one carries half.
  EXPECT_NEAR(table.packCoverage[0], 0.5, 1.0e-9);
}

// --- coefficient packs as a metric ------------------------------------------

// The columns of a table are the chain's, not the seed's, so a nuclide's index is not the order
// it was written in.
int columnOfLabel(const ResponseTable& table, const std::string& label) {
  for (int c = 0; c < table.contributorCount(); ++c) {
    if (table.labels[static_cast<std::size_t>(c)] == label) {
      return c;
    }
  }
  throw NusiftError("test: no column labeled " + label);
}

std::string packText(const char* progeny, const char* rows, const char* domain = "instant") {
  return std::string("# pack: test\n# version: 1\n# quantity: an index\n# unit: 1\n") +
         "# basis: activity\n# domain: " + domain + "\n# progeny: " + progeny +
         "\n# source: invented\n" + rows;
}

// A pack is a metric like any other once resolved: a fixed weight per nuclide, post-multiplied
// on the same solve, with the total the weighted sum and nothing more.
TEST(ResponsePack, IsAWeightedSumLikeEveryOtherMetric) {
  StoreArrays arrays = synth::linearChain({1.0e-3, 5.0e-4});
  const NuclearData data = NuclearData::fromArrays(std::move(arrays));
  Inventory inv;
  inv.add(Zai{50, 100, 0}, 1.0e20);

  const CoefficientPack pack =
      packFrom(packText("excluded", "nuclide,coefficient\nSn-100,2\nSb-100,4\n"));
  const ResolvedPack resolved = resolvePack(pack, data, inv);

  ResponseSpec spec;
  spec.metric = Metric::Pack;
  spec.unit = Unit::PackDefined;
  spec.pack = &resolved;

  const ResponseTable table = buildResponse(data, decay(data, inv, std::vector<double>{0.0}), spec);
  const int parent = data.indexOfKey(Zai{50, 100, 0}.key());
  EXPECT_NEAR(table.totals[0], 2.0 * data.decayConstant(parent) * 1.0e20, table.totals[0] * 1e-12);
  EXPECT_EQ(table.unitLabel, "1") << "the pack's own spelling, not the enum's placeholder";
  ASSERT_EQ(table.packCoverage.size(), 1u);
  EXPECT_DOUBLE_EQ(table.packCoverage[0], 1.0);
}

// The fold, end to end: a daughter whose parent is seeded contributes nothing of its own and
// says why, while the same daughter seeded alone uses its own row.
TEST(ResponsePack, AFoldedDaughterIsNotWeightedTwice) {
  StoreArrays arrays = synth::linearChain({1.0e-3, 5.0e-4});
  const NuclearData data = NuclearData::fromArrays(std::move(arrays));
  const CoefficientPack pack =
      packFrom(packText("folded", "nuclide,coefficient,folded\nSn-100,2,Sb-100\nSb-100,4,\n"));

  Inventory withParent;
  withParent.add(Zai{50, 100, 0}, 1.0e20);
  Inventory daughterOnly;
  daughterOnly.add(Zai{51, 100, 0}, 1.0e20);

  const ResolvedPack under = resolvePack(pack, data, withParent);
  const ResolvedPack alone = resolvePack(pack, data, daughterOnly);

  ResponseSpec spec;
  spec.metric = Metric::Pack;
  spec.unit = Unit::PackDefined;

  spec.pack = &under;
  const ResponseTable folded =
      buildResponse(data, decay(data, withParent, std::vector<double>{5000.0}), spec);
  spec.pack = &alone;
  const ResponseTable ownRow =
      buildResponse(data, decay(data, daughterOnly, std::vector<double>{5000.0}), spec);

  const int daughterColumn = columnOfLabel(folded, "Sb-100");
  EXPECT_DOUBLE_EQ(folded.valuesAt(0)[static_cast<std::size_t>(daughterColumn)], 0.0)
      << "its contribution is inside its parent's coefficient";
  EXPECT_NE(folded.flags[static_cast<std::size_t>(daughterColumn)] & kFlagFoldedInPack, 0);
  // Covered even though it carries no weight: it IS accounted for, by the row above it.
  EXPECT_DOUBLE_EQ(folded.packCoverage[0], 1.0);

  EXPECT_GT(ownRow.valuesAt(0)[static_cast<std::size_t>(columnOfLabel(ownRow, "Sb-100"))], 0.0)
      << "seeded alone, its own row applies";
}

// A total over a fifth of an inventory looks exactly like one over all of it. The coverage
// figure is the only thing that tells them apart.
TEST(ResponsePack, CoverageReportsWhatThePackCouldNotSpeakFor) {
  StoreArrays arrays = synth::linearChain({1.0e-3, 5.0e-4});
  const NuclearData data = NuclearData::fromArrays(std::move(arrays));
  Inventory inv;
  inv.add(Zai{50, 100, 0}, 1.0e20);

  // A pack that knows the daughter and not the parent, asked at t=0 where the parent holds
  // every becquerel.
  const CoefficientPack pack = packFrom(packText("excluded", "nuclide,coefficient\nSb-100,4\n"));
  const ResolvedPack resolved = resolvePack(pack, data, inv);
  ResponseSpec spec;
  spec.metric = Metric::Pack;
  spec.unit = Unit::PackDefined;
  spec.pack = &resolved;

  const ResponseTable table = buildResponse(data, decay(data, inv, std::vector<double>{0.0}), spec);
  EXPECT_DOUBLE_EQ(table.totals[0], 0.0);
  EXPECT_DOUBLE_EQ(table.packCoverage[0], 0.0);
  EXPECT_NE(table.flags[static_cast<std::size_t>(columnOfLabel(table, "Sn-100"))] & kFlagNotInPack,
            0);
}

TEST(ResponsePack, RefusesASpecItCannotAnswer) {
  StoreArrays arrays = synth::linearChain({1.0e-3, 5.0e-4});
  const NuclearData data = NuclearData::fromArrays(std::move(arrays));
  Inventory inv;
  inv.add(Zai{50, 100, 0}, 1.0e20);
  const DecayResult result = decay(data, inv, std::vector<double>{0.0});

  ResponseSpec missing;
  missing.metric = Metric::Pack;
  missing.unit = Unit::PackDefined;
  EXPECT_THROW(buildResponse(data, result, missing), InputError) << "a pack metric with no pack";

  const CoefficientPack instantOnly =
      packFrom(packText("excluded", "nuclide,coefficient\nSn-100,2\n", "instant"));
  const ResolvedPack resolved = resolvePack(instantOnly, data, inv);
  ResponseSpec spec;
  spec.metric = Metric::Pack;
  spec.unit = Unit::PackDefined;
  spec.pack = &resolved;

  // An instantaneous quantity against an integral is a different quantity, not the same one
  // summed, and the pack said which it is.
  std::vector<std::int64_t> keys;
  const std::vector<double> integral = intervalIntegral(data, inv, 0.0, 100.0, &keys);
  EXPECT_THROW(buildIntervalResponse(data, keys, integral, 0.0, 100.0, spec), InputError);

  ResponseSpec wrongUnit = spec;
  wrongUnit.unit = Unit::Becquerel;
  EXPECT_THROW(buildResponse(data, result, wrongUnit), InputError);
}

// The table carries the optical depth of the air path -- mu(E) rho d for a single line -- and
// the geometry it was built in, because both are part of what an exposure number means: the
// caveat about scattered photons is judged on the first, and a report states the second.
TEST(ResponseExposure, CarriesTheAirPathAndTheGeometryItWasBuiltIn) {
  const NuclearData data = chainWithOnePhotonEmitter();
  Inventory inv;
  inv.add(Zai{50, 100, 0}, 1.0e20);
  const DecayResult result = decay(data, inv, std::vector<double>{2000.0});

  ResponseSpec spec;
  spec.metric = Metric::Exposure;
  spec.unit = Unit::RoentgenPerHour;
  spec.geometry.distanceM = 100.0;
  const ResponseTable table = buildResponse(data, result, spec);

  ASSERT_EQ(table.meanOpticalDepth.size(), 1u);
  const double expected =
      exposure::airMassAttenuation(661657.0) * spec.geometry.airDensityKgM3 * 100.0;
  EXPECT_NEAR(table.meanOpticalDepth[0], expected, expected * 1e-12);
  EXPECT_EQ(table.geometry.distanceM, 100.0);
  EXPECT_EQ(rank(table, 0, RankRequest{}).meanOpticalDepth, table.meanOpticalDepth[0]);

  spec.geometry.airAttenuation = false;
  EXPECT_EQ(buildResponse(data, result, spec).meanOpticalDepth[0], 0.0);

  // Activity has no air path, and the table must not carry one that reads as a measurement.
  ResponseSpec activity;
  activity.metric = Metric::Activity;
  activity.unit = Unit::Becquerel;
  EXPECT_TRUE(buildResponse(data, result, activity).meanOpticalDepth.empty());
}

// Exposure is computed in R/h, and an interval weights atom-SECONDS, so the accrued total is
// the rate integrated over the window in HOURS. Leaving the window in seconds inflates every
// integrated exposure by 3600 -- invisible in a table of bare figures, and wrong by more than
// three orders of magnitude. The check is one anyone can do on paper: an hour at 1 R/h is 1 R.
TEST(ResponseExposure, IntegratedExposureAccruesPerHourNotPerSecond) {
  const double lambda = 1.0e-9;  // ~22 y, so almost nothing decays over the hour
  const double n0 = 1.0e20;
  StoreArrays arrays = synth::linearChain({lambda});
  synth::addLines(arrays, 0, {661657.0}, {0.9});  // the seeded nuclide is the emitter
  const NuclearData data = NuclearData::fromArrays(std::move(arrays));

  Inventory inv;
  inv.add(Zai{50, 100, 0}, n0);
  const double window = units::kSecondsPerHour;

  ResponseSpec rateSpec;
  rateSpec.metric = Metric::Exposure;
  rateSpec.unit = Unit::RoentgenPerHour;
  const double rate =
      buildResponse(data, decay(data, inv, std::vector<double>{0.0}), rateSpec).totals[0];

  std::vector<std::int64_t> keys;
  const std::vector<double> integral = intervalIntegral(data, inv, 0.0, window, &keys);
  ResponseSpec accruedSpec = rateSpec;
  accruedSpec.unit = Unit::Roentgen;
  const double accrued =
      buildIntervalResponse(data, keys, integral, 0.0, window, accruedSpec).totals[0];

  // The window in hours, per seeded atom: the exact integral rather than a flat hour, so the
  // comparison stays a strict equality instead of an approximation.
  const double hours = synth::batemanIntegralN0(1.0, lambda, window) / units::kSecondsPerHour;
  EXPECT_NEAR(accrued, rate * hours, rate * hours * 1e-9);
  EXPECT_LT(accrued, rate) << "an hour of decay accrues slightly less than the initial rate";
}

// A store with no photon lines cannot answer an exposure question. Returning zeros would be
// indistinguishable from "nothing here emits photons", which is a different claim entirely.
TEST(ResponseExposure, RefusesAStoreWithNoPhotonLines) {
  const NuclearData data = NuclearData::fromArrays(synth::linearChain({1.0e-3}));
  ASSERT_FALSE(data.hasPhotonLines());
  Inventory inv;
  inv.add(Zai{50, 100, 0}, 1.0e20);
  const DecayResult result = decay(data, inv, std::vector<double>{100.0});

  ResponseSpec spec;
  spec.metric = Metric::Exposure;
  spec.unit = Unit::RoentgenPerHour;
  try {
    buildResponse(data, result, spec);
    FAIL() << "expected InputError";
  } catch (const InputError& e) {
    const std::string what = e.what();
    EXPECT_NE(what.find("no photon lines"), std::string::npos) << what;
  }
}

// Becquerel does not measure exposure and roentgen does not measure activity. That is a
// category error rather than a rounding one, so it is refused rather than converted.
TEST(ResponseExposure, RefusesAUnitThatDoesNotMeasureTheMetric) {
  const NuclearData data = chainWithOnePhotonEmitter();
  Inventory inv;
  inv.add(Zai{50, 100, 0}, 1.0e20);
  const DecayResult result = decay(data, inv, std::vector<double>{100.0});

  ResponseSpec exposureInBq;
  exposureInBq.metric = Metric::Exposure;
  exposureInBq.unit = Unit::Becquerel;
  EXPECT_THROW(buildResponse(data, result, exposureInBq), InputError);

  ResponseSpec activityInRoentgen;
  activityInRoentgen.metric = Metric::Activity;
  activityInRoentgen.unit = Unit::RoentgenPerHour;
  EXPECT_THROW(buildResponse(data, result, activityInRoentgen), InputError);
}

// --- per-line aggregation ----------------------------------------------------

// An emitter with two strong lines and two below the column floor, so thresholding and
// ordering are both observable. The 500 keV line at 2e-7 is about 1e-7 of the emitter's
// exposure: far enough under the 1e-6 floor to be dropped as a column, and far enough above
// round-off that a total which forgot it would miss the tolerance below by four decades.
NuclearData chainWithStrongAndTraceLines() {
  StoreArrays arrays = synth::linearChain({1.0e-3, 5.0e-4});
  synth::addLines(arrays, 1, {1332492.0, 661657.0, 5.0e5, 100.0}, {0.5, 0.9, 2.0e-7, 1.0e-9});
  return NuclearData::fromArrays(std::move(arrays));
}

// Splitting a nuclide's exposure across its lines must not change the total. This is the
// invariant that says per-line aggregation is a finer view of the same quantity rather than a
// second, separately-derived calculation -- and it holds exactly, dropped columns included,
// because a line the table declines to carry still contributes to what it reports as the whole.
TEST(ResponseLines, LineTotalMatchesTheNuclideTotal) {
  const NuclearData data = chainWithStrongAndTraceLines();
  Inventory inv;
  inv.add(Zai{50, 100, 0}, 1.0e20);
  const DecayResult result = decay(data, inv, std::vector<double>{2000.0});

  ResponseSpec byNuclide;
  byNuclide.metric = Metric::Exposure;
  byNuclide.unit = Unit::RoentgenPerHour;

  ResponseSpec byLine = byNuclide;
  byLine.aggregate = Aggregate::GammaLine;

  const double nuclideTotal = buildResponse(data, result, byNuclide).totals[0];
  const double lineTotal = buildResponse(data, result, byLine).totals[0];
  EXPECT_NEAR(lineTotal, nuclideTotal, nuclideTotal * 1e-12);
}

// The strongest line by emitted exposure leads, and every column names its emitter and energy.
TEST(ResponseLines, RanksIndividualLinesAndLabelsThem) {
  const NuclearData data = chainWithStrongAndTraceLines();
  Inventory inv;
  inv.add(Zai{50, 100, 0}, 1.0e20);
  const DecayResult result = decay(data, inv, std::vector<double>{2000.0});

  ResponseSpec spec;
  spec.metric = Metric::Exposure;
  spec.unit = Unit::RoentgenPerHour;
  spec.aggregate = Aggregate::GammaLine;
  const Ranking ranking = rank(buildResponse(data, result, spec), 0, RankRequest{});

  ASSERT_GE(ranking.contributors.size(), 2u);
  // 1332 keV at 0.5 carries more energy than 662 keV at 0.9, so it leads.
  EXPECT_NE(ranking.contributors[0].label.find("1332.5 keV"), std::string::npos)
      << ranking.contributors[0].label;
  EXPECT_NE(ranking.contributors[0].label.find("Sb-100"), std::string::npos)
      << "a line column must name its emitter: " << ranking.contributors[0].label;
  EXPECT_NE(ranking.contributors[1].label.find("661.7 keV"), std::string::npos)
      << ranking.contributors[1].label;

  // The energy is carried numerically as well as in the label, for a consumer binning a
  // spectrum rather than reading one.
  EXPECT_NEAR(ranking.contributors[0].id.lineEnergyEv, 1332492.0, 1.0);
  EXPECT_EQ(ranking.contributors[0].id.dominantMemberKey, (Zai{51, 100, 0}.key()));
}

// A full evaluation carries 86000 lines and most contribute nothing. The threshold is relative
// to the emitter, so a minor nuclide keeps its own spectrum instead of vanishing wholesale.
TEST(ResponseLines, DropsLinesFarBelowTheirEmittersTotal) {
  const NuclearData data = chainWithStrongAndTraceLines();
  Inventory inv;
  inv.add(Zai{50, 100, 0}, 1.0e20);
  const DecayResult result = decay(data, inv, std::vector<double>{2000.0});

  ResponseSpec spec;
  spec.metric = Metric::Exposure;
  spec.unit = Unit::RoentgenPerHour;
  spec.aggregate = Aggregate::GammaLine;
  const ResponseTable table = buildResponse(data, result, spec);

  // Four lines were staged; the 500 keV line at 2e-7 and the 100 eV line at 1e-9 are both far
  // below the floor, and neither gets a column.
  EXPECT_EQ(table.contributorCount(), 2);
}

// A photon line has no activity of its own -- it is a way its emitter's decays get out. Asking
// for activity by line is a category error, not a rounding one.
TEST(ResponseLines, RefusesActivityRankedByLine) {
  const NuclearData data = chainWithStrongAndTraceLines();
  Inventory inv;
  inv.add(Zai{50, 100, 0}, 1.0e20);
  const DecayResult result = decay(data, inv, std::vector<double>{100.0});

  ResponseSpec spec;
  spec.metric = Metric::Activity;
  spec.aggregate = Aggregate::GammaLine;
  spec.unit = Unit::Becquerel;
  try {
    buildResponse(data, result, spec);
    FAIL() << "expected InputError";
  } catch (const InputError& e) {
    EXPECT_NE(std::string(e.what()).find("no activity of its own"), std::string::npos) << e.what();
  }
}

// Integrated over a window, per-line columns must still sum to the per-nuclide total -- the
// same invariant, in the other domain.
TEST(ResponseLines, IntervalLineTotalMatchesTheNuclideTotal) {
  const NuclearData data = chainWithStrongAndTraceLines();
  Inventory inv;
  inv.add(Zai{50, 100, 0}, 1.0e20);

  std::vector<std::int64_t> keys;
  const std::vector<double> integral = intervalIntegral(data, inv, 0.0, 2000.0, &keys);

  ResponseSpec byNuclide;
  byNuclide.metric = Metric::Exposure;
  byNuclide.unit = Unit::Roentgen;
  ResponseSpec byLine = byNuclide;
  byLine.aggregate = Aggregate::GammaLine;

  const double nuclideTotal =
      buildIntervalResponse(data, keys, integral, 0.0, 2000.0, byNuclide).totals[0];
  const double lineTotal =
      buildIntervalResponse(data, keys, integral, 0.0, 2000.0, byLine).totals[0];
  EXPECT_NEAR(lineTotal, nuclideTotal, nuclideTotal * 1e-12);
}

// Matching totals is not enough to say an interval was aggregated by line: per-nuclide columns
// sum to exactly the same number. The columns themselves have to be lines, carrying an energy
// and a labeled emitter, or `--by line` over a window silently answers a different question
// from the same flag over an instant.
TEST(ResponseLines, IntervalColumnsAreLinesRatherThanNuclides) {
  const NuclearData data = chainWithStrongAndTraceLines();
  Inventory inv;
  inv.add(Zai{50, 100, 0}, 1.0e20);

  std::vector<std::int64_t> keys;
  const std::vector<double> integral = intervalIntegral(data, inv, 0.0, 2000.0, &keys);

  ResponseSpec byLine;
  byLine.metric = Metric::Exposure;
  byLine.unit = Unit::Roentgen;
  byLine.aggregate = Aggregate::GammaLine;
  const ResponseTable table = buildIntervalResponse(data, keys, integral, 0.0, 2000.0, byLine);

  ResponseSpec instantByLine = byLine;
  instantByLine.unit = Unit::RoentgenPerHour;
  const ResponseTable instant =
      buildResponse(data, decay(data, inv, std::vector<double>{2000.0}), instantByLine);

  ASSERT_EQ(table.contributorCount(), instant.contributorCount());
  ASSERT_GT(table.contributorCount(), 0);
  for (int c = 0; c < table.contributorCount(); ++c) {
    const std::size_t i = static_cast<std::size_t>(c);
    EXPECT_GT(table.contributors[i].lineEnergyEv, 0.0) << "column " << c << " carries no energy";
    EXPECT_NE(table.labels[i].find("keV"), std::string::npos) << table.labels[i];
  }
  EXPECT_NE(table.labels[0].find("Sb-100"), std::string::npos) << table.labels[0];
}

// --- how much is missing, not just how many are flagged -----------------------

// A count of flagged nuclides says nothing about magnitude: a hundred negligible ones and
// three dominant ones look identical. The energy fraction is what separates them, and it is
// activity-weighted so a large unmodeled fraction on a negligible nuclide stays negligible.
TEST(ResponseExposure, ReportsHowMuchPhotonEnergyIsUnmodelled) {
  StoreArrays arrays = synth::linearChain({1.0e-3, 5.0e-4});
  // The emitter carries a 1 MeV line at unit intensity, plus an equal amount of continuum.
  // Half its photon energy is therefore outside the model.
  arrays.emEnergyEv = {0.0, 2.0e6, 0.0};
  arrays.continuumPhotonEv = {0.0, 1.0e6, 0.0};
  synth::addLines(arrays, 1, {1.0e6}, {1.0});
  const NuclearData data = NuclearData::fromArrays(std::move(arrays));

  Inventory inv;
  inv.add(Zai{50, 100, 0}, 1.0e20);
  const DecayResult result = decay(data, inv, std::vector<double>{2000.0});

  ResponseSpec spec;
  spec.metric = Metric::Exposure;
  spec.unit = Unit::RoentgenPerHour;
  const ResponseTable table = buildResponse(data, result, spec);

  ASSERT_EQ(table.unmodeledEnergyFraction.size(), 1u);
  EXPECT_NEAR(table.unmodeledEnergyFraction[0], 0.5, 1e-9);

  // And it reaches the ranking, which is what a report reads.
  EXPECT_NEAR(rank(table, 0, RankRequest{}).unmodeledEnergyFraction, 0.5, 1e-9);
}

// Nothing unmodeled means nothing to warn about, rather than a zero that still prints.
TEST(ResponseExposure, UnmodelledFractionIsZeroWhenEveryPhotonIsAccountedFor) {
  StoreArrays arrays = synth::linearChain({1.0e-3, 5.0e-4});
  arrays.emEnergyEv = {0.0, 1.0e6, 0.0};
  synth::addLines(arrays, 1, {1.0e6}, {1.0});
  const NuclearData data = NuclearData::fromArrays(std::move(arrays));

  Inventory inv;
  inv.add(Zai{50, 100, 0}, 1.0e20);
  const DecayResult result = decay(data, inv, std::vector<double>{2000.0});

  ResponseSpec spec;
  spec.metric = Metric::Exposure;
  spec.unit = Unit::RoentgenPerHour;
  EXPECT_DOUBLE_EQ(buildResponse(data, result, spec).unmodeledEnergyFraction[0], 0.0);
}

// Activity has no photon model, so there is nothing to be missing from it.
TEST(ResponseExposure, ActivityCarriesNoUnmodelledFraction) {
  const NuclearData data = chainWithOnePhotonEmitter();
  Inventory inv;
  inv.add(Zai{50, 100, 0}, 1.0e20);
  const DecayResult result = decay(data, inv, std::vector<double>{2000.0});
  EXPECT_TRUE(buildResponse(data, result, ResponseSpec{}).unmodeledEnergyFraction.empty());
}

// The per-contributor FLAG has to agree with the fraction about when it applies. It says a
// photon spectrum is incomplete, which understates an exposure and says nothing whatever about
// a count of decays -- and a report builds its footnote by scanning these flags, so setting
// them here put a paragraph about understated exposures under every activity ranking.
TEST(ResponseExposure, ActivityCarriesNoUnmodelledFlagEither) {
  StoreArrays arrays = synth::linearChain({1.0e-3, 5.0e-4});
  arrays.emEnergyEv = {0.0, 2.0e6, 0.0};
  arrays.continuumPhotonEv = {0.0, 1.0e6, 0.0};  // half the emitter's photon energy
  synth::addLines(arrays, 1, {1.0e6}, {1.0});
  const NuclearData data = NuclearData::fromArrays(std::move(arrays));

  Inventory inv;
  inv.add(Zai{50, 100, 0}, 1.0e20);
  const DecayResult result = decay(data, inv, std::vector<double>{2000.0});

  ResponseSpec exposureSpec;
  exposureSpec.metric = Metric::Exposure;
  exposureSpec.unit = Unit::RoentgenPerHour;
  int flagged = 0;
  for (const int flags : buildResponse(data, result, exposureSpec).flags) {
    flagged += (flags & kFlagUnmodeledContinuum) != 0 ? 1 : 0;
  }
  ASSERT_EQ(flagged, 1) << "the emitter should be flagged when the metric is exposure";

  for (const int flags : buildResponse(data, result, ResponseSpec{}).flags) {
    EXPECT_EQ(flags & kFlagUnmodeledContinuum, 0);
  }
}

// --- photon metric -----------------------------------------------------------
//
// The third metric exists because activity answers "how many decays" and exposure answers
// "how much dose", and neither says how many photons actually leave the source. Strength is
// the only geometry-free metric, which is what makes it the one to compare inventories with
// before a site has even been chosen; fluence is the same transport quoted at a point.

// Strength is activity weighted by the photons each decay emits. Checked against the
// Bateman solution at three times, including t = 0 where the daughter has not grown in:
// a photon table that is nonzero there has an initialisation bug that later times hide.
TEST(ResponsePhoton, StrengthEqualsActivityTimesTheDiscreteYield) {
  const NuclearData data = chainWithOnePhotonEmitter();
  const double lambda0 = 1.0e-3;
  const double lambda1 = 5.0e-4;
  const double n0 = 1.0e20;
  Inventory inv;
  inv.add(Zai{50, 100, 0}, n0);

  ResponseSpec spec;
  spec.metric = Metric::Photon;
  spec.unit = Unit::PhotonsPerSecond;

  for (const double t : {0.0, 100.0, 2000.0}) {
    // The daughter is the sole emitter and its only line carries 0.9 photons per decay, so
    // the whole table is lambda * yield * N1, no geometry anywhere in the expression.
    const double expected = lambda1 * 0.9 * synth::batemanN1(n0, lambda0, lambda1, t);
    const double total =
        buildResponse(data, decay(data, inv, std::vector<double>{t}), spec).totals[0];
    EXPECT_NEAR(total, expected, std::max(1.0, expected) * 1e-12) << "at t = " << t;
  }
}

// Nothing in the spec touches the distance, the buildup, or the air density, and the answer
// must not move: strength is a property of the inventory alone. The geometry record still
// has to be usable -- the refusal is about the geometry being nonsensical, not about it
// being used.
TEST(ResponsePhoton, StrengthDoesNotDependOnTheGeometry) {
  const NuclearData data = chainWithOnePhotonEmitter();
  Inventory inv;
  inv.add(Zai{50, 100, 0}, 1.0e20);
  const DecayResult result = decay(data, inv, std::vector<double>{2000.0});

  ResponseSpec near;
  near.metric = Metric::Photon;
  near.unit = Unit::PhotonsPerSecond;
  near.geometry.distanceM = 1.0;

  ResponseSpec far = near;
  far.geometry.distanceM = 100.0;
  far.geometry.buildup = 3.0;
  far.geometry.airDensityKgM3 = 0.5;
  far.geometry.airAttenuation = false;

  EXPECT_DOUBLE_EQ(buildResponse(data, result, near).totals[0],
                   buildResponse(data, result, far).totals[0]);
}

// Fluence is where the geometry is part of the answer. At one distance it must equal the
// hand calculation: yield-weighted activity over 4 pi r^2, times the exponential
// attenuation, times the buildup. The vacuum variant pins the geometric term separately, so
// a bug in the transport factor cannot hide behind it.
TEST(ResponsePhoton, FluenceMatchesSpreadingAttenuationAndBuildup) {
  const NuclearData data = chainWithOnePhotonEmitter();
  const double lambda0 = 1.0e-3;
  const double lambda1 = 5.0e-4;
  const double n0 = 1.0e20;
  Inventory inv;
  inv.add(Zai{50, 100, 0}, n0);
  const double t = 2000.0;
  const DecayResult result = decay(data, inv, std::vector<double>{t});

  const double d = 3.0;
  const double buildup = 2.5;

  ResponseSpec spec;
  spec.metric = Metric::Photon;
  spec.unit = Unit::PhotonsPerSquareMeterPerSecond;
  spec.geometry.distanceM = d;
  spec.geometry.buildup = buildup;

  const double activity = lambda1 * synth::batemanN1(n0, lambda0, lambda1, t);
  const double withAir =
      activity * 0.9 / (4.0 * M_PI * d * d) *
      std::exp(-exposure::airMassAttenuation(661657.0) * spec.geometry.airDensityKgM3 * d) *
      buildup;
  EXPECT_NEAR(buildResponse(data, result, spec).totals[0], withAir, withAir * 1e-12);

  spec.geometry.airAttenuation = false;
  const double inVacuum = activity * 0.9 / (4.0 * M_PI * d * d) * buildup;
  EXPECT_NEAR(buildResponse(data, result, spec).totals[0], inVacuum, inVacuum * 1e-12);
}

// Over a window, photons are counts: every atom that decayed emitted its yield, and the
// window is in seconds. Dividing by an hour (or failing to) moves the total by 3600 --
// invisible in a table of bare figures, and visible only if the units are read.
TEST(ResponsePhoton, IntervalAccruesPhotonsPerDecayNotPerHour) {
  const double lambda = 1.0e-9;  // ~5000 y, so almost nothing decays over the hour
  const double n0 = 1.0e20;
  StoreArrays arrays = synth::linearChain({lambda});
  synth::addLines(arrays, 0, {661657.0}, {0.9});  // the seeded nuclide is the emitter
  const NuclearData data = NuclearData::fromArrays(std::move(arrays));

  Inventory inv;
  inv.add(Zai{50, 100, 0}, n0);
  const double window = units::kSecondsPerHour;

  ResponseSpec rateSpec;
  rateSpec.metric = Metric::Photon;
  rateSpec.unit = Unit::PhotonsPerSecond;
  const double rate =
      buildResponse(data, decay(data, inv, std::vector<double>{0.0}), rateSpec).totals[0];

  std::vector<std::int64_t> keys;
  const std::vector<double> integral = intervalIntegral(data, inv, 0.0, window, &keys);
  ResponseSpec countSpec = rateSpec;
  countSpec.unit = Unit::Photons;
  const double count =
      buildIntervalResponse(data, keys, integral, 0.0, window, countSpec).totals[0];

  // The exact integral: lambda * yield * (atom-seconds over the window).
  const double expected = lambda * 0.9 * synth::batemanIntegralN0(n0, lambda, window);
  EXPECT_NEAR(count, expected, expected * 1e-9);
  EXPECT_LT(count, rate * window)
      << "an hour of decay emits slightly fewer photons than the initial rate";
}

// Splitting a nuclide's photons across its lines must not change the total, for the strength
// unit and the fluence unit alike -- the same invariant the per-nuclide exposure test pins,
// in the other two unit families, dropped columns included.
TEST(ResponsePhotonLines, LineTotalMatchesTheNuclideTotal) {
  const NuclearData data = chainWithStrongAndTraceLines();
  Inventory inv;
  inv.add(Zai{50, 100, 0}, 1.0e20);
  const DecayResult result = decay(data, inv, std::vector<double>{2000.0});

  for (const Unit unit : {Unit::PhotonsPerSecond, Unit::PhotonsPerSquareMeterPerSecond}) {
    ResponseSpec byNuclide;
    byNuclide.metric = Metric::Photon;
    byNuclide.unit = unit;
    ResponseSpec byLine = byNuclide;
    byLine.aggregate = Aggregate::GammaLine;

    const double nuclideTotal = buildResponse(data, result, byNuclide).totals[0];
    const double lineTotal = buildResponse(data, result, byLine).totals[0];
    EXPECT_NEAR(lineTotal, nuclideTotal, nuclideTotal * 1e-12) << unitName(unit);
  }
}

// Strength ranks lines by intensity alone, so the 662 keV line at 0.9 leads and the 1332 keV
// line at 0.5 follows -- the reverse of exposure's order, where the energy-weighted line
// leads. Same data, different question; a test that only checked the total would pass even
// if the lines were ranked by the wrong weight.
TEST(ResponsePhotonLines, RanksByIntensityNotEnergyWeightedIntensity) {
  const NuclearData data = chainWithStrongAndTraceLines();
  Inventory inv;
  inv.add(Zai{50, 100, 0}, 1.0e20);
  const DecayResult result = decay(data, inv, std::vector<double>{2000.0});

  ResponseSpec spec;
  spec.metric = Metric::Photon;
  spec.unit = Unit::PhotonsPerSecond;
  spec.aggregate = Aggregate::GammaLine;
  const Ranking ranking = rank(buildResponse(data, result, spec), 0, RankRequest{});

  ASSERT_GE(ranking.contributors.size(), 2u);
  EXPECT_NE(ranking.contributors[0].label.find("661.7 keV"), std::string::npos)
      << ranking.contributors[0].label;
  EXPECT_NE(ranking.contributors[1].label.find("1332.5 keV"), std::string::npos)
      << ranking.contributors[1].label;
  // Same emitter, same time: the values differ by the intensities alone.
  EXPECT_NEAR(ranking.contributors[0].value / ranking.contributors[1].value, 0.9 / 0.5, 1e-12);
}

// Fluence is the unit that carries the geometry, so an unusable geometry must fail here.
// Strength never touches the geometry, so the same spec with a strength unit still works:
// the refusal attaches to the fluence quantity, not to the geometry record.
TEST(ResponsePhoton, FluenceRefusesAnImpossibleGeometryButStrengthDoesNot) {
  const NuclearData data = chainWithOnePhotonEmitter();
  Inventory inv;
  inv.add(Zai{50, 100, 0}, 1.0e20);
  const DecayResult result = decay(data, inv, std::vector<double>{100.0});

  ResponseSpec fluence;
  fluence.metric = Metric::Photon;
  fluence.unit = Unit::PhotonsPerSquareMeterPerSecond;
  fluence.geometry.distanceM = 0.0;
  EXPECT_THROW(buildResponse(data, result, fluence), InputError);

  ResponseSpec strength = fluence;
  strength.unit = Unit::PhotonsPerSecond;
  const ResponseTable table = buildResponse(data, result, strength);
  EXPECT_GT(table.totals[0], 0.0);
}

// A store with no photon lines cannot answer a photon question any more than an exposure
// one: zeros would be indistinguishable from "nothing here emits photons".
TEST(ResponsePhoton, RefusesAStoreWithNoPhotonLines) {
  const NuclearData data = NuclearData::fromArrays(synth::linearChain({1.0e-3}));
  ASSERT_FALSE(data.hasPhotonLines());
  Inventory inv;
  inv.add(Zai{50, 100, 0}, 1.0e20);
  const DecayResult result = decay(data, inv, std::vector<double>{100.0});

  for (const Unit unit : {Unit::PhotonsPerSecond, Unit::PhotonsPerSquareMeterPerSecond}) {
    ResponseSpec spec;
    spec.metric = Metric::Photon;
    spec.unit = unit;
    try {
      buildResponse(data, result, spec);
      FAIL() << "expected InputError for " << unitName(unit);
    } catch (const InputError& e) {
      EXPECT_NE(std::string(e.what()).find("no photon lines"), std::string::npos) << e.what();
    }
  }
}

// Photons do not measure exposure and roentgen does not measure photons: category errors,
// refused at the boundary. The rate-versus-count line the other unit families get applies
// to the photon units the same way.
TEST(ResponsePhoton, RefusesAUnitThatDoesNotMeasureTheMetric) {
  const NuclearData data = chainWithOnePhotonEmitter();
  Inventory inv;
  inv.add(Zai{50, 100, 0}, 1.0e20);
  const DecayResult result = decay(data, inv, std::vector<double>{100.0});

  ResponseSpec photonsInRoentgen;
  photonsInRoentgen.metric = Metric::Photon;
  photonsInRoentgen.unit = Unit::RoentgenPerHour;
  EXPECT_THROW(buildResponse(data, result, photonsInRoentgen), InputError);

  ResponseSpec exposureInPhotons;
  exposureInPhotons.metric = Metric::Exposure;
  exposureInPhotons.unit = Unit::PhotonsPerSecond;
  EXPECT_THROW(buildResponse(data, result, exposureInPhotons), InputError);

  ResponseSpec activityInPhotons;
  activityInPhotons.metric = Metric::Activity;
  activityInPhotons.unit = Unit::Photons;
  EXPECT_THROW(buildResponse(data, result, activityInPhotons), InputError);

  std::vector<std::int64_t> keys;
  const std::vector<double> integral = intervalIntegral(data, inv, 0.0, 100.0, &keys);
  // The counts are interval-only and the rates instant-only, as for every other family.
  for (const Unit intervalOnly : {Unit::Photons, Unit::PhotonsPerSquareMeter}) {
    ResponseSpec countAtAnInstant;
    countAtAnInstant.metric = Metric::Photon;
    countAtAnInstant.unit = intervalOnly;
    EXPECT_THROW(buildResponse(data, result, countAtAnInstant), InputError);
  }
  for (const Unit instantOnly : {Unit::PhotonsPerSecond, Unit::PhotonsPerSquareMeterPerSecond}) {
    ResponseSpec rateOverAWindow;
    rateOverAWindow.metric = Metric::Photon;
    rateOverAWindow.unit = instantOnly;
    EXPECT_THROW(buildIntervalResponse(data, keys, integral, 0.0, 100.0, rateOverAWindow),
                 InputError);
  }
}

// An unmodeled continuum understates a photon count exactly as it understates an exposure:
// the flag, the energy fraction, and the ranking must all carry it. And the optical-depth
// caveat attaches only to the quantity that actually used the point geometry -- a strength
// number has no path to be thick, a fluence number carries the path it was built in.
TEST(ResponsePhoton, CarriesTheCaveatsTheMetricActuallyUsed) {
  StoreArrays arrays = synth::linearChain({1.0e-3, 5.0e-4});
  arrays.emEnergyEv = {0.0, 2.0e6, 0.0};
  arrays.continuumPhotonEv = {0.0, 1.0e6, 0.0};  // half the emitter's photon energy
  synth::addLines(arrays, 1, {1.0e6}, {1.0});
  const NuclearData data = NuclearData::fromArrays(std::move(arrays));

  Inventory inv;
  inv.add(Zai{50, 100, 0}, 1.0e20);
  const DecayResult result = decay(data, inv, std::vector<double>{2000.0});

  ResponseSpec strength;
  strength.metric = Metric::Photon;
  strength.unit = Unit::PhotonsPerSecond;
  const ResponseTable strengthTable = buildResponse(data, result, strength);

  int flagged = 0;
  for (const int flags : strengthTable.flags) {
    flagged += (flags & kFlagUnmodeledContinuum) != 0 ? 1 : 0;
  }
  ASSERT_EQ(flagged, 1) << "the emitter should be flagged when the metric is photon strength";
  ASSERT_EQ(strengthTable.unmodeledEnergyFraction.size(), 1u);
  EXPECT_NEAR(strengthTable.unmodeledEnergyFraction[0], 0.5, 1e-9);
  EXPECT_NEAR(rank(strengthTable, 0, RankRequest{}).unmodeledEnergyFraction, 0.5, 1e-9);
  EXPECT_TRUE(strengthTable.meanOpticalDepth.empty())
      << "a strength number at no distance has no air path to be thick";

  ResponseSpec fluence = strength;
  fluence.unit = Unit::PhotonsPerSquareMeterPerSecond;
  fluence.geometry.distanceM = 100.0;
  const ResponseTable fluenceTable = buildResponse(data, result, fluence);

  flagged = 0;
  for (const int flags : fluenceTable.flags) {
    flagged += (flags & kFlagUnmodeledContinuum) != 0 ? 1 : 0;
  }
  ASSERT_EQ(flagged, 1) << "the emitter should be flagged when the metric is photon fluence";
  ASSERT_EQ(fluenceTable.meanOpticalDepth.size(), 1u);
  const double expectedDepth =
      exposure::airMassAttenuation(1.0e6) * fluence.geometry.airDensityKgM3 * 100.0;
  EXPECT_NEAR(fluenceTable.meanOpticalDepth[0], expectedDepth, expectedDepth * 1e-12);

  // Activity remains caveat-free: the flags describe a photon model it never ran.
  for (const int flags : buildResponse(data, result, ResponseSpec{}).flags) {
    EXPECT_EQ(flags & kFlagUnmodeledContinuum, 0);
  }
}

// --- unit spellings -----------------------------------------------------------

// The CLI and the Python binding each carried their own table of spellings, and they had
// drifted: `--units bq` worked while `units="bq"` raised, for no reason a user could see. One
// table, derived from the names the reports print, is what makes "a notebook and a terminal
// never disagree" structural rather than a promise.
TEST(Units, AcceptEverySpellingTheyPrint) {
  const Unit all[] = {Unit::Becquerel,
                      Unit::Curie,
                      Unit::Decays,
                      Unit::RoentgenPerHour,
                      Unit::GrayPerHour,
                      Unit::SievertPerHour,
                      Unit::Roentgen,
                      Unit::Gray,
                      Unit::Sievert,
                      Unit::PhotonsPerSecond,
                      Unit::Photons,
                      Unit::PhotonsPerSquareMeterPerSecond,
                      Unit::PhotonsPerSquareMeter};
  for (const Unit unit : all) {
    Unit parsed = Unit::Becquerel;
    ASSERT_TRUE(parseUnit(unitName(unit), parsed)) << unitName(unit);
    EXPECT_EQ(parsed, unit) << unitName(unit);
  }
}

TEST(Units, AreSpelledCaseInsensitively) {
  Unit parsed = Unit::Becquerel;
  ASSERT_TRUE(parseUnit("bq", parsed));
  EXPECT_EQ(parsed, Unit::Becquerel);
  ASSERT_TRUE(parseUnit("gy/h", parsed));
  EXPECT_EQ(parsed, Unit::GrayPerHour);
  ASSERT_TRUE(parseUnit("sv", parsed));
  EXPECT_EQ(parsed, Unit::Sievert);
  ASSERT_TRUE(parseUnit("DECAYS", parsed));
  EXPECT_EQ(parsed, Unit::Decays);
}

TEST(Units, RejectWhatIsNotAUnit) {
  Unit parsed = Unit::Becquerel;
  EXPECT_FALSE(parseUnit("rem", parsed));
  EXPECT_FALSE(parseUnit("", parsed));
  EXPECT_FALSE(parseUnit("Bq/h", parsed));
}

// Naming no unit is the ordinary invocation, so the default has to satisfy both the metric and
// the domain -- otherwise the simplest command fails on a unit nobody chose.
TEST(Units, DefaultSatisfiesBothTheMetricAndTheDomain) {
  for (const Metric metric : {Metric::Activity, Metric::Exposure, Metric::Photon}) {
    for (const Domain domain : {Domain::Instant, Domain::Interval}) {
      const Unit unit = defaultUnit(metric, domain);
      EXPECT_TRUE(unitSuitsMetric(unit, metric)) << unitName(unit);
      EXPECT_TRUE(unitSuitsDomain(unit, domain)) << unitName(unit);
      EXPECT_EQ(requireUnit("", metric, domain), unit);
    }
  }
  // The photon default is the STRENGTH unit in both domains: a photon answer that does not
  // name a distance must not silently become one at some distance.
  EXPECT_EQ(defaultUnit(Metric::Photon, Domain::Instant), Unit::PhotonsPerSecond);
  EXPECT_EQ(defaultUnit(Metric::Photon, Domain::Interval), Unit::Photons);
  EXPECT_FALSE(isFluenceUnit(defaultUnit(Metric::Photon, Domain::Instant)));
  EXPECT_FALSE(isFluenceUnit(defaultUnit(Metric::Photon, Domain::Interval)));
  EXPECT_TRUE(isFluenceUnit(Unit::PhotonsPerSquareMeterPerSecond));
  EXPECT_TRUE(isFluenceUnit(Unit::PhotonsPerSquareMeter));
}

// The message has to name the spellings that would have worked; "not a unit" alone leaves the
// user guessing at the one thing the error knows.
TEST(Units, RequireThrowsNamingTheSpellingsThatWouldHaveWorked) {
  try {
    requireUnit("rem", Metric::Exposure, Domain::Instant);
    FAIL() << "expected InputError";
  } catch (const InputError& e) {
    const std::string what = e.what();
    EXPECT_NE(what.find("rem"), std::string::npos) << what;
    EXPECT_NE(what.find("Gy/h"), std::string::npos) << what;
    EXPECT_NE(what.find("decays"), std::string::npos) << what;
  }

  try {
    requireUnit("rem", Metric::Photon, Domain::Instant);
    FAIL() << "expected InputError";
  } catch (const InputError& e) {
    const std::string what = e.what();
    EXPECT_NE(what.find("rem"), std::string::npos) << what;
    EXPECT_NE(what.find("photons/s"), std::string::npos) << what;
    EXPECT_NE(what.find("photons/m2/s"), std::string::npos) << what;
  }
}

}  // namespace
}  // namespace nusift
