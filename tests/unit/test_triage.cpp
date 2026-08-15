#include <gtest/gtest.h>

#include <string>
#include <vector>

#include "nusift/core/error.hpp"
#include "nusift/engine/decay_engine.hpp"
#include "nusift/exposure/point_source.hpp"
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

// Sv/h is R/h times the air-kerma conversion, applied once at the very end.
TEST(ResponseExposure, SievertIsTheRoentgenValueConverted) {
  const NuclearData data = chainWithOnePhotonEmitter();
  Inventory inv;
  inv.add(Zai{50, 100, 0}, 1.0e20);
  const DecayResult result = decay(data, inv, std::vector<double>{2000.0});

  ResponseSpec roentgen;
  roentgen.metric = Metric::Exposure;
  roentgen.unit = Unit::RoentgenPerHour;
  ResponseSpec sievert = roentgen;
  sievert.unit = Unit::SievertPerHour;

  const double inR = buildResponse(data, result, roentgen).totals[0];
  const double inSv = buildResponse(data, result, sievert).totals[0];
  EXPECT_NEAR(inSv, inR * exposure::kGrayPerRoentgen, inSv * 1e-12);
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

// An emitter with three lines of very different strength, so thresholding and ordering are
// both observable.
NuclearData chainWithThreeLines() {
  StoreArrays arrays = synth::linearChain({1.0e-3, 5.0e-4});
  synth::addLines(arrays, 1, {1332492.0, 661657.0, 100.0}, {0.5, 0.9, 1.0e-9});
  return NuclearData::fromArrays(std::move(arrays));
}

// Splitting a nuclide's exposure across its lines must not change the total. This is the
// invariant that says per-line aggregation is a finer view of the same quantity rather than a
// second, separately-derived calculation.
TEST(ResponseLines, LineTotalMatchesTheNuclideTotal) {
  const NuclearData data = chainWithThreeLines();
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
  EXPECT_NEAR(lineTotal, nuclideTotal, nuclideTotal * 1e-9);
}

// The strongest line by emitted exposure leads, and every column names its emitter and energy.
TEST(ResponseLines, RanksIndividualLinesAndLabelsThem) {
  const NuclearData data = chainWithThreeLines();
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
  const NuclearData data = chainWithThreeLines();
  Inventory inv;
  inv.add(Zai{50, 100, 0}, 1.0e20);
  const DecayResult result = decay(data, inv, std::vector<double>{2000.0});

  ResponseSpec spec;
  spec.metric = Metric::Exposure;
  spec.unit = Unit::RoentgenPerHour;
  spec.aggregate = Aggregate::GammaLine;
  const ResponseTable table = buildResponse(data, result, spec);

  // Three lines were staged; the 100 eV line at 1e-9 intensity is far below the floor.
  EXPECT_EQ(table.contributorCount(), 2);
}

// A photon line has no activity of its own -- it is a way its emitter's decays get out. Asking
// for activity by line is a category error, not a rounding one.
TEST(ResponseLines, RefusesActivityRankedByLine) {
  const NuclearData data = chainWithThreeLines();
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
  const NuclearData data = chainWithThreeLines();
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
  EXPECT_NEAR(lineTotal, nuclideTotal, nuclideTotal * 1e-9);
}

// Matching totals is not enough to say an interval was aggregated by line: per-nuclide columns
// sum to exactly the same number. The columns themselves have to be lines, carrying an energy
// and a labelled emitter, or `--by line` over a window silently answers a different question
// from the same flag over an instant.
TEST(ResponseLines, IntervalColumnsAreLinesRatherThanNuclides) {
  const NuclearData data = chainWithThreeLines();
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
// activity-weighted so a large unmodelled fraction on a negligible nuclide stays negligible.
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

// Nothing unmodelled means nothing to warn about, rather than a zero that still prints.
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

// --- unit spellings -----------------------------------------------------------

// The CLI and the Python binding each carried their own table of spellings, and they had
// drifted: `--units bq` worked while `units="bq"` raised, for no reason a user could see. One
// table, derived from the names the reports print, is what makes "a notebook and a terminal
// never disagree" structural rather than a promise.
TEST(Units, AcceptEverySpellingTheyPrint) {
  const Unit all[] = {Unit::Becquerel,       Unit::Curie,       Unit::Decays,
                      Unit::RoentgenPerHour, Unit::GrayPerHour, Unit::SievertPerHour,
                      Unit::Roentgen,        Unit::Gray,        Unit::Sievert};
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
  for (const Metric metric : {Metric::Activity, Metric::Exposure}) {
    for (const Domain domain : {Domain::Instant, Domain::Interval}) {
      const Unit unit = defaultUnit(metric, domain);
      EXPECT_TRUE(unitSuitsMetric(unit, metric)) << unitName(unit);
      EXPECT_TRUE(unitSuitsDomain(unit, domain)) << unitName(unit);
      EXPECT_EQ(requireUnit("", metric, domain), unit);
    }
  }
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
}

}  // namespace
}  // namespace nusift
