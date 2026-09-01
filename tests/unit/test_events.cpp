#include <gtest/gtest.h>

#include <cmath>
#include <functional>
#include <vector>

#include "nusift/core/error.hpp"
#include "nusift/engine/decay_engine.hpp"
#include "nusift/io/time_spec.hpp"
#include "nusift/nucdata/nuclear_data.hpp"
#include "nusift/triage/events.hpp"
#include "synthetic_chain.hpp"

namespace nusift {
namespace {

// Sample a closed-form curve on a grid. `refine` decides whether the series can be evaluated
// between samples, which is the whole difference between a located event and an interpolated
// one -- so nearly every test here runs both ways.
EventSeries sampled(const std::function<double(double)>& curve, const std::vector<double>& times,
                    bool refine) {
  EventSeries series;
  series.times = times;
  for (const double t : times) {
    series.values.push_back(curve(t));
  }
  if (refine) {
    series.evaluate = curve;
  }
  return series;
}

std::vector<double> linspace(double from, double to, int count) {
  std::vector<double> times;
  times.reserve(static_cast<std::size_t>(count));
  for (int i = 0; i < count; ++i) {
    times.push_back(from + (to - from) * i / (count - 1));
  }
  return times;
}

// --- crossings ---------------------------------------------------------------

// N e^{-lambda t} = L at t = ln(N / L) / lambda. A closed form to check against, for the same
// reason the forecast boundary test uses one: a grid point would only prove the code returns a
// grid point.
TEST(TrajectoryEvents, CrossingMatchesTheAnalyticRoot) {
  const double lambda = 1.0e-3;
  const double n0 = 1.0e6;
  const double level = 1.0e3;
  const double expected = std::log(n0 / level) / lambda;

  const auto curve = [&](double t) { return n0 * std::exp(-lambda * t); };
  const std::vector<double> times = logspace(1.0, 1.0e5, 40);

  for (const bool refine : {false, true}) {
    const std::vector<TrajectoryEvent> found = crossings(sampled(curve, times, refine), level);
    ASSERT_EQ(found.size(), 1u) << "one monotone decay crosses one level once";
    EXPECT_EQ(found[0].kind, EventKind::Falling);
    // Log-linear interpolation is exact for a pure exponential, so both paths land on the
    // root; what separates them is how well each one can claim to know it.
    EXPECT_NEAR(found[0].timeSeconds, expected, expected * 1.0e-6);
    EXPECT_EQ(found[0].refined, refine);
  }
}

// The bracket is the honest error bar. Without an evaluator it is the grid interval; with one
// it is whatever the refinement achieved, which on a decade-wide log grid is orders better.
TEST(TrajectoryEvents, RefinementNarrowsTheBracketAndInterpolationDoesNot) {
  const auto curve = [](double t) { return 1.0e6 * std::exp(-1.0e-3 * t); };
  const std::vector<double> times = logspace(1.0, 1.0e5, 40);

  const TrajectoryEvent interpolated = crossings(sampled(curve, times, false), 1.0e3).front();
  const TrajectoryEvent refined = crossings(sampled(curve, times, true), 1.0e3).front();

  const double bracket = interpolated.bracketEndSeconds - interpolated.bracketStartSeconds;
  ASSERT_GT(bracket, 0.0);

  EXPECT_DOUBLE_EQ(interpolated.locatedToSeconds, bracket)
      << "an interpolated event is known no better than the interval it sits in";
  EXPECT_FALSE(interpolated.refined);

  EXPECT_TRUE(refined.refined);
  EXPECT_TRUE(refined.converged);
  EXPECT_LT(refined.locatedToSeconds, bracket * 1.0e-3) << "refinement has to actually narrow it";
  EXPECT_LE(refined.locatedToSeconds, refined.timeSeconds * 1.0e-6 + 1.0e-12)
      << "and reach the tolerance it was asked for";

  // Both stay inside the bracket the grid observed -- refinement never leaves it.
  EXPECT_GE(refined.timeSeconds, refined.bracketStartSeconds);
  EXPECT_LE(refined.timeSeconds, refined.bracketEndSeconds);
}

// The documented missed-event rule, asserted rather than trusted: an excursion that rises and
// falls back between two consecutive samples leaves no sign change and is NOT found. Supplying
// an evaluator does not change that -- detection is by sample, and only refinement evaluates.
TEST(TrajectoryEvents, AnExcursionBetweenTwoSamplesIsNotFound) {
  // Flat at 1, with a narrow spike to 100 centred at t = 15, invisible to samples 10 apart.
  const auto curve = [](double t) { return 1.0 + 99.0 * std::exp(-((t - 15.0) * (t - 15.0))); };
  const std::vector<double> times = {0.0, 10.0, 20.0, 30.0, 40.0};

  ASSERT_GT(curve(15.0), 50.0) << "the spike really does cross the level between the samples";
  ASSERT_LT(curve(10.0), 50.0);
  ASSERT_LT(curve(20.0), 50.0);

  for (const bool refine : {false, true}) {
    EXPECT_TRUE(crossings(sampled(curve, times, refine), 50.0).empty())
        << "the grid did not bracket it, so nothing may be invented between the samples";
  }

  // A grid that does sample the spike finds both of its crossings.
  const std::vector<double> dense = linspace(0.0, 40.0, 81);
  EXPECT_EQ(crossings(sampled(curve, dense, true), 50.0).size(), 2u);
}

TEST(TrajectoryEvents, RisingAndFallingAreDistinguishedAndOrderedInTime) {
  const auto curve = [](double t) {
    return 1.0 + 99.0 * std::exp(-((t - 15.0) * (t - 15.0)) / 9.0);
  };
  const std::vector<double> times = linspace(0.0, 40.0, 161);

  const std::vector<TrajectoryEvent> found = crossings(sampled(curve, times, true), 50.0);
  ASSERT_EQ(found.size(), 2u);
  EXPECT_EQ(found[0].kind, EventKind::Rising);
  EXPECT_EQ(found[1].kind, EventKind::Falling);
  EXPECT_LT(found[0].timeSeconds, found[1].timeSeconds);
  EXPECT_NEAR(found[0].value, 50.0, 1.0e-9) << "a crossing's value is the level by definition";
}

TEST(TrajectoryEvents, FirstAndLastCrossingAgreeWithTheFullList) {
  const auto curve = [](double t) {
    return 1.0 + 99.0 * std::exp(-((t - 15.0) * (t - 15.0)) / 9.0);
  };
  const EventSeries series = sampled(curve, linspace(0.0, 40.0, 161), true);

  const std::vector<TrajectoryEvent> all = crossings(series, 50.0);
  ASSERT_EQ(all.size(), 2u);

  const std::optional<TrajectoryEvent> first = firstCrossing(series, 50.0);
  const std::optional<TrajectoryEvent> last = lastCrossing(series, 50.0);
  ASSERT_TRUE(first.has_value());
  ASSERT_TRUE(last.has_value());
  EXPECT_DOUBLE_EQ(first->timeSeconds, all.front().timeSeconds);
  EXPECT_DOUBLE_EQ(last->timeSeconds, all.back().timeSeconds);
}

TEST(TrajectoryEvents, ALevelNeverReachedProducesNoEventAtAll) {
  const auto curve = [](double t) { return 1.0e6 * std::exp(-1.0e-3 * t); };
  const EventSeries series = sampled(curve, logspace(1.0, 1.0e3, 30), true);

  EXPECT_TRUE(crossings(series, 1.0e12).empty()) << "never that high";
  EXPECT_TRUE(crossings(series, 1.0e-12).empty()) << "never that low";
  EXPECT_FALSE(firstCrossing(series, 1.0e12).has_value());
}

// A curve that touches the level and retreats without passing through it flips the
// "is below" flag twice, so the touch is reported once on the way in and once on the way out
// rather than being silently dropped or double-counted within one interval.
TEST(TrajectoryEvents, ASampleExactlyOnTheLevelIsNotBelowIt) {
  EventSeries series;
  series.times = {0.0, 1.0, 2.0};
  series.values = {1.0, 5.0, 1.0};

  const std::vector<TrajectoryEvent> found = crossings(series, 5.0);
  ASSERT_EQ(found.size(), 2u);
  EXPECT_EQ(found[0].kind, EventKind::Rising);
  EXPECT_DOUBLE_EQ(found[0].timeSeconds, 1.0) << "the crossing is the sample itself";
  EXPECT_EQ(found[1].kind, EventKind::Falling);
}

// --- extrema -----------------------------------------------------------------

// e^{-a t} - e^{-b t} peaks at ln(b/a)/(b-a). This is the ingrowth shape -- a daughter growing
// in while its parent decays -- so the peak is exactly the "when does waiting stop helping"
// instant the capability exists to locate.
TEST(TrajectoryEvents, IngrowthPeakMatchesTheAnalyticMaximum) {
  const double a = 1.0e-4;
  const double b = 1.0e-3;
  const double expected = std::log(b / a) / (b - a);

  const auto curve = [&](double t) { return std::exp(-a * t) - std::exp(-b * t); };
  const std::vector<double> times = logspace(10.0, 1.0e6, 60);

  const std::vector<TrajectoryEvent> refined = extrema(sampled(curve, times, true), {});
  ASSERT_EQ(refined.size(), 1u);
  EXPECT_EQ(refined[0].kind, EventKind::Maximum);
  EXPECT_NEAR(refined[0].timeSeconds, expected, expected * 1.0e-4);
  EXPECT_NEAR(refined[0].value, curve(expected), curve(expected) * 1.0e-6);
  EXPECT_TRUE(refined[0].refined);
  EXPECT_TRUE(refined[0].converged);

  // The interpolated path is coarser -- a parabola through three log-spaced samples -- but
  // still lands inside the bracket and within a few percent of the true turn.
  const std::vector<TrajectoryEvent> guessed = extrema(sampled(curve, times, false), {});
  ASSERT_EQ(guessed.size(), 1u);
  EXPECT_NEAR(guessed[0].timeSeconds, expected, expected * 0.05);
  EXPECT_FALSE(guessed[0].refined);
  EXPECT_GE(guessed[0].timeSeconds, guessed[0].bracketStartSeconds);
  EXPECT_LE(guessed[0].timeSeconds, guessed[0].bracketEndSeconds);
}

TEST(TrajectoryEvents, AMinimumIsFoundAndNamedAsOne) {
  const auto curve = [](double t) { return (t - 5.0) * (t - 5.0) + 1.0; };
  const std::vector<TrajectoryEvent> found = extrema(sampled(curve, linspace(0.0, 10.0, 21), true));

  ASSERT_EQ(found.size(), 1u);
  EXPECT_EQ(found[0].kind, EventKind::Minimum);
  EXPECT_NEAR(found[0].timeSeconds, 5.0, 1.0e-4);
  EXPECT_NEAR(found[0].value, 1.0, 1.0e-6);
}

// A turn needs three samples to be distinguished from a monotone run, so a peak in the first
// or last interval is outside what the sampling can resolve and is not reported.
// With no evaluator an extremum is a parabola through three samples, so an exactly parabolic
// curve on an evenly spaced grid has to come back exactly -- no tolerance to hide behind. This
// pins the abscissa choice as much as the formula: fitted in log-time these samples are uneven
// and the vertex would be off by a percent or so.
TEST(TrajectoryEvents, AnExactParabolaIsInterpolatedExactly) {
  const auto curve = [](double t) { return -3.0 * (t - 5.0) * (t - 5.0) + 17.0; };
  const std::vector<TrajectoryEvent> found =
      extrema(sampled(curve, linspace(1.0, 11.0, 11), false));

  ASSERT_EQ(found.size(), 1u);
  EXPECT_EQ(found[0].kind, EventKind::Maximum);
  EXPECT_FALSE(found[0].refined);
  EXPECT_NEAR(found[0].timeSeconds, 5.0, 1.0e-9);
  EXPECT_NEAR(found[0].value, 17.0, 1.0e-9)
      << "the reported peak is the parabola's, not a sample's";
}

TEST(TrajectoryEvents, AMonotoneCurveAndAFlatOneHaveNoExtrema) {
  const auto falling = [](double t) { return std::exp(-t); };
  EXPECT_TRUE(extrema(sampled(falling, linspace(0.0, 10.0, 21), true)).empty());

  EventSeries flat;
  flat.times = {0.0, 1.0, 2.0, 3.0};
  flat.values = {7.0, 7.0, 7.0, 7.0};
  EXPECT_TRUE(extrema(flat).empty()) << "a flat stretch is not turning";
}

TEST(TrajectoryEvents, ExtremaNeedThreeSamples) {
  EventSeries two;
  two.times = {0.0, 1.0};
  two.values = {1.0, 2.0};
  EXPECT_TRUE(extrema(two).empty());
}

// --- windows -----------------------------------------------------------------

TEST(TrajectoryEvents, WindowAboveHasBothEdgesRefinedAndFlaggedObserved) {
  const auto curve = [](double t) {
    return 1.0 + 99.0 * std::exp(-((t - 15.0) * (t - 15.0)) / 9.0);
  };
  const EventSeries series = sampled(curve, linspace(0.0, 40.0, 161), true);

  const std::vector<LevelWindow> windows = windowsAbove(series, 50.0);
  ASSERT_EQ(windows.size(), 1u);
  EXPECT_TRUE(windows[0].entryObserved);
  EXPECT_TRUE(windows[0].exitObserved);
  EXPECT_LT(windows[0].startSeconds, windows[0].endSeconds);

  // The edges are the crossings, at the same instants and by the same code path.
  const std::vector<TrajectoryEvent> found = crossings(series, 50.0);
  ASSERT_EQ(found.size(), 2u);
  EXPECT_DOUBLE_EQ(windows[0].startSeconds, found[0].timeSeconds);
  EXPECT_DOUBLE_EQ(windows[0].endSeconds, found[1].timeSeconds);
}

// A grid that begins or ends inside the window never observed that edge. Reporting the grid's
// own endpoint as a crossing would invent one, so the flags say which ends are bounds.
TEST(TrajectoryEvents, WindowsOpenAtTheGridEdgeAreFlaggedNotClipped) {
  const auto decay = [](double t) { return 1.0e6 * std::exp(-1.0e-3 * t); };
  const EventSeries series = sampled(decay, logspace(1.0, 1.0e5, 40), true);

  const std::vector<LevelWindow> above = windowsAbove(series, 1.0e3);
  ASSERT_EQ(above.size(), 1u);
  EXPECT_FALSE(above[0].entryObserved) << "it was already above when the grid started";
  EXPECT_TRUE(above[0].exitObserved);
  EXPECT_DOUBLE_EQ(above[0].startSeconds, series.times.front());

  const std::vector<LevelWindow> below = windowsBelow(series, 1.0e3);
  ASSERT_EQ(below.size(), 1u);
  EXPECT_TRUE(below[0].entryObserved);
  EXPECT_FALSE(below[0].exitObserved) << "it is still below when the grid ends";
  EXPECT_DOUBLE_EQ(below[0].endSeconds, series.times.back());
}

// Above and below partition the grid: they meet at the same instant and neither overlaps.
TEST(TrajectoryEvents, AboveAndBelowPartitionTheGrid) {
  const auto decay = [](double t) { return 1.0e6 * std::exp(-1.0e-3 * t); };
  const EventSeries series = sampled(decay, logspace(1.0, 1.0e5, 40), true);

  const std::vector<LevelWindow> above = windowsAbove(series, 1.0e3);
  const std::vector<LevelWindow> below = windowsBelow(series, 1.0e3);
  ASSERT_EQ(above.size(), 1u);
  ASSERT_EQ(below.size(), 1u);
  EXPECT_DOUBLE_EQ(above[0].endSeconds, below[0].startSeconds);
}

TEST(TrajectoryEvents, ACurveWhollyBelowTheLevelHasNoWindowAbove) {
  const auto decay = [](double t) { return std::exp(-t); };
  const EventSeries series = sampled(decay, linspace(0.0, 10.0, 21), true);
  EXPECT_TRUE(windowsAbove(series, 1.0e6).empty());

  const std::vector<LevelWindow> below = windowsBelow(series, 1.0e6);
  ASSERT_EQ(below.size(), 1u);
  EXPECT_FALSE(below[0].entryObserved);
  EXPECT_FALSE(below[0].exitObserved) << "open at both ends: never observed to enter or leave";
}

// --- the grid has to be usable ------------------------------------------------

TEST(TrajectoryEvents, AnUnusableGridIsRefusedRatherThanSearched) {
  EventSeries mismatched;
  mismatched.times = {0.0, 1.0, 2.0};
  mismatched.values = {1.0, 2.0};
  EXPECT_THROW(crossings(mismatched, 1.0), InputError);

  EventSeries single;
  single.times = {1.0};
  single.values = {1.0};
  EXPECT_THROW(crossings(single, 1.0), InputError);

  EventSeries repeated;
  repeated.times = {0.0, 1.0, 1.0};
  repeated.values = {1.0, 2.0, 3.0};
  EXPECT_THROW(crossings(repeated, 1.0), InputError)
      << "a repeated time is a zero-width bracket nothing can be located in";

  EventSeries backwards;
  backwards.times = {0.0, 2.0, 1.0};
  backwards.values = {1.0, 2.0, 3.0};
  EXPECT_THROW(extrema(backwards), InputError);
}

// --- series drawn from a response table ---------------------------------------

NuclearData oneEmitter(double lambda) {
  StoreArrays a;
  a.provenance.version = 1;
  a.nuclideKey = {Zai{50, 100, 0}.key(), Zai{51, 100, 0}.key()};
  a.halfLife = {synth::halfLifeFor(lambda), 0.0};
  a.modeOffset = {0, 1, 1};
  a.modeRtyp = {synth::kBetaMinus};
  a.modeBranching = {1.0};
  a.modeFinalState = {0};
  a.modeIsFission = {0};
  return NuclearData::fromArrays(std::move(a));
}

// The end-to-end shape: decay a real inventory, build a response table, and ask when its total
// falls below a limit. Activity is lambda*N0*e^{-lambda t}, so the crossing has a closed form.
TEST(TrajectoryEvents, TotalSeriesLocatesALimitCrossingOnARealResponse) {
  const double lambda = 1.0e-4;
  const double atoms = 1.0e20;
  const double level = 1.0e13;
  const double expected = std::log(lambda * atoms / level) / lambda;
  ASSERT_GT(expected, 0.0);

  const NuclearData data = oneEmitter(lambda);
  Inventory inv;
  inv.add(Zai{50, 100, 0}, atoms);

  const std::vector<double> times = logspace(1.0, 1.0e6, 60);
  const ResponseTable table = buildResponse(data, decay(data, inv, times), ResponseSpec{});

  const EventSeries series = totalSeries(table);
  ASSERT_EQ(series.times.size(), times.size());

  const std::optional<TrajectoryEvent> crossing = firstCrossing(series, level);
  ASSERT_TRUE(crossing.has_value());
  EXPECT_EQ(crossing->kind, EventKind::Falling);
  EXPECT_NEAR(crossing->timeSeconds, expected, expected * 0.01);
  EXPECT_FALSE(crossing->refined) << "a table carries samples, not an evaluator";
}

TEST(TrajectoryEvents, ContributorSeriesTracksOneColumnAndRefusesAnIndexItDoesNotHave) {
  const NuclearData data = oneEmitter(1.0e-4);
  Inventory inv;
  inv.add(Zai{50, 100, 0}, 1.0e20);

  const std::vector<double> times = logspace(1.0, 1.0e6, 20);
  const ResponseTable table = buildResponse(data, decay(data, inv, times), ResponseSpec{});

  ASSERT_GT(table.contributorCount(), 0);
  const EventSeries series = contributorSeries(table, 0);
  EXPECT_EQ(series.values.size(), times.size());

  EXPECT_THROW(contributorSeries(table, table.contributorCount()), InputError);
  EXPECT_THROW(contributorSeries(table, -1), InputError);
}

// A ratio clock does not start until the nuclide it divides by exists. The samples before that
// are dropped rather than reported as an infinite ratio that would sit in front of every real
// crossing as a spurious one.
TEST(TrajectoryEvents, RatioSeriesDropsTheSamplesBeforeItsDenominatorExists) {
  // Parent -> daughter -> stable, so the daughter is absent at t=0 and grows in.
  StoreArrays a = synth::linearChain({1.0e-3, 1.0e-4});
  const NuclearData data = NuclearData::fromArrays(std::move(a));

  Inventory inv;
  inv.add(Zai{50, 100, 0}, 1.0e20);

  std::vector<double> times = {0.0};
  for (const double t : logspace(1.0, 1.0e5, 30)) {
    times.push_back(t);
  }
  const ResponseTable table = buildResponse(data, decay(data, inv, times), ResponseSpec{});
  ASSERT_GE(table.contributorCount(), 2);

  // Column order follows the table, so find the two by label rather than assuming.
  int parent = -1;
  int daughter = -1;
  for (int c = 0; c < table.contributorCount(); ++c) {
    if (table.labels[static_cast<std::size_t>(c)] == "Sn-100") {
      parent = c;
    }
    if (table.labels[static_cast<std::size_t>(c)] == "Sb-100") {
      daughter = c;
    }
  }
  ASSERT_GE(parent, 0);
  ASSERT_GE(daughter, 0);

  const EventSeries series = ratioSeries(table, parent, daughter);
  EXPECT_LT(series.times.size(), times.size())
      << "the sample where the daughter had not grown in yet is not part of the clock";
  EXPECT_GT(series.times.front(), 0.0);
  for (const double value : series.values) {
    EXPECT_TRUE(std::isfinite(value));
  }
}

// A denominator that is present, then absent again, is a pole inside the grid. Interpolating
// across it would produce a bracket whose endpoints straddle infinity, so this refuses instead.
// Built by hand because no physical chain does this -- which is the point: the guard exists for
// a table assembled some other way, not for one decay() produced.
TEST(TrajectoryEvents, RatioSeriesRefusesADenominatorThatReturnsToZero) {
  ResponseTable table;
  table.times = {1.0, 2.0, 3.0, 4.0};
  table.contributors = {ContributorId{1, 0, 0.0}, ContributorId{2, 0, 0.0}};
  table.labels = {"A", "B"};
  table.flags = {kFlagNone, kFlagNone};
  // Row-major by time, two columns: the second is absent, present, absent, present.
  table.values = {1.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0};
  table.totals = {1.0, 2.0, 1.0, 2.0};

  EXPECT_THROW(ratioSeries(table, 0, 1), InputError);

  // With the denominator absent only at the start, the same table is usable -- the leading
  // sample is dropped and the clock starts where the denominator does.
  table.values = {1.0, 0.0, 1.0, 1.0, 1.0, 2.0, 1.0, 4.0};
  const EventSeries series = ratioSeries(table, 0, 1);
  ASSERT_EQ(series.times.size(), 3u);
  EXPECT_DOUBLE_EQ(series.times.front(), 2.0);
  EXPECT_DOUBLE_EQ(series.values.front(), 1.0);
  EXPECT_DOUBLE_EQ(series.values.back(), 0.25);
}

// Fewer than two usable samples is not enough grid to locate anything in, and saying so beats
// returning a series that every search would then reject for a less obvious reason.
TEST(TrajectoryEvents, RatioSeriesRefusesADenominatorPresentTooBriefly) {
  ResponseTable table;
  table.times = {1.0, 2.0, 3.0};
  table.contributors = {ContributorId{1, 0, 0.0}, ContributorId{2, 0, 0.0}};
  table.labels = {"A", "B"};
  table.flags = {kFlagNone, kFlagNone};
  table.values = {1.0, 0.0, 1.0, 0.0, 1.0, 1.0};
  table.totals = {1.0, 1.0, 2.0};

  EXPECT_THROW(ratioSeries(table, 0, 1), InputError);
}

// --- refining against a real response -----------------------------------------

// The point of an evaluator: the same grid, the same crossing, placed by solving instead of by
// interpolating. Activity from one nuclide is a pure exponential, where log-linear
// interpolation is exact, so this is not a test that refinement is more ACCURATE here -- it is
// a test that it narrows the bracket it reports and says it did.
TEST(TrajectoryEvents, AnEvaluatorNarrowsACrossingTheTableAloneOnlyInterpolates) {
  const double lambda = 1.0e-4;
  const double atoms = 1.0e20;
  const double level = 1.0e13;
  const double expected = std::log(lambda * atoms / level) / lambda;

  const NuclearData data = oneEmitter(lambda);
  Inventory inv;
  inv.add(Zai{50, 100, 0}, atoms);

  const std::vector<double> times = logspace(1.0, 1.0e6, 20);
  const ResponseTable table = buildResponse(data, decay(data, inv, times), ResponseSpec{});

  const ResponseEvaluator evaluator(data, inv, ResponseSpec{});
  const EventSeries refined = totalSeries(table, evaluator);
  const EventSeries sampled = totalSeries(table);
  ASSERT_EQ(refined.values, sampled.values) << "refinement changes how an event is placed, "
                                               "not what the grid saw";

  const std::optional<TrajectoryEvent> narrow = firstCrossing(refined, level);
  const std::optional<TrajectoryEvent> wide = firstCrossing(sampled, level);
  ASSERT_TRUE(narrow.has_value());
  ASSERT_TRUE(wide.has_value());

  EXPECT_TRUE(narrow->refined);
  EXPECT_FALSE(wide->refined);
  EXPECT_TRUE(narrow->converged);
  EXPECT_LT(narrow->locatedToSeconds, 1.0e-3 * wide->locatedToSeconds)
      << "a refined event is placed far more tightly than the interval that found it";
  EXPECT_NEAR(narrow->timeSeconds, expected, expected * 1.0e-5);
  EXPECT_GT(evaluator.solves(), 0) << "refinement is solves, and they are countable";
}

// A curve where interpolation is NOT exact, so refinement moves the answer rather than only
// tightening the bracket it is reported with. The daughter of a decay chain peaks where its
// parabola through three log-spaced samples does not.
TEST(TrajectoryEvents, RefinementMovesAPeakInterpolationPlacesWrongly) {
  StoreArrays arrays = synth::linearChain({1.0e-3, 1.0e-4});
  const NuclearData data = NuclearData::fromArrays(std::move(arrays));
  Inventory inv;
  inv.add(Zai{50, 100, 0}, 1.0e20);

  // The Bateman maximum of the daughter: where its two exponentials balance.
  const double l0 = 1.0e-3;
  const double l1 = 1.0e-4;
  const double expected = std::log(l0 / l1) / (l0 - l1);

  const std::vector<double> times = logspace(10.0, 1.0e6, 25);
  const ResponseTable table = buildResponse(data, decay(data, inv, times), ResponseSpec{});
  const int daughter = 1;
  ASSERT_GT(table.contributorCount(), daughter);

  const ResponseEvaluator evaluator(data, inv, ResponseSpec{});
  const std::vector<TrajectoryEvent> interpolated = extrema(contributorSeries(table, daughter));
  const std::vector<TrajectoryEvent> refined =
      extrema(contributorSeries(table, daughter, evaluator));
  ASSERT_EQ(interpolated.size(), 1u);
  ASSERT_EQ(refined.size(), 1u);
  EXPECT_EQ(refined[0].kind, EventKind::Maximum);

  EXPECT_NEAR(refined[0].timeSeconds, expected, expected * 1.0e-4);
  EXPECT_LT(std::abs(refined[0].timeSeconds - expected),
            std::abs(interpolated[0].timeSeconds - expected))
      << "solving inside the bracket beats fitting a parabola across it";
}

// An evaluator built for another curve would narrow an event using values that are not on it.
// Refused at the point the two are joined, rather than producing a plausible wrong instant.
TEST(TrajectoryEvents, AnEvaluatorForADifferentCurveIsRefused) {
  const NuclearData data = oneEmitter(1.0e-4);
  Inventory inv;
  inv.add(Zai{50, 100, 0}, 1.0e20);

  const std::vector<double> times = logspace(1.0, 1.0e6, 10);
  const ResponseTable table = buildResponse(data, decay(data, inv, times), ResponseSpec{});

  ResponseSpec other;
  other.unit = Unit::Curie;
  const ResponseEvaluator wrongUnit(data, inv, other);
  EXPECT_THROW(totalSeries(table, wrongUnit), InputError);

  ResponseSpec elsewhere;
  elsewhere.geometry.distanceM = 2.0;
  const ResponseEvaluator wrongDistance(data, inv, elsewhere);
  EXPECT_THROW(totalSeries(table, wrongDistance), InputError);
}

// --- the fixed-duration task --------------------------------------------------

// Decays accrued by a task of length D starting at t are N0 e^{-lambda t} (1 - e^{-lambda D}),
// so the task curve is itself an exponential and the start time at which a budget binds has a
// closed form. That is the whole inversion the interval integral makes possible.
TEST(TrajectoryEvents, TaskCurveMatchesTheAnalyticAccrualAndInvertsToAStartTime) {
  const double lambda = 1.0e-4;
  const double atoms = 1.0e20;
  const double duration = 3600.0;
  const double perTask = atoms * (1.0 - std::exp(-lambda * duration));

  const NuclearData data = oneEmitter(lambda);
  Inventory inv;
  inv.add(Zai{50, 100, 0}, atoms);

  ResponseSpec spec;
  spec.unit = Unit::Decays;

  const std::vector<double> starts = logspace(1.0, 1.0e6, 40);
  const EventSeries series = taskSeries(data, inv, spec, starts, duration);
  ASSERT_EQ(series.values.size(), starts.size());
  EXPECT_NEAR(series.values.front(), perTask * std::exp(-lambda * starts.front()), perTask * 1.0e-6)
      << "every sample is the exact integral over its own window";

  const double budget = 0.1 * perTask;
  const double expected = std::log(perTask / budget) / lambda;
  const std::optional<TrajectoryEvent> fits = firstCrossing(series, budget);
  ASSERT_TRUE(fits.has_value());
  EXPECT_EQ(fits->kind, EventKind::Falling) << "waiting makes a job cheaper, not dearer";
  EXPECT_NEAR(fits->timeSeconds, expected, expected * 0.01);

  // And the windows say the same thing as intervals: from that start onwards, the job fits.
  const std::vector<LevelWindow> allowed = windowsBelow(series, budget);
  ASSERT_EQ(allowed.size(), 1u);
  EXPECT_NEAR(allowed.front().startSeconds, expected, expected * 0.01);
  EXPECT_FALSE(allowed.front().exitObserved) << "it never stops fitting";
}

// The question the task curve exists for. Where a daughter grows in faster than its parent
// decays, the accrued total RISES first, so the worst time to do a fixed-length job is neither
// as soon as possible nor as late as possible -- it is a turn on a curve that no ranking and no
// forecast shows.
//
// Sr-90 -> Y-90 in miniature: a long-lived parent feeding a fast daughter, where the total
// activity climbs to roughly twice the parent's own while the daughter fills in. Its maximum is
// closed form -- 2*lambda0*e^{-lambda0 t} = lambda1*e^{-lambda1 t} -- and for a task short
// against the daughter's half-life the accrued curve turns within a fraction of a percent of
// the rate curve it integrates.
TEST(TrajectoryEvents, TaskCurveFindsTheWorstTimeToStartAnIngrowthFedJob) {
  const double l0 = 1.0e-9;
  const double l1 = 1.0e-5;
  StoreArrays arrays = synth::linearChain({l0, l1});
  const NuclearData data = NuclearData::fromArrays(std::move(arrays));
  Inventory inv;
  inv.add(Zai{50, 100, 0}, 1.0e20);

  ResponseSpec spec;
  spec.unit = Unit::Decays;

  const double expected = std::log(l1 / (2.0 * l0)) / (l1 - l0);
  const std::vector<double> starts = logspace(1.0e3, 1.0e7, 30);
  const EventSeries series =
      taskSeries(data, inv, spec, starts, 3600.0, DecayOptions{}, /*refine=*/true);

  const std::vector<TrajectoryEvent> turns = extrema(series);
  ASSERT_EQ(turns.size(), 1u);
  EXPECT_EQ(turns[0].kind, EventKind::Maximum);
  EXPECT_TRUE(turns[0].refined) << "the curve was given an evaluator, so the turn was solved for";
  EXPECT_NEAR(turns[0].timeSeconds, expected, expected * 0.01);
  EXPECT_LT(turns[0].locatedToSeconds, turns[0].bracketEndSeconds - turns[0].bracketStartSeconds);
}

TEST(TrajectoryEvents, ATaskRefusesADurationOrAGridItCannotBeAskedAbout) {
  const NuclearData data = oneEmitter(1.0e-4);
  Inventory inv;
  inv.add(Zai{50, 100, 0}, 1.0e20);

  ResponseSpec spec;
  spec.unit = Unit::Decays;
  const std::vector<double> starts = {1.0, 2.0, 3.0};

  EXPECT_THROW(taskSeries(data, inv, spec, starts, 0.0), InputError);
  EXPECT_THROW(taskSeries(data, inv, spec, std::vector<double>{1.0}, 60.0), InputError);
  EXPECT_THROW(taskSeries(data, inv, spec, std::vector<double>{2.0, 1.0}, 60.0), InputError);

  // A rate cannot describe what a window accrued, and the unit gating says so here as it does
  // everywhere else.
  ResponseSpec rate;
  rate.unit = Unit::Becquerel;
  EXPECT_THROW(taskSeries(data, inv, rate, starts, 60.0), InputError);
}

// --- the stay time: the other inversion --------------------------------------

// A single decaying nuclide with a stable terminator, so the total decay rate is the parent's
// alone and the accrual over [t0, t0 + D] is
//
//     A(D) = N0 (e^{-lambda t0} - e^{-lambda (t0 + D)})
//
// which inverts in closed form. Sampling the task curve would never produce this number: it is
// a root in the LENGTH of the window, and every point of that curve has a fixed length.
TEST(StayTime, MatchesTheAnalyticInversionOfTheAccrual) {
  const double lambda = 1.0e-4;
  const double atoms = 1.0e20;
  const double start = 5.0e3;

  const NuclearData data = oneEmitter(lambda);
  Inventory inv;
  inv.add(Zai{50, 100, 0}, atoms);

  ResponseSpec spec;
  spec.unit = Unit::Decays;

  const double remaining = atoms * std::exp(-lambda * start);
  const double budget = 0.25 * remaining;
  const double expected = -std::log(1.0 - budget / remaining) / lambda;

  const StayTime stay = stayTime(data, inv, spec, start, budget, /*maxDurationSeconds=*/1.0e6);
  ASSERT_TRUE(stay.bounded);
  EXPECT_TRUE(stay.converged);
  EXPECT_NEAR(stay.durationSeconds, expected, expected * 1.0e-5);
  // The width the root was narrowed to is the honest error bar, and it has to be tighter than
  // the answer it qualifies or it says nothing.
  EXPECT_GT(stay.locatedToSeconds, 0.0);
  EXPECT_LT(stay.locatedToSeconds, expected * 1.0e-4);
  EXPECT_GT(stay.samples, 1);
}

// The cross-check between the two inversions. If the stay time says D, then a task of length D
// starting at the same instant must accrue exactly the budget -- the task curve and the stay
// time are one integral read two ways, and a disagreement would mean one of them is wrong about
// a quantity the other reports.
TEST(StayTime, AgreesWithTheTaskCurveAtTheDurationItFinds) {
  const double lambda = 3.0e-5;
  const double atoms = 1.0e20;
  const double start = 1.0e4;

  const NuclearData data = oneEmitter(lambda);
  Inventory inv;
  inv.add(Zai{50, 100, 0}, atoms);

  ResponseSpec spec;
  spec.unit = Unit::Decays;

  const double budget = 1.0e18;
  const StayTime stay = stayTime(data, inv, spec, start, budget, 1.0e6);
  ASSERT_TRUE(stay.bounded);

  const std::vector<double> starts = {start, start + 1.0};
  const EventSeries task = taskSeries(data, inv, spec, starts, stay.durationSeconds);
  EXPECT_NEAR(task.values.front(), budget, budget * 1.0e-5);
}

// The case that makes the flag worth having. A budget larger than everything the inventory has
// left to give is never spent, and reporting that as a very large duration -- or worse, as a
// zero -- would invert the answer.
TEST(StayTime, SaysWhenTheBudgetIsNeverSpentRatherThanInventingADuration) {
  const double lambda = 1.0e-4;
  const double atoms = 1.0e20;

  const NuclearData data = oneEmitter(lambda);
  Inventory inv;
  inv.add(Zai{50, 100, 0}, atoms);

  ResponseSpec spec;
  spec.unit = Unit::Decays;

  // Twice every decay the inventory can ever produce, so no ceiling would reach it.
  const StayTime stay = stayTime(data, inv, spec, 0.0, 2.0 * atoms, 1.0e6);
  EXPECT_FALSE(stay.bounded);
  EXPECT_DOUBLE_EQ(stay.durationSeconds, 0.0) << "unset, not 'leave immediately'";
  EXPECT_GT(stay.accruedAtMax, 0.0);
  EXPECT_LT(stay.accruedAtMax, stay.budget);
  EXPECT_EQ(stay.samples, 1) << "one integral establishes there is no bracket to search";

  // The same budget over a ceiling short enough to matter is the ordinary occupancy case: not
  // spent in an hour, and the report has a number to say how much of it an hour costs.
  const StayTime hour = stayTime(data, inv, spec, 0.0, 2.0 * atoms, 3600.0);
  EXPECT_FALSE(hour.bounded);
  EXPECT_LT(hour.accruedAtMax, stay.accruedAtMax);
}

// A decaying source is cheaper to stand next to later, so the same budget has to buy a longer
// stay from a later start. Monotone in the direction physics requires, which is the property
// that makes the root unique in the first place.
TEST(StayTime, BuysLongerTheLongerYouWait) {
  const double lambda = 1.0e-4;
  const double atoms = 1.0e20;

  const NuclearData data = oneEmitter(lambda);
  Inventory inv;
  inv.add(Zai{50, 100, 0}, atoms);

  ResponseSpec spec;
  spec.unit = Unit::Decays;
  const double budget = 1.0e17;

  const StayTime early = stayTime(data, inv, spec, 0.0, budget, 1.0e6);
  const StayTime late = stayTime(data, inv, spec, 2.0e4, budget, 1.0e6);
  ASSERT_TRUE(early.bounded);
  ASSERT_TRUE(late.bounded);
  EXPECT_GT(late.durationSeconds, early.durationSeconds);

  // And a larger budget buys longer from the same start, which is the other monotonicity the
  // uniqueness of the root rests on.
  const StayTime richer = stayTime(data, inv, spec, 0.0, 2.0 * budget, 1.0e6);
  ASSERT_TRUE(richer.bounded);
  EXPECT_GT(richer.durationSeconds, early.durationSeconds);
}

TEST(StayTime, RefusesAQuestionItCannotAnswer) {
  const NuclearData data = oneEmitter(1.0e-4);
  Inventory inv;
  inv.add(Zai{50, 100, 0}, 1.0e20);

  ResponseSpec spec;
  spec.unit = Unit::Decays;

  EXPECT_THROW(stayTime(data, inv, spec, -1.0, 1.0e17, 1.0e6), InputError);
  EXPECT_THROW(stayTime(data, inv, spec, 0.0, 0.0, 1.0e6), InputError);
  EXPECT_THROW(stayTime(data, inv, spec, 0.0, -1.0, 1.0e6), InputError);
  EXPECT_THROW(stayTime(data, inv, spec, 0.0, 1.0e17, 0.0), InputError);

  // A budget is an accrued total, so a rate unit is the wrong dimension for it -- the same
  // refusal a task curve makes, from the same place.
  ResponseSpec rate;
  rate.unit = Unit::Becquerel;
  EXPECT_THROW(stayTime(data, inv, rate, 0.0, 1.0e17, 1.0e6), InputError);
}

}  // namespace
}  // namespace nusift
