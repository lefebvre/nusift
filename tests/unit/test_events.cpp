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

}  // namespace
}  // namespace nusift
