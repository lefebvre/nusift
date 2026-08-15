#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <string>

#include "nusift/core/error.hpp"
#include "nusift/io/time_spec.hpp"
#include "nusift/units.hpp"

namespace nusift {
namespace {

TEST(TimeSpec, ParsesEveryUnitSuffix) {
  EXPECT_DOUBLE_EQ(parseDuration("30s"), 30.0);
  EXPECT_DOUBLE_EQ(parseDuration("5m"), 300.0);
  EXPECT_DOUBLE_EQ(parseDuration("2h"), 7200.0);
  EXPECT_DOUBLE_EQ(parseDuration("30d"), 30.0 * 86400.0);
  EXPECT_DOUBLE_EQ(parseDuration("1y"), units::kSecondsPerYear);
}

// A bare number is seconds. Anything else would silently reinterpret machine-generated
// input that was already in SI.
TEST(TimeSpec, BareNumberIsSeconds) {
  EXPECT_DOUBLE_EQ(parseDuration("3600"), 3600.0);
  EXPECT_DOUBLE_EQ(parseDuration("1.5e3"), 1500.0);
}

TEST(TimeSpec, AcceptsFractionsAndSurroundingSpace) {
  EXPECT_DOUBLE_EQ(parseDuration("1.5y"), 1.5 * units::kSecondsPerYear);
  EXPECT_DOUBLE_EQ(parseDuration("  30d  "), 30.0 * 86400.0);
  EXPECT_DOUBLE_EQ(parseDuration("0.5h"), 1800.0);
}

// The Julian year, stated in the header and in --help. Over a long decay the choice is
// worth more than a year of elapsed time, so it is pinned rather than left to drift.
TEST(TimeSpec, AYearIs365Point25Days) {
  EXPECT_NEAR(parseDuration("1y") / parseDuration("1d"), 365.25, 1e-12);
}

TEST(TimeSpec, CaseInsensitiveSuffixes) {
  EXPECT_DOUBLE_EQ(parseDuration("30D"), parseDuration("30d"));
  EXPECT_DOUBLE_EQ(parseDuration("1Y"), parseDuration("1y"));
}

TEST(TimeSpec, RejectsGarbageNamingTheToken) {
  EXPECT_THROW(parseDuration(""), InputError);
  EXPECT_THROW(parseDuration("d"), InputError);
  EXPECT_THROW(parseDuration("30x"), InputError);
  EXPECT_THROW(parseDuration("abc"), InputError);
  EXPECT_THROW(parseDuration("-5d"), InputError);

  try {
    parseDuration("30x");
    FAIL() << "expected InputError";
  } catch (const InputError& e) {
    const std::string what = e.what();
    EXPECT_NE(what.find("30x"), std::string::npos) << what;
  }
}

// A bare "inf" is caught by the unit check, since 'f' is not a unit -- but "infs" is infinity
// seconds, and from_chars reads it as a number like any other.
TEST(TimeSpec, RejectsInfinityAndNaNWrittenAsNumbers) {
  EXPECT_THROW(parseDuration("infs"), InputError);
  EXPECT_THROW(parseDuration("nans"), InputError);
  EXPECT_THROW(parseDuration("-infs"), InputError);

  try {
    parseDuration("infs");
    FAIL() << "expected InputError";
  } catch (const InputError& e) {
    const std::string what = e.what();
    EXPECT_NE(what.find("finite"), std::string::npos) << what;
  }
}

// Both endpoints must land exactly where asked. A reported time 1e-16 off the requested one
// is confusing in output, and can reorder against an interval endpoint meant to coincide.
TEST(TimeGrid, LogGridHitsBothEndpointsExactly) {
  const std::vector<double> times = logspace(3600.0, 3.15576e9, 40);
  ASSERT_EQ(times.size(), 40u);
  EXPECT_DOUBLE_EQ(times.front(), 3600.0);
  EXPECT_DOUBLE_EQ(times.back(), 3.15576e9);
}

TEST(TimeGrid, LogGridHasAConstantRatio) {
  const std::vector<double> times = logspace(1.0, 1000.0, 4);
  ASSERT_EQ(times.size(), 4u);
  for (std::size_t i = 1; i + 1 < times.size(); ++i) {
    EXPECT_NEAR(times[i + 1] / times[i], times[i] / times[i - 1], 1e-9);
  }
  EXPECT_NEAR(times[1], 10.0, 1e-9);
  EXPECT_NEAR(times[2], 100.0, 1e-9);
}

TEST(TimeGrid, LinearGridIsEvenlySpaced) {
  const std::vector<double> times = linspace(0.0, 100.0, 5);
  ASSERT_EQ(times.size(), 5u);
  EXPECT_DOUBLE_EQ(times.front(), 0.0);
  EXPECT_DOUBLE_EQ(times.back(), 100.0);
  EXPECT_DOUBLE_EQ(times[2], 50.0);
}

// Log spacing from zero is undefined. Saying so beats silently substituting a small number,
// which would put the first reported time somewhere the user never asked about.
TEST(TimeGrid, LogGridFromZeroIsRejectedWithAnExplanation) {
  try {
    logspace(0.0, 100.0, 10);
    FAIL() << "expected InputError";
  } catch (const InputError& e) {
    const std::string what = e.what();
    EXPECT_NE(what.find("positive"), std::string::npos) << what;
  }
}

// An infinite endpoint is not a time, and a grid built between one and a finite endpoint is
// infinite at every point rather than at one -- each of which reaches the solver as a time to
// evaluate at. Refused where the grid is built, since that is the last place the endpoint is
// still identifiable as the thing that was wrong.
TEST(TimeGrid, RejectsNonFiniteEndpoints) {
  const double inf = std::numeric_limits<double>::infinity();
  EXPECT_THROW(logspace(1.0, inf, 2), InputError);
  EXPECT_THROW(logspace(inf, 100.0, 2), InputError);
  EXPECT_THROW(logspace(1.0, std::nan(""), 10), InputError);
  EXPECT_THROW(linspace(0.0, inf, 5), InputError);
  EXPECT_THROW(linspace(std::nan(""), 100.0, 5), InputError);
  // A single-point grid takes the same endpoints and so the same rules.
  EXPECT_THROW(logspace(inf, inf, 1), InputError);

  // Log of a non-positive stop is no more defined than log of a non-positive start.
  EXPECT_THROW(logspace(1.0, -100.0, 10), InputError);
  EXPECT_THROW(logspace(1.0, 0.0, 10), InputError);
}

TEST(TimeGrid, ParsesTheGridSyntax) {
  const std::vector<double> times = parseTimeGrid("1h:100y:log:60");
  ASSERT_EQ(times.size(), 60u);
  EXPECT_DOUBLE_EQ(times.front(), 3600.0);
  EXPECT_DOUBLE_EQ(times.back(), 100.0 * units::kSecondsPerYear);

  const std::vector<double> linear = parseTimeGrid("0:1d:lin:25");
  ASSERT_EQ(linear.size(), 25u);
  EXPECT_DOUBLE_EQ(linear.front(), 0.0);
  EXPECT_DOUBLE_EQ(linear.back(), 86400.0);
}

TEST(TimeGrid, RejectsMalformedGrids) {
  EXPECT_THROW(parseTimeGrid("1h:100y:log"), InputError);       // too few fields
  EXPECT_THROW(parseTimeGrid("1h:100y:cubic:10"), InputError);  // unknown spacing
  EXPECT_THROW(parseTimeGrid("100y:1h:log:10"), InputError);    // stop before start
  EXPECT_THROW(parseTimeGrid("1h:100y:log:1"), InputError);     // needs at least 2
  EXPECT_THROW(parseTimeGrid("1h:100y:log:x"), InputError);     // count is not a number
}

// The engine rejects duplicate times outright, and two points differing in the last bit
// cost a full factorization while carrying no information.
TEST(MergeTimes, SortsAndCollapsesNearDuplicates) {
  const std::vector<double> merged = mergeTimes({100.0, 5.0, 100.0, 50.0, 100.0 + 1e-15});
  ASSERT_EQ(merged.size(), 3u);
  EXPECT_DOUBLE_EQ(merged[0], 5.0);
  EXPECT_DOUBLE_EQ(merged[1], 50.0);
  EXPECT_DOUBLE_EQ(merged[2], 100.0);
}

TEST(MergeTimes, DropsNegativeTimes) {
  const std::vector<double> merged = mergeTimes({-5.0, 10.0, 0.0});
  ASSERT_EQ(merged.size(), 2u);
  EXPECT_DOUBLE_EQ(merged[0], 0.0);
  EXPECT_DOUBLE_EQ(merged[1], 10.0);
}

// One format at every magnitude, so a column of times lines up rather than mixing "30 d"
// with "1.50 y".
TEST(FormatDuration, ChoosesTheLargestReadableUnit) {
  EXPECT_EQ(formatDuration(30.0), "30 s");
  EXPECT_EQ(formatDuration(300.0), "5 m");
  EXPECT_EQ(formatDuration(7200.0), "2 h");
  EXPECT_EQ(formatDuration(30.0 * 86400.0), "30 d");
  EXPECT_EQ(formatDuration(units::kSecondsPerYear), "1 y");
  EXPECT_EQ(formatDuration(1.5 * units::kSecondsPerYear), "1.5 y");
}

TEST(FormatDuration, RoundTripsThroughTheParser) {
  for (const double seconds : {30.0, 3600.0, 86400.0, units::kSecondsPerYear}) {
    const std::string text = formatDuration(seconds);
    // Strip the space the formatter inserts for readability; the parser tolerates the rest.
    std::string compact;
    for (const char c : text) {
      if (c != ' ') {
        compact += c;
      }
    }
    EXPECT_NEAR(parseDuration(compact), seconds, seconds * 1e-3) << text;
  }
}

}  // namespace
}  // namespace nusift
