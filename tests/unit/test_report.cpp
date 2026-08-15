#include <gtest/gtest.h>

#include <cstdlib>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "nusift/core/nuclide.hpp"
#include "nusift/io/report.hpp"
#include "nusift/triage/ranking.hpp"

namespace nusift {
namespace {

Ranking oneRow(const std::string& label) {
  Ranking ranking;
  ranking.total = 1.0;
  ranking.coveredFraction = 1.0;

  Contributor c;
  c.id = ContributorId{Zai{55, 137, 0}.key(), 0};
  c.label = label;
  c.value = 1.0;
  c.fraction = 1.0;
  c.cumulativeFraction = 1.0;
  c.rank = 1;
  ranking.contributors.push_back(c);
  return ranking;
}

// An exposure ranking whose photon model is slightly incomplete, but not incompletely enough
// for any single contributor to cross the 5% flag threshold. That combination -- a magnitude
// to report and no names to report it against -- is the one the footer got wrong.
Ranking exposureWithATraceUnmodelled() {
  Ranking ranking = oneRow("Cs-137");
  ranking.metric = Metric::Exposure;
  ranking.unit = Unit::RoentgenPerHour;
  ranking.unmodeledEnergyFraction = 7.9e-8;
  return ranking;
}

std::string asJson(const Ranking& ranking, const ReportContext& context) {
  std::ostringstream out;
  writeRanking(out, ranking, context, ReportFormat::Json);
  return out.str();
}

std::string asText(const Ranking& ranking, const ReportContext& context) {
  std::ostringstream out;
  writeRanking(out, ranking, context, ReportFormat::Text);
  return out.str();
}

std::string asCsv(const Ranking& ranking) {
  std::ostringstream out;
  writeRanking(out, ranking, ReportContext{}, ReportFormat::Csv);
  return out.str();
}

// The text of one JSON value, found by its key. Just enough of a parser to test a writer
// with -- the point is to read the digits back, not to validate the document.
std::string jsonValue(const std::string& text, const std::string& key) {
  const std::string needle = "\"" + key + "\": ";
  const std::size_t at = text.find(needle);
  if (at == std::string::npos) {
    return {};
  }
  const std::size_t start = at + needle.size();
  return text.substr(start, text.find_first_of(",}\n", start) - start);
}

double jsonNumberValue(const std::string& text, const std::string& key) {
  return std::strtod(jsonValue(text, key).c_str(), nullptr);
}

// Line `index` of a CSV split on commas. Deliberately naive: these tests write labels with
// commas in them, and a splitter that understood quoting would hide the very column shift the
// quoting exists to prevent.
std::vector<std::string> csvFields(const std::string& text, std::size_t index) {
  std::istringstream lines(text);
  std::string line;
  for (std::size_t k = 0; k <= index; ++k) {
    if (!std::getline(lines, line)) {
      return {};
    }
  }
  std::vector<std::string> fields;
  std::istringstream row(line);
  std::string field;
  while (std::getline(row, field, ',')) {
    fields.push_back(field);
  }
  return fields;
}

// Provenance is a path or a command line the user supplied, and a label can be anything the
// store spells a nuclide as -- so either can carry a quote, a newline, or a tab. Unescaped,
// each of them closes the JSON string early and the whole report becomes unparseable, in a
// consumer somewhere downstream rather than here where it was produced.
TEST(Report, JsonEscapesEveryCharacterThatWouldBreakAString) {
  ReportContext context;
  context.seedProvenance = "line one\nline two\ttabbed";
  context.storeLibrary = "back\\slash";
  const std::string text = asJson(oneRow("a \"quoted\" label"), context);

  // Named rather than written inline: MSVC's default preprocessor mangles a raw string that
  // contains backslashes when it appears as a macro argument, and every expectation here is
  // about backslashes.
  const std::string escapedSeed = "\"seed\": \"line one\\nline two\\ttabbed\"";
  const std::string escapedLibrary = "\"library\": \"back\\\\slash\"";
  const std::string escapedLabel = "\"label\": \"a \\\"quoted\\\" label\"";

  EXPECT_NE(text.find(escapedSeed), std::string::npos) << text;
  EXPECT_NE(text.find(escapedLibrary), std::string::npos) << text;
  EXPECT_NE(text.find(escapedLabel), std::string::npos) << text;
}

// Control characters with no short escape still have to leave as \u00XX. A vertical tab is
// not something anyone types, but it is something a file path copied out of a terminal can
// carry, and one of them is enough to make the document invalid.
TEST(Report, JsonEscapesControlCharactersWithNoShorterForm) {
  ReportContext context;
  context.seedProvenance = std::string("bell\x07vtab\x0b");
  const std::string text = asJson(oneRow("Cs-137"), context);

  const std::string escaped = "\"seed\": \"bell\\u0007vtab\\u000b\"";
  EXPECT_NE(text.find(escaped), std::string::npos) << text;

  // Nothing below 0x20 may survive into the document except the newlines the writer lays out
  // with -- pad is spaces, so every other control character in the output would be one that
  // escaped from a string.
  for (const char c : text) {
    const auto byte = static_cast<unsigned char>(c);
    EXPECT_TRUE(byte >= 0x20 || c == '\n') << "raw control character 0x" << std::hex
                                           << static_cast<unsigned>(byte) << " in the output";
  }
}

// Multi-byte UTF-8 is already valid inside a JSON string. Escaping it byte by byte would turn
// a correct name into mojibake, so it passes through untouched.
TEST(Report, JsonLeavesUtf8Alone) {
  ReportContext context;
  // Ends with U+00B0 DEGREE SIGN, written as the two UTF-8 bytes it is made of.
  context.geometry = "point source at 1 m, air at 20\xC2\xB0";
  const std::string text = asJson(oneRow("Cs-137"), context);
  EXPECT_NE(text.find("20\xC2\xB0"), std::string::npos) << text;
}

// --- the text footer ----------------------------------------------------------

// The energy-fraction paragraph used to be terminated only by the named list that sometimes
// follows it. With a fraction above zero and nothing crossing the per-nuclide flag threshold
// -- which is the ordinary case on a real store -- the report ended mid-line.
TEST(Report, TextFooterTerminatesItsLastLine) {
  const std::string text = asText(exposureWithATraceUnmodelled(), ReportContext{});
  ASSERT_FALSE(text.empty());
  EXPECT_EQ(text.back(), '\n') << text;
}

// And the consequence of not terminating it: writeRankings separates rankings with a newline,
// which an unterminated footer consumes finishing its own line. The blank line between two
// --at times disappears, and the two reports run together.
TEST(Report, RankingsStaySeparatedWhenTheFooterEndsOnTheEnergyFraction) {
  const Ranking ranking = exposureWithATraceUnmodelled();
  std::ostringstream out;
  writeRankings(out, {ranking, ranking}, ReportContext{}, ReportFormat::Text);

  const std::string text = out.str();
  EXPECT_NE(text.find("that much\n\nNuSIFT"), std::string::npos) << text;
}

// Nothing missing means no paragraph at all, rather than one reporting zero.
TEST(Report, TextFooterSaysNothingWhenNothingIsUnmodelled) {
  const std::string text = asText(oneRow("Cs-137"), ReportContext{});
  EXPECT_EQ(text.find("does not model"), std::string::npos) << text;
}

// Each ranking is footnoted with its OWN flagged emitters. Integrating several intervals gives
// each one a different set, and sharing one context named the last interval's nuclides under
// every ranking -- including ones where they do not appear.
TEST(Report, EachRankingCarriesItsOwnContext) {
  std::vector<ReportContext> contexts(2);
  contexts[0].unmodeledContinuum = {"Y-90"};
  contexts[1].unmodeledContinuum = {"Rb-90"};

  const Ranking ranking = exposureWithATraceUnmodelled();
  std::ostringstream out;
  writeRankings(out, {ranking, ranking}, contexts, ReportFormat::Text);

  const std::string text = out.str();
  const std::size_t first = text.find("Y-90");
  const std::size_t second = text.find("Rb-90");
  ASSERT_NE(first, std::string::npos) << text;
  ASSERT_NE(second, std::string::npos) << text;
  EXPECT_LT(first, second) << text;
  EXPECT_EQ(text.find("Y-90", first + 1), std::string::npos) << text;
}

// --- pinned rows ---------------------------------------------------------------

// A ranking cut at one row, with a contributor pinned from well below it.
Ranking withAPinnedTail(int trueRank, double value) {
  Ranking ranking = oneRow("Ba-140");
  ranking.total = 100.0;
  ranking.contributors[0].id = ContributorId{Zai{56, 140, 0}.key(), 0};
  ranking.contributors[0].value = 80.0;
  ranking.contributors[0].fraction = 0.80;
  ranking.contributors[0].cumulativeFraction = 0.80;
  ranking.coveredFraction = 0.80 + value / 100.0;
  ranking.omittedCount = 12;

  Contributor pinned;
  pinned.id = ContributorId{Zai{55, 137, 0}.key(), 0};
  pinned.label = "Cs-137";
  pinned.value = value;
  pinned.fraction = value / 100.0;
  pinned.cumulativeFraction = trueRank > 0 ? 0.998 : 0.0;
  pinned.rank = trueRank;
  pinned.pinned = true;
  ranking.contributors.push_back(pinned);
  return ranking;
}

// Run together with the prefix above them, pinned rows read as one list with numbers missing
// out of it. The heading is what says the rows below it were asked for rather than reached.
TEST(Report, TextSeparatesPinnedRowsFromTheRankingAboveThem) {
  const std::string text = asText(withAPinnedTail(37, 0.2), ReportContext{});

  const std::size_t heading = text.find("  pinned:\n");
  const std::size_t leader = text.find("Ba-140");
  const std::size_t pinned = text.find("Cs-137");
  ASSERT_NE(heading, std::string::npos) << text;
  EXPECT_LT(leader, heading) << text;
  EXPECT_LT(heading, pinned) << text;
  // The place it actually holds, which is the whole reason the row is worth printing.
  EXPECT_NE(text.find("  37  Cs-137"), std::string::npos) << text;
}

// A pinned contributor that contributes nothing holds no place in the ordering. A 0 in the rank
// column would look like one, and a cumulative is meaningless where there is nothing above.
TEST(Report, TextPrintsNoRankRatherThanZeroForAContributorWithNoPlace) {
  const std::string text = asText(withAPinnedTail(0, 0.0), ReportContext{});

  EXPECT_NE(text.find("   -  Cs-137"), std::string::npos) << text;
  EXPECT_NE(text.find("contributes nothing to this activity at this time"), std::string::npos)
      << text;
}

TEST(Report, TextSaysNothingAboutPinsWhenNoneWereGiven) {
  const std::string text = asText(oneRow("Cs-137"), ReportContext{});
  EXPECT_EQ(text.find("pinned"), std::string::npos) << text;
}

// Without this a loaded table cannot tell a row that placed from one fetched from below the
// cut, which is the difference between a top-N and a top-N plus an aside.
TEST(Report, MachineReadableFormatsMarkWhichRowsWerePinned) {
  const Ranking ranking = withAPinnedTail(37, 0.2);

  const std::vector<std::string> ranked = csvFields(asCsv(ranking), 1);
  const std::vector<std::string> pinned = csvFields(asCsv(ranking), 2);
  ASSERT_EQ(ranked.size(), 11u);
  ASSERT_EQ(pinned.size(), 11u);
  EXPECT_EQ(ranked.back(), "0");
  EXPECT_EQ(pinned.back(), "1");

  const std::string json = asJson(ranking, ReportContext{});
  EXPECT_NE(json.find("\"rank\": 1, \"label\": \"Ba-140\""), std::string::npos) << json;
  EXPECT_NE(json.find("\"rank\": 37, \"label\": \"Cs-137\""), std::string::npos) << json;
  EXPECT_NE(json.find("\"pinned\": true"), std::string::npos) << json;
  EXPECT_NE(json.find("\"pinned\": false"), std::string::npos) << json;
}

// --- machine-readable precision -----------------------------------------------

// CSV and JSON exist to be parsed again. At the stream default of six significant digits
// `--at 1.23456789y` comes back as a different time than the one the report describes, and a
// cumulative fraction of 0.999990 is indistinguishable from 0.99999.
TEST(Report, JsonNumbersReadBackAsTheValuesTheyCameFrom) {
  Ranking ranking = oneRow("Cs-137");
  ranking.domain = Domain::Interval;
  ranking.time = 38955600.123456789;
  ranking.timeEnd = 1.0 / 3.0;
  ranking.total = 1049.6234567890123;
  ranking.coveredFraction = 0.99999000000000005;
  ranking.contributors[0].value = 6.02214076e23 / 7.0;
  ranking.contributors[0].cumulativeFraction = ranking.coveredFraction;

  const std::string text = asJson(ranking, ReportContext{});
  EXPECT_EQ(jsonNumberValue(text, "time_s"), ranking.time) << text;
  EXPECT_EQ(jsonNumberValue(text, "time_end_s"), ranking.timeEnd) << text;
  EXPECT_EQ(jsonNumberValue(text, "total"), ranking.total) << text;
  EXPECT_EQ(jsonNumberValue(text, "covered_fraction"), ranking.coveredFraction) << text;
  EXPECT_EQ(jsonNumberValue(text, "value"), ranking.contributors[0].value) << text;
}

// A value that needs no extra digits must not grow any. Round-tripping at a fixed 17 digits
// would render 0.99999 as 0.99999000000000005, which is correct and unreadable.
TEST(Report, JsonPrintsNoMoreDigitsThanTheValueHas) {
  Ranking ranking = oneRow("Cs-137");
  ranking.coveredFraction = 0.5;
  ranking.total = 1049.5;

  const std::string text = asJson(ranking, ReportContext{});
  EXPECT_EQ(jsonValue(text, "covered_fraction"), "0.5") << text;
  EXPECT_EQ(jsonValue(text, "total"), "1049.5") << text;
}

// JSON has no infinity and no NaN. A bare `inf` makes the whole document unparseable, and it
// fails in whatever consumes the report rather than here where it was written.
TEST(Report, JsonWritesNullRatherThanAnUnparseableInfinity) {
  Ranking ranking = oneRow("Cs-137");
  ranking.total = std::numeric_limits<double>::infinity();
  ranking.contributors[0].value = std::numeric_limits<double>::quiet_NaN();

  const std::string text = asJson(ranking, ReportContext{});
  EXPECT_EQ(jsonValue(text, "total"), "null") << text;
  EXPECT_EQ(jsonValue(text, "value"), "null") << text;
  EXPECT_EQ(text.find("inf"), std::string::npos) << text;
  EXPECT_EQ(text.find("nan"), std::string::npos) << text;
}

TEST(Report, CsvNumbersReadBackAsTheValuesTheyCameFrom) {
  Ranking ranking = oneRow("Cs-137");
  ranking.time = 38955600.123456789;
  ranking.contributors[0].value = 1049.6234567890123;
  ranking.contributors[0].cumulativeFraction = 0.99999000000000005;

  // Columns of the single data row: time_s, time_end_s, rank, contributor, key, value, unit,
  // fraction, cumulative_fraction, flags, pinned. The assertion is that the digits parse back to
  // the same double, not that they are spelled the way the literal above was -- the shortest
  // form of 38955600.123456789 is 38955600.12345679, and both name the same value.
  const std::vector<std::string> row = csvFields(asCsv(ranking), 1);
  ASSERT_EQ(row.size(), 11u);
  EXPECT_EQ(std::strtod(row[0].c_str(), nullptr), ranking.time);
  EXPECT_EQ(std::strtod(row[5].c_str(), nullptr), ranking.contributors[0].value);
  EXPECT_EQ(std::strtod(row[8].c_str(), nullptr), ranking.contributors[0].cumulativeFraction);
}

// --- CSV quoting ---------------------------------------------------------------

// No label carries a comma today. That is a property of the labels, not of the format: one
// that did would shift every column right of it by one, silently, in a file nobody re-reads
// by eye.
TEST(Report, CsvQuotesALabelThatWouldOtherwiseShiftTheColumns) {
  const std::string text = asCsv(oneRow("A=140 (La-140, Ba-140)"));
  EXPECT_NE(text.find("\"A=140 (La-140, Ba-140)\""), std::string::npos) << text;
}

TEST(Report, CsvDoublesAQuoteInsideALabel) {
  const std::string text = asCsv(oneRow("odd \"name\", quoted"));
  EXPECT_NE(text.find("\"odd \"\"name\"\", quoted\""), std::string::npos) << text;
}

// A label with nothing to escape stays bare, so the ordinary file is unchanged.
TEST(Report, CsvLeavesAnOrdinaryLabelUnquoted) {
  const std::string text = asCsv(oneRow("Cs-137"));
  EXPECT_NE(text.find(",Cs-137,"), std::string::npos) << text;
}

}  // namespace
}  // namespace nusift
