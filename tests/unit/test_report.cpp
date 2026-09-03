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

// Two seed shares whose importances differ far below what %.4e can hold -- the Xe-140/Cs-140
// case from a real fission seed, where the two serialize identically at five digits.
SeedAttribution twoNearlyIdenticalShares() {
  SeedAttribution a;
  a.time = 2592000.0;
  a.unit = Unit::Becquerel;
  a.total = 4.0;
  a.coveredFraction = 1.0;

  SeedShare first;
  first.key = Zai{54, 140, 0}.key();
  first.label = "Xe-140";
  first.seedAtoms = 1.2345678901234567e+22;
  first.importance = 2.6505432109876543e-07;
  first.value = first.seedAtoms * first.importance;
  first.fraction = 0.5;
  first.cumulativeFraction = 0.5;
  first.rank = 1;

  SeedShare second = first;
  second.key = Zai{55, 140, 0}.key();
  second.label = "Cs-140";
  second.importance = 2.6505105987654321e-07;  // differs from the first at the seventh digit
  second.value = second.seedAtoms * second.importance;
  second.cumulativeFraction = 1.0;
  second.rank = 2;

  a.shares = {first, second};
  return a;
}

std::string attributionAsCsv(const SeedAttribution& a) {
  std::ostringstream out;
  writeAttribution(out, a, ReportContext{}, ReportFormat::Csv);
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
  EXPECT_NE(text.find("usually is\n\nNuSIFT"), std::string::npos) << text;
}

// The fraction is one of emitted energy, and exposure per unit energy climbs steeply below
// 100 keV, so a soft continuum costs more exposure than its share of the energy. The footer
// has to say "of that order" and which way the error leans, not claim an equality it cannot.
TEST(Report, TextFooterStatesTheEnergyFractionAsAnOrderNotAnEquality) {
  const std::string text = asText(exposureWithATraceUnmodelled(), ReportContext{});
  EXPECT_NE(text.find("of that order"), std::string::npos) << text;
  EXPECT_NE(text.find("softer than the lines"), std::string::npos) << text;
  EXPECT_EQ(text.find("by roughly that much"), std::string::npos) << text;
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

// --- the buildup caveat --------------------------------------------------------

Ranking exposureThroughAir(double opticalDepth, double buildup) {
  Ranking ranking = oneRow("Cs-137");
  ranking.metric = Metric::Exposure;
  ranking.unit = Unit::RoentgenPerHour;
  ranking.meanOpticalDepth = opticalDepth;
  ranking.buildup = buildup;
  return ranking;
}

// Past half a mean free path the scattered photons an uncollided calculation leaves out are the
// largest thing the number is missing, and a reader who never set --buildup has to be told so
// beside the number rather than in a document they did not open.
TEST(Report, TextFootnotesAThickAirPathLeftUncorrected) {
  const std::string text = asText(exposureThroughAir(1.9, 1.0), ReportContext{});
  EXPECT_NE(text.find("1.9 mean free paths"), std::string::npos) << text;
  EXPECT_NE(text.find("--buildup"), std::string::npos) << text;
  EXPECT_EQ(text.back(), '\n') << text;
}

// A meter of air is a few hundredths of a mean free path. A paragraph about scatter there would
// be noise beside a percent-level effect, and would teach readers to skip the paragraph.
TEST(Report, TextSaysNothingAboutScatterOverAThinPath) {
  const std::string text = asText(exposureThroughAir(0.02, 1.0), ReportContext{});
  EXPECT_EQ(text.find("mean free path"), std::string::npos) << text;
}

// A caller who set a buildup factor has made their own assumption about scatter, and telling
// them they left it out would be false.
TEST(Report, TextTrustsAnExplicitBuildupFactor) {
  const std::string text = asText(exposureThroughAir(1.9, 2.5), ReportContext{});
  EXPECT_EQ(text.find("mean free path"), std::string::npos) << text;
}

// Activity has no air path. The fields exist on every ranking, so the writer has to go by the
// metric rather than by whether the numbers happen to be zero.
TEST(Report, TextSaysNothingAboutAirForAnActivityRanking) {
  Ranking ranking = oneRow("Cs-137");
  ranking.meanOpticalDepth = 3.0;
  const std::string text = asText(ranking, ReportContext{});
  EXPECT_EQ(text.find("mean free path"), std::string::npos) << text;
}

// The same caveats, as numbers a script can act on -- and absent from an activity ranking,
// where a zero would read as a measurement of a path that does not exist.
TEST(Report, JsonCarriesTheExposureCaveatsAsNumbers) {
  const std::string exposure = asJson(exposureThroughAir(1.9, 1.0), ReportContext{});
  EXPECT_EQ(jsonNumberValue(exposure, "mean_optical_depth"), 1.9) << exposure;
  EXPECT_EQ(jsonNumberValue(exposure, "buildup"), 1.0) << exposure;
  EXPECT_NE(exposure.find("\"unmodeled_energy_fraction\""), std::string::npos) << exposure;

  const std::string activity = asJson(oneRow("Cs-137"), ReportContext{});
  EXPECT_EQ(activity.find("mean_optical_depth"), std::string::npos) << activity;
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

// --- attribution CSV -----------------------------------------------------------

// The attribution table is written to be loaded again, so it is held to the same standard the
// ranking CSV is: every number reads back as the double it came from. Written with %.4e, the
// two importances below serialize to the same five digits, and `value` stops equalling
// seed_atoms x importance in the reloaded frame.
TEST(Report, AttributionCsvNumbersReadBackAsTheValuesTheyCameFrom) {
  const SeedAttribution a = twoNearlyIdenticalShares();
  const std::string text = attributionAsCsv(a);

  // Columns: time_s, rank, seed, key, seed_atoms, importance, value, unit, fraction,
  // cumulative, pinned.
  const std::vector<std::string> first = csvFields(text, 1);
  const std::vector<std::string> second = csvFields(text, 2);
  ASSERT_EQ(first.size(), 11u);
  ASSERT_EQ(second.size(), 11u);

  EXPECT_EQ(std::strtod(first[0].c_str(), nullptr), a.time);
  EXPECT_EQ(std::strtod(first[4].c_str(), nullptr), a.shares[0].seedAtoms);
  EXPECT_EQ(std::strtod(first[5].c_str(), nullptr), a.shares[0].importance);
  EXPECT_EQ(std::strtod(first[6].c_str(), nullptr), a.shares[0].value);
  EXPECT_EQ(std::strtod(second[5].c_str(), nullptr), a.shares[1].importance);

  // The distinction survives the round trip, which is the whole point: two seeds a fraction of
  // a percent apart must not reload as the same number.
  EXPECT_NE(first[5], second[5]) << text;

  // And the identity a consumer will recompute still holds.
  const double atoms = std::strtod(first[4].c_str(), nullptr);
  const double importance = std::strtod(first[5].c_str(), nullptr);
  EXPECT_EQ(atoms * importance, a.shares[0].value);
}

// A share is a number of becquerel or of R/h at one instant. A table stating neither is
// uninterpretable once it leaves the terminal, however exactly its digits round-trip.
TEST(Report, AttributionCsvCarriesItsTimeAndUnit) {
  SeedAttribution a = twoNearlyIdenticalShares();
  a.metric = Metric::Exposure;
  a.unit = Unit::RoentgenPerHour;

  const std::string text = attributionAsCsv(a);
  const std::vector<std::string> header = csvFields(text, 0);
  ASSERT_FALSE(header.empty());
  EXPECT_EQ(header[0], "time_s");
  EXPECT_EQ(header[7], "unit");

  const std::vector<std::string> row = csvFields(text, 1);
  ASSERT_EQ(row.size(), 11u);
  EXPECT_EQ(std::strtod(row[0].c_str(), nullptr), a.time);
  EXPECT_EQ(row[7], unitName(Unit::RoentgenPerHour));
}

// An attributed exposure is the ranking's number seen from the other side, so the JSON says
// how far the model was stretched to get it in the same words and under the same keys.
TEST(Report, AttributionJsonCarriesTheExposureCaveats) {
  SeedAttribution a = twoNearlyIdenticalShares();
  a.metric = Metric::Exposure;
  a.unit = Unit::RoentgenPerHour;
  a.unmodeledEnergyFraction = 0.031;
  a.meanOpticalDepth = 0.78;
  a.buildup = 1.0;

  std::ostringstream out;
  ReportContext context;
  context.storeLibrary = "ENDF/B-VIII.1";
  context.geometry = "point source at 100 m";
  writeAttribution(out, a, context, ReportFormat::Json);
  const std::string text = out.str();

  EXPECT_NE(text.find("\"unmodeled_energy_fraction\":0.031"), std::string::npos) << text;
  EXPECT_NE(text.find("\"mean_optical_depth\":0.78"), std::string::npos) << text;
  EXPECT_NE(text.find("\"buildup\":1"), std::string::npos) << text;
  EXPECT_NE(text.find("\"library\":\"ENDF/B-VIII.1\""), std::string::npos) << text;
  EXPECT_NE(text.find("\"model\":\"point source at 100 m\""), std::string::npos) << text;
}

// An activity attribution never computed an exposure, so it carries no reservation about one.
TEST(Report, AttributionJsonOmitsExposureCaveatsForActivity) {
  std::ostringstream out;
  writeAttribution(out, twoNearlyIdenticalShares(), ReportContext{}, ReportFormat::Json);
  const std::string text = out.str();
  EXPECT_EQ(text.find("unmodeled_energy_fraction"), std::string::npos) << text;
  EXPECT_EQ(text.find("mean_optical_depth"), std::string::npos) << text;
}

// The warnings the forward ranking prints, printed here too. Reported without them, the same
// exposure looked better characterised for having been decomposed.
TEST(Report, AttributionTextFootnotesTheSameCaveatsTheRankingDoes) {
  SeedAttribution a = twoNearlyIdenticalShares();
  a.metric = Metric::Exposure;
  a.unit = Unit::RoentgenPerHour;
  a.unmodeledEnergyFraction = 0.031;
  a.meanOpticalDepth = 0.78;
  a.buildup = 1.0;
  a.unmodeledContinuum = {"Y-90"};

  std::ostringstream out;
  writeAttribution(out, a, ReportContext{}, ReportFormat::Text);
  const std::string text = out.str();

  EXPECT_NE(text.find("NuSIFT does not model"), std::string::npos) << text;
  EXPECT_NE(text.find("Y-90"), std::string::npos) << text;
  EXPECT_NE(text.find("mean free paths"), std::string::npos) << text;
}

// --- located events -----------------------------------------------------------

EventReport oneCrossing() {
  EventReport report;
  report.metric = "activity";
  report.curve = "the total";
  report.unit = "Bq";
  report.gridStartSeconds = 3600.0;
  report.gridEndSeconds = 3.0e9;
  report.gridPoints = 60;
  report.hasLevel = true;
  report.level = 5.0e13;

  TrajectoryEvent event;
  event.kind = EventKind::Falling;
  event.timeSeconds = 2.0e9;
  event.value = 5.0e13;
  event.bracketStartSeconds = 1.9e9;
  event.bracketEndSeconds = 2.2e9;
  event.locatedToSeconds = 3.0e8;
  event.refined = false;
  report.events.push_back(event);

  LevelWindow window;
  window.startSeconds = 3600.0;
  window.endSeconds = 2.0e9;
  window.entryObserved = false;
  window.exitObserved = true;
  report.windows.push_back(window);
  return report;
}

std::string eventsAs(const EventReport& report, ReportFormat format) {
  std::ostringstream out;
  writeEvents(out, report, ReportContext{}, format);
  return out.str();
}

// The bracket is the honest error bar on the instant. A report that printed the crossing time
// alone would claim a precision the sampling does not support, so every format carries it.
TEST(Report, EveryFormatCarriesTheBracketAnEventWasFoundIn) {
  const EventReport report = oneCrossing();

  const std::string text = eventsAs(report, ReportFormat::Text);
  EXPECT_NE(text.find("bracket"), std::string::npos) << text;
  EXPECT_NE(text.find("interpolated"), std::string::npos)
      << "an interpolated event has to say so: " << text;

  const std::string json = eventsAs(report, ReportFormat::Json);
  EXPECT_DOUBLE_EQ(jsonNumberValue(json, "bracket_start_s"), 1.9e9);
  EXPECT_DOUBLE_EQ(jsonNumberValue(json, "located_to_s"), 3.0e8);
  EXPECT_EQ(jsonValue(json, "refined"), "false");

  const std::string csv = eventsAs(report, ReportFormat::Csv);
  EXPECT_NE(csv.find("bracket_start_s"), std::string::npos) << csv;
  EXPECT_NE(csv.find("located_to_s"), std::string::npos) << csv;
}

// A reader who does not know the grid decides what is findable will read an empty list as "it
// never happens" rather than "the sampling did not resolve it".
TEST(Report, TheTextReportSaysWhatTheGridCouldNotHaveSeen) {
  EventReport report = oneCrossing();
  report.events.clear();
  report.windows.clear();

  const std::string text = eventsAs(report, ReportFormat::Text);
  EXPECT_NE(text.find("crossings: none"), std::string::npos) << text;
  EXPECT_NE(text.find("rises and falls back between two samples"), std::string::npos) << text;
}

// An edge outside the grid is a bound, not a crossing, and the report has to say which it is.
TEST(Report, AnUnobservedWindowEdgeIsNamedRatherThanPrintedAsACrossing) {
  const std::string text = eventsAs(oneCrossing(), ReportFormat::Text);
  EXPECT_NE(text.find("already above when the grid started"), std::string::npos) << text;

  const std::string json = eventsAs(oneCrossing(), ReportFormat::Json);
  EXPECT_NE(json.find("\"entry_observed\": false"), std::string::npos) << json;
}

// Two tables would not be a CSV, and a reader given only the crossings would lose which grid
// edges were never observed. One table, with `kind` saying which row is which.
TEST(Report, EventsAndWindowsShareOneCsvUnderAKindColumn) {
  const std::string csv = eventsAs(oneCrossing(), ReportFormat::Csv);
  EXPECT_EQ(csv.substr(0, 4), "kind");
  EXPECT_NE(csv.find("\nfalling,"), std::string::npos) << csv;
  EXPECT_NE(csv.find("\nwindow,"), std::string::npos) << csv;
}

// The same pair of instants means opposite things on the two sides of a level -- a stay time
// above it, a waiting period below it -- so every format has to say which side it is reporting
// rather than leaving the reader to assume the usual one.
TEST(Report, WindowsSayWhichSideOfTheLevelTheyHold) {
  EventReport report = oneCrossing();
  report.windowsBelowLevel = true;

  const std::string text = eventsAs(report, ReportFormat::Text);
  EXPECT_NE(text.find("below the level:"), std::string::npos) << text;
  EXPECT_EQ(text.find("above the level:"), std::string::npos) << text;
  EXPECT_NE(text.find("already below when the grid started"), std::string::npos)
      << "an unobserved edge is described on the side it is open on: " << text;

  const std::string json = eventsAs(report, ReportFormat::Json);
  EXPECT_NE(json.find("\"windows_side\": \"below\""), std::string::npos) << json;

  const std::string csv = eventsAs(report, ReportFormat::Csv);
  EXPECT_NE(csv.find(",side"), std::string::npos) << csv;
  EXPECT_NE(csv.find(",below\n"), std::string::npos) << csv;

  // And the default is the other side, said just as explicitly.
  const std::string above = eventsAs(oneCrossing(), ReportFormat::Json);
  EXPECT_NE(above.find("\"windows_side\": \"above\""), std::string::npos) << above;
}

TEST(Report, ARatioCurveIsNotGivenAUnit) {
  EventReport report = oneCrossing();
  report.curve = "Zr-95 / Nb-95";
  report.unit.clear();

  const std::string json = eventsAs(report, ReportFormat::Json);
  EXPECT_NE(json.find("\"unit\": null"), std::string::npos)
      << "a ratio is dimensionless, and naming a unit would say it is a count of becquerel: "
      << json;
}

// --- maximum allowable scale ---------------------------------------------------

std::vector<Criterion> oneCriterion() {
  Criterion criterion;
  criterion.name = "A2 transport";
  criterion.limit = 3.7e13;
  // ResponseSpec defaults to instantaneous activity in becquerel, which is what the limit above
  // is quoted in.
  return std::vector<Criterion>{criterion};
}

std::vector<AllowableScale> boundedThenNot() {
  std::vector<AllowableScale> scaled(2);

  scaled[0].timeSeconds = 3600.0;
  scaled[0].bounded = true;
  scaled[0].scale = 0.25;
  scaled[0].bindingIndex = 0;
  CriterionHeadroom bound;
  bound.name = "A2 transport";
  bound.limit = 3.7e13;
  bound.response = 1.48e14;
  bound.fraction = 4.0;
  bound.scale = 0.25;
  bound.binding = true;
  scaled[0].criteria.push_back(bound);
  scaled[0].limiting.push_back(LimitingContributor{ContributorId{}, "Cs-137", 0.62});

  scaled[1].timeSeconds = 3.0e9;
  scaled[1].bounded = false;
  CriterionHeadroom free;
  free.name = "A2 transport";
  free.limit = 3.7e13;
  free.unbounded = true;
  scaled[1].criteria.push_back(free);
  return scaled;
}

std::string allowableAs(ReportFormat format) {
  const std::vector<Criterion> criteria = oneCriterion();
  std::ostringstream out;
  writeAllowable(out, boundedThenNot(), criteria, ReportContext{}, format);
  return out.str();
}

// "Nothing here is limited by this" and "this allows an enormous multiple" are different
// statements, and a report that printed the second when it meant the first would be lying in
// the direction that matters.
TEST(Report, AnUnconstrainedTimeIsNamedRatherThanGivenAHugeNumber) {
  const std::string text = allowableAs(ReportFormat::Text);
  EXPECT_NE(text.find("unbounded"), std::string::npos) << text;
  EXPECT_NE(text.find("not\n  that a very large multiple is permitted"), std::string::npos)
      << "the footnote has to explain what unbounded means: " << text;

  // null, not a number: the only encoding a parser cannot mistake for a bound of zero.
  const std::string json = allowableAs(ReportFormat::Json);
  EXPECT_NE(json.find("\"scale\": null"), std::string::npos) << json;
  EXPECT_NE(json.find("\"binding\": null"), std::string::npos) << json;
  EXPECT_NE(json.find("\"unbounded\": true"), std::string::npos) << json;
}

TEST(Report, TheAllowableTextNamesTheBindingCriterionAndWhatDrivesIt) {
  const std::string text = allowableAs(ReportFormat::Text);
  EXPECT_NE(text.find("A2 transport"), std::string::npos) << text;
  EXPECT_NE(text.find("Cs-137"), std::string::npos) << text;
}

// Long format, one row per time per criterion: every criterion's headroom is in the file, not
// only the binding one's, which is what makes "how close was the runner-up" answerable.
TEST(Report, TheAllowableCsvCarriesEveryCriterionAtEveryTime) {
  const std::string csv = allowableAs(ReportFormat::Csv);
  EXPECT_EQ(csv.substr(0, 6), "time_s");
  EXPECT_NE(csv.find("criterion_scale"), std::string::npos) << csv;
  EXPECT_NE(csv.find("allowed_scale"), std::string::npos) << csv;

  // Two times, one criterion, so two data rows follow the header.
  int rows = 0;
  for (const char c : csv) {
    if (c == '\n') {
      ++rows;
    }
  }
  EXPECT_EQ(rows, 3) << csv;
}

// --- counterfactual interventions -----------------------------------------------

InterventionStudy twoAlternatives() {
  InterventionStudy study;
  study.metric = Metric::Exposure;
  study.unit = Unit::SievertPerHour;
  study.interventionTimeSeconds = 2.592e6;
  study.responseTimeSeconds = 9.4672e8;
  study.baseline = 4.6102;

  InterventionEffect big;
  big.name = "Cs separation";
  big.removed = 4.61;
  big.response = study.baseline - big.removed;
  big.removedFraction = big.removed / study.baseline;
  big.contributors.push_back(RemovedContributor{551370, "Cs-137", 1.6e23, 4.61, 1.0});
  study.effects.push_back(big);

  InterventionEffect small;
  small.name = "Sr separation";
  small.removed = 3.1e-5;
  small.response = study.baseline - small.removed;
  small.removedFraction = small.removed / study.baseline;
  small.contributors.push_back(RemovedContributor{380900, "Sr-90", 2.0e22, 3.1e-5, 1.0});
  study.effects.push_back(small);
  return study;
}

std::string interventionsAs(ReportFormat format) {
  std::ostringstream out;
  writeInterventions(out, twoAlternatives(), ReportContext{}, format);
  return out.str();
}

// The rows are alternatives against one baseline, and they look perfectly addable. A reader who
// summed two of them would be describing a schedule nobody computed, so the report says so.
TEST(Report, TheInterventionReportSaysItsRowsAreAlternativesNotASequence) {
  const std::string text = interventionsAs(ReportFormat::Text);
  EXPECT_NE(text.find("two rows do not add"), std::string::npos) << text;
  EXPECT_NE(text.find("baseline"), std::string::npos) << text;
  EXPECT_NE(text.find("Cs separation"), std::string::npos) << text;
}

// A benefit is only interpretable against the baseline it was measured from, so every format
// carries it rather than leaving the reader to reconstruct it from two other columns.
TEST(Report, EveryInterventionFormatCarriesTheBaseline) {
  const std::string json = interventionsAs(ReportFormat::Json);
  EXPECT_DOUBLE_EQ(jsonNumberValue(json, "baseline"), 4.6102);
  EXPECT_NE(json.find("\"removed_fraction\""), std::string::npos) << json;

  const std::string csv = interventionsAs(ReportFormat::Csv);
  EXPECT_EQ(csv.substr(0, 12), "intervention");
  EXPECT_NE(csv.find("removed_fraction"), std::string::npos) << csv;
}

// Long format, one row per intervention per nuclide: a summary-only table would drop what the
// benefit was actually made of, which is the half that says whether it is worth doing.
TEST(Report, TheInterventionCsvCarriesWhatEachBenefitWasMadeOf) {
  const std::string csv = interventionsAs(ReportFormat::Csv);
  EXPECT_NE(csv.find("Cs-137"), std::string::npos) << csv;
  EXPECT_NE(csv.find("Sr-90"), std::string::npos) << csv;
  EXPECT_NE(csv.find("atoms_removed"), std::string::npos) << csv;
}

}  // namespace
}  // namespace nusift
