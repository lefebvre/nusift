#include <gtest/gtest.h>

#include <sstream>
#include <string>

#include "nusift/core/error.hpp"
#include "nusift/io/inventory_io.hpp"
#include "nusift/nucdata/nuclear_data.hpp"
#include "synthetic_chain.hpp"

namespace nusift {
namespace {

// Sn-100 (unstable) and Sb-100 (its stable daughter), both with atomic weights so mass and
// activity units are exercisable.
NuclearData chain() {
  StoreArrays arrays = synth::linearChain({1.0e-3});
  arrays.awr = {99.1, 99.2};
  return NuclearData::fromArrays(std::move(arrays));
}

Inventory parse(const std::string& text, const NuclearData& data,
                const InventoryReadOptions& options = {}) {
  std::istringstream in(text);
  return readInventoryCsv(in, data, "test.csv", options);
}

const Zai kUnstable{50, 100, 0};
const Zai kStable{51, 100, 0};

TEST(InventoryIo, ReadsAPlainThreeColumnFile) {
  const NuclearData data = chain();
  const Inventory inventory = parse("Sn-100, 1.0e20, atoms\nSb-100, 5.0e19, atoms\n", data);

  EXPECT_EQ(inventory.size(), 2);
  EXPECT_DOUBLE_EQ(inventory.atomsOf(kUnstable), 1.0e20);
  EXPECT_DOUBLE_EQ(inventory.atomsOf(kStable), 5.0e19);
}

// Every one of these turns up in a file that came out of a spreadsheet, and any one of them
// breaking the read would look like a corrupt inventory rather than a formatting quirk.
TEST(InventoryIo, ToleratesWhatSpreadsheetsActuallyProduce) {
  const NuclearData data = chain();
  const std::string text =
      "\xEF\xBB\xBF"  // UTF-8 BOM, which Excel writes by default
      "# a leading comment\r\n"
      "\r\n"                         // blank line
      "nuclide, quantity, unit\r\n"  // header, detected rather than required
      "  Sn-100 ,  1.0e20 , atoms \r\n"
      "\r\n"
      "Sb-100, 5.0e19, atoms   # a trailing comment\r\n";
  const Inventory inventory = parse(text, data);

  EXPECT_EQ(inventory.size(), 2);
  EXPECT_DOUBLE_EQ(inventory.atomsOf(kUnstable), 1.0e20);
  EXPECT_DOUBLE_EQ(inventory.atomsOf(kStable), 5.0e19);
}

// A comment before the header must not consume the "first row" state. When it did, the
// header parsed as data and the error reported that "nuclide" is not a nuclide name --
// technically true, and useless.
TEST(InventoryIo, HeaderIsStillDetectedAfterLeadingComments) {
  const NuclearData data = chain();
  const Inventory inventory =
      parse("# comment\n# another\nnuclide,quantity,unit\nSn-100,1.0e20,atoms\n", data);
  EXPECT_EQ(inventory.size(), 1);
  EXPECT_DOUBLE_EQ(inventory.atomsOf(kUnstable), 1.0e20);
}

// A file with no header is equally valid, and must not lose its first row to header
// detection.
TEST(InventoryIo, FirstRowIsKeptWhenThereIsNoHeader) {
  const NuclearData data = chain();
  const Inventory inventory = parse("Sn-100,1.0e20,atoms\nSb-100,2.0e20,atoms\n", data);
  EXPECT_EQ(inventory.size(), 2);
  EXPECT_DOUBLE_EQ(inventory.atomsOf(kUnstable), 1.0e20);
}

// A header that names the nuclide and the quantity but spells another column some way the
// reader does not know used to leave that column unmapped. For the unit, every row then
// defaulted to atoms -- 1e6 Bq read as 1e6 atoms -- and for the date, dated assays were summed
// without being reconciled. Both are refused now, naming the column.
TEST(InventoryIo, AHeaderColumnTheReaderDoesNotKnowIsRefused) {
  const NuclearData data = chain();
  for (const char* header : {"nuclide,quantity,unit_of_measure", "nuclide,quantity,unit,assay date",
                             "nuclide,quantity,unit,notes"}) {
    try {
      parse(std::string(header) + "\nSn-100,1.0e6,Bq,2024-01-01\n", data);
      FAIL() << "expected InputError for " << header;
    } catch (const InputError& e) {
      EXPECT_NE(std::string(e.what()).find("line 1"), std::string::npos) << e.what();
      EXPECT_NE(std::string(e.what()).find("does not know"), std::string::npos) << e.what();
    }
  }
  EXPECT_THROW(parse("nuclide,unit\nSn-100,Bq\n", data), InputError);  // no quantity named
  EXPECT_THROW(parse("nuclide,quantity,unit,units\nSn-100,1,Bq,Bq\n", data), InputError);
}

// The synonyms the reader documents still map, in any order, and a blank trailing column that a
// spreadsheet exported is not a column anyone named.
TEST(InventoryIo, AHeaderMapsColumnsByNameInAnyOrder) {
  const NuclearData data = chain();
  const Inventory inventory = parse("Units,Isotope,Amount,\nBq,Sn-100,1.0e6,\n", data);
  EXPECT_NEAR(inventory.atomsOf(kUnstable), 1.0e6 / data.decayConstant(data.indexOf(kUnstable)),
              1.0);
}

// A first data row whose quantity is malformed used to pass for a header -- its second field was
// not a number -- and was dropped without a word. A row that opens with a nuclide is data.
TEST(InventoryIo, AMalformedFirstRowIsRefusedRatherThanTakenForAHeader) {
  const NuclearData data = chain();
  for (const char* row : {"Sn-100,1.0e6Bq\n", "Sn-100,,Bq\n", "Sn-100,nan,Bq\n"}) {
    try {
      parse(std::string(row) + "Sb-100,1.0e20,atoms\n", data);
      FAIL() << "expected InputError for " << row;
    } catch (const InputError& e) {
      EXPECT_NE(std::string(e.what()).find("line 1"), std::string::npos) << e.what();
    }
  }
}

// A header naming none of the columns is not describing them, so it is skipped and the columns
// read positionally, as they always were. Nothing in it could be misread.
TEST(InventoryIo, AHeaderNamingNoKnownColumnIsReadPositionally) {
  const NuclearData data = chain();
  const Inventory inventory = parse("name,qty,uom\nSn-100,1.0e20,atoms\n", data);
  EXPECT_DOUBLE_EQ(inventory.atomsOf(kUnstable), 1.0e20);
}

// Mixed units in one file is the normal case, not an edge case: an inventory is assembled
// from sources that each report in their own units.
TEST(InventoryIo, EachRowCarriesItsOwnUnit) {
  const NuclearData data = chain();
  const Inventory inventory = parse("Sn-100, 1.0e6, Bq\nSb-100, 1.0, g\n", data);

  const int unstableIndex = data.indexOf(kUnstable);
  const int stableIndex = data.indexOf(kStable);
  EXPECT_NEAR(inventory.atomsOf(kUnstable), 1.0e6 / data.decayConstant(unstableIndex), 1e6);
  EXPECT_NEAR(inventory.atomsOf(kStable), units::kAvogadro / data.molarMassGPerMol(stableIndex),
              1e12);
}

// Atoms is the only unit that needs no nuclear data, so it is the only safe default.
TEST(InventoryIo, MissingUnitColumnMeansAtoms) {
  const NuclearData data = chain();
  const Inventory inventory = parse("Sn-100, 1.0e20\n", data);
  EXPECT_DOUBLE_EQ(inventory.atomsOf(kUnstable), 1.0e20);
}

TEST(InventoryIo, RepeatedNuclideAccumulates) {
  const NuclearData data = chain();
  const Inventory inventory = parse("Sn-100, 1.0e20\nSn-100, 5.0e19\n", data);
  EXPECT_EQ(inventory.size(), 1);
  EXPECT_DOUBLE_EQ(inventory.atomsOf(kUnstable), 1.5e20);
}

// Errors have to name the line, or a bad row in a long file is unfindable.
TEST(InventoryIo, ErrorsNameTheOffendingLine) {
  const NuclearData data = chain();
  try {
    parse("Sn-100, 1.0e20\nXx-137, 5.0\n", data);
    FAIL() << "expected InputError";
  } catch (const InputError& e) {
    const std::string what = e.what();
    EXPECT_NE(what.find("line 2"), std::string::npos) << what;
    EXPECT_NE(what.find("Xx-137"), std::string::npos) << what;
  }
}

TEST(InventoryIo, RejectsMalformedRows) {
  const NuclearData data = chain();
  EXPECT_THROW(parse("Sn-100\n", data), InputError);                 // no quantity
  EXPECT_THROW(parse("Sn-100, banana, atoms\n", data), InputError);  // not a number
  EXPECT_THROW(parse("Sn-100, -1.0, atoms\n", data), InputError);    // negative
  EXPECT_THROW(parse("Sn-100, 1.0, furlongs\n", data), InputError);  // not a unit
  EXPECT_THROW(parse("# only comments\n", data), InputError);        // nothing usable
}

// A nuclide with no evaluated data cannot be decayed. Failing by default matters: a silently
// dropped row understates every ranking that follows, with nothing in the output to say so.
TEST(InventoryIo, UnknownNuclideIsAnErrorByDefault) {
  const NuclearData data = chain();
  try {
    parse("Cm-244, 1.0e20, atoms\n", data);
    FAIL() << "expected InputError";
  } catch (const InputError& e) {
    const std::string what = e.what();
    EXPECT_NE(what.find("Cm-244"), std::string::npos) << what;
    EXPECT_NE(what.find("--ignore-unknown"), std::string::npos)
        << "the message should say how to proceed: " << what;
  }
}

TEST(InventoryIo, IgnoreUnknownSkipsAndReports) {
  const NuclearData data = chain();
  std::ostringstream warnings;
  InventoryReadOptions options;
  options.ignoreUnknown = true;
  options.warnings = &warnings;

  const Inventory inventory = parse("Cm-244, 1.0e20\nSn-100, 5.0e19\n", data, options);
  EXPECT_EQ(inventory.size(), 1);
  EXPECT_DOUBLE_EQ(inventory.atomsOf(kUnstable), 5.0e19);
  EXPECT_NE(warnings.str().find("Cm-244"), std::string::npos)
      << "a skipped row must still be reported: " << warnings.str();
}

// A conversion the store cannot support is reported against the row that asked for it.
TEST(InventoryIo, ImpossibleConversionNamesTheRow) {
  const NuclearData data = chain();
  try {
    parse("Sn-100, 1.0e20, atoms\nSb-100, 1.0e6, Bq\n", data);  // Sb-100 is stable
    FAIL() << "expected InputError";
  } catch (const InputError& e) {
    const std::string what = e.what();
    EXPECT_NE(what.find("line 2"), std::string::npos) << what;
    EXPECT_NE(what.find("stable"), std::string::npos) << what;
  }
}

TEST(InventoryIo, CsvRoundTrips) {
  const NuclearData data = chain();
  const Inventory original = parse("Sn-100, 1.0e20\nSb-100, 5.0e19\n", data);

  std::ostringstream written;
  writeInventoryCsv(written, original, data, Quantity::Atoms);

  std::istringstream reread(written.str());
  const Inventory back = readInventoryCsv(reread, data, "roundtrip.csv");

  ASSERT_EQ(back.size(), original.size());
  EXPECT_NEAR(back.atomsOf(kUnstable), original.atomsOf(kUnstable), 1e10);
  EXPECT_NEAR(back.atomsOf(kStable), original.atomsOf(kStable), 1e10);
}

// A snapshot written in a MASS unit and read back is the same inventory. This is the path
// `decay --write --write-units kg` takes, and it runs the molar mass in both directions: out
// through atoms -> grams, back in through grams -> atoms. A wrong or missing atomic weight
// cancels itself in a round trip only if it is the same wrong weight both ways, so this pins
// the plumbing rather than the value -- which is what the validation suite is for.
TEST(InventoryIo, MassUnitsRoundTrip) {
  const NuclearData data = chain();
  const Inventory original = parse("Sn-100, 1.0e20\nSb-100, 5.0e19\n", data);

  std::ostringstream written;
  writeInventoryCsv(written, original, data, Quantity::Kilograms);
  EXPECT_NE(written.str().find("kg"), std::string::npos) << written.str();

  std::istringstream reread(written.str());
  const Inventory back = readInventoryCsv(reread, data, "roundtrip-kg.csv");

  ASSERT_EQ(back.size(), original.size());
  EXPECT_NEAR(back.atomsOf(kUnstable), original.atomsOf(kUnstable), 1.0e8);
  EXPECT_NEAR(back.atomsOf(kStable), original.atomsOf(kStable), 1.0e8);
}

// Writing in a unit one nuclide cannot express must not abort the whole file or emit a wrong
// number: that row falls back to atoms and says so in its unit column.
TEST(InventoryIo, WriteFallsBackPerRowRatherThanFailing) {
  const NuclearData data = chain();
  const Inventory inventory = parse("Sn-100, 1.0e20\nSb-100, 5.0e19\n", data);

  std::ostringstream written;
  writeInventoryCsv(written, inventory, data, Quantity::Becquerel);
  const std::string text = written.str();

  EXPECT_NE(text.find("Sn-100"), std::string::npos);
  EXPECT_NE(text.find("Bq"), std::string::npos) << text;
  // Sb-100 is stable, so it cannot be written in becquerel and falls back.
  EXPECT_NE(text.find("Sb-100"), std::string::npos);
  EXPECT_NE(text.find("atoms"), std::string::npos) << text;
}

TEST(InventoryIo, JsonRoundTrips) {
  const NuclearData data = chain();
  const Inventory original = parse("Sn-100, 1.0e20\nSb-100, 5.0e19\n", data);

  std::ostringstream written;
  writeInventoryJson(written, original, data, Quantity::Atoms);

  std::istringstream reread(written.str());
  const Inventory back = readInventoryJson(reread, data, "roundtrip.json");

  ASSERT_EQ(back.size(), original.size());
  EXPECT_NEAR(back.atomsOf(kUnstable), original.atomsOf(kUnstable), 1e10);
  EXPECT_NEAR(back.atomsOf(kStable), original.atomsOf(kStable), 1e10);
}

// An inventory written here is the seed of a later run. Ten significant digits looked like
// plenty and cost 1e-10 of every count on the way out, so a solve re-run from its own written
// inventory started from a different one; every digit the value has must survive the trip.
TEST(InventoryIo, WrittenQuantitiesReadBackExactly) {
  const NuclearData data = chain();
  Inventory original;
  original.add(kUnstable, 1.2345678901234567e20);
  original.add(kStable, 6.02214076e23 / 7.0);

  std::ostringstream csv;
  writeInventoryCsv(csv, original, data, Quantity::Atoms);
  std::istringstream rereadCsv(csv.str());
  const Inventory fromCsv = readInventoryCsv(rereadCsv, data, "roundtrip.csv");
  EXPECT_EQ(fromCsv.atomsOf(kUnstable), original.atomsOf(kUnstable)) << csv.str();
  EXPECT_EQ(fromCsv.atomsOf(kStable), original.atomsOf(kStable)) << csv.str();

  std::ostringstream json;
  writeInventoryJson(json, original, data, Quantity::Atoms);
  std::istringstream rereadJson(json.str());
  const Inventory fromJson = readInventoryJson(rereadJson, data, "roundtrip.json");
  EXPECT_EQ(fromJson.atomsOf(kUnstable), original.atomsOf(kUnstable)) << json.str();
  EXPECT_EQ(fromJson.atomsOf(kStable), original.atomsOf(kStable)) << json.str();
}

TEST(InventoryIo, JsonAcceptsAQuotedOrBareQuantity) {
  const NuclearData data = chain();
  std::istringstream in(
      R"([{"nuclide": "Sn-100", "quantity": 1.0e20, "unit": "atoms"},
          {"nuclide": "Sb-100", "quantity": "5.0e19"}])");
  const Inventory inventory = readInventoryJson(in, data, "test.json");

  EXPECT_EQ(inventory.size(), 2);
  EXPECT_DOUBLE_EQ(inventory.atomsOf(kUnstable), 1.0e20);
  EXPECT_DOUBLE_EQ(inventory.atomsOf(kStable), 5.0e19);
}

TEST(InventoryIo, JsonRejectsMalformedInput) {
  const NuclearData data = chain();
  const auto readIt = [&](const std::string& text) {
    std::istringstream in(text);
    return readInventoryJson(in, data, "test.json");
  };
  EXPECT_THROW(readIt("not json at all"), InputError);
  EXPECT_THROW(readIt(R"([{"quantity": 1.0}])"), InputError);      // no nuclide
  EXPECT_THROW(readIt(R"([{"nuclide": "Sn-100"}])"), InputError);  // no quantity
  EXPECT_THROW(readIt(R"([{"nuclide": "Sn-100", "quantity": -1}])"), InputError);
  // from_chars reads "inf" and "nan" as numbers. An inventory quantity is a physical amount.
  EXPECT_THROW(readIt(R"([{"nuclide": "Sn-100", "quantity": inf}])"), InputError);
  EXPECT_THROW(readIt(R"([{"nuclide": "Sn-100", "quantity": nan}])"), InputError);
}

// Keys are read one by one, by the names the CSV header accepts. A substring search ignored
// "units" -- so 1e6 Bq read as 1e6 atoms -- found keys inside values, and ended an object at a
// brace inside a string.
TEST(InventoryIo, JsonReadsKeysByNameAndRefusesOnesItDoesNotKnow) {
  const NuclearData data = chain();
  const auto readIt = [&](const std::string& text) {
    std::istringstream in(text);
    return readInventoryJson(in, data, "test.json");
  };
  const double perBq = 1.0 / data.decayConstant(data.indexOf(kUnstable));

  EXPECT_NEAR(
      readIt(R"([{"isotope": "Sn-100", "amount": 1.0e6, "units": "Bq"}])").atomsOf(kUnstable),
      1.0e6 * perBq, 1.0);
  EXPECT_NEAR(
      readIt(R"([{"unit": "Bq", "quantity": 1.0e6, "nuclide": "Sn-100"}])").atomsOf(kUnstable),
      1.0e6 * perBq, 1.0);
  // A key that only CONTAINS a known name, a brace inside a string, and an escaped quote.
  EXPECT_THROW(readIt(R"([{"nuclide": "Sn-100", "quantity": 1, "unit_of_measure": "Bq"}])"),
               InputError);
  EXPECT_THROW(readIt(R"([{"nuclide": "Sn-100", "quantity": 1, "note": "a } b"}])"), InputError);
  // Outside the macro, because MSVC misreads a backslash in a raw string passed as an argument.
  const std::string escapedQuote = R"([{"nuclide": "Sn-100", "quantity": 1, "unit": "B\"q"}])";
  EXPECT_THROW(readIt(escapedQuote), InputError);
  // Twice, under two spellings of one column.
  EXPECT_THROW(readIt(R"([{"nuclide": "Sn-100", "quantity": 1, "unit": "Bq", "units": "Ci"}])"),
               InputError);
  EXPECT_THROW(readIt(R"([{"nuclide": "Sn-100", "quantity": {"v": 1}}])"), InputError);
  // null says nothing, as an absent key does.
  EXPECT_DOUBLE_EQ(
      readIt(R"([{"nuclide": "Sn-100", "quantity": 1e20, "unit": null}])").atomsOf(kUnstable),
      1.0e20);
}

// Being liberal about the punctuation is not kindness. A reader that skips any run of commas
// accepts a file no other JSON tool will read, and one that stops at the first ] accepts the
// first half of a file that is two documents long -- and then reports confidently on an
// inventory nobody wrote.
TEST(InventoryIo, JsonRejectsWhatIsNotActuallyJson) {
  const NuclearData data = chain();
  const auto readIt = [&](const std::string& text) {
    std::istringstream in(text);
    return readInventoryJson(in, data, "test.json");
  };
  const std::string entry = R"({"nuclide": "Sn-100", "quantity": 1.0e20})";

  EXPECT_NO_THROW(readIt("[" + entry + "," + entry + "]"));
  EXPECT_THROW(readIt("[" + entry + ",," + entry + "]"), InputError);  // doubled separator
  EXPECT_THROW(readIt("[" + entry + " " + entry + "]"), InputError);   // no separator
  EXPECT_THROW(readIt("[," + entry + "]"), InputError);                // leading separator
  EXPECT_THROW(readIt("[" + entry + ",]"), InputError);                // trailing comma
  EXPECT_THROW(readIt("[" + entry + "] " + entry), InputError);        // content after the array
  EXPECT_THROW(readIt("[" + entry + "][" + entry + "]"), InputError);  // two documents
}

}  // namespace
}  // namespace nusift
