#include "nusift/io/inventory_io.hpp"

#include <algorithm>
#include <cctype>
#include <charconv>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <istream>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

#include "nusift/core/error.hpp"
#include "nusift/core/nuclide_name.hpp"
#include "nusift/engine/reconcile.hpp"
#include "nusift/io/number_format.hpp"
#include "nusift/io/time_spec.hpp"
#include "nusift/nucdata/nuclear_data.hpp"

namespace nusift {
namespace {

constexpr const char* kModule = "inventory file";

std::string_view trim(std::string_view s) {
  while (!s.empty() && std::isspace(static_cast<unsigned char>(s.front())) != 0) {
    s.remove_prefix(1);
  }
  while (!s.empty() && std::isspace(static_cast<unsigned char>(s.back())) != 0) {
    s.remove_suffix(1);
  }
  return s;
}

// A UTF-8 BOM at the start of the file. Excel writes one by default, and without stripping
// it the first nuclide name silently becomes unparseable in a way that is invisible in a
// text editor.
void stripBom(std::string& line) {
  if (line.size() >= 3 && static_cast<unsigned char>(line[0]) == 0xEF &&
      static_cast<unsigned char>(line[1]) == 0xBB && static_cast<unsigned char>(line[2]) == 0xBF) {
    line.erase(0, 3);
  }
}

std::vector<std::string_view> splitFields(std::string_view line) {
  std::vector<std::string_view> fields;
  std::size_t start = 0;
  while (true) {
    const std::size_t at = line.find(',', start);
    if (at == std::string_view::npos) {
      fields.push_back(trim(line.substr(start)));
      break;
    }
    fields.push_back(trim(line.substr(start, at - start)));
    start = at + 1;
  }
  return fields;
}

bool parseNumber(std::string_view text, double& out) {
  if (text.empty()) {
    return false;
  }
  const auto result = std::from_chars(text.data(), text.data() + text.size(), out);
  // from_chars spells "inf" and "nan" as numbers. An inventory quantity is a physical amount,
  // so neither is one -- and refusing them here is what keeps the complaint on the row that
  // carries the value rather than on the inventory it eventually corrupts.
  return result.ec == std::errc{} && result.ptr == text.data() + text.size() && std::isfinite(out);
}

[[noreturn]] void failAt(const std::string& source, int line, const std::string& what) {
  throw InputError(tagged(kModule, source + " line " + std::to_string(line) + ": " + what));
}

std::string toLower(std::string_view s) {
  std::string out(s);
  std::transform(out.begin(), out.end(), out.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return out;
}

// Which column holds what. Positional by default -- nuclide, quantity, unit, assayed,
// uncertainty -- and taken from the header row when one names the columns.
//
// The mapping exists because the positional order cannot serve both optional columns at once.
// `assayed` shipped as the fourth field, so a file wanting an uncertainty and no dates would
// have to write "Cs-137,1e14,Bq,,5%" and hope the empty field is noticed. A header that says
// `nuclide,quantity,unit,uncertainty` says it instead.
struct ColumnMap {
  int nuclide = 0;
  int quantity = 1;
  int unit = 2;
  int assayed = 3;
  int uncertainty = 4;
};

enum class Column { Nuclide, Quantity, Unit, Assayed, Uncertainty, Unknown };

// What a column name means, with the obvious synonyms a spreadsheet is likely to carry. Shared by
// the CSV header and the JSON keys, so the two formats accept exactly the same names.
Column columnNamed(std::string_view text) {
  const std::string name = toLower(trim(text));
  if (name == "nuclide" || name == "isotope" || name == "radionuclide") {
    return Column::Nuclide;
  }
  if (name == "quantity" || name == "amount" || name == "value") {
    return Column::Quantity;
  }
  if (name == "unit" || name == "units") {
    return Column::Unit;
  }
  if (name == "assayed" || name == "date" || name == "assay_date") {
    return Column::Assayed;
  }
  if (name == "uncertainty" || name == "sigma" || name == "error") {
    return Column::Uncertainty;
  }
  return Column::Unknown;
}

constexpr const char* kColumnNames =
    "nuclide, quantity, unit, assayed, uncertainty (or isotope, amount, units, date, sigma)";

// Classify a first row whose quantity field is not a number, and map the columns if it is a
// header. Returns false when the row is data after all -- one whose quantity is simply bad, which
// the caller then refuses on its own line.
//
// Each refusal here closes a way a file was read wrongly with no error. A header that names the
// nuclide and the quantity but spells the unit column some other way used to leave the unit
// unmapped, and every row then defaulted to ATOMS: 1e14 Bq of Cs-137 became 1e14 atoms, seven
// decades low. A misnamed date column did the same to the dates, so dated assays were summed
// without being reconciled. And a first data row with a malformed quantity used to pass for a
// header and vanish.
//
// What is still accepted is a header naming none of the columns at all, which is read
// positionally as it always was: such a header is not describing the columns, so ignoring it
// cannot misread one.
bool takeHeader(const std::vector<std::string_view>& fields, const std::string& source, int line,
                ColumnMap& out) {
  if (parseNuclideName(fields.front()).has_value()) {
    return false;
  }
  ColumnMap mapped{-1, -1, -1, -1, -1};
  std::vector<std::string> unknown;
  int recognised = 0;
  for (int i = 0; i < static_cast<int>(fields.size()); ++i) {
    const std::string_view name = trim(fields[static_cast<std::size_t>(i)]);
    // A blank name is an unnamed trailing column, which spreadsheets export freely. Nothing reads
    // it, and nothing could be misread through it.
    if (name.empty()) {
      continue;
    }
    int* slot = nullptr;
    switch (columnNamed(name)) {
      case Column::Nuclide:
        slot = &mapped.nuclide;
        break;
      case Column::Quantity:
        slot = &mapped.quantity;
        break;
      case Column::Unit:
        slot = &mapped.unit;
        break;
      case Column::Assayed:
        slot = &mapped.assayed;
        break;
      case Column::Uncertainty:
        slot = &mapped.uncertainty;
        break;
      case Column::Unknown:
        unknown.emplace_back(name);
        continue;
    }
    if (*slot >= 0) {
      failAt(source, line, "the header names the column \"" + std::string(name) + "\" twice");
    }
    *slot = i;
    ++recognised;
  }
  if (recognised == 0) {
    return true;
  }
  if (!unknown.empty()) {
    failAt(source, line,
           "the header names a column \"" + unknown.front() +
               "\" this reader does not know, so it cannot say what the column holds -- and a "
               "unit, date or uncertainty it missed would be read as absent. The columns are " +
               kColumnNames);
  }
  if (mapped.nuclide < 0 || mapped.quantity < 0) {
    failAt(source, line,
           std::string("the header names some columns but not ") +
               (mapped.nuclide < 0 ? "the nuclide" : "the quantity") + "; the columns are " +
               kColumnNames);
  }
  out = mapped;
  return true;
}

// The field at `column`, or empty when the row is short or the column absent.
std::string_view fieldAt(const std::vector<std::string_view>& fields, int column) {
  if (column < 0 || column >= static_cast<int>(fields.size())) {
    return {};
  }
  return fields[static_cast<std::size_t>(column)];
}

// An uncertainty as a sheet writes one: absolute in the row's own unit, or relative with a
// trailing percent. Both are common and they are not distinguishable by magnitude, so the
// spelling has to carry it -- "5" beside a quantity of 100 Bq could be either.
double parseUncertainty(std::string_view text, double quantity, const std::string& source,
                        int line) {
  std::string_view field = trim(text);
  const bool relative = !field.empty() && field.back() == '%';
  if (relative) {
    field.remove_suffix(1);
    field = trim(field);
  }
  double value = 0.0;
  if (!parseNumber(field, value)) {
    failAt(source, line,
           "\"" + std::string(trim(text)) +
               "\" is not an uncertainty; give it in the row's own unit, or as a percentage "
               "like 5%");
  }
  if (value < 0.0) {
    failAt(source, line, "an uncertainty cannot be negative");
  }
  return relative ? 0.01 * value * quantity : value;
}

// Convert one row to atoms, or report why it cannot be. Returns false when the row should be
// skipped under ignoreUnknown.
bool rowToAtoms(const Zai& zai, double value, Quantity unit, const NuclearData& data,
                const std::string& source, int line, const InventoryReadOptions& options,
                double& atomsOut) {
  const int index = data.indexOf(zai);
  if (index < 0) {
    if (!options.ignoreUnknown) {
      failAt(source, line,
             "the data store has no nuclide " + formatNuclideName(zai) +
                 " (pass --ignore-unknown to skip rows like this)");
    }
    if (options.warnings != nullptr) {
      *options.warnings << "  skipping " << formatNuclideName(zai) << " (" << source << " line "
                        << line << "): not in the data store\n";
    }
    return false;
  }
  // toAtoms throws for a conversion the data cannot support; re-raise with the row attached
  // so a bad line in a 500-row file is findable.
  try {
    atomsOut = toAtoms(value, unit, index, data);
  } catch (const InputError& e) {
    failAt(source, line, e.what());
  }
  return true;
}

// One row as read, before anything is grouped. Kept flat rather than accumulated straight into
// an Inventory because rows measured on different dates are not addable: merging them first and
// discovering the dates afterwards would already have lost the thing that makes them separable.
struct DatedRow {
  Zai zai;
  double atoms = 0.0;
  double atomsFromActivity = 0.0;
  double sigmaAtoms = 0.0;
  double dateSeconds = 0.0;
  bool dated = false;
  int line = 0;
};

// Read the optional date field, and hold the file to one convention. `firstDated` and
// `firstUndated` remember where each kind was first seen, so the refusal can name both lines --
// "you dated line 4 and not line 9" is actionable where "mixed dates" is not.
void takeDate(std::string_view text, const std::string& source, int line, DatedRow& row,
              int& firstDated, int& firstUndated) {
  const std::string_view field = trim(text);
  if (!field.empty()) {
    if (!looksLikeCalendarDate(field)) {
      failAt(source, line,
             "\"" + std::string(field) +
                 "\" is not an assay date; dates are ISO-8601, as 2024-03-15 or "
                 "2024-03-15T09:30:00Z");
    }
    try {
      row.dateSeconds = parseCalendarDate(field);
    } catch (const InputError& e) {
      failAt(source, line, e.what());
    }
    row.dated = true;
    if (firstDated == 0) {
      firstDated = line;
    }
  } else if (firstUndated == 0) {
    firstUndated = line;
  }
}

// The refusal the header argues for. Checked once, after the whole file is read, because a file
// may date its later rows and not its earlier ones and the complaint should name both.
void requireOneConvention(const std::string& source, int firstDated, int firstUndated) {
  if (firstDated != 0 && firstUndated != 0) {
    throw InputError(
        tagged(kModule, source + ": line " + std::to_string(firstDated) +
                            " carries an assay date and line " + std::to_string(firstUndated) +
                            " does not. A file dates every row or none: an undated row among "
                            "dated ones has no reading that is not a guess about when it was "
                            "measured, and the guess would not be visible in any answer"));
  }
}

// Rows to assays. Rows sharing a date are one assay; an undated file is one assay whose date is
// never used.
DatedInventory groupRows(std::vector<DatedRow> rows, const std::string& sourceName, bool dated) {
  DatedInventory out;
  out.dated = dated;
  std::stable_sort(rows.begin(), rows.end(), [](const DatedRow& a, const DatedRow& b) {
    return a.dateSeconds < b.dateSeconds;
  });
  for (const DatedRow& row : rows) {
    if (out.groups.empty() || out.groups.back().dateSeconds != row.dateSeconds) {
      AssayGroup group;
      group.dateSeconds = row.dateSeconds;
      group.dated = dated;
      group.label = dated ? sourceName + " @ " + formatCalendarDate(row.dateSeconds) : sourceName;
      group.inventory.setProvenance(group.label);
      out.groups.push_back(std::move(group));
    }
    out.groups.back().inventory.add(row.zai, row.atoms, row.sigmaAtoms, row.atomsFromActivity);
  }
  return out;
}

}  // namespace

DatedInventory readInventoryDatedCsv(std::istream& in, const NuclearData& data,
                                     const std::string& sourceName,
                                     const InventoryReadOptions& options) {
  std::vector<DatedRow> rows;
  ColumnMap columns;
  int firstDated = 0;
  int firstUndated = 0;
  std::string line;
  int lineNumber = 0;
  // "No data row seen yet", which is NOT the same as "on the first line": a file that opens
  // with comments still has its header on the first row that carries content. Letting a
  // comment consume this flag makes the header parse as data, and the resulting error names
  // the word "nuclide" as a bad nuclide name -- confusing, and pointing at the wrong line.
  bool beforeFirstRow = true;
  int accepted = 0;

  while (std::getline(in, line)) {
    ++lineNumber;
    if (lineNumber == 1) {
      stripBom(line);
    }

    // Trailing CR from a file written on Windows and read on a POSIX runner.
    if (!line.empty() && line.back() == '\r') {
      line.pop_back();
    }
    // Comments may be whole-line or trailing.
    if (const std::size_t hash = line.find('#'); hash != std::string::npos) {
      line.erase(hash);
    }
    const std::string_view content = trim(line);
    if (content.empty()) {
      continue;
    }

    const std::vector<std::string_view> fields = splitFields(content);
    double ignored = 0.0;
    if (beforeFirstRow && fields.size() >= 2 && !parseNumber(fields[1], ignored) &&
        takeHeader(fields, sourceName, lineNumber, columns)) {
      beforeFirstRow = false;
      continue;
    }
    beforeFirstRow = false;

    if (fields.size() < 2) {
      failAt(sourceName, lineNumber,
             "expected at least \"nuclide, quantity\", got \"" + std::string(content) + "\"");
    }

    const Zai zai = requireNuclideName(fieldAt(fields, columns.nuclide),
                                       sourceName + " line " + std::to_string(lineNumber));
    double value = 0.0;
    if (!parseNumber(fieldAt(fields, columns.quantity), value)) {
      failAt(sourceName, lineNumber,
             "\"" + std::string(fieldAt(fields, columns.quantity)) + "\" is not a number");
    }
    if (!(value >= 0.0)) {
      failAt(sourceName, lineNumber, "a quantity cannot be negative");
    }

    // A missing unit column means atoms, which is the only unit that needs no nuclear data
    // and so the only safe default.
    Quantity unit = Quantity::Atoms;
    const std::string_view unitText = fieldAt(fields, columns.unit);
    if (!unitText.empty()) {
      if (!parseQuantity(unitText, unit)) {
        failAt(sourceName, lineNumber,
               "\"" + std::string(unitText) +
                   "\" is not a unit (try atoms, mol, g, kg, mg, Bq, kBq, MBq, GBq, TBq, "
                   "Ci, mCi, uCi)");
      }
    }

    DatedRow row;
    row.zai = zai;
    row.line = lineNumber;
    takeDate(fieldAt(fields, columns.assayed), sourceName, lineNumber, row, firstDated,
             firstUndated);

    if (!rowToAtoms(zai, value, unit, data, sourceName, lineNumber, options, row.atoms)) {
      continue;
    }
    row.atomsFromActivity = isActivityQuantity(unit) ? row.atoms : 0.0;
    const std::string_view sigmaText = fieldAt(fields, columns.uncertainty);
    if (!sigmaText.empty()) {
      // Converted through the same call the quantity was, which is exactly right: toAtoms is
      // value * k in every branch, so a sigma in the row's unit becomes a sigma in atoms under
      // the identical factor. The measurement basis needs no separate bookkeeping.
      const double sigma = parseUncertainty(sigmaText, value, sourceName, lineNumber);
      if (!rowToAtoms(zai, sigma, unit, data, sourceName, lineNumber, options, row.sigmaAtoms)) {
        continue;
      }
    }
    rows.push_back(row);
    ++accepted;
  }

  if (accepted == 0) {
    throw InputError(tagged(kModule, sourceName + " contains no usable inventory rows"));
  }
  requireOneConvention(sourceName, firstDated, firstUndated);
  return groupRows(std::move(rows), sourceName, firstDated != 0);
}

DatedInventory readInventoryDatedJson(std::istream& in, const NuclearData& data,
                                      const std::string& sourceName,
                                      const InventoryReadOptions& options) {
  // Deliberately a small hand-rolled reader rather than a JSON dependency: the accepted
  // shape is one flat array of flat objects with five known keys, and the error messages a
  // purpose-built parser can give ("line 12: ...") are better than a generic one's.
  std::ostringstream buffer;
  buffer << in.rdbuf();
  const std::string text = buffer.str();

  std::vector<DatedRow> rows;
  int firstDated = 0;
  int firstUndated = 0;
  int accepted = 0;
  std::size_t pos = 0;

  const auto skipSpace = [&]() {
    while (pos < text.size() && std::isspace(static_cast<unsigned char>(text[pos])) != 0) {
      ++pos;
    }
  };
  const auto lineOf = [&](std::size_t at) {
    return 1 +
           static_cast<int>(std::count(text.begin(), text.begin() + static_cast<long>(at), '\n'));
  };

  skipSpace();
  if (pos >= text.size() || text[pos] != '[') {
    throw InputError(tagged(kModule, sourceName +
                                         ": expected a JSON array of {nuclide, quantity, unit} "
                                         "objects"));
  }
  ++pos;

  // Entries are comma-separated, and the separator is checked rather than skipped over. A
  // reader that swallows any run of commas accepts "[{...},,{...},]" -- which no other JSON
  // tool will read -- and a file that only NuSIFT can parse is a file that will eventually be
  // handed to something else and rejected, long after it was written.
  bool afterEntry = false;
  while (true) {
    skipSpace();
    if (pos >= text.size()) {
      throw InputError(tagged(kModule, sourceName + ": unterminated JSON array"));
    }
    if (text[pos] == ']') {
      ++pos;
      break;
    }
    if (afterEntry) {
      if (text[pos] != ',') {
        failAt(sourceName, lineOf(pos), "expected a comma between entries");
      }
      ++pos;
      skipSpace();
      if (pos >= text.size()) {
        throw InputError(tagged(kModule, sourceName + ": unterminated JSON array"));
      }
      if (text[pos] == ']') {
        failAt(sourceName, lineOf(pos), "a trailing comma is not valid JSON");
      }
    }
    if (text[pos] != '{') {
      failAt(sourceName, lineOf(pos), "expected an object");
    }
    const std::size_t objectStart = pos;
    ++pos;

    // One object, read key by key rather than searched for substrings. Searching matched a key
    // inside another key or inside a value, ended the object at a '}' that sat in a string, and
    // ignored a key it was not looking for -- so {"units": "Bq"} read as atoms, seven decades low
    // for 1e14 Bq of Cs-137. Every key now has to name a column, by the same names the CSV
    // header accepts, and none may appear twice.
    const auto readString = [&]() {
      // At the opening quote. Escapes are decoded as JSON defines them, and one JSON does not
      // define is refused. Taking the character after a backslash literally read "\tBq" -- a tab,
      // then Bq -- as "tBq", terabecquerels, twelve decades high.
      ++pos;
      std::string value;
      while (pos < text.size() && text[pos] != '"') {
        const char c = text[pos];
        if (static_cast<unsigned char>(c) < 0x20) {
          failAt(sourceName, lineOf(pos), "a control character inside a string must be escaped");
        }
        if (c != '\\') {
          value += c;
          ++pos;
          continue;
        }
        ++pos;
        if (pos >= text.size()) {
          break;
        }
        switch (text[pos]) {
          case '"':
            value += '"';
            break;
          case '\\':
            value += '\\';
            break;
          case '/':
            value += '/';
            break;
          case 'b':
            value += '\b';
            break;
          case 'f':
            value += '\f';
            break;
          case 'n':
            value += '\n';
            break;
          case 'r':
            value += '\r';
            break;
          case 't':
            value += '\t';
            break;
          case 'u':
            failAt(sourceName, lineOf(pos), "\\u escapes are not supported in an inventory");
          default:
            failAt(sourceName, lineOf(pos),
                   std::string("\\") + text[pos] + " is not a JSON escape");
        }
        ++pos;
      }
      if (pos >= text.size()) {
        failAt(sourceName, lineOf(objectStart), "unterminated string");
      }
      ++pos;  // the closing quote
      return value;
    };

    std::string values[5];
    bool seen[5] = {false, false, false, false, false};
    bool afterPair = false;
    while (true) {
      skipSpace();
      if (pos >= text.size()) {
        failAt(sourceName, lineOf(objectStart), "unterminated object");
      }
      if (text[pos] == '}') {
        ++pos;
        break;
      }
      if (afterPair) {
        if (text[pos] != ',') {
          failAt(sourceName, lineOf(pos), "expected a comma between keys");
        }
        ++pos;
        skipSpace();
      }
      if (pos >= text.size() || text[pos] != '"') {
        failAt(sourceName, lineOf(pos), "expected a quoted key");
      }
      const std::string key = readString();
      skipSpace();
      if (pos >= text.size() || text[pos] != ':') {
        failAt(sourceName, lineOf(pos), "expected ':' after \"" + key + "\"");
      }
      ++pos;
      skipSpace();
      if (pos >= text.size()) {
        failAt(sourceName, lineOf(objectStart), "unterminated object");
      }
      std::string value;
      if (text[pos] == '"') {
        value = readString();
      } else if (text[pos] == '{' || text[pos] == '[') {
        failAt(sourceName, lineOf(pos),
               "\"" + key + "\" holds a nested value; give a number or a string");
      } else {
        const std::size_t end = text.find_first_of(",}", pos);
        if (end == std::string::npos) {
          failAt(sourceName, lineOf(objectStart), "unterminated object");
        }
        value = std::string(trim(std::string_view(text).substr(pos, end - pos)));
        pos = end;
      }

      const Column column = columnNamed(key);
      if (column == Column::Unknown) {
        failAt(sourceName, lineOf(objectStart),
               "\"" + key +
                   "\" is not a key this reader knows, so it cannot say what it holds. "
                   "The keys are " +
                   kColumnNames);
      }
      const std::size_t slot = static_cast<std::size_t>(column);
      if (seen[slot]) {
        failAt(sourceName, lineOf(objectStart), "\"" + key + "\" names a column given twice");
      }
      seen[slot] = true;
      // JSON's null says nothing, which is what an absent key says too.
      values[slot] = value == "null" ? std::string() : std::string(trim(value));
      afterPair = true;
    }

    if (!seen[static_cast<std::size_t>(Column::Nuclide)] ||
        !seen[static_cast<std::size_t>(Column::Quantity)]) {
      failAt(sourceName, lineOf(objectStart), "an entry needs \"nuclide\" and \"quantity\"");
    }
    const std::string& nuclideText = values[static_cast<std::size_t>(Column::Nuclide)];
    const std::string& quantityText = values[static_cast<std::size_t>(Column::Quantity)];
    const std::string& unitText = values[static_cast<std::size_t>(Column::Unit)];
    const std::string& assayedText = values[static_cast<std::size_t>(Column::Assayed)];
    const std::string& uncertaintyText = values[static_cast<std::size_t>(Column::Uncertainty)];

    const Zai zai = requireNuclideName(nuclideText,
                                       sourceName + " line " + std::to_string(lineOf(objectStart)));
    double value = 0.0;
    if (!parseNumber(quantityText, value) || !(value >= 0.0)) {
      failAt(sourceName, lineOf(objectStart),
             "\"" + quantityText + "\" is not a non-negative number");
    }
    Quantity unit = Quantity::Atoms;
    if (!unitText.empty() && !parseQuantity(unitText, unit)) {
      failAt(sourceName, lineOf(objectStart), "\"" + unitText + "\" is not a unit");
    }

    DatedRow row;
    row.zai = zai;
    row.line = lineOf(objectStart);
    takeDate(assayedText, sourceName, row.line, row, firstDated, firstUndated);

    if (rowToAtoms(zai, value, unit, data, sourceName, row.line, options, row.atoms)) {
      row.atomsFromActivity = isActivityQuantity(unit) ? row.atoms : 0.0;
      bool usable = true;
      if (!uncertaintyText.empty()) {
        const double sigma = parseUncertainty(uncertaintyText, value, sourceName, row.line);
        usable = rowToAtoms(zai, sigma, unit, data, sourceName, row.line, options, row.sigmaAtoms);
      }
      if (usable) {
        rows.push_back(row);
        ++accepted;
      }
    }
    afterEntry = true;
  }

  // Nothing may follow the array. Trailing content is the signature of a file that is not what
  // it appears to be -- two documents concatenated, or an array pasted over a longer one --
  // and quietly reporting on the first half would be worse than refusing the file.
  skipSpace();
  if (pos < text.size()) {
    failAt(sourceName, lineOf(pos), "unexpected content after the closing ]");
  }

  if (accepted == 0) {
    throw InputError(tagged(kModule, sourceName + " contains no usable inventory entries"));
  }
  requireOneConvention(sourceName, firstDated, firstUndated);
  return groupRows(std::move(rows), sourceName, firstDated != 0);
}

DatedInventory readInventoryDated(const std::string& path, const NuclearData& data,
                                  const InventoryReadOptions& options) {
  std::ifstream in(path, std::ios::binary);
  if (!in) {
    throw InputError(tagged(kModule, "cannot open \"" + path + "\""));
  }
  const std::string lower = toLower(path);
  if (lower.size() >= 5 && lower.compare(lower.size() - 5, 5, ".json") == 0) {
    return readInventoryDatedJson(in, data, path, options);
  }
  return readInventoryDatedCsv(in, data, path, options);
}

namespace {

// Assays to one inventory. Undated files never reach the solver at all, which keeps the
// overwhelmingly common case exactly as cheap as it was; a dated file pays one solve per assay
// older than the epoch, which is the price of the rows being combinable in the first place.
Inventory collapse(const NuclearData& data, DatedInventory dated) {
  if (!dated.dated) {
    Inventory inventory = std::move(dated.groups.front().inventory);
    return inventory;
  }
  const double epoch = latestAssayDate(dated.groups);
  Reconciliation reconciled = reconcile(data, dated.groups, epoch);
  // The epoch is part of what the inventory IS, so it travels in the provenance every report
  // header prints. A reconciled seed whose header said only the file name would look identical
  // to the same file read on a different date.
  std::string where = dated.groups.front().inventory.provenance();
  const std::size_t at = where.find(" @ ");
  if (at != std::string::npos) {
    where.erase(at);
  }
  reconciled.inventory.setProvenance(where + " (" + std::to_string(dated.groups.size()) +
                                     " assays reconciled to " + formatCalendarDate(epoch) + ")");
  return std::move(reconciled.inventory);
}

}  // namespace

Inventory readInventory(const std::string& path, const NuclearData& data,
                        const InventoryReadOptions& options) {
  return collapse(data, readInventoryDated(path, data, options));
}

Inventory readInventoryCsv(std::istream& in, const NuclearData& data, const std::string& sourceName,
                           const InventoryReadOptions& options) {
  return collapse(data, readInventoryDatedCsv(in, data, sourceName, options));
}

Inventory readInventoryJson(std::istream& in, const NuclearData& data,
                            const std::string& sourceName, const InventoryReadOptions& options) {
  return collapse(data, readInventoryDatedJson(in, data, sourceName, options));
}

namespace {

// Express one entry in `unit`, falling back to atoms when the conversion is impossible for
// that particular nuclide. Reporting a fallback beats aborting a whole report over one row,
// and beats emitting a zero that would read as a real value.
double emitValue(double atoms, Quantity unit, int index, const NuclearData& data,
                 Quantity& usedUnit) {
  usedUnit = unit;
  if (unit == Quantity::Atoms || unit == Quantity::Moles) {
    return fromAtoms(atoms, unit, index, data);
  }
  const bool massWithoutWeight =
      (unit == Quantity::Grams || unit == Quantity::Kilograms || unit == Quantity::Milligrams) &&
      data.molarMassGPerMol(index) <= 0.0;
  const bool activityOfStable = data.decayConstant(index) <= 0.0 && unit != Quantity::Grams &&
                                unit != Quantity::Kilograms && unit != Quantity::Milligrams;
  if (massWithoutWeight || activityOfStable) {
    usedUnit = Quantity::Atoms;
    return atoms;
  }
  return fromAtoms(atoms, unit, index, data);
}

}  // namespace

void writeInventoryCsv(std::ostream& out, const Inventory& inventory, const NuclearData& data,
                       Quantity unit) {
  out << "# nusift inventory\n";
  if (!inventory.provenance().empty()) {
    out << "# source: " << inventory.provenance() << "\n";
  }
  out << "nuclide,quantity,unit\n";
  for (const InventoryEntry& entry : inventory.entries()) {
    const Zai zai = Zai::fromKey(entry.zaiKey);
    const int index = data.indexOf(zai);
    // A nuclide the store does not carry is written in atoms, the one unit that needs no data
    // about it -- and labelled atoms, rather than with the unit its number is not in.
    Quantity used = index >= 0 ? unit : Quantity::Atoms;
    const double value = index >= 0 ? emitValue(entry.atoms, unit, index, data, used) : entry.atoms;
    // Every digit the value has, and no more: an inventory written here is read back as the
    // seed of a later run, and a count rounded on the way out is a different inventory.
    out << formatNuclideName(zai) << ',' << shortestRoundTrip(value) << ',' << quantityName(used)
        << '\n';
  }
}

void writeInventoryJson(std::ostream& out, const Inventory& inventory, const NuclearData& data,
                        Quantity unit) {
  out << "[\n";
  bool firstEntry = true;
  for (const InventoryEntry& entry : inventory.entries()) {
    const Zai zai = Zai::fromKey(entry.zaiKey);
    const int index = data.indexOf(zai);
    // A nuclide the store does not carry is written in atoms, the one unit that needs no data
    // about it -- and labelled atoms, rather than with the unit its number is not in.
    Quantity used = index >= 0 ? unit : Quantity::Atoms;
    const double value = index >= 0 ? emitValue(entry.atoms, unit, index, data, used) : entry.atoms;
    if (!firstEntry) {
      out << ",\n";
    }
    firstEntry = false;
    out << "  {\"nuclide\": \"" << formatNuclideName(zai)
        << "\", \"quantity\": " << shortestRoundTrip(value) << ", \"unit\": \""
        << quantityName(used) << "\"}";
  }
  out << "\n]\n";
}

}  // namespace nusift
