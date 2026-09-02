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

// A first row whose second field is not a number is a header naming the columns. Detecting
// it beats requiring one, since half the spreadsheets in the world have one and half do not.
bool looksLikeHeader(const std::vector<std::string_view>& fields) {
  if (fields.size() < 2) {
    return false;
  }
  double ignored = 0.0;
  return !parseNumber(fields[1], ignored);
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

// Recognised spellings, with the obvious synonyms a spreadsheet is likely to carry. A header
// that names neither `nuclide` nor `quantity` is not describing these columns at all, and the
// reader falls back to positional rather than erroring -- that is what it did before this
// existed, and a file that worked yesterday must not stop working today.
bool mapColumns(const std::vector<std::string_view>& header, ColumnMap& out) {
  ColumnMap mapped{-1, -1, -1, -1, -1};
  for (int i = 0; i < static_cast<int>(header.size()); ++i) {
    const std::string name = toLower(header[static_cast<std::size_t>(i)]);
    if (name == "nuclide" || name == "isotope" || name == "radionuclide") {
      mapped.nuclide = i;
    } else if (name == "quantity" || name == "amount" || name == "value") {
      mapped.quantity = i;
    } else if (name == "unit" || name == "units") {
      mapped.unit = i;
    } else if (name == "assayed" || name == "date" || name == "assay_date") {
      mapped.assayed = i;
    } else if (name == "uncertainty" || name == "sigma" || name == "error") {
      mapped.uncertainty = i;
    }
  }
  if (mapped.nuclide < 0 || mapped.quantity < 0) {
    return false;
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
    if (beforeFirstRow && looksLikeHeader(fields)) {
      beforeFirstRow = false;
      mapColumns(fields, columns);
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
  // shape is one flat array of objects with four known keys, and the error messages a
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
    const std::size_t objectEnd = text.find('}', pos);
    if (objectEnd == std::string::npos) {
      failAt(sourceName, lineOf(pos), "unterminated object");
    }
    const std::string object = text.substr(objectStart, objectEnd - objectStart + 1);
    pos = objectEnd + 1;

    const auto field = [&](const char* name, std::string& out) {
      const std::string needle = std::string("\"") + name + "\"";
      const std::size_t at = object.find(needle);
      if (at == std::string::npos) {
        return false;
      }
      std::size_t colon = object.find(':', at + needle.size());
      if (colon == std::string::npos) {
        return false;
      }
      ++colon;
      while (colon < object.size() &&
             std::isspace(static_cast<unsigned char>(object[colon])) != 0) {
        ++colon;
      }
      const bool quoted = colon < object.size() && object[colon] == '"';
      if (quoted) {
        ++colon;
      }
      const std::size_t end = quoted ? object.find('"', colon) : object.find_first_of(",}", colon);
      if (end == std::string::npos) {
        return false;
      }
      out = object.substr(colon, end - colon);
      out = std::string(trim(out));
      return true;
    };

    std::string nuclideText;
    std::string quantityText;
    std::string unitText;
    std::string assayedText;
    std::string uncertaintyText;
    if (!field("nuclide", nuclideText) || !field("quantity", quantityText)) {
      failAt(sourceName, lineOf(objectStart), "an entry needs \"nuclide\" and \"quantity\"");
    }
    field("unit", unitText);
    field("assayed", assayedText);
    field("uncertainty", uncertaintyText);

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
    Quantity used = unit;
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
    Quantity used = unit;
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
