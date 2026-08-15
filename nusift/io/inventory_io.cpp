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

// A first row whose second field is not a number is a header naming the columns. Detecting
// it beats requiring one, since half the spreadsheets in the world have one and half do not.
bool looksLikeHeader(const std::vector<std::string_view>& fields) {
  if (fields.size() < 2) {
    return false;
  }
  double ignored = 0.0;
  return !parseNumber(fields[1], ignored);
}

std::string toLower(std::string_view s) {
  std::string out(s);
  std::transform(out.begin(), out.end(), out.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return out;
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

}  // namespace

Inventory readInventoryCsv(std::istream& in, const NuclearData& data, const std::string& sourceName,
                           const InventoryReadOptions& options) {
  Inventory inventory;
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
      continue;
    }
    beforeFirstRow = false;

    if (fields.size() < 2) {
      failAt(sourceName, lineNumber,
             "expected at least \"nuclide, quantity\", got \"" + std::string(content) + "\"");
    }

    const Zai zai =
        requireNuclideName(fields[0], sourceName + " line " + std::to_string(lineNumber));
    double value = 0.0;
    if (!parseNumber(fields[1], value)) {
      failAt(sourceName, lineNumber, "\"" + std::string(fields[1]) + "\" is not a number");
    }
    if (!(value >= 0.0)) {
      failAt(sourceName, lineNumber, "a quantity cannot be negative");
    }

    // A missing unit column means atoms, which is the only unit that needs no nuclear data
    // and so the only safe default.
    Quantity unit = Quantity::Atoms;
    if (fields.size() >= 3 && !fields[2].empty()) {
      if (!parseQuantity(fields[2], unit)) {
        failAt(sourceName, lineNumber,
               "\"" + std::string(fields[2]) +
                   "\" is not a unit (try atoms, mol, g, kg, mg, Bq, kBq, MBq, GBq, TBq, "
                   "Ci, mCi, uCi)");
      }
    }

    double atoms = 0.0;
    if (rowToAtoms(zai, value, unit, data, sourceName, lineNumber, options, atoms)) {
      inventory.add(zai, atoms);
      ++accepted;
    }
  }

  if (accepted == 0) {
    throw InputError(tagged(kModule, sourceName + " contains no usable inventory rows"));
  }
  inventory.setProvenance(sourceName);
  return inventory;
}

Inventory readInventoryJson(std::istream& in, const NuclearData& data,
                            const std::string& sourceName, const InventoryReadOptions& options) {
  // Deliberately a small hand-rolled reader rather than a JSON dependency: the accepted
  // shape is one flat array of objects with three known keys, and the error messages a
  // purpose-built parser can give ("line 12: ...") are better than a generic one's.
  std::ostringstream buffer;
  buffer << in.rdbuf();
  const std::string text = buffer.str();

  Inventory inventory;
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
    if (!field("nuclide", nuclideText) || !field("quantity", quantityText)) {
      failAt(sourceName, lineOf(objectStart), "an entry needs \"nuclide\" and \"quantity\"");
    }
    field("unit", unitText);

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

    double atoms = 0.0;
    if (rowToAtoms(zai, value, unit, data, sourceName, lineOf(objectStart), options, atoms)) {
      inventory.add(zai, atoms);
      ++accepted;
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
  inventory.setProvenance(sourceName);
  return inventory;
}

Inventory readInventory(const std::string& path, const NuclearData& data,
                        const InventoryReadOptions& options) {
  std::ifstream in(path, std::ios::binary);
  if (!in) {
    throw InputError(tagged(kModule, "cannot open \"" + path + "\""));
  }
  const std::string lower = toLower(path);
  if (lower.size() >= 5 && lower.compare(lower.size() - 5, 5, ".json") == 0) {
    return readInventoryJson(in, data, path, options);
  }
  return readInventoryCsv(in, data, path, options);
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
    char buffer[64];
    std::snprintf(buffer, sizeof(buffer), "%.10g", value);
    out << formatNuclideName(zai) << ',' << buffer << ',' << quantityName(used) << '\n';
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
    char buffer[64];
    std::snprintf(buffer, sizeof(buffer), "%.10g", value);
    out << "  {\"nuclide\": \"" << formatNuclideName(zai) << "\", \"quantity\": " << buffer
        << ", \"unit\": \"" << quantityName(used) << "\"}";
  }
  out << "\n]\n";
}

}  // namespace nusift
