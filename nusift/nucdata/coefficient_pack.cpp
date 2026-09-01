#include "nusift/nucdata/coefficient_pack.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <istream>
#include <optional>
#include <sstream>
#include <stdexcept>

#include "nusift/core/error.hpp"
#include "nusift/core/nuclide.hpp"
#include "nusift/core/nuclide_name.hpp"
#include "nusift/nucdata/nuclear_data.hpp"
#include "nusift/units.hpp"

namespace nusift {
namespace {

constexpr const char* kModule = "pack";

std::string trimmed(std::string_view text) {
  const std::size_t first = text.find_first_not_of(" \t\r\n");
  if (first == std::string_view::npos) {
    return {};
  }
  const std::size_t last = text.find_last_not_of(" \t\r\n");
  return std::string(text.substr(first, last - first + 1));
}

std::string lowered(std::string text) {
  std::transform(text.begin(), text.end(), text.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return text;
}

PackBasis basisFrom(const std::string& text, const std::string& where) {
  const std::string value = lowered(text);
  if (value == "activity") {
    return PackBasis::Activity;
  }
  if (value == "atoms") {
    return PackBasis::Atoms;
  }
  if (value == "mass") {
    return PackBasis::Mass;
  }
  throw InputError(tagged(kModule, where + ": basis \"" + text +
                                       "\" is not one of activity, atoms, mass. It says what "
                                       "the coefficient multiplies, and guessing it would "
                                       "scale every answer by a decay constant"));
}

PackDomains domainsFrom(const std::string& text, const std::string& where) {
  const std::string value = lowered(text);
  if (value == "instant") {
    return PackDomains::InstantOnly;
  }
  if (value == "interval") {
    return PackDomains::IntervalOnly;
  }
  if (value == "both") {
    return PackDomains::Both;
  }
  throw InputError(
      tagged(kModule, where + ": domain \"" + text + "\" is not one of instant, interval, both"));
}

// A header line is `# key: value`, and a line whose value continues is appended to the last key
// -- which is how a citation runs to several lines without becoming several fields.
bool splitHeader(const std::string& line, std::string& key, std::string& value) {
  const std::string body = trimmed(std::string_view(line).substr(1));
  const std::size_t colon = body.find(':');
  if (colon == std::string::npos) {
    return false;
  }
  key = lowered(trimmed(std::string_view(body).substr(0, colon)));
  value = trimmed(std::string_view(body).substr(colon + 1));
  // A key is one word: `# source: ... http://x` has a colon inside its value, and the text
  // before the first colon there is a sentence rather than a field name.
  return !key.empty() && key.find(' ') == std::string::npos;
}

}  // namespace

const char* packBasisName(PackBasis basis) {
  switch (basis) {
    case PackBasis::Activity:
      return "activity";
    case PackBasis::Atoms:
      return "atoms";
    case PackBasis::Mass:
      return "mass";
  }
  return "?";
}

CoefficientPack CoefficientPack::read(std::istream& in, const std::string& sourceName) {
  CoefficientPack pack;
  pack.provenance_.path = sourceName;

  std::unordered_map<std::string, std::string> fields;
  std::string lastKey;
  std::string line;
  int lineNumber = 0;
  bool sawColumns = false;
  int foldedColumn = -1;
  int noteColumn = -1;
  int rows = 0;

  while (std::getline(in, line)) {
    ++lineNumber;
    const std::string text = trimmed(line);
    if (text.empty()) {
      continue;
    }
    const std::string where = sourceName + ":" + std::to_string(lineNumber);

    if (text[0] == '#') {
      std::string key;
      std::string value;
      if (splitHeader(text, key, value)) {
        fields[key] = value;
        lastKey = key;
      } else if (!lastKey.empty()) {
        // A continuation of the field above, joined with a space so a citation reads as one.
        fields[lastKey] += " " + trimmed(std::string_view(text).substr(1));
      }
      continue;
    }

    if (!sawColumns) {
      // The column line, checked rather than skipped: a pack whose columns are in another order
      // would load silently and weight everything by the wrong number.
      const std::string columns = lowered(text);
      if (columns.rfind("nuclide,coefficient", 0) != 0) {
        throw InputError(tagged(kModule, where +
                                             ": expected the column line "
                                             "`nuclide,coefficient[,folded][,note]`, got \"" +
                                             text + "\""));
      }
      // Which optional columns are present, in order, so a row is read against the header it
      // declared rather than by counting commas and hoping.
      foldedColumn = columns.find(",folded") != std::string::npos ? 2 : -1;
      noteColumn = columns.find(",note") != std::string::npos ? (foldedColumn > 0 ? 3 : 2) : -1;
      sawColumns = true;
      continue;
    }

    std::vector<std::string> cells;
    std::size_t at = 0;
    while (at <= text.size()) {
      const std::size_t comma = text.find(',', at);
      cells.push_back(trimmed(comma == std::string::npos
                                  ? std::string_view(text).substr(at)
                                  : std::string_view(text).substr(at, comma - at)));
      if (comma == std::string::npos) {
        break;
      }
      at = comma + 1;
    }
    if (cells.size() < 2) {
      throw InputError(tagged(kModule, where + ": a row is `nuclide,coefficient[,folded][,note]`"));
    }
    const std::string name = cells[0];
    const std::string number = cells[1];
    const std::string folded =
        foldedColumn > 0 && cells.size() > static_cast<std::size_t>(foldedColumn)
            ? cells[static_cast<std::size_t>(foldedColumn)]
            : std::string();
    const std::string note = noteColumn > 0 && cells.size() > static_cast<std::size_t>(noteColumn)
                                 ? cells[static_cast<std::size_t>(noteColumn)]
                                 : std::string();

    const std::optional<Zai> zai = parseNuclideName(name);
    if (!zai) {
      throw InputError(tagged(kModule, where + ": \"" + name + "\" is not a nuclide name"));
    }
    double value = 0.0;
    try {
      std::size_t consumed = 0;
      value = std::stod(number, &consumed);
      if (consumed != number.size()) {
        throw std::invalid_argument("trailing");
      }
    } catch (const std::exception&) {
      throw InputError(tagged(kModule, where + ": \"" + number + "\" is not a number"));
    }
    if (!(value >= 0.0)) {
      throw InputError(tagged(kModule, where + ": a coefficient cannot be negative"));
    }
    if (!pack.coefficients_.emplace(zai->key(), value).second) {
      throw InputError(tagged(kModule, where + ": " + name +
                                           " appears twice. Two coefficients for one nuclide is "
                                           "not a stated quantity, and choosing between them "
                                           "would be a guess"));
    }
    if (!note.empty()) {
      pack.notes_[zai->key()] = note;
    }
    if (!folded.empty()) {
      std::istringstream daughters(folded);
      std::string daughter;
      while (daughters >> daughter) {
        const std::optional<Zai> child = parseNuclideName(daughter);
        if (!child) {
          throw InputError(tagged(kModule, where + ": \"" + daughter + "\" is not a nuclide name"));
        }
        // Several parents may fold the same daughter, and that is chains nesting rather than
        // a contradiction: each parent's coefficient covers everything below it, so any one of
        // them being present accounts for the daughter.
        std::vector<std::int64_t>& parents = pack.foldedInto_[child->key()];
        if (std::find(parents.begin(), parents.end(), zai->key()) == parents.end()) {
          parents.push_back(zai->key());
        }
      }
    }
    ++rows;
  }

  const auto required = [&](const char* key) {
    const auto it = fields.find(key);
    if (it == fields.end() || it->second.empty()) {
      throw InputError(tagged(kModule, sourceName + ": the header has no `" + key +
                                           "`. Every field of it is part of what an answer from "
                                           "this pack means"));
    }
    return it->second;
  };

  pack.provenance_.name = required("pack");
  pack.provenance_.version = required("version");
  pack.provenance_.quantity = required("quantity");
  pack.provenance_.unit = required("unit");
  pack.provenance_.source = required("source");
  pack.provenance_.basis = basisFrom(required("basis"), sourceName);
  pack.provenance_.domains = domainsFrom(required("domain"), sourceName);

  // Progeny handling is not a preference, it is whether the numbers can be used at all: a table
  // that folds daughters into the parent double-counts against a chain that tracks them.
  const std::string progeny = lowered(required("progeny"));
  if (progeny == "excluded") {
    if (!pack.foldedInto_.empty()) {
      throw InputError(tagged(
          kModule, sourceName + ": progeny is `excluded` but rows name folded daughters. One "
                                "of the two is wrong, and which cannot be guessed"));
    }
  } else if (progeny == "folded") {
    pack.provenance_.foldsProgeny = true;
    if (pack.foldedInto_.empty()) {
      throw InputError(
          tagged(kModule, sourceName + ": progeny is `folded` but no row says what it folds. The "
                                       "list is what makes the declaration usable rather than a "
                                       "warning"));
    }
    // A daughter holding its own row as well as a parent's fold is NOT an error, and refusing
    // it here was the first thing this loader got wrong. SSR-6 does it 36 times deliberately:
    // Y-90 has an A2 for being shipped alone and is also inside Sr-90's. Which one applies is a
    // question about an inventory, so it is answered in weights() and not here.
  } else {
    throw InputError(
        tagged(kModule, sourceName + ": progeny is \"" + progeny +
                            "\", and only `excluded` or `folded` can be used. Whether a "
                            "coefficient already accounts for its daughters decides "
                            "whether a chain that tracks them double-counts"));
  }

  const auto scenario = fields.find("scenario");
  if (scenario != fields.end()) {
    pack.provenance_.scenario = scenario->second;
  }

  if (rows == 0) {
    throw InputError(tagged(kModule, sourceName + ": no coefficients"));
  }
  return pack;
}

CoefficientPack CoefficientPack::open(const std::string& path) {
  std::ifstream in(path);
  if (!in) {
    throw InputError(tagged(kModule, "cannot read \"" + path + "\""));
  }
  return read(in, path);
}

bool CoefficientPack::covers(std::int64_t zaiKey) const {
  return coefficients_.find(zaiKey) != coefficients_.end();
}

double CoefficientPack::coefficient(std::int64_t zaiKey) const {
  const auto it = coefficients_.find(zaiKey);
  return it == coefficients_.end() ? 0.0 : it->second;
}

const std::string& CoefficientPack::note(std::int64_t zaiKey) const {
  static const std::string kNone;
  const auto it = notes_.find(zaiKey);
  return it == notes_.end() ? kNone : it->second;
}

namespace {

// Whether the parent that folds `key`, if any, is present in the seed. Presence is judged on
// the seed rather than on the decayed inventory so that the weight vector stays fixed down the
// time axis, which is the property every metric in NuSIFT is built on.
bool parentIsSeeded(const std::unordered_map<std::int64_t, std::vector<std::int64_t>>& foldedInto,
                    std::int64_t key, const Inventory& seed) {
  const auto it = foldedInto.find(key);
  if (it == foldedInto.end()) {
    return false;
  }
  for (const std::int64_t parent : it->second) {
    if (seed.atomsOf(Zai::fromKey(parent)) > 0.0) {
      return true;
    }
  }
  return false;
}

}  // namespace

std::vector<double> CoefficientPack::weights(const NuclearData& data, const Inventory& seed) const {
  std::vector<double> weights(static_cast<std::size_t>(data.size()), 0.0);
  for (int i = 0; i < data.size(); ++i) {
    const std::int64_t key = data.zaiAt(i).key();
    const auto it = coefficients_.find(key);
    if (it == coefficients_.end()) {
      continue;
    }
    // Its parent is here, so this nuclide's contribution is already inside that parent's
    // coefficient. Weighting it again is the double count the fold list exists to prevent.
    if (parentIsSeeded(foldedInto_, key, seed)) {
      continue;
    }
    switch (provenance_.basis) {
      case PackBasis::Activity:
        // Per becquerel against atoms is per becquerel times lambda.
        weights[static_cast<std::size_t>(i)] = it->second * data.decayConstant(i);
        break;
      case PackBasis::Atoms:
        weights[static_cast<std::size_t>(i)] = it->second;
        break;
      case PackBasis::Mass: {
        const double molarMass = data.molarMassGPerMol(i);
        if (!(molarMass > 0.0)) {
          throw InputError(tagged(kModule,
                                  "this pack is per gram, and the store has no atomic weight "
                                  "for " +
                                      formatNuclideName(data.zaiAt(i)) +
                                      ", so its mass cannot be formed"));
        }
        weights[static_cast<std::size_t>(i)] = it->second * molarMass / units::kAvogadro;
        break;
      }
    }
  }
  return weights;
}

const std::vector<std::int64_t>& CoefficientPack::foldedInto(std::int64_t zaiKey) const {
  static const std::vector<std::int64_t> kNone;
  const auto it = foldedInto_.find(zaiKey);
  return it == foldedInto_.end() ? kNone : it->second;
}

PackCoverage CoefficientPack::coverageOf(std::int64_t zaiKey, const Inventory& seed) const {
  // Folded first: a daughter whose parent is present is accounted for by that parent even when
  // it has a row of its own, and reporting it as Own would say the wrong thing about which
  // number was used.
  if (parentIsSeeded(foldedInto_, zaiKey, seed)) {
    return PackCoverage::Folded;
  }
  return covers(zaiKey) ? PackCoverage::Own : PackCoverage::None;
}

std::vector<PackCoverage> CoefficientPack::covered(const NuclearData& data,
                                                   const Inventory& seed) const {
  std::vector<PackCoverage> covered(static_cast<std::size_t>(data.size()), PackCoverage::None);
  for (int i = 0; i < data.size(); ++i) {
    const std::int64_t key = data.zaiAt(i).key();
    covered[static_cast<std::size_t>(i)] = coverageOf(key, seed);
  }
  return covered;
}

}  // namespace nusift
