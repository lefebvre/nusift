#include "nusift/engine/inventory.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <string>

#include "nusift/core/error.hpp"
#include "nusift/core/nuclide_name.hpp"
#include "nusift/nucdata/nuclear_data.hpp"
#include "nusift/units.hpp"

namespace nusift {
namespace {

constexpr const char* kModule = "inventory";

std::string lower(std::string_view s) {
  std::string out(s);
  std::transform(out.begin(), out.end(), out.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return out;
}

// Multiplier from a mass unit to grams; 0 for anything that is not a mass.
double gramsPer(Quantity q) {
  switch (q) {
    case Quantity::Grams:
      return 1.0;
    case Quantity::Kilograms:
      return 1000.0;
    case Quantity::Milligrams:
      return 0.001;
    default:
      return 0.0;
  }
}

// Multiplier from an activity unit to becquerel; 0 for anything that is not an activity.
double becquerelPer(Quantity q) {
  switch (q) {
    case Quantity::Becquerel:
      return 1.0;
    case Quantity::Kilobecquerel:
      return 1.0e3;
    case Quantity::Megabecquerel:
      return 1.0e6;
    case Quantity::Gigabecquerel:
      return 1.0e9;
    case Quantity::Terabecquerel:
      return 1.0e12;
    case Quantity::Curie:
      return units::kBqPerCi;
    case Quantity::Millicurie:
      return units::kBqPerCi * 1.0e-3;
    case Quantity::Microcurie:
      return units::kBqPerCi * 1.0e-6;
    default:
      return 0.0;
  }
}

}  // namespace

bool parseQuantity(std::string_view text, Quantity& out) {
  const std::string key = lower(text);
  struct Entry {
    const char* name;
    Quantity quantity;
  };
  static constexpr Entry kTable[] = {
      {"atoms", Quantity::Atoms},       {"atom", Quantity::Atoms},
      {"mol", Quantity::Moles},         {"mole", Quantity::Moles},
      {"moles", Quantity::Moles},       {"g", Quantity::Grams},
      {"gram", Quantity::Grams},        {"grams", Quantity::Grams},
      {"kg", Quantity::Kilograms},      {"mg", Quantity::Milligrams},
      {"bq", Quantity::Becquerel},      {"kbq", Quantity::Kilobecquerel},
      {"mbq", Quantity::Megabecquerel}, {"gbq", Quantity::Gigabecquerel},
      {"tbq", Quantity::Terabecquerel}, {"ci", Quantity::Curie},
      {"mci", Quantity::Millicurie},    {"uci", Quantity::Microcurie},
      {"µci", Quantity::Microcurie},
  };
  for (const Entry& entry : kTable) {
    if (key == entry.name) {
      out = entry.quantity;
      return true;
    }
  }
  return false;
}

const char* quantityName(Quantity quantity) {
  switch (quantity) {
    case Quantity::Atoms:
      return "atoms";
    case Quantity::Moles:
      return "mol";
    case Quantity::Grams:
      return "g";
    case Quantity::Kilograms:
      return "kg";
    case Quantity::Milligrams:
      return "mg";
    case Quantity::Becquerel:
      return "Bq";
    case Quantity::Kilobecquerel:
      return "kBq";
    case Quantity::Megabecquerel:
      return "MBq";
    case Quantity::Gigabecquerel:
      return "GBq";
    case Quantity::Terabecquerel:
      return "TBq";
    case Quantity::Curie:
      return "Ci";
    case Quantity::Millicurie:
      return "mCi";
    case Quantity::Microcurie:
      return "uCi";
  }
  return "?";
}

void Inventory::addKey(std::int64_t zaiKey, double atoms) {
  // A count of atoms is non-negative and finite by definition. Refused here rather than
  // downstream because the failure a bad count causes is silent: a negative seed decays into
  // negative activities that rank as the smallest contributors and vanish off the bottom of
  // every table, and a NaN or infinity propagates through CRAM into a column of NaNs whose
  // origin is no longer visible.
  if (!std::isfinite(atoms) || atoms < 0.0) {
    throw InputError(tagged(kModule, "atoms of " + formatNuclideName(Zai::fromKey(zaiKey)) +
                                         " must be non-negative and finite; got " +
                                         std::to_string(atoms)));
  }
  const auto it = std::lower_bound(
      entries_.begin(), entries_.end(), zaiKey,
      [](const InventoryEntry& entry, std::int64_t key) { return entry.zaiKey < key; });
  if (it != entries_.end() && it->zaiKey == zaiKey) {
    it->atoms += atoms;
    return;
  }
  entries_.insert(it, InventoryEntry{zaiKey, atoms});
}

void Inventory::add(const Zai& zai, double atoms) {
  addKey(zai.key(), atoms);
}

double Inventory::atomsOf(const Zai& zai) const {
  const std::int64_t key = zai.key();
  const auto it = std::lower_bound(
      entries_.begin(), entries_.end(), key,
      [](const InventoryEntry& entry, std::int64_t k) { return entry.zaiKey < k; });
  return (it != entries_.end() && it->zaiKey == key) ? it->atoms : 0.0;
}

double Inventory::totalAtoms() const {
  double total = 0.0;
  for (const InventoryEntry& entry : entries_) {
    total += entry.atoms;
  }
  return total;
}

double toAtoms(double value, Quantity quantity, int index, const NuclearData& data) {
  if (quantity == Quantity::Atoms) {
    return value;
  }
  if (quantity == Quantity::Moles) {
    return value * units::kAvogadro;
  }

  const Zai zai = data.zaiAt(index);

  if (const double grams = gramsPer(quantity); grams > 0.0) {
    const double molar = data.molarMassGPerMol(index);
    if (molar <= 0.0) {
      // Deliberately not falling back on A as the molar mass. That approximation is wrong by
      // ~0.1% for mid-A nuclides and worse for light ones, and silently applying it would put
      // an unattributable error into every mass-specified inventory. A store without atomic
      // weights genuinely cannot answer this.
      throw InputError(tagged(kModule, "cannot convert mass for " + formatNuclideName(zai) +
                                           ": the data store carries no atomic weight ratio "
                                           "(stage from ENDF, or give this nuclide in atoms, "
                                           "moles, or an activity unit)"));
    }
    return value * grams / molar * units::kAvogadro;
  }

  if (const double bq = becquerelPer(quantity); bq > 0.0) {
    const double lambda = data.decayConstant(index);
    if (lambda <= 0.0) {
      throw InputError(tagged(kModule, formatNuclideName(zai) +
                                           " is stable, so an activity cannot be converted to "
                                           "an atom count (N = A/lambda is undefined); give it "
                                           "in atoms, moles, or a mass unit"));
    }
    return value * bq / lambda;
  }

  throw InputError(tagged(kModule, "unhandled quantity for " + formatNuclideName(zai)));
}

double fromAtoms(double atoms, Quantity quantity, int index, const NuclearData& data) {
  if (quantity == Quantity::Atoms) {
    return atoms;
  }
  if (quantity == Quantity::Moles) {
    return atoms / units::kAvogadro;
  }
  if (const double grams = gramsPer(quantity); grams > 0.0) {
    const double molar = data.molarMassGPerMol(index);
    return molar > 0.0 ? atoms / units::kAvogadro * molar / grams : 0.0;
  }
  if (const double bq = becquerelPer(quantity); bq > 0.0) {
    return atoms * data.decayConstant(index) / bq;
  }
  return 0.0;
}

}  // namespace nusift
