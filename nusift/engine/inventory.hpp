#pragma once
/**
 * @file
 * @brief An isotopic inventory in atoms, and the unit conversions that produce one.
 * @ingroup engine
 */
//
// The inventory is what NuSIFT is handed and what it decays forward. Internally it is always
// atoms: every input unit converts on the way in, so nothing downstream ever asks what units
// a number is in. Conversions that cannot be performed -- grams without a staged atomic
// weight, becquerel for a stable nuclide -- throw naming the nuclide, rather than producing
// a zero or an infinity that would silently distort every ranking that follows.
//
#include <cstdint>
#include <span>
#include <string>
#include <vector>

#include "nusift/core/nuclide.hpp"

namespace nusift {

class NuclearData;

// Units an inventory quantity may be written in. Mass and activity both appear in real
// inventories -- a spreadsheet from a burnup calculation is usually in grams, while a
// measured source is in curies -- so both convert here rather than at each call site.
enum class Quantity {
  Atoms,
  Moles,
  Grams,
  Kilograms,
  Milligrams,
  Becquerel,
  Kilobecquerel,
  Megabecquerel,
  Gigabecquerel,
  Terabecquerel,
  Curie,
  Millicurie,
  Microcurie,
};

// Parse a unit spelling from an inventory file or CLI argument, case-insensitively:
// "atoms", "g", "kg", "mg", "mol", "bq", "kbq", "mbq", "gbq", "tbq", "ci", "mci", "uci".
// Returns false if unrecognised, leaving `out` untouched.
bool parseQuantity(std::string_view text, Quantity& out);

// Canonical spelling, for round-tripping an inventory back out.
const char* quantityName(Quantity quantity);

// Whether a quantity is an ACTIVITY -- becquerel or curie and their multiples. Named because
// what it decides is not a conversion but a derivative: a row given this way fixes A0 rather
// than n0, so its atom count moves when a decay constant does.
bool isActivityQuantity(Quantity quantity);

struct InventoryEntry {
  std::int64_t zaiKey = 0;
  double atoms = 0.0;

  // 1-sigma uncertainty on `atoms`, or zero when the row stated none. Zero and "known exactly"
  // are deliberately the same encoding here, because no inventory row is ever known exactly and
  // a reader who supplied no sigma has said nothing rather than said zero -- what consumes this
  // reports how much of the answer rests on rows carrying no sigma, which is the honest form of
  // the distinction.
  //
  // MEANINGFUL ONLY ON AN AS-MEASURED INVENTORY. Carrying an assay forward turns a diagonal
  // covariance into D Sigma D^T, and D is not diagonal, so a per-row sigma on a RECONCILED
  // inventory would be a lie -- reconcile() drops it rather than propagating it, and
  // triage/uncertainty.hpp reaches back to the assays instead.
  double sigmaAtoms = 0.0;

  // How many of `atoms` came from a row written as an ACTIVITY. Zero for a row given in atoms,
  // moles or a mass unit, and equal to `atoms` for one given in becquerel or curie.
  //
  // Carried because a derivative with respect to a decay constant needs it and nothing else
  // does. A row measured in atoms fixes n0, and lambda may move without it; a row measured in
  // becquerel fixes A0, so n0 = A0/lambda MOVES with lambda and dR/dlambda gains a third term
  // dR/dn0 * (-n0/lambda). Held as an amount rather than as a basis enum so that two rows of one
  // nuclide given in different units -- a drum weighed and a source counted -- carry the right
  // split instead of one overwriting the other's basis.
  double atomsFromActivity = 0.0;
};

// A set of nuclides and their atom counts, merged by nuclide and kept sorted by key so that
// two inventories built in different orders compare and serialize identically.
class Inventory {
public:
  // Adds to any existing entry for this nuclide rather than replacing it, so a file listing
  // a nuclide twice accumulates instead of silently keeping only the last row.
  //
  // Uncertainties on two rows of one nuclide combine IN QUADRATURE, which is the reading that
  // matches the accumulation above: two rows are two measurements -- two containers, two
  // aliquots -- and two independent measurements add in variance. Two rows that are the same
  // measurement written twice would not, but a file cannot say that and neither can this.
  //
  // Throws InputError on a negative, infinite, or NaN count: an atom count is non-negative
  // and finite by definition, and none of those survive as anything but a corrupted ranking.
  // A sigma is held to the same standard and to no more than that -- in particular one LARGER
  // than its own quantity is accepted, because that is what a measurement near a detection
  // limit honestly reports.
  void add(const Zai& zai, double atoms, double sigmaAtoms = 0.0, double atomsFromActivity = 0.0);
  void addKey(std::int64_t zaiKey, double atoms, double sigmaAtoms = 0.0,
              double atomsFromActivity = 0.0);

  // Whether any row carries an uncertainty. What decides if an error bar can be offered at all.
  bool hasUncertainties() const;

  std::span<const InventoryEntry> entries() const { return entries_; }
  bool empty() const { return entries_.empty(); }
  int size() const { return static_cast<int>(entries_.size()); }

  // Atoms of one nuclide, or 0 if absent.
  double atomsOf(const Zai& zai) const;

  double totalAtoms() const;

  // How this inventory came to be -- a file path, or a description of the seeding that
  // produced it. Reported in every output header.
  const std::string& provenance() const { return provenance_; }
  void setProvenance(std::string text) { provenance_ = std::move(text); }

private:
  std::vector<InventoryEntry> entries_;  // sorted by zaiKey, unique
  std::string provenance_;
};

// Convert `value` in `quantity` to atoms of the nuclide at `index` in `data`.
//
// Throws InputError naming the nuclide when the conversion is impossible:
//   - a mass unit when the store carries no atomic weight ratio for it
//   - an activity unit for a stable nuclide, where N = A/lambda is undefined
double toAtoms(double value, Quantity quantity, int index, const NuclearData& data);

// The inverse, for reporting an inventory back in the units it was written in. Returns 0
// rather than throwing for an impossible conversion, since this is used in output paths
// where a missing value is better than an aborted report.
double fromAtoms(double atoms, Quantity quantity, int index, const NuclearData& data);

}  // namespace nusift
