#pragma once
/**
 * @file
 * @brief A versioned table of per-nuclide coefficients, loaded from a file rather than compiled in.
 * @ingroup nucdata
 */
//
// Every metric NuSIFT reports is a fixed weight vector post-multiplied on a solve. For activity
// the weights are decay constants; for exposure and effective dose they are photon sums through
// a kernel. For most of the rest -- transport indices, intake dose, gross alpha, skin dose,
// transfer factors -- they are simply a published number per nuclide, and there is no reason
// each one should require a new enumerator, a new switch arm, a new unit, a CLI flag, a report
// field and a binding. A pack is that table as DATA: one file, one quantity, loaded at runtime.
//
// The format is CSV with a machine-readable header, chosen so that a pack diffs line by line
// when it is revised. That matters more here than elsewhere: the version is part of the answer,
// and "what changed between editions" has to be a question a reader can put to `git diff`.
//
//     # pack: iaea-ssr6-a2
//     # version: 2012 edition
//     # quantity: A2 sum of fractions
//     # unit: 1
//     # basis: activity
//     # domain: instant
//     # progeny: folded            (or `excluded`)
//     # scenario: slow lung absorption      (optional)
//     # source: ... (free text, may span lines)
//     nuclide,coefficient,folded,note
//     Co-60,2.7027e-15,,
//     Cs-137,1.666667e-12,Ba-137m,
//     Ca-41,0,,unlimited: SSR-6 sets no A2 for this nuclide
//
// THREE RULES, and they are the reason this is a type rather than a std::map.
//
//   A nuclide the pack does not carry is FLAGGED and counted against a coverage figure, never
//   silently weighted zero. "This pack covers 97% of the activity" is the difference between an
//   index a reader can act on and one that is quietly missing its largest term. A zero written
//   IN the file is a different statement -- SSR-6 setting no limit at all -- and carries its
//   reason in the note column so the two can never be confused.
//
//   Progeny handling is DECLARED, and a pack that folds daughters into their parents has to say
//   which ones. Tables written in the "+D" convention are common -- SSR-6's A2 values include
//   every progeny with a half-life under ten days, which is 129 daughters across 77 parents,
//   Ba-137m into Cs-137 and Y-90 into Sr-90 among them -- and against a chain that tracks those
//   daughters explicitly they are a trap in both directions. Weighting the daughter as well as
//   the parent double-counts it; weighting neither, and reporting the daughter as uncovered,
//   understates the coverage figure by exactly the nuclides that usually dominate.
//
//   So `progeny: folded` requires a per-row list of what each coefficient already accounts for.
//   `progeny: excluded` is the other legal value and forbids the column entirely. Isomers are
//   distinct entries throughout, as Ba-137m is distinct from Cs-137.
//
//   WHICH COEFFICIENT APPLIES TO A FOLDED DAUGHTER DEPENDS ON THE INVENTORY, and that is not a
//   complication invented here -- it is what the regulation says. SSR-6 gives Y-90 an A2 of its
//   own AND folds Y-90 into Sr-90's; 36 of its 129 folded daughters are like this. Yttrium-90
//   shipped alone is limited by its own value, and yttrium-90 accompanying strontium-90 is
//   already inside its parent's. So the fold is resolved against the SEED: a daughter whose
//   parent is in the inventory takes no weight of its own, and one seeded alone keeps its row.
//
//   Resolving it against the seed rather than against the inventory at each time is deliberate
//   and is what keeps a metric a FIXED weight vector -- the whole architectural bet -- rather
//   than one that moves down the time axis. It is also the physical answer: a folded daughter
//   is short-lived by construction (SSR-6's cut is ten days), so it exists at a later time only
//   because its parent is feeding it, and a seed with the parent stays a seed with the parent.
//
//   The version and the scenario travel with every answer. A1 and A2 changed between SSR-6
//   editions; inhalation coefficients depend on absorption type and particle size. A number
//   from this file is meaningless without knowing which edition and which scenario produced it,
//   so the provenance is carried on the pack, into the response table, and into the report.
//
#include <cstdint>
#include <iosfwd>
#include <string>
#include <unordered_map>
#include <vector>

#include "nusift/engine/inventory.hpp"

namespace nusift {

class NuclearData;

// What the coefficient multiplies. A published table states this, usually in its units: an A2
// fraction is per becquerel, an intake dose coefficient is per becquerel, a gross alpha weight
// is per decay, a nCi/g index is per gram. Getting it wrong scales every answer by a decay
// constant, so it is declared rather than inferred.
enum class PackBasis {
  Activity,  // coefficient x activity   [per Bq]
  Atoms,     // coefficient x atoms      [per atom]
  Mass,      // coefficient x mass       [per gram]
  // coefficient x activity CONCENTRATION [per Bq/m2, Bq/m3 or Bq/kg]
  //
  // The basis a distributed source needs, and the one an Inventory cannot supply on its own: it
  // carries atoms, and a concentration is atoms divided by the extent they are spread through.
  // That extent is a fact about the situation rather than about the material, so the caller
  // states it at resolution and `per` says which unit it must be in.
  Concentration,
};

const char* packBasisName(PackBasis basis);

// Which domains the quantity is meaningful in. A rate-like coefficient against an interval
// integral is a different quantity, and one a pack has no business implying.
enum class PackDomains {
  InstantOnly,
  IntervalOnly,
  Both,
};

// How a nuclide is accounted for by a pack. The three cases are genuinely different answers to
// "is this nuclide covered", and collapsing them into a bool is what makes a coverage figure
// lie: a folded daughter is covered, by its parent's row, and has no coefficient of its own.
enum class PackCoverage : char {
  None,    // the pack does not carry this nuclide at all
  Own,     // it has its own row
  Folded,  // its contribution is inside a parent's coefficient
};

struct PackProvenance {
  std::string name;
  std::string version;
  std::string quantity;
  std::string unit;      // printed as the column header; "1" for a dimensionless index
  std::string scenario;  // may be empty; part of the metric's identity when it is not
  std::string source;
  std::string path;  // where it was read from, for a report that has to say
  // Concentration packs only: the denominator of the concentration, "m2", "m3" or "kg". Ground
  // deposition is per square metre and a cloud is per cubic metre, and confusing the two is not
  // a units slip but a different question answered.
  std::string per;
  PackBasis basis = PackBasis::Activity;
  PackDomains domains = PackDomains::Both;
  // True when coefficients absorb their progeny, in which case foldedInto() answers for the
  // daughters. Reported alongside the version, because it changes what the number means.
  bool foldsProgeny = false;
};

// The extent an inventory is spread through, for a concentration pack. Zero means none was
// given, which is what every other basis requires.
struct PackExtent {
  double value = 0.0;
  std::string unit;  // must match the pack's `per`
};

class CoefficientPack {
public:
  static CoefficientPack read(std::istream& in, const std::string& sourceName);
  static CoefficientPack open(const std::string& path);

  const PackProvenance& provenance() const { return provenance_; }
  int size() const { return static_cast<int>(coefficients_.size()); }

  // Whether the pack carries this nuclide at all -- which is NOT the same question as whether
  // its coefficient is zero, and the two must never be answered by one call.
  bool covers(std::int64_t zaiKey) const;

  // The coefficient, or zero for a nuclide the pack does not carry. Always ask covers() first
  // unless a zero is genuinely what you want for both cases.
  double coefficient(std::int64_t zaiKey) const;

  // The note against a row, empty when it has none. Carried because a zero that means "no limit
  // applies" has to be able to say so.
  const std::string& note(std::int64_t zaiKey) const;

  // The parents whose coefficients already account for this nuclide, empty when none does.
  // This is what lets a report say "Ba-137m is inside Cs-137's limit" rather than "Ba-137m is
  // not covered", which is the difference between a true coverage figure and a misleading one.
  //
  // SEVERAL parents, because decay chains nest and each parent's coefficient covers everything
  // below it: SSR-6 folds Tl-208 into Bi-212, and Bi-212 with Tl-208 into Pb-212, and all of
  // them into Ra-224 and Th-228. Any one of them being seeded is enough to account for it.
  const std::vector<std::int64_t>& foldedInto(std::int64_t zaiKey) const;

  // How this nuclide is accounted for against a seed: its own row, a parent's, or not at all.
  PackCoverage coverageOf(std::int64_t zaiKey, const Inventory& seed) const;

  // The per-nuclide weight vector over `data`'s whole index space, already carrying the basis:
  // lambda for an activity-basis pack, one for an atom-basis pack, the molar mass over
  // Avogadro's number for a mass-basis pack. This is what the response layer post-multiplies,
  // and it is the same shape responseWeights() returns for a built-in metric.
  //
  // `seed` decides the folds: a daughter whose parent is seeded is accounted for by that
  // parent's row and is weighted zero here, while one seeded alone keeps its own coefficient.
  // An empty inventory weights every row on its own terms, which is what a pack with no folds
  // does in any case.
  //
  // `extent` divides, for a concentration pack: the coefficient is per unit concentration, so
  // the weight is the coefficient over the area, volume or mass the inventory occupies. It is
  // required for such a pack and refused for every other, because a number per becquerel does
  // not become a different number when told how much space the becquerels are in.
  std::vector<double> weights(const NuclearData& data, const Inventory& seed,
                              const PackExtent& extent = {}) const;

  // How each of `data`'s nuclides is accounted for, parallel to weights(). A response uses this
  // to report coverage and to flag the contributors that fall outside the pack.
  std::vector<PackCoverage> covered(const NuclearData& data, const Inventory& seed) const;

private:
  PackProvenance provenance_;
  std::unordered_map<std::int64_t, double> coefficients_;
  std::unordered_map<std::int64_t, std::string> notes_;
  // Daughter -> the parents that fold it, built once at load. The reverse direction is what
  // every question at runtime actually asks: given this nuclide, is anything already
  // accounting for it.
  std::unordered_map<std::int64_t, std::vector<std::int64_t>> foldedInto_;
};

}  // namespace nusift
