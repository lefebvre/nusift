#include "nusift/nucdata/store_arrays.hpp"

#include <algorithm>
#include <string>

#include "nusift/core/error.hpp"

namespace nusift {
namespace {

constexpr const char* kModule = "nucdata store";

[[noreturn]] void fail(const std::string& what) {
  throw NusiftError(tagged(kModule, what));
}

// Every per-nuclide array is either exactly N long or empty. Empty means "this field was
// never staged", which is a legitimate state for the optional fields (a store built from a
// depletion-chain XML has no photon lines or AWR) and is distinguishable from "staged as
// zero" precisely because the array is absent rather than full of zeros.
void requireNuclideArray(const std::vector<double>& v, int n, const char* name, bool optional) {
  const int size = static_cast<int>(v.size());
  if (size == n) {
    return;
  }
  if (optional && size == 0) {
    return;
  }
  fail(std::string(name) + " has " + std::to_string(size) + " entries, expected " +
       std::to_string(n));
}

// A CSR offset array must be length n+1, start at 0, be non-decreasing, and end exactly at
// the length of the value arrays it indexes. Checking all four is what makes a malformed
// store fail at load with a precise message instead of reading out of bounds later.
void requireCsr(const std::vector<int>& offset, int n, std::size_t valueCount, const char* name) {
  if (static_cast<int>(offset.size()) != n + 1) {
    fail(std::string(name) + " must have " + std::to_string(n + 1) + " entries, has " +
         std::to_string(offset.size()));
  }
  if (offset.front() != 0) {
    fail(std::string(name) + " must start at 0, starts at " + std::to_string(offset.front()));
  }
  for (int i = 0; i < n; ++i) {
    if (offset[i + 1] < offset[i]) {
      fail(std::string(name) + " decreases at index " + std::to_string(i) + " (" +
           std::to_string(offset[i]) + " -> " + std::to_string(offset[i + 1]) + ")");
    }
  }
  if (static_cast<std::size_t>(offset.back()) != valueCount) {
    fail(std::string(name) + " ends at " + std::to_string(offset.back()) + " but indexes " +
         std::to_string(valueCount) + " values");
  }
}

}  // namespace

const char* dataSourceName(DataSource source) {
  switch (source) {
    case DataSource::Endf:
      return "endf";
    case DataSource::OpenmcChainXml:
      return "openmc-chain-xml";
    case DataSource::None:
      break;
  }
  return "none";
}

void validateStoreArrays(const StoreArrays& a) {
  const int n = a.nuclideCount();

  // Half-life is the one genuinely mandatory per-nuclide field: without it a nuclide has no
  // decay constant and cannot participate in the matrix at all.
  requireNuclideArray(a.halfLife, n, "nuclide_half_life", /*optional=*/false);
  requireNuclideArray(a.awr, n, "nuclide_awr", /*optional=*/true);
  requireNuclideArray(a.emEnergyEv, n, "nuclide_em_energy_ev", /*optional=*/true);
  requireNuclideArray(a.lpEnergyEv, n, "nuclide_lp_energy_ev", /*optional=*/true);
  requireNuclideArray(a.hpEnergyEv, n, "nuclide_hp_energy_ev", /*optional=*/true);
  requireNuclideArray(a.continuumPhotonEv, n, "nuclide_continuum_photon_ev", /*optional=*/true);

  // The sort order is a contract, not a convenience: NuclearData binary-searches the key
  // array, and a store whose axis is unsorted would silently fail those lookups. Duplicates
  // are rejected for the same reason -- two rows for one nuclide means one of them is
  // unreachable and the inventory it holds is silently dropped.
  for (int i = 1; i < n; ++i) {
    if (a.nuclideKey[i] <= a.nuclideKey[i - 1]) {
      fail("nuclide_key must be sorted ascending and unique; index " + std::to_string(i) +
           " has key " + std::to_string(a.nuclideKey[i]) + " after " +
           std::to_string(a.nuclideKey[i - 1]));
    }
  }

  requireCsr(a.modeOffset, n, a.modeRtyp.size(), "mode_offset");
  if (a.modeBranching.size() != a.modeRtyp.size() || a.modeFinalState.size() != a.modeRtyp.size() ||
      a.modeIsFission.size() != a.modeRtyp.size()) {
    fail("mode_* arrays disagree in length");
  }

  // Lines are optional wholesale -- a chain-XML store has none -- but if the offsets are
  // present they must be consistent.
  if (!a.lineOffset.empty()) {
    requireCsr(a.lineOffset, n, a.lineEnergyEv.size(), "line_offset");
    if (a.lineIntensity.size() != a.lineEnergyEv.size() ||
        a.lineStyp.size() != a.lineEnergyEv.size()) {
      fail("line_* arrays disagree in length");
    }
  }

  const int sets = a.yieldSetCount();
  if (static_cast<int>(a.nfyEnergyEv.size()) != sets) {
    fail("nfy_energy_ev has " + std::to_string(a.nfyEnergyEv.size()) + " entries, expected " +
         std::to_string(sets));
  }
  if (sets > 0 || !a.nfySetOffset.empty()) {
    requireCsr(a.nfySetOffset, sets, a.nfyProductKey.size(), "nfy_set_offset");
    if (a.nfyProductYield.size() != a.nfyProductKey.size()) {
      fail("nfy_product_* arrays disagree in length");
    }
  }

  const int targets = static_cast<int>(a.xsTargetKey.size());
  if (targets > 0 || !a.xsOffset.empty()) {
    requireCsr(a.xsOffset, targets, a.xsReactionType.size(), "xs_offset");
    if (a.xsProductKey.size() != a.xsReactionType.size() ||
        a.xsSigmaBarn.size() != a.xsReactionType.size()) {
      fail("xs_* arrays disagree in length");
    }
  }
}

}  // namespace nusift
