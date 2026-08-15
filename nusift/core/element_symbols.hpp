#pragma once
/**
 * @file
 * @brief IUPAC element-symbol lookup by atomic number Z (1..118), and its inverse.
 * @ingroup core
 */
//
// cram exposes (Z, A, I) but no chemical symbol, so every human-readable nuclide label
// NuSIFT prints goes through this table. The reverse lookup is here because NuSIFT parses
// nuclide names out of user-authored inventory CSVs, so it has to go from "Cs" back to
// Z=55 as well as the other way.
//
// Header-only, no link dependency. cram grows its own elementSymbol() on the burnup branch;
// the signature here is deliberately identical so that adopting it is a deletion rather than
// a refactor.
//
#include <cstddef>
#include <string_view>

namespace nusift {
namespace detail {

// Index 0 is a sentinel so the array can be indexed directly by Z.
inline constexpr const char* kElementSymbols[] = {
    "?",                                                         // 0 (sentinel)
    "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne",  // 1-10
    "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca",  // 11-20
    "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",  // 21-30
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr",  // 31-40
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",  // 41-50
    "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",  // 51-60
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",  // 61-70
    "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg",  // 71-80
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",  // 81-90
    "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",  // 91-100
    "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",  // 101-110
    "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",              // 111-118
};

}  // namespace detail

inline constexpr int kMaxAtomicNumber = 118;

// Symbol for atomic number `z` (1..118). Returns "?" out of range rather than throwing:
// this is called from label formatting, where a malformed nuclide should produce a visibly
// wrong string in the output instead of aborting a long run.
inline constexpr const char* elementSymbol(int z) {
  return (z >= 1 && z <= kMaxAtomicNumber) ? detail::kElementSymbols[z] : "?";
}

// Atomic number for an element symbol, case-insensitively ("Cs", "cs", "CS" all give 55).
// Returns 0 if the symbol is not an element, which is the caller's cue to raise an
// InputError naming the offending token.
inline constexpr int atomicNumber(std::string_view symbol) {
  if (symbol.empty() || symbol.size() > 2) {
    return 0;
  }
  const auto lower = [](char c) constexpr {
    return (c >= 'A' && c <= 'Z') ? static_cast<char>(c - 'A' + 'a') : c;
  };
  for (int z = 1; z <= kMaxAtomicNumber; ++z) {
    const std::string_view candidate(detail::kElementSymbols[z]);
    if (candidate.size() != symbol.size()) {
      continue;
    }
    bool match = true;
    for (std::size_t k = 0; k < symbol.size(); ++k) {
      if (lower(candidate[k]) != lower(symbol[k])) {
        match = false;
        break;
      }
    }
    if (match) {
      return z;
    }
  }
  return 0;
}

}  // namespace nusift
