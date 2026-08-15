#pragma once
/**
 * @file
 * @brief Cross-PIMPL access to the cram depletion chain. NOT installed.
 * @ingroup nucdata
 */
//
// This header exposes cram (and through it Eigen) types, so it is deliberately excluded from
// the install rules and must be included only by translation units that genuinely assemble or
// solve a matrix -- today nucdata/nuclear_data.cpp and engine/decay_engine.cpp. A CI check
// enforces that; if the list grows, the growth should be a decision rather than a drift.
//
#include "cram/chain.hpp"
#include "nusift/core/nuclide.hpp"
#include "nusift/nucdata/nuclear_data.hpp"

namespace nusift {

// The conversion between NuSIFT's public Zai and cram's. They are layout-compatible, but
// converting explicitly rather than reinterpreting keeps the two types free to diverge and
// costs nothing at any call site that matters.
inline cram::Zai toCram(const Zai& zai) {
  return cram::Zai{zai.z, zai.a, zai.i};
}
inline Zai fromCram(const cram::Zai& zai) {
  return Zai{zai.z, zai.a, zai.i};
}

struct NuclearDataAccess {
  static const cram::DepletionChain& chain(const NuclearData& data);
};

// Convenience for the two call sites that want it inline.
inline const cram::DepletionChain& chainOf(const NuclearData& data) {
  return NuclearDataAccess::chain(data);
}

}  // namespace nusift
