#pragma once
/**
 * @file
 * @brief The runtime nuclear-data model: nuclides, decay data, photon spectra, and the
 *        depletion chain the engine solves on.
 * @ingroup nucdata
 */
//
// NuclearData owns the cram::DepletionChain and everything hanging off it, and hides both
// behind a PIMPL. That is what allows cram, Eigen, and HDF5 to be PRIVATE link dependencies:
// nothing a consumer includes pulls in a sparse-matrix template. Only the handful of .cpp
// files that genuinely need the chain include nuclear_data_internal.hpp to reach it.
//
// Two construction paths, deliberately sharing everything downstream of StoreArrays:
//
//   open(path)         reads a versioned HDF5 store  -- the production path
//   fromArrays(arrays) takes the flat arrays directly -- the test path
//
// The second is why the whole engine can be tested on synthetic chains with no HDF5 file and
// no ENDF tape: a Bateman chain is twenty lines of StoreArrays.
//
// INDEX SPACE. Indices run over the CLOSED chain, which may be larger than the staged
// nuclide axis: closing registers any decay daughter that was reachable but absent, so the
// matrix can never silently drop production into a nuclide it does not know. Staged nuclides
// keep their store index, and closure-added daughters follow. A closure-added daughter is a
// stable terminator with no decay data and no photon lines, which is exactly how it should
// behave in every metric.
//
#include <cstdint>
#include <memory>
#include <span>
#include <string>

#include "nusift/core/nuclide.hpp"
#include "nusift/nucdata/fission_yield.hpp"
#include "nusift/nucdata/photon_lines.hpp"
#include "nusift/nucdata/store_arrays.hpp"

namespace nusift {

class NuclearData {
public:
  // Read a versioned HDF5 store. Throws NusiftError on a missing file, a file that is not a
  // NuSIFT store, or a store written by a newer NuSIFT than this one.
  static NuclearData open(const std::string& storePath);

  // Build directly from flat arrays. Validates them first, so a malformed synthetic chain in
  // a test fails the same way a malformed file does.
  static NuclearData fromArrays(StoreArrays arrays);

  ~NuclearData();
  NuclearData(NuclearData&&) noexcept;
  NuclearData& operator=(NuclearData&&) noexcept;
  NuclearData(const NuclearData&) = delete;
  NuclearData& operator=(const NuclearData&) = delete;

  // Number of nuclides in the closed chain: everything the matrix can address.
  int size() const;

  // Number of nuclides the STORE actually carries data for, which is generally far fewer.
  // The chain grows past the staged axis for two reasons -- closure registers reachable decay
  // daughters, and loading fission yields registers every product so a fission source has
  // somewhere to deposit. Both are necessary, and both are bookkeeping rather than evaluated
  // data. Reporting size() as though it were the store's coverage overstates it by orders of
  // magnitude: a store of three nuclides plus one fission-yield set loads as a chain of
  // twelve hundred.
  int stagedCount() const;

  // Chain index for a nuclide, or -1 if absent. Absence is a normal, reportable condition --
  // an inventory naming a nuclide with no evaluated data -- not an error at this level.
  int indexOf(const Zai& zai) const;
  int indexOfKey(std::int64_t zaiKey) const;

  Zai zaiAt(int index) const;
  std::span<const std::int64_t> nuclideKeys() const;

  // Seconds. <= 0 means stable, which is the encoding used throughout: a stable nuclide is a
  // terminator with no removal rate, not an error and not an infinity.
  double halfLifeSeconds(int index) const;

  // ln(2)/halfLife, or 0 for a stable nuclide. This is the per-nuclide weight the activity
  // metric multiplies by, so it is precomputed rather than derived at every use.
  double decayConstant(int index) const;

  // g/mol, from the staged ENDF atomic weight ratio. Returns 0 when AWR was never staged,
  // which the inventory reader turns into a specific error rather than a wrong mass: a store
  // built from a depletion-chain XML carries no AWR, so gram input is genuinely unavailable
  // rather than merely imprecise.
  double molarMassGPerMol(int index) const;

  // The nuclide's discrete photon lines, as a view into the store's flat arrays. Empty for a
  // stable nuclide, for a closure-added daughter, and for any store staged without lines.
  LineSpectrum lines(int index) const;

  // Average electromagnetic energy per decay [eV] from ENDF MT457, covering both the discrete
  // lines and any continuum.
  double emEnergyEv(int index) const;

  // Photon energy per decay carried by a CONTINUOUS spectrum, which NuSIFT does not model.
  // Non-zero here means the exposure this nuclide contributes is understated. Reporting it is
  // what keeps "no photons" distinguishable from "photons we ignore" -- without it, both look
  // like a zero contribution and the second is a silent error.
  double continuumPhotonEv(int index) const;

  // Fraction of the nuclide's photon energy that NuSIFT does not model, in [0, 1]. Zero when
  // there is no continuum or nothing to compare against.
  double unmodeledPhotonFraction(int index) const;

  // Independent fission yields, empty when the store carries none. Exposed as NuSIFT's own
  // value type rather than the chain's, so seeding never needs to reach across the PIMPL.
  const FissionYieldTable& fissionYields() const;

  const StoreProvenance& provenance() const;

  // True when the store carries discrete photon lines at all. False makes every exposure
  // metric unavailable, and the CLI turns that into an error naming the missing source
  // rather than reporting zeros.
  bool hasPhotonLines() const;

  // True when atomic weight ratios were staged, so gram <-> atom conversion is possible.
  bool hasAtomicWeights() const;

private:
  NuclearData();
  struct Impl;
  std::unique_ptr<Impl> impl_;

  friend struct NuclearDataAccess;
};

}  // namespace nusift
