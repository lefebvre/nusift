#pragma once
/**
 * @file
 * @brief The nuclear-data store as flat CSR arrays — the one representation the HDF5 reader
 *        and every in-memory consumer share.
 * @ingroup nucdata
 */
//
// StoreArrays is the seam between the file format and the runtime model. The HDF5 reader
// fills it, NuclearData is built from it, and tests construct it directly -- which is why
// the whole engine can be exercised on synthetic chains with no HDF5 file and no ENDF tape
// anywhere in sight.
//
// Everything variable-length is CSR-packed: an offset array of length N+1 indexes into a
// flat value array, so nuclide i owns [offset[i], offset[i+1]). Flat arrays mean the HDF5
// datasets are plain 1-D, which keeps the file readable by h5py and any other tool without
// knowing NuSIFT's types.
//
// The nuclide axis is sorted ascending by ZAI key. That is a deliberate property of the
// store rather than an accident of how the chain was built: it makes a restaged store
// byte-comparable against its predecessor, keeps golden tests stable, and means NuSIFT
// never depends on the iteration order of whatever built the chain.
//
#include <cstdint>
#include <string>
#include <vector>

namespace nusift {

// Where a field in the store came from. Recorded per field rather than per file because the
// two ingestion paths cover different fields: an OpenMC depletion-chain XML carries decay
// data, branchings, and fission yields but no photon lines and no atomic weight ratios. A
// store built only from XML must therefore report that exposure and gram conversions are
// unavailable, rather than silently returning zeros.
enum class DataSource : int {
  None = 0,
  Endf = 1,
  OpenmcChainXml = 2,
};

const char* dataSourceName(DataSource source);

// Provenance stamped into the store's root attributes. For a triage tool the evaluation that
// produced an answer is part of the answer, so this travels with the data and is reported by
// `nusift data info` and in every report header.
struct StoreProvenance {
  int version = 0;            // nusift_store_version
  std::string library;        // e.g. "ENDF/B-VIII.1"
  std::string createdUtc;     // ISO-8601, when the store was staged
  std::string nusiftVersion;  // the tool version that wrote it
  int stagedTapeCount = 0;
  DataSource decaySource = DataSource::None;
  DataSource linesSource = DataSource::None;
  DataSource yieldsSource = DataSource::None;
};

// The store as flat arrays. Field order and names mirror the HDF5 datasets exactly, so the
// reader is a straight loop and a mismatch is obvious on inspection.
struct StoreArrays {
  StoreProvenance provenance;

  // --- nuclide axis, length N, sorted ascending by key --------------------
  std::vector<std::int64_t> nuclideKey;   // Z*10000 + A*10 + I
  std::vector<double> halfLife;           // s; <= 0 means stable (a chain terminator)
  std::vector<double> awr;                // ENDF atomic weight ratio; 0 if unstaged
  std::vector<double> emEnergyEv;         // MT457 average electromagnetic energy per decay
  std::vector<double> lpEnergyEv;         // MT457 average light-particle energy (decay heat)
  std::vector<double> hpEnergyEv;         // MT457 average heavy-particle energy (decay heat)
  std::vector<double> continuumPhotonEv;  // photon energy per decay in a continuum, not lines

  // --- decay modes, CSR by nuclide ----------------------------------------
  std::vector<int> modeOffset;        // length N+1
  std::vector<double> modeRtyp;       // ENDF RTYP; may be multi-step, e.g. 1.5
  std::vector<double> modeBranching;  // branch fraction
  std::vector<int> modeFinalState;    // RFS: isomeric state of the daughter
  std::vector<int> modeIsFission;     // spontaneous fission flag

  // --- discrete photon lines, CSR by nuclide ------------------------------
  std::vector<int> lineOffset;  // length N+1
  std::vector<double> lineEnergyEv;
  std::vector<double> lineIntensity;  // ABSOLUTE, photons per decay
  std::vector<int> lineStyp;          // ENDF STYP: 0 gamma, 9 X-ray/annihilation

  // --- independent fission yields, CSR by set -----------------------------
  std::vector<std::int64_t> nfyParentKey;
  std::vector<double> nfyEnergyEv;  // incident neutron energy; 0 for spontaneous
  std::vector<int> nfySetOffset;    // length S+1
  std::vector<std::int64_t> nfyProductKey;
  std::vector<double> nfyProductYield;  // atoms per fission

  // --- one-group activation cross sections, CSR by target -----------------
  // Reserved in schema v1 and written empty, so adding activation later forces neither a
  // schema bump nor a restage of everything else. The fields mirror cram's ReactionXS rather
  // than a bare ENDF MT number, so the loader can feed setReactions() with no translation
  // table to keep in sync.
  std::vector<std::int64_t> xsTargetKey;
  std::vector<int> xsOffset;  // length T+1
  std::vector<int> xsReactionType;
  std::vector<std::int64_t> xsProductKey;
  std::vector<double> xsSigmaBarn;
  std::vector<double> xsQEv;
  std::vector<double> xsEnergyEv;
  std::vector<int> xsSpectrumId;
  std::vector<std::string> xsSpectrumLabel;

  int nuclideCount() const { return static_cast<int>(nuclideKey.size()); }
  int yieldSetCount() const { return static_cast<int>(nfyParentKey.size()); }
};

// Check every CSR invariant and array-length agreement, throwing NusiftError naming the
// offending field on the first violation. Called by both the HDF5 reader and the in-memory
// constructor, so a hand-built StoreArrays in a test is held to exactly the same contract as
// a file -- which is the point of having one validator rather than two.
void validateStoreArrays(const StoreArrays& arrays);

}  // namespace nusift
