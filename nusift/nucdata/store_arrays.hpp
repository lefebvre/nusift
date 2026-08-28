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

// The wire encoding of xs_reaction_type, and a STORE FORMAT rather than a program detail:
// these integers are written into HDF5 and must never be renumbered. They are spelled out
// here instead of being a static_cast of cram::ReactionType precisely because cram's enum is
// an implementation detail of the solver whose ordering carries no stability promise -- a
// reordering there would silently reinterpret every cross section in every store staged
// before it, with no version bump to notice and no error to catch it. The translation to
// cram is one explicit switch, written where the chain is built.
//
// Unknown is 0 so that a zero-filled column is not a valid channel: an unset reaction type
// should fail to resolve, not decode as fission.
enum class ReactionChannel : int {
  Unknown = 0,
  Fission = 1,
  NGamma = 2,
  N2n = 3,
  N3n = 4,
  N4n = 5,
  NAlpha = 6,
  NProton = 7,
};

// The channel's conventional spelling ("(n,gamma)", "fission", ...), which is also the name
// OpenMC's depletion chain uses for it. "unknown" for an unrecognized code.
const char* reactionChannelName(ReactionChannel channel);

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
  std::vector<std::int64_t> nuclideKey;  // Z*10000 + A*10 + I
  std::vector<double> halfLife;          // s; <= 0 means stable (a chain terminator)
  // 1-sigma absolute uncertainty on halfLife [s]; <= 0 means the evaluation states none.
  // Reserved and written empty, like the xs_* block below: an adjoint sensitivity dR/dlambda
  // is only half an answer without the sigma_lambda it multiplies, and staging reads MT457
  // for the half-life already, so the tape data is in hand the moment this is wanted.
  std::vector<double> halfLifeUncertainty;
  std::vector<double> awr;                // ENDF atomic weight ratio; 0 if unstaged
  std::vector<double> emEnergyEv;         // MT457 average electromagnetic energy per decay
  std::vector<double> lpEnergyEv;         // MT457 average light-particle energy (decay heat)
  std::vector<double> hpEnergyEv;         // MT457 average heavy-particle energy (decay heat)
  std::vector<double> continuumPhotonEv;  // photon energy per decay in a continuum, not lines

  // --- decay modes, CSR by nuclide ----------------------------------------
  std::vector<int> modeOffset;                   // length N+1
  std::vector<double> modeRtyp;                  // ENDF RTYP; may be multi-step, e.g. 1.5
  std::vector<double> modeBranching;             // branch fraction
  std::vector<double> modeBranchingUncertainty;  // 1-sigma absolute; <= 0 means none stated
  std::vector<int> modeFinalState;               // RFS: isomeric state of the daughter
  std::vector<int> modeIsFission;                // spontaneous fission flag

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
  // 1-sigma absolute uncertainty on the yield; <= 0 means the evaluation states none. Reserved
  // and written empty like the other two, and reserved at the same time as them because the
  // three are wanted together: a seed is n0_i = fissions * Y_i, so dR/dY and dR/dn0 are the
  // same elasticity, and yields carry an uncertainty comparable to the half-lives' -- an error
  // budget that staged sigma_lambda and sigma_branching but not sigma_Y would omit a parameter
  // class of the same order as the two it kept.
  std::vector<double> nfyProductYieldUncertainty;

  // --- one-group activation cross sections, CSR by irradiated parent ------
  // Reserved and written empty, so adding activation later forces neither a schema bump nor a
  // restage of everything else. The fields mirror cram's ReactionXS so the loader can feed
  // setReactions() with no translation table to keep in sync -- with one deliberate exception:
  // the channel is stored as this file's own ReactionChannel encoding rather than a cast of
  // cram::ReactionType, for the reason given on that enum.
  //
  // Parent and product are named the way cram names them, which is the opposite of what an
  // earlier draft of this block called them. The CSR axis is the IRRADIATED nuclide -- cram's
  // `parent`, the key of setReactions() -- and xsProductKey is what the channel produces,
  // which is cram's ReactionXS::target. Getting those two the wrong way round transposes every
  // reaction in the store, so they are spelled out rather than left to "target".
  std::vector<std::int64_t> xsParentKey;   // length P: the irradiated nuclide
  std::vector<int> xsOffset;               // length P+1
  std::vector<int> xsReactionType;         // ReactionChannel, not cram::ReactionType
  std::vector<std::int64_t> xsProductKey;  // cram's ReactionXS::target
  std::vector<double> xsSigmaBarn;
  std::vector<double> xsQEv;
  std::vector<double> xsEnergyEv;
  // Which one-group spectrum the sigmas were collapsed over. A cross section is meaningless
  // without it, so it travels per entry. The spectrum's human-readable name is a root
  // attribute when spectra land, not a column: a per-entry string is neither writable by the
  // numeric writeArray() nor useful repeated once per channel.
  std::vector<int> xsSpectrumId;

  int nuclideCount() const { return static_cast<int>(nuclideKey.size()); }
  int yieldSetCount() const { return static_cast<int>(nfyParentKey.size()); }
};

// Check every CSR invariant and array-length agreement, throwing NusiftError naming the
// offending field on the first violation. Called by both the HDF5 reader and the in-memory
// constructor, so a hand-built StoreArrays in a test is held to exactly the same contract as
// a file -- which is the point of having one validator rather than two.
void validateStoreArrays(const StoreArrays& arrays);

}  // namespace nusift
