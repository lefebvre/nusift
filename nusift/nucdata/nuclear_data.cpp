#include "nusift/nucdata/nuclear_data.hpp"

#include <algorithm>
#include <string>
#include <utility>
#include <vector>

#include "nusift/core/error.hpp"
#include "nusift/core/nuclide_name.hpp"
#include "nusift/nucdata/data_store.hpp"
#include "nusift/nucdata/nuclear_data_internal.hpp"
#include "nusift/units.hpp"

namespace nusift {
namespace {

constexpr const char* kModule = "nucdata";

[[noreturn]] void fail(const std::string& what) {
  throw NusiftError(tagged(kModule, what));
}

}  // namespace

// Everything cram- and Eigen-flavoured lives here and nowhere a consumer can see it.
struct NuclearData::Impl {
  cram::DepletionChain chain;
  StoreProvenance provenance;

  // Sized to the CLOSED chain, so every index in [0, chain.size()) is addressable. Entries
  // past the staged axis belong to closure-added daughters and stay at their defaults, which
  // is precisely the stable-terminator behavior they should have.
  std::vector<std::int64_t> keys;
  std::vector<double> halfLife;
  std::vector<double> halfLifeUncertainty;
  std::vector<int> modeOffset;
  std::vector<double> modeBranchingUncertainty;
  std::vector<double> lambda;
  std::vector<double> molarMass;
  std::vector<MassSource> massSource;
  std::vector<double> emEnergy;
  std::vector<double> lpEnergy;
  std::vector<double> hpEnergy;
  std::vector<double> continuumPhoton;

  // Photon lines, CSR over the closed chain. Held as GammaLine rather than as parallel
  // arrays so lines() can hand out a span with no per-call assembly.
  std::vector<int> lineOffset;
  std::vector<GammaLine> lines;

  FissionYieldTable yields;
  std::vector<Zai> sfWithoutYields;

  bool hasLines = false;
  bool hasAwr = false;
  int stagedCount = 0;
};

const cram::DepletionChain& NuclearDataAccess::chain(const NuclearData& data) {
  return data.impl_->chain;
}

NuclearData::NuclearData() : impl_(std::make_unique<Impl>()) {}
NuclearData::~NuclearData() = default;
NuclearData::NuclearData(NuclearData&&) noexcept = default;
NuclearData& NuclearData::operator=(NuclearData&&) noexcept = default;

NuclearData NuclearData::fromArrays(StoreArrays a) {
  validateStoreArrays(a);

  NuclearData data;
  Impl& impl = *data.impl_;
  impl.provenance = a.provenance;

  const int staged = a.nuclideCount();
  impl.stagedCount = staged;

  // Register in store order. add() is idempotent and appends, and the store axis is validated
  // unique, so staged nuclide i lands at chain index i exactly -- which is what lets the
  // per-nuclide arrays below be indexed by chain index without a translation table.
  for (int i = 0; i < staged; ++i) {
    const Zai zai = Zai::fromKey(a.nuclideKey[i]);
    const int index = impl.chain.add(toCram(zai));
    if (index != i) {
      fail("internal: nuclide " + formatNuclideName(zai) + " landed at chain index " +
           std::to_string(index) + ", expected " + std::to_string(i));
    }
  }

  // Decay data. A non-positive half-life is a stable terminator: it gets no DecayData at all,
  // so cram leaves it off the diagonal and nothing decays out of it.
  for (int i = 0; i < staged; ++i) {
    if (a.halfLife[i] <= 0.0) {
      continue;
    }
    // Braced, not member-assigned after a bare declaration: cram removed the default
    // initializer on halfLife precisely so -Wmissing-field-initializers catches an omitted
    // one, and a bare `cram::DecayData decay;` opts out of that check by leaving the field
    // indeterminate instead. modes is appended to below rather than supplied here, which is
    // why it keeps its initializer upstream and needs no entry.
    cram::DecayData decay{.halfLife = a.halfLife[i],
                          .decayConstant = units::decayConstant(a.halfLife[i]),
                          .gammaEnergyPerDecay = a.emEnergyEv.empty() ? 0.0 : a.emEnergyEv[i]};
    for (int m = a.modeOffset[i]; m < a.modeOffset[i + 1]; ++m) {
      // Designated, and stopping short of cram's trailing `daughter`: the store carries no
      // explicit product, so leaving it unset is what tells cram to derive the daughter from
      // RTYP and the final state -- which is exactly what ENDF decay data implies. Positional
      // initialization would express the same thing today and silently absorb whatever member
      // cram appends next.
      decay.modes.push_back(cram::DecayMode{.rtyp = a.modeRtyp[m],
                                            .branching = a.modeBranching[m],
                                            .finalState = a.modeFinalState[m],
                                            .isFission = a.modeIsFission[m] != 0});
    }
    impl.chain.setDecay(toCram(Zai::fromKey(a.nuclideKey[i])), std::move(decay));
  }

  // Independent fission yields. Loaded now even though seeding from fission lands later, so
  // that a store round-trip is lossless and the chain is complete the moment it is needed.
  for (int s = 0; s < a.yieldSetCount(); ++s) {
    // Braced for the same reason as DecayData above, and it matters more here: an omitted
    // energy is not a missing value to cram but a meaningful one -- 0 eV classifies the set
    // as spontaneous fission and sends every nearestYields() lookup to the wrong table.
    cram::FissionYields yields{.energy = a.nfyEnergyEv[s]};
    for (int p = a.nfySetOffset[s]; p < a.nfySetOffset[s + 1]; ++p) {
      yields.products.emplace_back(toCram(Zai::fromKey(a.nfyProductKey[p])), a.nfyProductYield[p]);
    }
    impl.chain.addFissionYields(toCram(Zai::fromKey(a.nfyParentKey[s])), std::move(yields));

    // A second, cram-free copy so seeding can reach the yields without crossing the PIMPL.
    // Duplicated deliberately: the chain needs them to assemble a fission source, and the
    // public API needs them to build a seed inventory, and neither should have to know about
    // the other's representation.
    FissionYieldSet set;
    set.parent = Zai::fromKey(a.nfyParentKey[s]);
    set.energyEv = a.nfyEnergyEv[s];
    for (int p = a.nfySetOffset[s]; p < a.nfySetOffset[s + 1]; ++p) {
      // The uncertainty column is optional wholesale -- a store staged before uncertainties
      // were carried has none at all -- so it is indexed only once it is known to be there.
      const double sigma =
          a.nfyProductYieldUncertainty.empty() ? 0.0 : a.nfyProductYieldUncertainty[p];
      set.products.push_back(
          FissionProduct{Zai::fromKey(a.nfyProductKey[p]), a.nfyProductYield[p], sigma});
    }
    impl.yields.add(std::move(set));
  }

  // A fission branch with no yields is the one production the matrix drops without a trace:
  // the parent's whole decay constant stays on the diagonal, and cram skips the products when
  // nearestEntry() finds no set at any energy. Recorded here, from the same arrays, so the
  // count `data info` prints and the one staging warns about are the same count.
  for (int i = 0; i < staged; ++i) {
    if (a.halfLife[i] <= 0.0) {
      continue;
    }
    bool fissions = false;
    for (int m = a.modeOffset[i]; m < a.modeOffset[i + 1]; ++m) {
      fissions = fissions || (a.modeIsFission[m] != 0 && a.modeBranching[m] > 0.0);
    }
    const Zai zai = Zai::fromKey(a.nuclideKey[i]);
    if (fissions && impl.yields.nearest(zai, kSpontaneousEv) == nullptr) {
      impl.sfWithoutYields.push_back(zai);
    }
  }

  // Register every reachable daughter that was not staged. Without this the matrix would
  // silently drop production into an unknown daughter, and atoms would vanish.
  impl.chain.close();

  const int total = impl.chain.size();
  impl.keys.resize(total);
  for (int i = 0; i < total; ++i) {
    impl.keys[i] = fromCram(impl.chain.nuclides()[i]).key();
  }

  impl.halfLife.assign(total, 0.0);
  impl.halfLifeUncertainty.assign(total, 0.0);
  // CSR over the STAGED nuclides only; closure-added daughters carry no modes and so no
  // offsets. Copied wholesale because the mode ordering has to stay the tape's -- a sigma
  // matched to the wrong mode is the failure the staging tool keys by (RTYP, RFS) to avoid,
  // and re-sorting here would undo it.
  impl.modeOffset.assign(static_cast<std::size_t>(total) + 1, 0);
  impl.modeBranchingUncertainty = a.modeBranchingUncertainty;
  impl.lambda.assign(total, 0.0);
  impl.molarMass.assign(total, 0.0);
  impl.massSource.assign(total, MassSource::None);
  impl.emEnergy.assign(total, 0.0);
  impl.lpEnergy.assign(total, 0.0);
  impl.hpEnergy.assign(total, 0.0);
  impl.continuumPhoton.assign(total, 0.0);

  impl.hasAwr = !a.awr.empty();
  for (int i = 0; i < staged; ++i) {
    impl.halfLife[i] = a.halfLife[i];
    // Guarded because the column is optional: a store staged before it existed carries an empty
    // vector rather than a short one, and reading past it would be the same silent corruption
    // the loader's length check exists to prevent.
    if (i < static_cast<int>(a.halfLifeUncertainty.size())) {
      impl.halfLifeUncertainty[i] = a.halfLifeUncertainty[static_cast<std::size_t>(i)];
    }
    if (i + 1 < static_cast<int>(a.modeOffset.size())) {
      impl.modeOffset[static_cast<std::size_t>(i)] = a.modeOffset[static_cast<std::size_t>(i)];
      impl.modeOffset[static_cast<std::size_t>(i) + 1] =
          a.modeOffset[static_cast<std::size_t>(i) + 1];
    }
    impl.lambda[i] = units::decayConstant(a.halfLife[i]);
    if (impl.hasAwr) {
      impl.molarMass[i] = units::molarMassFromAwr(a.awr[i]);
      // The source column is optional, and a store staged before masses could be mixed did
      // not carry it. Every AWR such a store holds came from an ENDF head record, so that is
      // what an absent column means -- not "unknown", which would make a perfectly good
      // legacy store unable to say where its masses came from.
      if (impl.molarMass[i] > 0.0) {
        impl.massSource[i] = i < static_cast<int>(a.awrSource.size())
                                 ? static_cast<MassSource>(a.awrSource[static_cast<std::size_t>(i)])
                                 : MassSource::Endf;
      }
    }
    if (!a.emEnergyEv.empty()) {
      impl.emEnergy[i] = a.emEnergyEv[i];
    }
    // Each guarded on its own emptiness: a store staged before these columns existed carries
    // neither, and one staged by a tool that read only some of them carries what it read.
    if (!a.lpEnergyEv.empty()) {
      impl.lpEnergy[i] = a.lpEnergyEv[i];
    }
    if (!a.hpEnergyEv.empty()) {
      impl.hpEnergy[i] = a.hpEnergyEv[i];
    }
    if (!a.continuumPhotonEv.empty()) {
      impl.continuumPhoton[i] = a.continuumPhotonEv[i];
    }
  }

  // Masses for the closure-added members, matched by key rather than by position because
  // that is the only thing the two sides share: their chain indices depend on what closure
  // added, which the staging tool cannot know when it writes the axis. Binary search over a
  // table the loader has already validated as sorted and unique.
  for (int i = staged; i < total; ++i) {
    const auto found = std::lower_bound(a.closureMassKey.begin(), a.closureMassKey.end(),
                                        impl.keys[static_cast<std::size_t>(i)]);
    if (found == a.closureMassKey.end() || *found != impl.keys[static_cast<std::size_t>(i)]) {
      continue;
    }
    const std::size_t at = static_cast<std::size_t>(found - a.closureMassKey.begin());
    impl.molarMass[i] = units::molarMassFromAwr(a.closureAwr[at]);
    if (impl.molarMass[i] > 0.0) {
      impl.massSource[i] = static_cast<MassSource>(a.closureAwrSource[at]);
    }
  }

  // Photon lines, re-CSR'd over the closed chain. The offset array is extended past the
  // staged axis with repeats of the final offset, so a closure-added daughter reports an
  // empty span rather than reading out of bounds.
  impl.hasLines = !a.lineOffset.empty() && !a.lineEnergyEv.empty();
  impl.lineOffset.assign(total + 1, 0);
  if (impl.hasLines) {
    impl.lines.reserve(a.lineEnergyEv.size());
    for (std::size_t k = 0; k < a.lineEnergyEv.size(); ++k) {
      impl.lines.push_back(GammaLine{a.lineEnergyEv[k], a.lineIntensity[k],
                                     static_cast<SpectrumType>(a.lineStyp[k])});
    }
    for (int i = 0; i <= staged; ++i) {
      impl.lineOffset[i] = a.lineOffset[i];
    }
    for (int i = staged + 1; i <= total; ++i) {
      impl.lineOffset[i] = a.lineOffset[staged];
    }
  }

  return data;
}

NuclearData NuclearData::open(const std::string& storePath) {
  return fromArrays(readStore(storePath));
}

int NuclearData::size() const {
  return impl_->chain.size();
}

int NuclearData::stagedCount() const {
  return impl_->stagedCount;
}

int NuclearData::indexOf(const Zai& zai) const {
  return impl_->chain.indexOf(toCram(zai));
}

int NuclearData::indexOfKey(std::int64_t zaiKey) const {
  return indexOf(Zai::fromKey(zaiKey));
}

Zai NuclearData::zaiAt(int index) const {
  if (index < 0 || index >= size()) {
    fail("nuclide index " + std::to_string(index) + " out of range [0, " + std::to_string(size()) +
         ")");
  }
  return fromCram(impl_->chain.nuclides()[index]);
}

std::span<const std::int64_t> NuclearData::nuclideKeys() const {
  return impl_->keys;
}

int NuclearData::modeCount(int index) const {
  const std::size_t i = static_cast<std::size_t>(index);
  if (i + 1 >= impl_->modeOffset.size()) {
    return 0;
  }
  return impl_->modeOffset[i + 1] - impl_->modeOffset[i];
}

double NuclearData::modeBranchingUncertainty(int index, int mode) const {
  const std::size_t i = static_cast<std::size_t>(index);
  if (i + 1 >= impl_->modeOffset.size() || mode < 0 || mode >= modeCount(index)) {
    return 0.0;
  }
  const std::size_t at = static_cast<std::size_t>(impl_->modeOffset[i] + mode);
  // Guarded because the column is optional: a store staged before it existed carries an empty
  // vector rather than a short one.
  return at < impl_->modeBranchingUncertainty.size() ? impl_->modeBranchingUncertainty[at] : 0.0;
}

double NuclearData::halfLifeUncertainty(int index) const {
  return impl_->halfLifeUncertainty[static_cast<std::size_t>(index)];
}

double NuclearData::halfLifeSeconds(int index) const {
  return impl_->halfLife[static_cast<std::size_t>(index)];
}

double NuclearData::decayConstant(int index) const {
  return impl_->lambda[static_cast<std::size_t>(index)];
}

double NuclearData::molarMassGPerMol(int index) const {
  return impl_->molarMass[static_cast<std::size_t>(index)];
}

MassSource NuclearData::massSource(int index) const {
  return impl_->massSource[static_cast<std::size_t>(index)];
}

LineSpectrum NuclearData::lines(int index) const {
  if (!impl_->hasLines) {
    return {};
  }
  const int begin = impl_->lineOffset[static_cast<std::size_t>(index)];
  const int end = impl_->lineOffset[static_cast<std::size_t>(index) + 1];
  return LineSpectrum(impl_->lines.data() + begin, static_cast<std::size_t>(end - begin));
}

double NuclearData::lpEnergyEv(int index) const {
  return impl_->lpEnergy[static_cast<std::size_t>(index)];
}

double NuclearData::hpEnergyEv(int index) const {
  return impl_->hpEnergy[static_cast<std::size_t>(index)];
}

double NuclearData::decayEnergyEv(int index) const {
  const std::size_t i = static_cast<std::size_t>(index);
  return impl_->emEnergy[i] + impl_->lpEnergy[i] + impl_->hpEnergy[i];
}

double NuclearData::emEnergyEv(int index) const {
  return impl_->emEnergy[static_cast<std::size_t>(index)];
}

double NuclearData::continuumPhotonEv(int index) const {
  return impl_->continuumPhoton[static_cast<std::size_t>(index)];
}

double NuclearData::unmodeledPhotonFraction(int index) const {
  const double continuum = continuumPhotonEv(index);
  if (continuum <= 0.0) {
    return 0.0;
  }
  // Against the discrete total rather than the staged EM average: the EM average and the
  // line sum come from different parts of the evaluation and disagree by a few percent, so
  // dividing by it can produce a fraction slightly outside [0, 1] for no physical reason.
  const double discrete = discretePhotonEnergyEv(lines(index));
  const double total = discrete + continuum;
  return total > 0.0 ? continuum / total : 0.0;
}

const FissionYieldTable& NuclearData::fissionYields() const {
  return impl_->yields;
}

std::vector<Zai> NuclearData::spontaneousFissionWithoutYields() const {
  return impl_->sfWithoutYields;
}

const StoreProvenance& NuclearData::provenance() const {
  return impl_->provenance;
}

bool NuclearData::hasPhotonLines() const {
  return impl_->hasLines;
}

bool NuclearData::hasAtomicWeights() const {
  return impl_->hasAwr;
}

}  // namespace nusift
