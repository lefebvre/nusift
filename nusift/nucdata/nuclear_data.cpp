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
  // is precisely the stable-terminator behaviour they should have.
  std::vector<std::int64_t> keys;
  std::vector<double> halfLife;
  std::vector<double> lambda;
  std::vector<double> molarMass;
  std::vector<double> emEnergy;
  std::vector<double> continuumPhoton;

  // Photon lines, CSR over the closed chain. Held as GammaLine rather than as parallel
  // arrays so lines() can hand out a span with no per-call assembly.
  std::vector<int> lineOffset;
  std::vector<GammaLine> lines;

  FissionYieldTable yields;

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
      decay.modes.push_back(cram::DecayMode{a.modeRtyp[m], a.modeBranching[m], a.modeFinalState[m],
                                            a.modeIsFission[m] != 0});
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
      set.products.push_back(
          FissionProduct{Zai::fromKey(a.nfyProductKey[p]), a.nfyProductYield[p]});
    }
    impl.yields.add(std::move(set));
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
  impl.lambda.assign(total, 0.0);
  impl.molarMass.assign(total, 0.0);
  impl.emEnergy.assign(total, 0.0);
  impl.continuumPhoton.assign(total, 0.0);

  impl.hasAwr = !a.awr.empty();
  for (int i = 0; i < staged; ++i) {
    impl.halfLife[i] = a.halfLife[i];
    impl.lambda[i] = units::decayConstant(a.halfLife[i]);
    if (impl.hasAwr) {
      impl.molarMass[i] = units::molarMassFromAwr(a.awr[i]);
    }
    if (!a.emEnergyEv.empty()) {
      impl.emEnergy[i] = a.emEnergyEv[i];
    }
    if (!a.continuumPhotonEv.empty()) {
      impl.continuumPhoton[i] = a.continuumPhotonEv[i];
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

double NuclearData::halfLifeSeconds(int index) const {
  return impl_->halfLife[static_cast<std::size_t>(index)];
}

double NuclearData::decayConstant(int index) const {
  return impl_->lambda[static_cast<std::size_t>(index)];
}

double NuclearData::molarMassGPerMol(int index) const {
  return impl_->molarMass[static_cast<std::size_t>(index)];
}

LineSpectrum NuclearData::lines(int index) const {
  if (!impl_->hasLines) {
    return {};
  }
  const int begin = impl_->lineOffset[static_cast<std::size_t>(index)];
  const int end = impl_->lineOffset[static_cast<std::size_t>(index) + 1];
  return LineSpectrum(impl_->lines.data() + begin, static_cast<std::size_t>(end - begin));
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
