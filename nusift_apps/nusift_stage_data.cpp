// nusift_stage_data -- turns ENDF tapes into a NuSIFT data store.
//
//   nusift_stage_data --decay-dir <dir>... [--nfy-dir <dir>...] [--sfy-dir <dir>...] -o out.h5
//   nusift_stage_data decay.endf... [--nfy <tape> <Z> <A>]... -o out.h5
//
// This is an offline step, run once per evaluation. Its output is a versioned HDF5 store that
// production runs read in milliseconds with no ENDFtk anywhere in the build. That separation
// is the point: ENDFtk is a heavy dependency, and a build carrying it can never be installed
// (see the install policy in the top-level CMakeLists), so it lives here and nowhere else.
//
// Decay data, branching ratios, and fission yields come through cram's readers. Photon LINE
// SPECTRA do not -- cram is deliberately a pure depletion library and carries no photon data
// -- so they are read directly from the same MF8/MT457 sections cram parses for the chain.
//
// Per-tape failures are reported to stderr and skipped rather than aborting the run. An
// evaluation is hundreds of files and one unparseable tape should cost that nuclide, not the
// whole store.
#include <CLI/CLI.hpp>
#include <ENDFtk.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <exception>
#include <filesystem>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include "cram/chain.hpp"
#include "cram/endf_reader.hpp"
#include "cram/nuclide.hpp"
#include "nusift/core/nuclide.hpp"
#include "nusift/core/nuclide_name.hpp"
#include "nusift/nucdata/data_store.hpp"
#include "nusift/units.hpp"
#include "nusift/version.hpp"

namespace fs = std::filesystem;

namespace {

using namespace nusift;

// Everything read from a tape that cram's chain does not carry.
struct ExtraData {
  double awr = 0.0;
  double emEnergyEv = 0.0;
  double continuumPhotonEv = 0.0;
  std::vector<double> lineEnergyEv;
  std::vector<double> lineIntensity;
  std::vector<int> lineStyp;
};

std::int64_t keyFromZa(int za, int liso) {
  return Zai{za / 1000, za % 1000, liso}.key();
}

// Read one decay tape's MF8/MT457 sections for the fields cram does not expose: the atomic
// weight ratio, the average electromagnetic decay energy, and the discrete photon lines.
//
// Intensities are made ABSOLUTE here (FD * RI) so that everything downstream is a plain sum
// over lines with no normalisation left to remember.
void readExtras(const std::string& path, std::map<std::int64_t, ExtraData>& extras) {
  using namespace njoy::ENDFtk;

  auto tape = tree::fromFile(path);
  for (const auto& material : tape.materials()) {
    if (!material.hasSection(8, 457)) {
      continue;
    }
    const auto section = material.section(8, 457).parse<8, 457>();
    if (section.isStable()) {
      continue;
    }

    const std::int64_t key =
        keyFromZa(static_cast<int>(section.ZA()), static_cast<int>(section.LISO()));
    ExtraData extra;
    extra.awr = static_cast<double>(section.atomicWeightRatio());

    const auto& energies = section.averageDecayEnergies();
    if (energies.numberDecayEnergies() >= 2) {
      extra.emEnergyEv = static_cast<double>(*energies.electromagneticDecayEnergy().begin());
    }

    double discreteEnergy = 0.0;
    for (const auto& spectrum : section.decaySpectra()) {
      const int styp = static_cast<int>(std::lround(spectrum.STYP()));
      // Photons only: gamma (STYP 0) and X-ray / annihilation radiation (STYP 9). Beta and
      // alpha spectra are in the same record and are not photons.
      if (styp != 0 && styp != 9) {
        continue;
      }
      const double normalisation = static_cast<double>(spectrum.discreteNormalisationFactor()[0]);
      for (const auto& line : spectrum.discreteSpectra()) {
        const double energyEv = static_cast<double>(line.discreteEnergy()[0]);
        const double intensity = normalisation * static_cast<double>(line.relativeIntensity()[0]);
        if (energyEv > 0.0 && intensity > 0.0) {
          extra.lineEnergyEv.push_back(energyEv);
          extra.lineIntensity.push_back(intensity);
          extra.lineStyp.push_back(styp);
          discreteEnergy += energyEv * intensity;
        }
      }
    }

    // The 511 keV double-counting check. Annihilation radiation is reported under STYP 9, but
    // nothing stops an evaluation from also carrying a 511 keV gamma under STYP 0, and the
    // exposure sum downstream is a plain sum over every line staged -- so an evaluation that
    // lists it twice would double the strongest line of every positron emitter.
    //
    // Reported rather than silently resolved: which of the two entries is the duplicate is a
    // judgement about that evaluation, and dropping the wrong one loses real data. ENDF/B-VIII.1
    // makes the choice in neither direction, so this is quiet on the shipped tapes -- which is
    // the point of checking rather than assuming.
    constexpr double kAnnihilationEv = 511.0e3;
    constexpr double kAnnihilationWindowEv = 100.0;  // a line NAMED 511 keV, not one near it
    bool annihilationAsGamma = false;
    bool annihilationAsXray = false;
    for (std::size_t k = 0; k < extra.lineEnergyEv.size(); ++k) {
      if (std::abs(extra.lineEnergyEv[k] - kAnnihilationEv) > kAnnihilationWindowEv) {
        continue;
      }
      if (extra.lineStyp[k] == 9) {
        annihilationAsXray = true;
      } else {
        annihilationAsGamma = true;
      }
    }
    if (annihilationAsGamma && annihilationAsXray) {
      std::fprintf(stderr,
                   "  warning: %s lists 511 keV under both STYP 0 and STYP 9; both are staged "
                   "and summed, so annihilation is counted twice for this nuclide\n",
                   formatNuclideName(Zai::fromKey(key)).c_str());
    }

    // Photon energy the discrete lines do not account for -- a continuous spectrum, most
    // often bremsstrahlung. Taken as the shortfall against the evaluated average rather than
    // by integrating the continuum record, which makes it robust to how the evaluation chose
    // to represent it and captures anything else the lines miss. The consumer needs to know
    // only that this much photon energy exists and is not modelled.
    //
    // Clamped at zero: the average and the line sum come from different parts of an
    // evaluation and disagree by a few percent, so a small negative shortfall means the lines
    // account for everything, not that there is negative continuum.
    extra.continuumPhotonEv = std::max(0.0, extra.emEnergyEv - discreteEnergy);

    extras[key] = std::move(extra);
  }
}

// The ENDF NFY incident-energy grid: thermal, fast, and 14 MeV, plus 0 for spontaneous
// fission. nearestYields snaps each probe to the closest tabulated set, and de-duplicating by
// the returned energy recovers exactly the distinct sets a tape provides.
constexpr double kProbeEnergiesEv[] = {0.0, 0.0253, 5.0e5, 1.4e7};

// Parse the fissile parent from an IAEA fission-yield filename, e.g. "nfpy_092-U-235_9228.dat"
// or "sfpy_098-Cf-252_9861.dat". The tape itself identifies its parent, but the directory
// staging path needs to know it before parsing to ask nearestYields for the right sets.
std::optional<cram::Zai> parseFissileFromName(const fs::path& path) {
  const std::string stem = path.stem().string();
  const auto first = stem.find('_');
  if (first == std::string::npos) {
    return std::nullopt;
  }
  const auto second = stem.find('_', first + 1);
  if (second == std::string::npos) {
    return std::nullopt;
  }
  const std::string middle = stem.substr(first + 1, second - first - 1);  // 092-U-235
  const auto firstDash = middle.find('-');
  const auto lastDash = middle.rfind('-');
  if (firstDash == std::string::npos || firstDash == lastDash) {
    return std::nullopt;
  }
  std::string massField = middle.substr(lastDash + 1);
  int isomer = 0;
  if (!massField.empty() && (massField.back() == 'M' || massField.back() == 'm')) {
    isomer = 1;
    massField.pop_back();
  }
  try {
    const int z = std::stoi(middle.substr(0, firstDash));
    const int a = std::stoi(massField);
    return cram::Zai{z, a, isomer};
  } catch (const std::exception&) {
    return std::nullopt;
  }
}

struct YieldSet {
  std::int64_t parentKey = 0;
  double energyEv = 0.0;
  std::vector<std::int64_t> productKey;
  std::vector<double> productYield;
};

// Load one fission-yield tape and append its distinct energy sets. Returns the number kept,
// or -1 if the tape could not be read.
int appendYields(std::vector<YieldSet>& sets, const std::string& path, const cram::Zai& parent) {
  cram::DepletionChain scratch;
  try {
    // Independent yields (MT454), never cumulative (MT459). The chain feeds precursors into
    // their daughters explicitly, so seeding with cumulative yields would count every
    // precursor decay twice.
    cram::loadFissionYields(scratch, path, /*useCumulative=*/false);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "  skip %s: %s\n", path.c_str(), e.what());
    return -1;
  }

  std::set<double> seen;
  int kept = 0;
  for (const double probe : kProbeEnergiesEv) {
    const cram::FissionYields* yields = scratch.nearestYields(parent, probe);
    if (yields == nullptr || seen.count(yields->energy) != 0) {
      continue;
    }
    seen.insert(yields->energy);

    YieldSet set;
    set.parentKey = Zai{parent.z, parent.a, parent.i}.key();
    set.energyEv = yields->energy;
    set.productKey.reserve(yields->products.size());
    set.productYield.reserve(yields->products.size());
    for (const auto& [zai, value] : yields->products) {
      set.productKey.push_back(Zai{zai.z, zai.a, zai.i}.key());
      set.productYield.push_back(value);
    }
    sets.push_back(std::move(set));
    ++kept;
  }
  return kept;
}

bool isEndfFile(const fs::path& path) {
  const std::string extension = path.extension().string();
  return extension == ".dat" || extension == ".endf" || extension == ".txt";
}

std::string utcNow() {
  const auto now = std::chrono::system_clock::now();
  const auto seconds = std::chrono::time_point_cast<std::chrono::seconds>(now);
  return std::format("{:%Y-%m-%dT%H:%M:%SZ}", seconds);
}

}  // namespace

int main(int argc, char** argv) {
  CLI::App app{"nusift_stage_data -- stage ENDF nuclear data into a NuSIFT store"};

  std::vector<std::string> decayFiles;
  std::vector<std::string> decayDirs;
  std::vector<std::string> yieldDirs;
  std::vector<std::string> sfyDirs;
  std::vector<std::tuple<std::string, int, int>> namedYieldTapes;
  std::string output;
  std::string library = "ENDF/B-VIII.1";

  app.add_option("decay.endf", decayFiles, "Decay tapes (or use --decay-dir)")
      ->check(CLI::ExistingFile);
  app.add_option("--decay-dir", decayDirs, "Directory of ENDF decay tapes")
      ->expected(1, -1)
      ->check(CLI::ExistingDirectory);
  // The validator is scoped to element 0 so it checks the tape path without rejecting the
  // Z and A integers that follow it.
  app.add_option("--nfy", namedYieldTapes, "Fission-yield tape plus parent Z A (repeatable)")
      ->check(CLI::ExistingFile.application_index(0));
  app.add_option("--nfy-dir", yieldDirs, "Directory of neutron-induced fission-yield tapes")
      ->expected(1, -1)
      ->check(CLI::ExistingDirectory);
  app.add_option("--sfy-dir", sfyDirs, "Directory of spontaneous fission-yield tapes")
      ->expected(1, -1)
      ->check(CLI::ExistingDirectory);
  app.add_option("-o,--output", output, "Output store path (.h5)")->required();
  app.add_option("--library", library, "Evaluation name recorded in the store's provenance");

  try {
    app.parse(argc, argv);
  } catch (const CLI::ParseError& e) {
    return app.exit(e);
  }

  if (decayFiles.empty() && decayDirs.empty()) {
    std::fprintf(stderr, "error: give at least one decay source (--decay-dir or a tape path)\n");
    return 2;
  }
  yieldDirs.insert(yieldDirs.end(), sfyDirs.begin(), sfyDirs.end());

  try {
    cram::DepletionChain chain;
    std::map<std::int64_t, ExtraData> extras;
    int tapesRead = 0;
    int tapesFailed = 0;

    const auto stageDecayTape = [&](const std::string& path) {
      try {
        cram::loadDecayData(chain, path);
        readExtras(path, extras);
        ++tapesRead;
      } catch (const std::exception& e) {
        std::fprintf(stderr, "  skip %s: %s\n", path.c_str(), e.what());
        ++tapesFailed;
      }
    };

    for (const std::string& path : decayFiles) {
      stageDecayTape(path);
    }
    for (const std::string& directory : decayDirs) {
      for (const fs::directory_entry& entry : fs::directory_iterator(directory)) {
        if (entry.is_regular_file() && isEndfFile(entry.path())) {
          stageDecayTape(entry.path().string());
        }
      }
      std::printf("staged decay directory %s\n", directory.c_str());
    }
    std::printf("read %d decay tape(s), %d skipped\n", tapesRead, tapesFailed);

    std::vector<YieldSet> yieldSets;
    for (const auto& [tape, z, a] : namedYieldTapes) {
      const int kept = appendYields(yieldSets, tape, cram::Zai{z, a, 0});
      std::printf("  %s: %d energy set(s) for Z=%d A=%d\n", tape.c_str(), kept, z, a);
    }
    for (const std::string& directory : yieldDirs) {
      int staged = 0;
      for (const fs::directory_entry& entry : fs::directory_iterator(directory)) {
        if (!entry.is_regular_file() || !isEndfFile(entry.path())) {
          continue;
        }
        const std::optional<cram::Zai> parent = parseFissileFromName(entry.path());
        if (!parent) {
          std::fprintf(stderr, "  skip %s: cannot read the fissile parent from the filename\n",
                       entry.path().string().c_str());
          continue;
        }
        if (appendYields(yieldSets, entry.path().string(), *parent) > 0) {
          ++staged;
        }
      }
      std::printf("staged yield directory %s: %d tape(s)\n", directory.c_str(), staged);
    }

    // Register any decay daughter that was reachable but not staged, so the matrix can never
    // drop production into a nuclide the chain does not know.
    if (const int added = chain.close(); added > 0) {
      std::printf("chain closure added %d reachable daughter(s) absent from the tapes\n", added);
    }

    // Sorted ascending by key. This is a property of the store rather than an accident of the
    // order tapes happened to be read in, which is what makes a restaged store comparable
    // against its predecessor.
    std::vector<std::int64_t> keys;
    keys.reserve(static_cast<std::size_t>(chain.size()));
    for (const cram::Zai& zai : chain.nuclides()) {
      keys.push_back(Zai{zai.z, zai.a, zai.i}.key());
    }
    std::sort(keys.begin(), keys.end());

    StoreArrays arrays;
    arrays.provenance.library = library;
    arrays.provenance.createdUtc = utcNow();
    arrays.provenance.nusiftVersion = kVersion;
    arrays.provenance.stagedTapeCount = tapesRead;
    arrays.provenance.decaySource = DataSource::Endf;
    arrays.provenance.linesSource = DataSource::Endf;
    arrays.provenance.yieldsSource = yieldSets.empty() ? DataSource::None : DataSource::Endf;

    arrays.modeOffset.push_back(0);
    arrays.lineOffset.push_back(0);
    int withLines = 0;
    for (const std::int64_t key : keys) {
      const Zai zai = Zai::fromKey(key);
      const cram::Zai cramZai{zai.z, zai.a, zai.i};
      const cram::DecayData* decay = chain.decay(cramZai);

      arrays.nuclideKey.push_back(key);
      arrays.halfLife.push_back(decay != nullptr ? decay->halfLife : 0.0);

      const auto extra = extras.find(key);
      const bool haveExtra = extra != extras.end();
      arrays.awr.push_back(haveExtra ? extra->second.awr : 0.0);
      arrays.emEnergyEv.push_back(haveExtra ? extra->second.emEnergyEv : 0.0);
      arrays.continuumPhotonEv.push_back(haveExtra ? extra->second.continuumPhotonEv : 0.0);

      if (decay != nullptr) {
        for (const cram::DecayMode& mode : decay->modes) {
          arrays.modeRtyp.push_back(mode.rtyp);
          arrays.modeBranching.push_back(mode.branching);
          arrays.modeFinalState.push_back(mode.finalState);
          arrays.modeIsFission.push_back(mode.isFission ? 1 : 0);
        }
      }
      arrays.modeOffset.push_back(static_cast<int>(arrays.modeRtyp.size()));

      if (haveExtra && !extra->second.lineEnergyEv.empty()) {
        for (std::size_t k = 0; k < extra->second.lineEnergyEv.size(); ++k) {
          arrays.lineEnergyEv.push_back(extra->second.lineEnergyEv[k]);
          arrays.lineIntensity.push_back(extra->second.lineIntensity[k]);
          arrays.lineStyp.push_back(extra->second.lineStyp[k]);
        }
        ++withLines;
      }
      arrays.lineOffset.push_back(static_cast<int>(arrays.lineEnergyEv.size()));
    }

    arrays.nfySetOffset.push_back(0);
    for (const YieldSet& set : yieldSets) {
      arrays.nfyParentKey.push_back(set.parentKey);
      arrays.nfyEnergyEv.push_back(set.energyEv);
      for (std::size_t k = 0; k < set.productKey.size(); ++k) {
        arrays.nfyProductKey.push_back(set.productKey[k]);
        arrays.nfyProductYield.push_back(set.productYield[k]);
      }
      arrays.nfySetOffset.push_back(static_cast<int>(arrays.nfyProductKey.size()));
    }

    // Average light- and heavy-particle decay energies are left unstaged. The store reserves
    // them for decay heat, which is not implemented, and an empty field is honestly "never
    // staged" whereas a column of zeros would read as "no energy".

    writeStore(output, arrays);
    std::printf("\nwrote %s\n", output.c_str());
    std::printf("  %d nuclides, %d with photon lines (%zu lines total)\n", arrays.nuclideCount(),
                withLines, arrays.lineEnergyEv.size());
    std::printf("  %zu fission-yield set(s)\n", yieldSets.size());
  } catch (const std::exception& e) {
    std::fprintf(stderr, "nusift_stage_data: %s\n", e.what());
    return 1;
  }
  return 0;
}
