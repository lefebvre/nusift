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
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <exception>
#include <filesystem>
#include <fstream>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "cram/chain.hpp"
#include "cram/endf_reader.hpp"
#include "cram/nuclide.hpp"
#include "nusift/core/nuclide.hpp"
#include "nusift/core/nuclide_name.hpp"
#include "nusift/nucdata/data_store.hpp"
#include "nusift/nucdata/nuclear_data.hpp"
#include "nusift/units.hpp"
#include "nusift/version.hpp"

namespace fs = std::filesystem;

namespace {

using namespace nusift;

// Everything read from a tape that cram's chain does not carry.
struct ExtraData {
  double awr = 0.0;
  // The three average decay energies of MT457, which together are the decay heat: every joule a
  // decay releases leaves as electromagnetic radiation, as light particles (betas, positrons,
  // Auger and conversion electrons, neutrinos excluded by the evaluation), or as heavy ones
  // (alphas, recoils, fission fragments). ENDF partitions them exactly so that the sum is the
  // total, which is why decay heat is a weight rather than a model.
  double emEnergyEv = 0.0;
  double lpEnergyEv = 0.0;
  double hpEnergyEv = 0.0;
  double continuumPhotonEv = 0.0;
  std::vector<double> lineEnergyEv;
  std::vector<double> lineIntensity;
  std::vector<int> lineStyp;

  // MT457's evaluated half-life, and the uncertainty on it. The VALUE is read here only to
  // check it against the one cram returns: the two come from the same record through different
  // parsers, so a disagreement means the two readers are not looking at the same nuclide and
  // every sigma attached on that basis would be attached to the wrong one.
  double halfLifeSeconds = 0.0;
  double halfLifeUncertainty = 0.0;

  // Branching uncertainties keyed by (RTYP, RFS) rather than held in tape order. cram builds
  // its mode list from this same record and almost certainly preserves the order, but "almost
  // certainly" is not a basis on which to attach an uncertainty to a decay mode -- a silent
  // off-by-one here would put beta-minus's sigma on the alpha branch and nothing downstream
  // could tell. Keyed lookup makes the assumption checkable, and unmatched modes stay zero.
  std::map<std::pair<double, double>, double> branchingUncertainty;
};

std::int64_t keyFromZa(int za, int liso) {
  return Zai{za / 1000, za % 1000, liso}.key();
}

// Read one decay tape's MF8/MT457 sections for the fields cram does not expose: the atomic
// weight ratio, the average electromagnetic decay energy, and the discrete photon lines.
//
// Intensities are made ABSOLUTE here (FD * RI) so that everything downstream is a plain sum
// over lines with no normalization left to remember.
void readExtras(const std::string& path, std::map<std::int64_t, ExtraData>& extras) {
  using namespace njoy::ENDFtk;

  auto tape = tree::fromFile(path);
  for (const auto& material : tape.materials()) {
    if (!material.hasSection(8, 457)) {
      continue;
    }
    const auto section = material.section(8, 457).parse<8, 457>();

    const std::int64_t key =
        keyFromZa(static_cast<int>(section.ZA()), static_cast<int>(section.LISO()));
    ExtraData extra;
    extra.awr = static_cast<double>(section.atomicWeightRatio());

    // A stable material has no decay data to read, but its head record still carries the
    // atomic weight ratio -- and that is the ONLY place a stable nuclide's molar mass can come
    // from, because it appears in none of the other sublibraries staged here. Skipping the
    // whole section left every stable chain end point unweighable, which is most of an
    // inventory's MASS once the short-lived activity has gone: at a year a fission source is
    // almost entirely stable end points, so a gram conversion was refusing the bulk of what it
    // was asked to weigh. So the AWR is taken and only the decay data below is skipped.
    if (section.isStable()) {
      extras[key] = std::move(extra);
      continue;
    }

    // [T, dT]: the value and its 1-sigma, in that order, as every uncertain quantity in MT457
    // is written.
    const auto halfLife = section.halfLife();
    if (halfLife.size() >= 2) {
      extra.halfLifeSeconds = static_cast<double>(halfLife[0]);
      extra.halfLifeUncertainty = static_cast<double>(halfLife[1]);
    }

    for (const auto& mode : section.decayModes().decayModes()) {
      const auto branching = mode.branchingRatio();
      if (branching.size() >= 2) {
        extra.branchingUncertainty[{static_cast<double>(mode.RTYP()),
                                    static_cast<double>(mode.RFS())}] =
            static_cast<double>(branching[1]);
      }
    }

    // All three, read together because they are one partition and are only meaningful summed.
    //
    // numberDecayEnergies() counts [value, uncertainty] PAIRS, not values, and the three sit at
    // pair indices 0, 1 and 2 in that order -- light particle, electromagnetic, heavy particle.
    // So each guard is its own index plus one, and the electromagnetic guard of 2 that was here
    // before is the middle case of the same rule rather than a different convention. A tape
    // carrying fewer pairs stopped early rather than disagreeing, so what is there is taken and
    // the rest stay zero.
    const auto& energies = section.averageDecayEnergies();
    if (energies.numberDecayEnergies() >= 1) {
      extra.lpEnergyEv = static_cast<double>(*energies.lightParticleDecayEnergy().begin());
    }
    if (energies.numberDecayEnergies() >= 2) {
      extra.emEnergyEv = static_cast<double>(*energies.electromagneticDecayEnergy().begin());
    }
    if (energies.numberDecayEnergies() >= 3) {
      extra.hpEnergyEv = static_cast<double>(*energies.heavyParticleDecayEnergy().begin());
    }

    double discreteEnergy = 0.0;
    for (const auto& spectrum : section.decaySpectra()) {
      const int styp = static_cast<int>(std::lround(spectrum.STYP()));
      // Photons only: gamma (STYP 0) and X-ray / annihilation radiation (STYP 9). Beta and
      // alpha spectra are in the same record and are not photons.
      if (styp != 0 && styp != 9) {
        continue;
      }
      const double normalization = static_cast<double>(spectrum.discreteNormalisationFactor()[0]);
      for (const auto& line : spectrum.discreteSpectra()) {
        const double energyEv = static_cast<double>(line.discreteEnergy()[0]);
        const double intensity = normalization * static_cast<double>(line.relativeIntensity()[0]);
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
    // only that this much photon energy exists and is not modeled.
    //
    // Clamped at zero: the average and the line sum come from different parts of an
    // evaluation and disagree by a few percent, so a small negative shortfall means the lines
    // account for everything, not that there is negative continuum.
    extra.continuumPhotonEv = std::max(0.0, extra.emEnergyEv - discreteEnergy);

    extras[key] = std::move(extra);
  }
}

// Surrounding blanks off a fixed-width field. The AME table right-justifies its numbers into
// columns, so every field arrives padded.
std::string trimmed(std::string text) {
  const auto notSpace = [](unsigned char c) { return std::isspace(c) == 0; };
  text.erase(text.begin(), std::find_if(text.begin(), text.end(), notSpace));
  text.erase(std::find_if(text.rbegin(), text.rend(), notSpace).base(), text.end());
  return text;
}

// One nuclide's mass as AME2020 states it.
struct AmeMass {
  double massAmu = 0.0;
  bool estimated = false;  // extrapolated from systematics rather than measured
};

// Read the AME2020 atomic mass table (mass_1.mas20.txt) into a map keyed by (Z, A).
//
// PARSED BY COLUMN, never by splitting on whitespace. The table is Fortran-formatted
//
//   a1,i3,i5,i5,i5,1x,a3,a4,1x,f14.6,f12.6,f13.5,1x,f10.5,1x,a2,f13.5,f11.5,1x,i3,1x,f13.6,f12.6
//
// and its fields run together: a wide enough mass excess touches the field beside it, and an
// element symbol can carry an origin flag like "-n" or "-pp" in the column after it. Splitting
// on spaces reads a different quantity for some rows than for others, and silently.
//
// The atomic mass is the last value pair, and is written in micro-u SPLIT ACROSS TWO FIELDS:
// an integer count of whole u (i3) and the remainder in micro-u (f13.6), so the neutron is
// "  1 008664.91590". The two are recombined here rather than stored apart.
//
// ESTIMATED VALUES carry '#' in place of the decimal point. They are taken -- see the note in
// data/ame/fetch_ame.sh on why an extrapolated mass is still far more precision than a gram
// conversion needs -- but the flag is kept so the store can report which is which.
std::map<std::pair<int, int>, AmeMass> readAmeMasses(const std::string& path) {
  // Field extents, 0-based and half-open, straight off the format above.
  constexpr std::size_t kZBegin = 9;
  constexpr std::size_t kZEnd = 14;
  constexpr std::size_t kABegin = 14;
  constexpr std::size_t kAEnd = 19;
  constexpr std::size_t kMassWholeBegin = 106;
  constexpr std::size_t kMassWholeEnd = 109;
  constexpr std::size_t kMassMicroBegin = 110;
  constexpr std::size_t kMassMicroEnd = 123;

  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("cannot open AME mass table " + path);
  }

  // The table opens with ~35 lines of prose and column headings. Rather than counting them --
  // a count that differs between AME releases -- every line is offered to the same parse and
  // one that does not yield an integer Z, an integer A, and a mass is skipped. A heading line
  // cannot pass all three.
  std::map<std::pair<int, int>, AmeMass> masses;
  std::string line;
  while (std::getline(in, line)) {
    if (line.size() < kMassMicroEnd) {
      continue;
    }
    const auto integerField = [&line](std::size_t begin, std::size_t end) -> std::optional<int> {
      const std::string text = trimmed(line.substr(begin, end - begin));
      if (text.empty()) {
        return std::nullopt;
      }
      try {
        std::size_t consumed = 0;
        const int value = std::stoi(text, &consumed);
        return consumed == text.size() ? std::optional<int>(value) : std::nullopt;
      } catch (const std::exception&) {
        return std::nullopt;
      }
    };

    const std::optional<int> z = integerField(kZBegin, kZEnd);
    const std::optional<int> a = integerField(kABegin, kAEnd);
    const std::optional<int> whole = integerField(kMassWholeBegin, kMassWholeEnd);
    if (!z || !a || !whole || *a <= 0) {
      continue;
    }

    std::string micro = trimmed(line.substr(kMassMicroBegin, kMassMicroEnd - kMassMicroBegin));
    const bool estimated = micro.find('#') != std::string::npos;
    std::replace(micro.begin(), micro.end(), '#', '.');
    double microU = 0.0;
    try {
      std::size_t consumed = 0;
      microU = std::stod(micro, &consumed);
      if (consumed != micro.size()) {
        continue;
      }
    } catch (const std::exception&) {
      continue;
    }

    constexpr double kMicroUPerU = 1.0e6;
    masses[{*z, *a}] =
        AmeMass{.massAmu = (static_cast<double>(*whole) * kMicroUPerU + microU) / kMicroUPerU,
                .estimated = estimated};
  }

  if (masses.empty()) {
    throw std::runtime_error("no masses parsed from " + path + "; is it the AME mass table?");
  }
  return masses;
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
  std::vector<double> productYieldUncertainty;
};

// MT454's yield uncertainties, by incident energy and then by product. Read directly for the
// same reason the photon lines are: cram's reader carries the yields and not the DFY column
// beside them.
//
// Keyed by product rather than held in tape order, on the same argument the branching sigmas
// are: cram de-duplicates and re-sorts the product list, so a positional pairing would attach
// each uncertainty to whichever product happened to land in that slot.
using YieldUncertainties = std::map<double, std::map<std::int64_t, double>>;

YieldUncertainties readYieldUncertainties(const std::string& path) {
  using namespace njoy::ENDFtk;

  YieldUncertainties out;
  auto tape = tree::fromFile(path);
  for (const auto& material : tape.materials()) {
    if (!material.hasSection(8, 454)) {
      continue;
    }
    const auto section = material.section(8, 454).parse<8, 454>();
    for (const auto& set : section.yields()) {
      auto& byProduct = out[static_cast<double>(set.E())];
      const auto identifiers = set.ZAFP();
      const auto states = set.FPS();
      const auto sigmas = set.DFY();
      const std::size_t n =
          static_cast<std::size_t>(std::distance(identifiers.begin(), identifiers.end()));
      auto id = identifiers.begin();
      auto state = states.begin();
      auto sigma = sigmas.begin();
      for (std::size_t k = 0; k < n; ++k, ++id, ++state, ++sigma) {
        const double value = static_cast<double>(*sigma);
        if (!(value > 0.0)) {
          continue;
        }
        const std::int64_t key =
            keyFromZa(static_cast<int>(std::lround(*id)), static_cast<int>(std::lround(*state)));
        // Summed rather than overwritten: a tape listing a product twice in one set is listing
        // two contributions to one yield, and cram sums the values, so the uncertainties have
        // to combine the same way or they would describe a different number than the yield does.
        double& held = byProduct[key];
        held = std::hypot(held, value);
      }
    }
  }
  return out;
}

// Load one fission-yield tape and append its distinct energy sets. Returns the number kept,
// or -1 if the tape could not be read.
int appendYields(std::vector<YieldSet>& sets, const std::string& path, const cram::Zai& parent,
                 int& withSigma, int& withoutSigma) {
  cram::DepletionChain scratch;
  YieldUncertainties uncertainties;
  try {
    uncertainties = readYieldUncertainties(path);
  } catch (const std::exception& e) {
    // A tape whose uncertainties cannot be read still has usable yields. Reported, and the
    // yields staged without sigmas rather than the whole tape dropped.
    std::fprintf(stderr, "  warning: no yield uncertainties from %s: %s\n", path.c_str(), e.what());
  }
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
    set.productYieldUncertainty.reserve(yields->products.size());

    // The energy cram reports came from this tape, so the exact match is the expected path;
    // the tolerance is for a reader that rounded on the way through rather than for a guess
    // about which set was meant.
    const std::map<std::int64_t, double>* forEnergy = nullptr;
    for (const auto& [energy, byProduct] : uncertainties) {
      const double scale = std::max(std::abs(energy), std::abs(yields->energy));
      if (std::abs(energy - yields->energy) <= std::max(1.0e-9, 1.0e-9 * scale)) {
        forEnergy = &byProduct;
        break;
      }
    }

    for (const auto& [zai, value] : yields->products) {
      const std::int64_t productKey = Zai{zai.z, zai.a, zai.i}.key();
      set.productKey.push_back(productKey);
      set.productYield.push_back(value);
      double sigma = 0.0;
      if (forEnergy != nullptr) {
        const auto found = forEnergy->find(productKey);
        if (found != forEnergy->end()) {
          sigma = found->second;
        }
      }
      set.productYieldUncertainty.push_back(sigma);
      if (sigma > 0.0) {
        ++withSigma;
      } else {
        ++withoutSigma;
      }
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
  std::string amePath;
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
  app.add_option("--ame", amePath,
                 "AME2020 atomic mass table (mass_1.mas20.txt), for the nuclides ENDF states "
                 "no atomic weight ratio for")
      ->check(CLI::ExistingFile);
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
    std::map<std::pair<int, int>, AmeMass> ameMasses;
    if (!amePath.empty()) {
      ameMasses = readAmeMasses(amePath);
      std::printf("staged AME mass table %s: %zu nuclide(s)\n", amePath.c_str(), ameMasses.size());
    }
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
    int yieldsWithSigma = 0;
    int yieldsWithout = 0;
    for (const auto& [tape, z, a] : namedYieldTapes) {
      const int kept =
          appendYields(yieldSets, tape, cram::Zai{z, a, 0}, yieldsWithSigma, yieldsWithout);
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
        if (appendYields(yieldSets, entry.path().string(), *parent, yieldsWithSigma,
                         yieldsWithout) > 0) {
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
    int withHalfLifeSigma = 0;
    int withBranchingSigma = 0;
    int unmatchedModes = 0;
    int halfLifeMismatches = 0;
    int massFromEndf = 0;
    int massFromAme = 0;
    int massFromAmeEstimated = 0;
    int massMismatches = 0;
    std::vector<Zai> withoutMass;
    for (const std::int64_t key : keys) {
      const Zai zai = Zai::fromKey(key);
      const cram::Zai cramZai{zai.z, zai.a, zai.i};
      const cram::DecayData* decay = chain.decay(cramZai);

      arrays.nuclideKey.push_back(key);
      arrays.halfLife.push_back(decay != nullptr ? decay->halfLife : 0.0);

      const auto extra = extras.find(key);
      const bool haveExtra = extra != extras.end();

      // The cross-check the sigma rests on. cram and the direct reader parse the same MT457
      // record, so their half-lives must agree; if they do not, the two are not describing the
      // same nuclide and every uncertainty attached on that basis would be attached to the
      // wrong one. Reported per nuclide and the sigma withheld, rather than staged on a
      // pairing that has just been shown to be wrong.
      double halfLifeSigma = 0.0;
      if (haveExtra && decay != nullptr && extra->second.halfLifeUncertainty > 0.0) {
        const double fromCram = decay->halfLife;
        const double fromTape = extra->second.halfLifeSeconds;
        const double scale = std::max(std::abs(fromCram), std::abs(fromTape));
        if (scale > 0.0 && std::abs(fromCram - fromTape) > 1.0e-6 * scale) {
          std::fprintf(stderr,
                       "  warning: %s half-life differs between readers (%g s vs %g s); its "
                       "uncertainty is not staged\n",
                       formatNuclideName(zai).c_str(), fromCram, fromTape);
          ++halfLifeMismatches;
        } else {
          halfLifeSigma = extra->second.halfLifeUncertainty;
          ++withHalfLifeSigma;
        }
      }
      arrays.halfLifeUncertainty.push_back(halfLifeSigma);
      // AME2020 FIRST, ENDF only where AME has no entry.
      //
      // Both are atomic weight ratios against the neutron mass by the time they land in this
      // column -- AME states the atomic mass the ratio is formed from, ENDF states the ratio
      // directly -- so nothing downstream has to know which of the two a nuclide's mass came
      // from. The source column records it anyway, because "the evaluation measured this" and
      // "a mass model extrapolates this" are different claims about the same number.
      //
      // AME wins because ENDF/B-VIII.1's decay sublibrary is demonstrably wrong for some
      // nuclides and AME is what its correct values are derived from in the first place. Seven
      // tapes -- Cu-81, Zr-110, Rh-123, Pd-125, Pd-126, I-145, Ba-153 -- put the atomic mass
      // in u in the AWR field instead of the ratio to the neutron mass, which is a uniform
      // +0.87% error in the molar mass, exactly the neutron-mass factor. The cross-check below
      // is what found them, and preferring ENDF would mean knowingly staging a value the check
      // has just reported as wrong.
      //
      // AME is keyed by (Z, A) and carries no isomeric states, so an isomer takes its GROUND
      // STATE mass. The excitation energy that omits is at most a few MeV, a part in 1e5 of a
      // fission product's mass and orders below anything a gram conversion resolves.
      const double endfAwr = haveExtra ? extra->second.awr : 0.0;
      const auto ame = ameMasses.find({zai.z, zai.a});
      const bool haveAme = ame != ameMasses.end();
      const double ameAwr = haveAme ? ame->second.massAmu / units::kNeutronMassAmu : 0.0;

      double awr = 0.0;
      MassSource massSource = MassSource::None;
      if (haveAme) {
        awr = ameAwr;
        massSource = ame->second.estimated ? MassSource::AmeEstimated : MassSource::Ame;
        if (ame->second.estimated) {
          ++massFromAmeEstimated;
        } else {
          ++massFromAme;
        }
        // The same kind of cross-check the half-life sigma rests on, and for the same reason:
        // two readers, two evaluations, one nuclide. ENDF's ratios derive from AME's masses,
        // so they agree to rounding wherever both exist, and a real disagreement means one of
        // the two tapes is stating something other than what it claims to. Reported per
        // nuclide rather than reconciled, because which one is wrong is a judgement about that
        // evaluation and the count belongs in front of whoever restages.
        constexpr double kMassAgreement = 1.0e-4;  // relative; ENDF rounds AWR, AME does not
        if (endfAwr > 0.0 && std::abs(endfAwr - ameAwr) > kMassAgreement * ameAwr) {
          std::fprintf(stderr,
                       "  warning: %s atomic weight differs between ENDF and AME2020 "
                       "(%.6f vs %.6f); the AME value is staged\n",
                       formatNuclideName(zai).c_str(), endfAwr, ameAwr);
          ++massMismatches;
        }
      } else if (endfAwr > 0.0) {
        awr = endfAwr;
        massSource = MassSource::Endf;
        ++massFromEndf;
      } else {
        withoutMass.push_back(zai);
      }
      arrays.awr.push_back(awr);
      arrays.awrSource.push_back(static_cast<int>(massSource));
      arrays.emEnergyEv.push_back(haveExtra ? extra->second.emEnergyEv : 0.0);
      arrays.lpEnergyEv.push_back(haveExtra ? extra->second.lpEnergyEv : 0.0);
      arrays.hpEnergyEv.push_back(haveExtra ? extra->second.hpEnergyEv : 0.0);
      arrays.continuumPhotonEv.push_back(haveExtra ? extra->second.continuumPhotonEv : 0.0);

      if (decay != nullptr) {
        for (const cram::DecayMode& mode : decay->modes) {
          arrays.modeRtyp.push_back(mode.rtyp);
          arrays.modeBranching.push_back(mode.branching);
          arrays.modeFinalState.push_back(mode.finalState);
          arrays.modeIsFission.push_back(mode.isFission ? 1 : 0);
          // Matched by (RTYP, RFS), never by position. A mode the tape's own list does not
          // carry under that pair gets no sigma rather than its neighbour's.
          double branchingSigma = 0.0;
          if (haveExtra) {
            const auto found = extra->second.branchingUncertainty.find(
                {mode.rtyp, static_cast<double>(mode.finalState)});
            if (found != extra->second.branchingUncertainty.end()) {
              branchingSigma = found->second;
              ++withBranchingSigma;
            } else {
              ++unmatchedModes;
            }
          }
          arrays.modeBranchingUncertainty.push_back(branchingSigma);
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
        arrays.nfyProductYieldUncertainty.push_back(set.productYieldUncertainty[k]);
      }
      arrays.nfySetOffset.push_back(static_cast<int>(arrays.nfyProductKey.size()));
    }

    // Average light- and heavy-particle decay energies are left unstaged. The store reserves
    // them for decay heat, which is not implemented, and an empty field is honestly "never
    // staged" whereas a column of zeros would read as "no energy".

    // Masses for the nuclides the LOADER will add but the axis does not carry. Which ones
    // those are is not something this tool can work out from its own chain: closure happens
    // when the store is read, from the fission-yield products and decay daughters the arrays
    // name. So the arrays are loaded here exactly as a run will load them, and every chain
    // member past the staged axis is looked up in AME -- the only evaluation that has a mass
    // for a nuclide nothing else evaluated at all.
    {
      const NuclearData closed = NuclearData::fromArrays(arrays);
      for (int i = closed.stagedCount(); i < closed.size(); ++i) {
        const Zai zai = closed.zaiAt(i);
        const auto ame = ameMasses.find({zai.z, zai.a});
        if (ame == ameMasses.end()) {
          withoutMass.push_back(zai);
          continue;
        }
        arrays.closureMassKey.push_back(zai.key());
        arrays.closureAwr.push_back(ame->second.massAmu / units::kNeutronMassAmu);
        arrays.closureAwrSource.push_back(
            static_cast<int>(ame->second.estimated ? MassSource::AmeEstimated : MassSource::Ame));
        if (ame->second.estimated) {
          ++massFromAmeEstimated;
        } else {
          ++massFromAme;
        }
      }
      // The axis is sorted ascending by key and closure appends in chain order, which is not
      // the same order. Sorted here because the loader binary-searches this table.
      std::vector<std::size_t> order(arrays.closureMassKey.size());
      for (std::size_t k = 0; k < order.size(); ++k) {
        order[k] = k;
      }
      std::sort(order.begin(), order.end(), [&](std::size_t l, std::size_t r) {
        return arrays.closureMassKey[l] < arrays.closureMassKey[r];
      });
      std::vector<std::int64_t> keys;
      std::vector<double> awrs;
      std::vector<int> sources;
      keys.reserve(order.size());
      awrs.reserve(order.size());
      sources.reserve(order.size());
      for (const std::size_t k : order) {
        keys.push_back(arrays.closureMassKey[k]);
        awrs.push_back(arrays.closureAwr[k]);
        sources.push_back(arrays.closureAwrSource[k]);
      }
      arrays.closureMassKey = std::move(keys);
      arrays.closureAwr = std::move(awrs);
      arrays.closureAwrSource = std::move(sources);
    }

    // The one coverage gap the decay matrix cannot report for itself: a spontaneous-fission
    // branch with no yield set of any energy removes the parent's atoms and produces nothing
    // in their place. Counted by loading the arrays exactly as a run will, so the number here
    // is the number `nusift data info` prints for the store this writes.
    const std::vector<Zai> leaking =
        NuclearData::fromArrays(arrays).spontaneousFissionWithoutYields();
    if (!leaking.empty()) {
      std::fprintf(stderr,
                   "  warning: %zu nuclide(s) have a spontaneous-fission branch but no "
                   "fission-yield set, so the decay matrix will lose their fissioning atoms "
                   "(first: %s). Stage the SFY sublibrary with --sfy-dir to close the gap.\n",
                   leaking.size(), formatNuclideName(leaking.front()).c_str());
    }

    writeStore(output, arrays);
    std::printf("\nwrote %s\n", output.c_str());
    std::printf("  %d nuclides, %d with photon lines (%zu lines total)\n", arrays.nuclideCount(),
                withLines, arrays.lineEnergyEv.size());
    std::printf("  %zu fission-yield set(s)\n", yieldSets.size());
    // The uncertainty census. Printed because how MUCH of an evaluation carries a sigma is the
    // thing that decides whether an error budget built on it means anything, and a store that
    // reported only that the columns exist would not say.
    std::printf("  uncertainties: %d/%d half-lives, %d/%zu branchings, %d/%zu yields\n",
                withHalfLifeSigma, arrays.nuclideCount(), withBranchingSigma,
                arrays.modeBranching.size(), yieldsWithSigma, arrays.nfyProductYield.size());
    if (unmatchedModes > 0) {
      std::fprintf(stderr,
                   "  warning: %d decay mode(s) had no (RTYP, RFS) match in the tape's own mode "
                   "list, so they carry no branching uncertainty\n",
                   unmatchedModes);
    }
    if (halfLifeMismatches > 0) {
      std::fprintf(stderr,
                   "  warning: %d nuclide(s) had disagreeing half-lives between the two readers "
                   "and carry no half-life uncertainty\n",
                   halfLifeMismatches);
    }
    // The mass census, by where each one came from. Reported like the uncertainty census
    // above and for the same reason: which nuclides can be given or reported in grams is a
    // property of the store a user has to be able to ask about, and one total would hide that
    // the answer differs per nuclide.
    std::printf(
        "  atomic weights: %d from AME2020, %d AME2020-estimated, %d from ENDF, "
        "%zu with none\n",
        massFromAme, massFromAmeEstimated, massFromEndf, withoutMass.size());
    if (!withoutMass.empty()) {
      std::fprintf(stderr,
                   "  warning: %zu nuclide(s) have no atomic weight from any source, so they "
                   "cannot be given or reported in grams (first: %s).%s\n",
                   withoutMass.size(), formatNuclideName(withoutMass.front()).c_str(),
                   ameMasses.empty() ? " Pass --ame with the AME2020 mass table to cover the"
                                       " nuclides ENDF states no AWR for."
                                     : "");
    }
    if (massMismatches > 0) {
      std::fprintf(stderr,
                   "  warning: %d nuclide(s) have ENDF and AME2020 atomic weights that "
                   "disagree by more than rounding; the ENDF value is staged for each\n",
                   massMismatches);
    }
  } catch (const std::exception& e) {
    std::fprintf(stderr, "nusift_stage_data: %s\n", e.what());
    return 1;
  }
  return 0;
}
