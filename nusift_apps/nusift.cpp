// nusift -- the command-line driver.
//
// Exit codes:
//   0  success
//   1  runtime failure (a solve diverged, a store is unreadable)
//   2  bad input (an unparseable argument, a malformed inventory row)
//
#include <CLI/CLI.hpp>

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <exception>
#include <fstream>
#include <iostream>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "nusift/core/error.hpp"
#include "nusift/core/nuclide_name.hpp"
#include "nusift/engine/decay_engine.hpp"
#include "nusift/engine/inventory.hpp"
#include "nusift/exposure/air_coefficients.hpp"
#include "nusift/exposure/point_source.hpp"
#include "nusift/io/inventory_io.hpp"
#include "nusift/io/report.hpp"
#include "nusift/io/time_spec.hpp"
#include "nusift/nucdata/nuclear_data.hpp"
#include "nusift/nucdata/store_locator.hpp"
#include "nusift/seed/seed_fission.hpp"
#include "nusift/triage/forecast.hpp"
#include "nusift/triage/ranking.hpp"
#include "nusift/triage/response.hpp"
#include "nusift/version.hpp"

namespace {

using namespace nusift;

// Options shared by every subcommand that decays something. Added by one function so the
// spelling, defaults, and help text cannot drift between subcommands.
struct CommonOptions {
  std::string storePath;
  std::string inventoryPath;
  bool ignoreUnknown = false;

  // Fission seeding, as an alternative source to --inventory.
  std::string seedFissile;
  std::string seedEnergy = "thermal";
  double seedFissions = 0.0;
  double seedYieldKt = 0.0;
  double seedEnergyJ = 0.0;
  std::string seedMeVPerFission = "explosive";

  std::vector<std::string> atTimes;  // repeatable --at
  std::string gridSpec;              // --times
  std::vector<std::string> intervals;

  std::string metric = "activity";
  std::string aggregate = "nuclide";
  std::string unit;
  double distanceM = 1.0;
  double airDensity = 1.205;
  bool noAirAttenuation = false;
  double buildup = 1.0;
  int topN = 10;
  double coverage = 0.0;
  double minFraction = 0.0;
  std::vector<std::string> pins;

  std::string output;
  std::string format = "text";
  int cramOrder = 48;
  bool noPrune = false;
  int threads = 0;
};

void addCommonOptions(CLI::App* app, CommonOptions& options, bool wantsTimes, bool wantsIntervals) {
  app->add_option("--store", options.storePath, "Nuclear-data store (.h5)");
  app->add_flag("--ignore-unknown", options.ignoreUnknown,
                "Skip inventory rows naming a nuclide the store does not carry");

  // Exactly one source. Enforced here rather than left to the first confusing failure
  // downstream, and the two are genuinely alternatives: an inventory read from a file, or one
  // generated from fission.
  CLI::Option* inventoryOption =
      app->add_option("-i,--inventory", options.inventoryPath, "Inventory CSV or JSON");
  CLI::Option* fissileOption =
      app->add_option("--seed-fission", options.seedFissile,
                      "Build the inventory from fission of this nuclide, e.g. U-235");
  inventoryOption->excludes(fissileOption);
  fissileOption->excludes(inventoryOption);

  app->add_option("--energy", options.seedEnergy,
                  "Incident energy: thermal, fast, 14mev, spontaneous, or a value in eV");
  app->add_option("--fissions", options.seedFissions, "Number of fissions")
      ->check(CLI::PositiveNumber);
  app->add_option("--yield-kt", options.seedYieldKt, "Fission yield in kilotons TNT")
      ->check(CLI::PositiveNumber);
  app->add_option("--energy-j", options.seedEnergyJ, "Fission energy release in joules")
      ->check(CLI::PositiveNumber);
  app->add_option("--mev-per-fission", options.seedMeVPerFission,
                  "explosive (180), recoverable (200), or a value in MeV");

  if (wantsTimes) {
    app->add_option("--at", options.atTimes, "Cooling time, repeatable (e.g. 30d, 1.5y)");
    app->add_option("--times", options.gridSpec,
                    "Time grid start:stop:log|lin:count (e.g. 1h:100y:log:60)");
  }
  if (wantsIntervals) {
    app->add_option("--interval", options.intervals,
                    "Integration window T1,T2, repeatable (e.g. 1h,30d)");
  }

  app->add_option("--metric", options.metric, "activity or exposure")
      ->check(CLI::IsMember({"activity", "exposure"}));
  app->add_option("--by", options.aggregate,
                  "Aggregate: nuclide, mass-chain, element, line (line is exposure only)")
      ->check(CLI::IsMember({"nuclide", "mass-chain", "element", "line"}));
  app->add_option("--units", options.unit,
                  "activity: Bq, Ci, decays;  exposure: R/h, Gy/h, Sv/h, R, Gy, Sv");

  // Exposure geometry. An exposure number is uninterpretable without the distance it was
  // computed at, so these are reported in the output header alongside the values.
  app->add_option("--distance", options.distanceM, "Point-source distance in metres")
      ->check(CLI::PositiveNumber);
  app->add_option("--air-density", options.airDensity, "Air density in kg/m^3 (lower at elevation)")
      ->check(CLI::PositiveNumber);
  app->add_flag("--no-air-attenuation", options.noAirAttenuation,
                "Pure inverse-square, no air path attenuation (matches published constants)");
  app->add_option("--buildup", options.buildup,
                  "Scatter buildup factor; 1.0 counts uncollided photons only")
      ->check(CLI::PositiveNumber);
  app->add_option("--top", options.topN, "Show this many contributors; 0 for all");
  app->add_option("--coverage", options.coverage,
                  "Show the fewest contributors reaching this fraction (e.g. 0.95)")
      ->check(CLI::Range(0.0, 1.0));
  app->add_option("--min-fraction", options.minFraction, "Drop contributors below this fraction")
      ->check(CLI::Range(0.0, 1.0));
  // The one option here that adds a row rather than removing one. Spelled the way the aggregate
  // in force names its contributors, and a nuclide name works for all of them.
  app->add_option("--pin", options.pins,
                  "Always show this contributor whatever it ranks, repeatable "
                  "(e.g. Cs-137, A=140, Cs)");

  app->add_option("-o,--output", options.output, "Write here instead of stdout");
  app->add_option("--format", options.format, "text, csv, or json")
      ->check(CLI::IsMember({"text", "csv", "json"}));
  app->add_option("--cram-order", options.cramOrder, "CRAM order: 16 for screening, 48 default")
      ->check(CLI::IsMember({16, 48}));
  app->add_flag("--no-prune", options.noPrune,
                "Solve the whole chain instead of the seed's forward closure");
  app->add_option("--threads", options.threads,
                  "Worker threads for the per-time solves; 0 uses every core")
      ->check(CLI::NonNegativeNumber);
}

Metric metricFrom(const std::string& text) {
  return text == "exposure" ? Metric::Exposure : Metric::Activity;
}

exposure::PointSourceGeometry geometryFrom(const CommonOptions& options) {
  exposure::PointSourceGeometry geometry;
  geometry.distanceM = options.distanceM;
  geometry.airDensityKgM3 = options.airDensity;
  geometry.airAttenuation = !options.noAirAttenuation;
  geometry.buildup = options.buildup;
  return geometry;
}

Aggregate aggregateFrom(const std::string& text) {
  if (text == "mass-chain") {
    return Aggregate::MassChain;
  }
  if (text == "element") {
    return Aggregate::Element;
  }
  if (text == "line") {
    return Aggregate::GammaLine;
  }
  return Aggregate::Nuclide;
}

DecayOptions decayOptionsFrom(const CommonOptions& options) {
  DecayOptions decayOptions;
  decayOptions.order = options.cramOrder == 16 ? CramOrder::Order16 : CramOrder::Order48;
  decayOptions.prune = !options.noPrune;
  decayOptions.threads = options.threads;
  return decayOptions;
}

// Pins are resolved against the table they will be applied to, not against the store: what a
// pin has to name is a CONTRIBUTOR, and which contributors exist depends on the aggregate and
// on how far the inventory's chain reaches. Resolving here is what lets `--pin Cs-137` be
// refused with a reason instead of quietly matching nothing.
std::vector<std::int64_t> pinsFrom(const CommonOptions& options, const ResponseTable& table) {
  std::vector<std::int64_t> keys;
  keys.reserve(options.pins.size());
  for (const std::string& text : options.pins) {
    keys.push_back(requirePin(table, text));
  }
  return keys;
}

RankRequest rankRequestFrom(const CommonOptions& options, const ResponseTable& table) {
  RankRequest request;
  request.topN = options.topN;
  request.coverage = options.coverage;
  request.minFraction = options.minFraction;
  request.pinned = pinsFrom(options, table);
  return request;
}

NuclearData openStore(const CommonOptions& options, const char* argv0, std::string& resolved) {
  StoreSearch search;
  search.explicitPath = options.storePath;
  search.executablePath = argv0 != nullptr ? argv0 : "";
  resolved = locateStore(search);
  return NuclearData::open(resolved);
}

// One line describing how an exposure was computed. Empty for activity, which needs no model
// beyond the decay constants.
std::string describeGeometry(const CommonOptions& options, Metric metric) {
  if (metric != Metric::Exposure) {
    return {};
  }
  char distance[32];
  std::snprintf(distance, sizeof(distance), "%.3g", options.distanceM);

  std::string text = "point source at ";
  text += distance;
  text += " m, ";
  text += options.noAirAttenuation ? "no air attenuation" : "air attenuation on";
  // Buildup is only worth naming when it is not the default, but when it IS set the number
  // matters more than anything else in the line -- it scales every value in the table.
  if (options.buildup != 1.0) {
    char factor[32];
    std::snprintf(factor, sizeof(factor), "%.3g", options.buildup);
    text += ", buildup x";
    text += factor;
  } else {
    text += ", uncollided only";
  }
  return text;
}

ReportContext contextFor(const NuclearData& data, const std::string& storePath,
                         const Inventory& inventory, const ResponseTable& table) {
  ReportContext context;
  context.storePath = storePath;
  context.storeLibrary = data.provenance().library;
  context.storeCreatedUtc = data.provenance().createdUtc;
  context.storeNuclideCount = data.stagedCount();
  context.seedProvenance = inventory.provenance();

  // Scanned over the whole response table, NOT over the ranked rows.
  //
  // This distinction is the difference between the warning working and not. A nuclide whose
  // photon output is entirely continuum -- Y-90's bremsstrahlung, say -- has zero MODELLED
  // exposure, so it never places in an exposure ranking and would carry its warning off the
  // page with it. That is exactly the case the reader most needs to be told about: the
  // contributor is absent from the table precisely because the part NuSIFT cannot model is
  // the only part it has.
  //
  // Named by EMITTER, not by column. Ranking by gamma line gives a flagged emitter one column
  // per line, so a column-wise list reads "Rb-90 196.8 keV, Rb-90 314.5 keV, Rb-90 543.6 keV"
  // and counts one nuclide eight times. What is understated is the nuclide; which of its lines
  // you are looking at does not change that.
  std::set<std::string> emitters;
  for (int c = 0; c < table.contributorCount(); ++c) {
    if ((table.flags[static_cast<std::size_t>(c)] & kFlagUnmodeledContinuum) == 0) {
      continue;
    }
    const ContributorId& id = table.contributors[static_cast<std::size_t>(c)];
    emitters.insert(table.aggregate == Aggregate::GammaLine
                        ? formatNuclideName(Zai::fromKey(id.dominantMemberKey))
                        : table.labels[static_cast<std::size_t>(c)]);
  }
  context.unmodeledContinuum.assign(emitters.begin(), emitters.end());
  return context;
}

// Resolve --at and --times into one ascending, de-duplicated set.
std::vector<double> timesFrom(const CommonOptions& options) {
  std::vector<double> times;
  for (const std::string& spec : options.atTimes) {
    times.push_back(parseDuration(spec));
  }
  if (!options.gridSpec.empty()) {
    for (const double t : parseTimeGrid(options.gridSpec)) {
      times.push_back(t);
    }
  }
  if (times.empty()) {
    throw InputError("time: give at least one --at or a --times grid");
  }
  return mergeTimes(std::move(times));
}

std::pair<double, double> parseInterval(const std::string& spec) {
  const std::size_t comma = spec.find(',');
  if (comma == std::string::npos) {
    throw InputError("time: an interval is T1,T2 (e.g. 1h,30d), got \"" + spec + "\"");
  }
  const double t1 = parseDuration(spec.substr(0, comma));
  const double t2 = parseDuration(spec.substr(comma + 1));
  if (t2 <= t1) {
    throw InputError("time: interval \"" + spec + "\" does not end after it starts");
  }
  return {t1, t2};
}

// Resolves -o, defaulting to stdout. Held by the caller for exactly as long as it writes.
class OutputStream {
public:
  explicit OutputStream(const std::string& path) {
    if (!path.empty()) {
      file_ = std::make_unique<std::ofstream>(path);
      if (!*file_) {
        throw InputError("output: cannot write to \"" + path + "\"");
      }
    }
  }
  std::ostream& get() { return file_ ? static_cast<std::ostream&>(*file_) : std::cout; }

private:
  std::unique_ptr<std::ofstream> file_;
};

// Resolve however the user described the fission source into a fission count.
seed::FissionSeed fissionSeedFrom(const CommonOptions& options) {
  seed::FissionSeed fissionSeed;
  fissionSeed.fissile = requireNuclideName(options.seedFissile);
  if (!parseIncidentEnergy(options.seedEnergy, fissionSeed.incidentEnergyEv)) {
    throw InputError("fission seed: \"" + options.seedEnergy +
                     "\" is not an incident energy (thermal, fast, 14mev, spontaneous, or eV)");
  }
  if (!seed::parseMeVPerFission(options.seedMeVPerFission, fissionSeed.meVPerFission)) {
    throw InputError("fission seed: \"" + options.seedMeVPerFission +
                     "\" is not an energy per fission (explosive, recoverable, or MeV)");
  }

  const int given = (options.seedFissions > 0.0 ? 1 : 0) + (options.seedYieldKt > 0.0 ? 1 : 0) +
                    (options.seedEnergyJ > 0.0 ? 1 : 0);
  if (given == 0) {
    throw InputError(
        "fission seed: give the source size as --fissions, --yield-kt, or "
        "--energy-j");
  }
  if (given > 1) {
    throw InputError(
        "fission seed: --fissions, --yield-kt and --energy-j are alternatives; "
        "give exactly one");
  }

  if (options.seedFissions > 0.0) {
    fissionSeed.fissions = options.seedFissions;
  } else if (options.seedYieldKt > 0.0) {
    fissionSeed.fissions = seed::fissionsFromKt(options.seedYieldKt, fissionSeed.meVPerFission);
  } else {
    fissionSeed.fissions =
        seed::fissionsFromEnergyJ(options.seedEnergyJ, fissionSeed.meVPerFission);
  }
  return fissionSeed;
}

Inventory loadInventory(const CommonOptions& options, const NuclearData& data) {
  if (!options.seedFissile.empty()) {
    return seed::seedFromFission(data, fissionSeedFrom(options));
  }
  if (options.inventoryPath.empty()) {
    throw InputError("give an inventory with --inventory, or generate one with --seed-fission");
  }
  InventoryReadOptions readOptions;
  readOptions.ignoreUnknown = options.ignoreUnknown;
  readOptions.warnings = &std::cerr;
  return readInventory(options.inventoryPath, data, readOptions);
}

int runRank(const CommonOptions& options, const char* argv0) {
  std::string storePath;
  const NuclearData data = openStore(options, argv0, storePath);
  const Inventory inventory = loadInventory(options, data);

  const std::vector<double> times = timesFrom(options);
  const DecayResult result = decay(data, inventory, times, decayOptionsFrom(options));

  ResponseSpec spec;
  spec.metric = metricFrom(options.metric);
  spec.aggregate = aggregateFrom(options.aggregate);
  spec.unit = requireUnit(options.unit, spec.metric, Domain::Instant);
  spec.geometry = geometryFrom(options);
  const ResponseTable table = buildResponse(data, result, spec);

  const std::vector<Ranking> rankings = rankAll(table, rankRequestFrom(options, table));

  ReportFormat format = ReportFormat::Text;
  parseReportFormat(options.format, format);
  OutputStream out(options.output);
  ReportContext context = contextFor(data, storePath, inventory, table);
  context.geometry = describeGeometry(options, spec.metric);
  writeRankings(out.get(), rankings, context, format);
  return 0;
}

int runIntegrate(const CommonOptions& options, const char* argv0) {
  if (options.intervals.empty()) {
    throw InputError("time: integrate needs at least one --interval T1,T2");
  }

  std::string storePath;
  const NuclearData data = openStore(options, argv0, storePath);
  const Inventory inventory = loadInventory(options, data);

  ResponseSpec spec;
  spec.metric = metricFrom(options.metric);
  spec.aggregate = aggregateFrom(options.aggregate);
  spec.unit = requireUnit(options.unit, spec.metric, Domain::Interval);
  spec.geometry = geometryFrom(options);

  // A context per interval. Each interval is solved separately over its own index space, so
  // the set of contributors carrying unmodelled continuum is a property of that interval --
  // building one context from the last table footnoted every ranking with the last interval's
  // emitters, which need not appear in the ranking they annotate.
  const std::string geometry = describeGeometry(options, spec.metric);
  std::vector<Ranking> rankings;
  std::vector<ReportContext> contexts;
  for (const std::string& text : options.intervals) {
    const auto [t1, t2] = parseInterval(text);
    std::vector<std::int64_t> keys;
    const std::vector<double> integral =
        intervalIntegral(data, inventory, t1, t2, &keys, decayOptionsFrom(options));
    const ResponseTable table = buildIntervalResponse(data, keys, integral, t1, t2, spec);
    // Each interval is solved over its own index space, so its pins are resolved against its own
    // table -- a pin naming something this interval's chain does not reach is refused for that
    // interval rather than silently dropped from one report out of several.
    rankings.push_back(rank(table, 0, rankRequestFrom(options, table)));
    contexts.push_back(contextFor(data, storePath, inventory, table));
    contexts.back().geometry = geometry;
  }

  ReportFormat format = ReportFormat::Text;
  parseReportFormat(options.format, format);
  OutputStream out(options.output);
  writeRankings(out.get(), rankings, contexts, format);
  return 0;
}

// Raw inventory versus time, with no ranking. The escape hatch for anyone who wants the
// numbers rather than the triage, and the first thing to reach for when a ranking looks
// wrong.
// Who leads, and when it changes. The same response table `rank` builds, read down the time
// axis instead of across it.
int runForecast(const CommonOptions& options, const char* argv0) {
  std::string storePath;
  const NuclearData data = openStore(options, argv0, storePath);
  const Inventory inventory = loadInventory(options, data);

  const std::vector<double> times = timesFrom(options);
  if (times.size() < 2) {
    throw InputError("time: a forecast needs a range of times -- try --times 1h:100y:log:60");
  }
  const DecayResult result = decay(data, inventory, times, decayOptionsFrom(options));

  ResponseSpec spec;
  spec.metric = metricFrom(options.metric);
  spec.aggregate = aggregateFrom(options.aggregate);
  spec.unit = requireUnit(options.unit, spec.metric, Domain::Instant);
  spec.geometry = geometryFrom(options);
  const ResponseTable table = buildResponse(data, result, spec);

  const std::vector<DominanceWindow> windows = dominanceWindows(table);
  // topN doubles as the depth of "ever near the top"; 0 means no limit, which for a forecast
  // would print every nuclide in the chain, so it falls back to a readable default.
  const std::vector<std::int64_t> pins = pinsFrom(options, table);
  const std::vector<RankTrack> tracks =
      unionTopN(table, options.topN > 0 ? options.topN : 10, pins);

  ReportFormat format = ReportFormat::Text;
  parseReportFormat(options.format, format);
  ReportContext context = contextFor(data, storePath, inventory, table);
  context.geometry = describeGeometry(options, spec.metric);

  OutputStream out(options.output);
  writeForecast(out.get(), windows, tracks, table, context, format);
  return 0;
}

int runDecay(const CommonOptions& options, const char* argv0) {
  std::string storePath;
  const NuclearData data = openStore(options, argv0, storePath);
  const Inventory inventory = loadInventory(options, data);

  const std::vector<double> times = timesFrom(options);
  const DecayResult result = decay(data, inventory, times, decayOptionsFrom(options));

  OutputStream out(options.output);
  std::ostream& os = out.get();
  os << "nuclide";
  for (const double t : result.times) {
    os << ',' << t;
  }
  os << '\n';
  for (int i = 0; i < result.nuclideCount(); ++i) {
    const Zai zai = Zai::fromKey(result.nuclideKeys[static_cast<std::size_t>(i)]);
    os << formatNuclideName(zai);
    for (int k = 0; k < result.timeCount(); ++k) {
      char buffer[32];
      std::snprintf(buffer, sizeof(buffer), "%.6e", result.atomsAt(k)[i]);
      os << ',' << buffer;
    }
    os << '\n';
  }
  return 0;
}

int runDataInfo(const std::string& storePath, const char* argv0) {
  CommonOptions options;
  options.storePath = storePath;
  std::string resolved;
  const NuclearData data = openStore(options, argv0, resolved);
  const StoreProvenance& provenance = data.provenance();

  std::printf("store:            %s\n", resolved.c_str());
  std::printf("schema version:   %d\n", provenance.version);
  std::printf("library:          %s\n",
              provenance.library.empty() ? "(unrecorded)" : provenance.library.c_str());
  std::printf("staged:           %s\n",
              provenance.createdUtc.empty() ? "(unrecorded)" : provenance.createdUtc.c_str());
  std::printf("staged by:        %s\n",
              provenance.nusiftVersion.empty() ? "(unrecorded)" : provenance.nusiftVersion.c_str());
  std::printf("tapes staged:     %d\n", provenance.stagedTapeCount);
  // Two counts, because they differ by orders of magnitude and only the first describes what
  // the store actually knows. Closure and fission-yield products inflate the chain with
  // nuclides that carry no evaluated data; reporting only the chain size would claim coverage
  // the store does not have.
  std::printf("nuclides staged:  %d\n", data.stagedCount());
  std::printf(
      "chain size:       %d  (staged, plus decay daughters and fission products\n"
      "                       registered so nothing is produced into a gap)\n",
      data.size());
  std::printf("decay data from:  %s\n", dataSourceName(provenance.decaySource));
  std::printf("photon lines:     %s\n", dataSourceName(provenance.linesSource));
  std::printf("fission yields:   %s\n", dataSourceName(provenance.yieldsSource));

  // Coverage diagnostics. What a store cannot answer is as much a part of its description as
  // what it can, and discovering a gap at report time is far worse than discovering it here.
  //
  // The two gaps are reported separately because they are not the same problem. A nuclide with
  // NO evaluated spectrum contributes exactly zero to an exposure ranking while genuinely
  // emitting photons -- it is invisible, not merely understated. One that has lines plus a
  // continuum tail is present but low. Lumping them into a single count hid the first behind
  // the second.
  int unstable = 0;
  int withLines = 0;
  int noSpectrumAtAll = 0;   // emits EM energy, but no discrete lines are evaluated
  int partialContinuum = 0;  // has lines, and a continuum tail above 5%
  int withWeights = 0;
  // Lines whose air coefficients are clamped end points rather than interpolations, and the
  // nuclides carrying them. Third gap, and a different one again: these lines are staged,
  // ranked, and counted, but the exposure model has nothing tabulated to evaluate them with.
  int clampedLines = 0;
  int withClampedLines = 0;
  for (int i = 0; i < data.size(); ++i) {
    if (data.decayConstant(i) > 0.0) {
      ++unstable;
      const LineSpectrum lines = data.lines(i);
      const bool hasLines = !lines.empty();
      if (hasLines) {
        ++withLines;
        if (data.unmodeledPhotonFraction(i) > 0.05) {
          ++partialContinuum;
        }
        int clampedHere = 0;
        for (const GammaLine& line : lines) {
          if (exposure::isOutsideTabulatedRange(line.energyEv)) {
            ++clampedHere;
          }
        }
        if (clampedHere > 0) {
          clampedLines += clampedHere;
          ++withClampedLines;
        }
      } else if (data.emEnergyEv(i) > 0.0) {
        ++noSpectrumAtAll;
      }
    }
    if (data.molarMassGPerMol(i) > 0.0) {
      ++withWeights;
    }
  }
  std::printf("\nunstable nuclides:              %d\n", unstable);
  std::printf("  with discrete photon lines:   %d\n", withLines);
  std::printf("    of those, >5%% continuum:    %d\n", partialContinuum);
  std::printf("    of those, clamped lines:    %d\n", withClampedLines);
  std::printf("  emit photons, no spectrum:    %d\n", noSpectrumAtAll);
  std::printf("nuclides with atomic weights:   %d\n", withWeights);

  if (noSpectrumAtAll > 0) {
    std::printf(
        "\n%d unstable nuclides have an evaluated average photon energy but no discrete\n"
        "spectrum in this evaluation. They contribute ZERO to an exposure ranking while\n"
        "genuinely emitting photons, so an exposure answer dominated by short-lived\n"
        "exotic species is understated in a way the ranking cannot show.\n",
        noSpectrumAtAll);
  }
  if (clampedLines > 0) {
    std::printf(
        "\n%d discrete lines across %d nuclides fall outside the tabulated air-coefficient\n"
        "range (%.0f keV to %.0f MeV). Their attenuation and absorption coefficients are the\n"
        "clamped end values rather than interpolations, so their contribution to an exposure\n"
        "is an order-of-magnitude figure. Nearly all are soft X-rays, which any real source\n"
        "encapsulation absorbs before they reach air.\n",
        clampedLines, withClampedLines, exposure::kMinTabulatedEv / 1.0e3,
        exposure::kMaxTabulatedEv / 1.0e6);
  }
  if (!data.hasPhotonLines()) {
    std::printf(
        "\nNote: this store carries no photon lines, so exposure metrics are unavailable.\n"
        "Stage from ENDF decay tapes to add them.\n");
  }
  if (!data.hasAtomicWeights()) {
    std::printf(
        "\nNote: this store carries no atomic weight ratios, so inventories cannot be given\n"
        "in mass units. Use atoms, moles, or an activity unit.\n");
  }
  return 0;
}

// `nusift data nuclide <name>...` -- everything the store knows about one nuclide.
//
// A triage answer is only as good as the evaluated data behind it, so being able to look at
// that data directly is part of the tool rather than a debugging aid. It is also the fastest
// way to understand why a nuclide ranks where it does -- or why it does not appear at all.
int runDataNuclide(const std::string& storePath, const std::vector<std::string>& names,
                   const char* argv0) {
  CommonOptions options;
  options.storePath = storePath;
  std::string resolved;
  const NuclearData data = openStore(options, argv0, resolved);

  for (const std::string& name : names) {
    const Zai zai = requireNuclideName(name);
    const int index = data.indexOf(zai);
    if (index < 0) {
      std::printf("%s: not in %s\n", formatNuclideName(zai).c_str(), resolved.c_str());
      continue;
    }

    const double halfLife = data.halfLifeSeconds(index);
    std::printf("%s   Z=%d A=%d I=%d   key=%lld\n", formatNuclideName(zai).c_str(), zai.z, zai.a,
                zai.i, static_cast<long long>(zai.key()));
    if (halfLife > 0.0) {
      std::printf("  half-life        %s\n", formatDuration(halfLife).c_str());
      std::printf("  decay constant   %.6g /s\n", data.decayConstant(index));
    } else {
      std::printf("  half-life        stable\n");
    }
    if (data.molarMassGPerMol(index) > 0.0) {
      std::printf("  molar mass       %.6g g/mol\n", data.molarMassGPerMol(index));
    }

    const LineSpectrum lines = data.lines(index);
    const double discrete = discretePhotonEnergyEv(lines);
    const double continuum = data.continuumPhotonEv(index);
    std::printf("  avg EM energy    %.6g eV/decay\n", data.emEnergyEv(index));
    std::printf("  discrete lines   %zu  (%.6g eV/decay)\n", lines.size(), discrete);
    if (continuum > 0.0) {
      std::printf("  continuum        %.6g eV/decay  (%.1f%% of photon energy, NOT modelled)\n",
                  continuum, data.unmodeledPhotonFraction(index) * 100.0);
    }
    if (!lines.empty()) {
      // The vacuum constant, which is what published tables quote, so a reader can check this
      // nuclide against a reference without running a whole ranking.
      const double gammaCgs = exposure::gammaConstant(lines) * 1.0e4 * 3.7e7;
      std::printf("  gamma constant   %.4g R.cm2/(h.mCi)  (vacuum, at 1 m)\n", gammaCgs);

      // Strongest first: exposure is usually driven by two or three lines out of dozens.
      std::vector<GammaLine> sorted(lines.begin(), lines.end());
      std::sort(sorted.begin(), sorted.end(), [](const GammaLine& a, const GammaLine& b) {
        return a.energyEv * a.intensity > b.energyEv * b.intensity;
      });
      const std::size_t shown = std::min<std::size_t>(sorted.size(), 10);
      std::printf("  strongest lines  (of %zu, by emitted energy)\n", sorted.size());
      for (std::size_t k = 0; k < shown; ++k) {
        // A line outside the tabulated range is evaluated with clamped air coefficients, so
        // what it contributes to an exposure is indicative rather than computed. Marked where
        // the line itself is shown, which is where anyone checking a number would look.
        const char* type = sorted[k].type == SpectrumType::XrayOrAnnih ? "X-ray/annih" : "gamma";
        if (exposure::isOutsideTabulatedRange(sorted[k].energyEv)) {
          std::printf("    %12.6g eV  x %-9.5g %-11s  ! air coefficients clamped\n",
                      sorted[k].energyEv, sorted[k].intensity, type);
        } else {
          std::printf("    %12.6g eV  x %-9.5g %s\n", sorted[k].energyEv, sorted[k].intensity,
                      type);
        }
      }
      if (sorted.size() > shown) {
        std::printf("    ... and %zu more\n", sorted.size() - shown);
      }
    }
    std::printf("\n");
  }
  return 0;
}

int runInventoryConvert(const std::string& storePath, const std::string& inputPath,
                        const std::string& outputPath, const std::string& unitText,
                        bool ignoreUnknown, const char* argv0) {
  CommonOptions options;
  options.storePath = storePath;
  std::string resolved;
  const NuclearData data = openStore(options, argv0, resolved);

  InventoryReadOptions readOptions;
  readOptions.ignoreUnknown = ignoreUnknown;
  readOptions.warnings = &std::cerr;
  const Inventory inventory = readInventory(inputPath, data, readOptions);

  Quantity unit = Quantity::Atoms;
  if (!unitText.empty() && !parseQuantity(unitText, unit)) {
    throw InputError("inventory: \"" + unitText + "\" is not a unit");
  }

  OutputStream out(outputPath);
  if (outputPath.size() >= 5 && outputPath.compare(outputPath.size() - 5, 5, ".json") == 0) {
    writeInventoryJson(out.get(), inventory, data, unit);
  } else {
    writeInventoryCsv(out.get(), inventory, data, unit);
  }
  return 0;
}

// Canonicalize nuclide names. Small, but it is the one path that needs no data store, and it
// answers the question "what does NuSIFT think this spelling in my spreadsheet means".
int runNuclide(const std::vector<std::string>& names) {
  for (const std::string& name : names) {
    const Zai zai = requireNuclideName(name);
    std::printf("%-12s Z=%-3d A=%-3d I=%d  key=%lld\n", formatNuclideName(zai).c_str(), zai.z,
                zai.a, zai.i, static_cast<long long>(zai.key()));
  }
  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  CLI::App app{"nusift -- nuclear source-term isotope forecasting and triage"};
  app.set_version_flag("--version", std::string(nusift::kVersion));
  app.require_subcommand(1);
  app.footer(
      "Times accept s, m, h, d, y suffixes; a bare number is seconds and a year is 365.25 d.");

  CommonOptions rankOptions;
  CLI::App* rankCmd = app.add_subcommand("rank", "Top contributors at one or more times");
  addCommonOptions(rankCmd, rankOptions, /*wantsTimes=*/true, /*wantsIntervals=*/false);

  CommonOptions integrateOptions;
  CLI::App* integrateCmd =
      app.add_subcommand("integrate", "Top contributors integrated over a time window");
  addCommonOptions(integrateCmd, integrateOptions, /*wantsTimes=*/false, /*wantsIntervals=*/true);

  CommonOptions spectrumOptions;
  spectrumOptions.metric = "exposure";
  spectrumOptions.aggregate = "line";
  CLI::App* spectrumCmd =
      app.add_subcommand("spectrum", "Top contributing photon lines (exposure, by line)");
  addCommonOptions(spectrumCmd, spectrumOptions, /*wantsTimes=*/true, /*wantsIntervals=*/false);

  CommonOptions forecastOptions;
  CLI::App* forecastCmd =
      app.add_subcommand("forecast", "Who dominates, and over which time windows");
  addCommonOptions(forecastCmd, forecastOptions, /*wantsTimes=*/true, /*wantsIntervals=*/false);

  CommonOptions decayCmdOptions;
  CLI::App* decayCmd = app.add_subcommand("decay", "Raw inventory versus time, unranked");
  addCommonOptions(decayCmd, decayCmdOptions, /*wantsTimes=*/true, /*wantsIntervals=*/false);

  std::string dataStorePath;
  CLI::App* dataCmd = app.add_subcommand("data", "Inspect the nuclear-data store");
  CLI::App* dataInfoCmd = dataCmd->add_subcommand("info", "Provenance and coverage");
  dataInfoCmd->add_option("--store", dataStorePath, "Nuclear-data store (.h5)");

  std::vector<std::string> inspectNames;
  std::string inspectStore;
  CLI::App* dataNuclideCmd =
      dataCmd->add_subcommand("nuclide", "What the store knows about a nuclide");
  dataNuclideCmd->add_option("name", inspectNames, "Nuclide names, e.g. Co-60 Cs-137")->required();
  dataNuclideCmd->add_option("--store", inspectStore, "Nuclear-data store (.h5)");

  std::string convertStore;
  std::string convertIn;
  std::string convertOut;
  std::string convertUnit;
  bool convertIgnoreUnknown = false;
  CLI::App* inventoryCmd = app.add_subcommand("inventory", "Inventory utilities");
  CLI::App* convertCmd =
      inventoryCmd->add_subcommand("convert", "Validate and convert an inventory's units");
  convertCmd->add_option("-i,--inventory", convertIn, "Inventory CSV or JSON")->required();
  convertCmd->add_option("-o,--output", convertOut, "Write here instead of stdout");
  convertCmd->add_option("--units", convertUnit, "Target unit (atoms, g, Bq, Ci, ...)");
  convertCmd->add_option("--store", convertStore, "Nuclear-data store (.h5)");
  convertCmd->add_flag("--ignore-unknown", convertIgnoreUnknown,
                       "Skip rows naming a nuclide the store does not carry");

  std::vector<std::string> nuclideNames;
  CLI::App* nuclideCmd = app.add_subcommand("nuclide", "Canonicalize nuclide names and show keys");
  nuclideCmd->add_option("name", nuclideNames, "Nuclide names, e.g. Cs-137 am242m 922350")
      ->required();

  try {
    app.parse(argc, argv);
  } catch (const CLI::ParseError& e) {
    return app.exit(e);
  }

  // The single catch site. Libraries throw; nothing below main handles an error, and nothing
  // above it needs to know how errors are reported.
  try {
    const char* argv0 = argc > 0 ? argv[0] : nullptr;
    if (rankCmd->parsed()) {
      return runRank(rankOptions, argv0);
    }
    if (integrateCmd->parsed()) {
      return runIntegrate(integrateOptions, argv0);
    }
    if (spectrumCmd->parsed()) {
      // A thin front on `rank --metric exposure --by line`: same code path, defaults set to
      // what someone asking about a spectrum means.
      return runRank(spectrumOptions, argv0);
    }
    if (forecastCmd->parsed()) {
      return runForecast(forecastOptions, argv0);
    }
    if (decayCmd->parsed()) {
      return runDecay(decayCmdOptions, argv0);
    }
    if (dataInfoCmd->parsed()) {
      return runDataInfo(dataStorePath, argv0);
    }
    if (dataNuclideCmd->parsed()) {
      return runDataNuclide(inspectStore, inspectNames, argv0);
    }
    if (convertCmd->parsed()) {
      return runInventoryConvert(convertStore, convertIn, convertOut, convertUnit,
                                 convertIgnoreUnknown, argv0);
    }
    if (nuclideCmd->parsed()) {
      return runNuclide(nuclideNames);
    }
    // A parent subcommand given with no leaf, e.g. `nusift data`.
    std::fprintf(stderr, "nusift: incomplete command; try --help\n");
    return 2;
  } catch (const nusift::InputError& e) {
    std::fprintf(stderr, "nusift: %s\n", e.what());
    return 2;
  } catch (const std::exception& e) {
    std::fprintf(stderr, "nusift: %s\n", e.what());
    return 1;
  }
}
