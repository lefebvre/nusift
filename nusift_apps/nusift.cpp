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
#include <optional>
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
#include "nusift/nucdata/coefficient_pack.hpp"
#include "nusift/nucdata/nuclear_data.hpp"
#include "nusift/nucdata/store_locator.hpp"
#include "nusift/seed/seed_fission.hpp"
#include "nusift/triage/allowable.hpp"
#include "nusift/triage/attribution.hpp"
#include "nusift/triage/events.hpp"
#include "nusift/triage/forecast.hpp"
#include "nusift/triage/intervention.hpp"
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

  // A coefficient pack, which becomes the metric when given. Named by path rather than by a
  // registry name: a pack is a file someone chose, and resolving a name would need a search
  // order that could quietly pick a different edition than the one intended.
  std::string packPath;

  std::string metric = "activity";
  std::string aggregate = "nuclide";
  std::string unit;
  std::string irradiation = "ap";
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

  app->add_option("--metric", options.metric, "activity, exposure, or photon")
      ->check(CLI::IsMember({"activity", "exposure", "photon"}));
  app->add_option("--pack", options.packPath,
                  "Rank by a coefficient pack instead: a CSV of per-nuclide coefficients with "
                  "its quantity, units, basis and version in its header");
  app->add_option("--by", options.aggregate,
                  "Aggregate: nuclide, mass-chain, element, line (line is exposure or photon "
                  "only)")
      ->check(CLI::IsMember({"nuclide", "mass-chain", "element", "line"}));
  app->add_option("--units", options.unit,
                  "activity: Bq, Ci, decays;  exposure: R/h, Gy/h, Sv/h, R, Gy, Sv;  photon: "
                  "photons/s, photons, photons/m2/s, photons/m2");

  // Exposure geometry. An exposure number is uninterpretable without the distance it was
  // computed at, so these are reported in the output header alongside the values.
  app->add_option("--distance", options.distanceM, "Point-source distance in metres")
      ->check(CLI::PositiveNumber);
  // Effective dose is defined for a body in a field, so how the body stands in it is half of
  // what a sievert means. Ignored by every unit that is not a sievert, which is why it sits
  // beside the distance rather than in front of it.
  app->add_option("--irradiation", options.irradiation,
                  "Irradiation geometry for sievert units: ap, pa, llat, rlat, rot, iso")
      ->check(CLI::IsMember(
          {"ap", "pa", "llat", "rlat", "rot", "iso", "AP", "PA", "LLAT", "RLAT", "ROT", "ISO"}));
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

// The metric a command is actually computing. --pack wins over --metric, and says so rather
// than silently ignoring the other: two ways to name one thing is how a report ends up
// describing a quantity it did not compute.
Metric metricFrom(const CommonOptions& options);

Metric metricFrom(const std::string& text) {
  if (text == "exposure") {
    return Metric::Exposure;
  }
  if (text == "photon") {
    return Metric::Photon;
  }
  return Metric::Activity;
}

// Owns a pack and its resolution against the seed for as long as a spec points at them. The
// spec holds a pointer rather than a copy because the resolved weight vector is the length of
// the store, and because there is exactly one right answer per (pack, seed) pair -- resolving
// it twice in one command would be two chances to resolve it differently.
class PackHolder {
public:
  PackHolder(const CommonOptions& options, const NuclearData& data, const Inventory& seed) {
    if (options.packPath.empty()) {
      return;
    }
    pack_ = CoefficientPack::open(options.packPath);
    resolved_ = resolvePack(*pack_, data, seed);
  }

  const ResolvedPack* resolved() const { return resolved_ ? &*resolved_ : nullptr; }

  // What the report prints in its header: the pack, the edition its numbers are from, and the
  // scenario they were tabulated under. All three, because any of them alone would leave a
  // number nobody could reproduce.
  std::string describe() const {
    if (!pack_) {
      return {};
    }
    const PackProvenance& p = pack_->provenance();
    std::string text = p.name + " " + p.version + ", " + p.quantity;
    if (!p.scenario.empty()) {
      text += " (" + p.scenario + ")";
    }
    return text;
  }

private:
  std::optional<CoefficientPack> pack_;
  std::optional<ResolvedPack> resolved_;
};

Metric metricFrom(const CommonOptions& options) {
  if (options.packPath.empty()) {
    return metricFrom(options.metric);
  }
  if (options.metric != "activity") {
    throw InputError("metric: --pack IS the metric, so --metric " + options.metric +
                     " cannot also apply. Drop one of them");
  }
  return Metric::Pack;
}

exposure::PointSourceGeometry geometryFrom(const CommonOptions& options) {
  exposure::PointSourceGeometry geometry;
  geometry.distanceM = options.distanceM;
  geometry.airDensityKgM3 = options.airDensity;
  geometry.airAttenuation = !options.noAirAttenuation;
  geometry.buildup = options.buildup;
  if (!exposure::parseIrradiation(options.irradiation, geometry.irradiation)) {
    throw InputError("exposure: \"" + options.irradiation +
                     "\" is not an irradiation geometry (ap, pa, llat, rlat, rot, iso)");
  }
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
  search.warnings = &std::cerr;
  resolved = locateStore(search);
  return NuclearData::open(resolved);
}

// One line describing how an exposure -- or a photon fluence -- was computed. Empty for
// activity, which needs no model beyond the decay constants, and for photon STRENGTH, which
// ignores the geometry by definition: a line under a number that does not use it would imply
// a model the number was not computed with.
std::string describeGeometry(const CommonOptions& options, Metric metric, Unit unit) {
  const bool usesGeometry =
      metric == Metric::Exposure || (metric == Metric::Photon && isFluenceUnit(unit));
  if (!usesGeometry) {
    return {};
  }
  char distance[32];
  std::snprintf(distance, sizeof(distance), "%.3g", options.distanceM);

  std::string text;
  // Which QUANTITY a sievert is cannot be left to the unit alone: the same metric in gray is
  // air kerma and in sievert is effective dose to a person, computed for a stated orientation.
  // A report that named neither would leave the reader to assume, and the assumption people
  // arrive with is the one this replaced.
  if (isEffectiveDoseUnit(unit)) {
    exposure::Irradiation irradiation = exposure::Irradiation::AP;
    exposure::parseIrradiation(options.irradiation, irradiation);
    text += "ICRP 116 effective dose, ";
    text += exposure::irradiationName(irradiation);
    text += ", ";
  }
  text += "point source at ";
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
                         const Inventory& inventory, const ResponseTable& table,
                         const PackHolder& packs) {
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
  context.pack = packs.describe();
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

int runAttribute(const CommonOptions& options, const char* argv0) {
  std::string storePath;
  const NuclearData data = openStore(options, argv0, storePath);
  const Inventory inventory = loadInventory(options, data);
  const PackHolder packs(options, data, inventory);

  const std::vector<double> times = timesFrom(options);
  if (times.size() != 1) {
    throw InputError(
        "time: attribute reports one time; pass a single --at (the shares partition the "
        "response at that instant, and a grid would be a different table per point)");
  }

  ResponseSpec spec;
  spec.metric = metricFrom(options);
  // Passed through rather than forced to Nuclide. attributeToSeed() refuses every other
  // aggregate with a message saying why -- an inventory row names a nuclide, so that is what a
  // share can name -- and pinning it here instead answered `--by element` with a nuclide
  // attribution, which is a different table than the one that was asked for.
  spec.aggregate = aggregateFrom(options.aggregate);
  spec.unit = requireUnit(options.unit, spec.metric, Domain::Instant);
  spec.pack = packs.resolved();
  spec.geometry = geometryFrom(options);

  RankRequest request;
  request.topN = options.topN;
  request.coverage = options.coverage;
  request.minFraction = options.minFraction;
  // Pins resolve against the SEED rather than against a response table: what can be followed
  // here is a nuclide that was seeded, which is a different set from what a ranking carries.
  for (const std::string& pin : options.pins) {
    request.pinned.push_back(requireSeedPin(data, inventory, pin));
  }

  const SeedAttribution attribution =
      attributeToSeed(data, inventory, times.front(), spec, request, decayOptionsFrom(options));

  ReportFormat format = ReportFormat::Text;
  parseReportFormat(options.format, format);
  OutputStream out(options.output);
  ReportContext context;
  context.storePath = storePath;
  context.storeLibrary = data.provenance().library;
  context.storeCreatedUtc = data.provenance().createdUtc;
  context.storeNuclideCount = data.size();
  context.seedProvenance = attribution.seedProvenance;
  context.geometry = describeGeometry(options, spec.metric, spec.unit);
  writeAttribution(out.get(), attribution, context, format);
  return 0;
}

int runRank(const CommonOptions& options, const char* argv0) {
  std::string storePath;
  const NuclearData data = openStore(options, argv0, storePath);
  const Inventory inventory = loadInventory(options, data);
  const PackHolder packs(options, data, inventory);

  const std::vector<double> times = timesFrom(options);
  const DecayResult result = decay(data, inventory, times, decayOptionsFrom(options));

  ResponseSpec spec;
  spec.metric = metricFrom(options);
  spec.aggregate = aggregateFrom(options.aggregate);
  spec.unit = requireUnit(options.unit, spec.metric, Domain::Instant);
  spec.pack = packs.resolved();
  spec.geometry = geometryFrom(options);
  const ResponseTable table = buildResponse(data, result, spec);

  const std::vector<Ranking> rankings = rankAll(table, rankRequestFrom(options, table));

  ReportFormat format = ReportFormat::Text;
  parseReportFormat(options.format, format);
  OutputStream out(options.output);
  ReportContext context = contextFor(data, storePath, inventory, table, packs);
  context.geometry = describeGeometry(options, spec.metric, spec.unit);
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
  const PackHolder packs(options, data, inventory);

  ResponseSpec spec;
  spec.metric = metricFrom(options);
  spec.aggregate = aggregateFrom(options.aggregate);
  spec.unit = requireUnit(options.unit, spec.metric, Domain::Interval);
  spec.pack = packs.resolved();
  spec.geometry = geometryFrom(options);

  // A context per interval. Each interval is solved separately over its own index space, so
  // the set of contributors carrying unmodelled continuum is a property of that interval --
  // building one context from the last table footnoted every ranking with the last interval's
  // emitters, which need not appear in the ranking they annotate.
  const std::string geometry = describeGeometry(options, spec.metric, spec.unit);
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
    contexts.push_back(contextFor(data, storePath, inventory, table, packs));
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
// Options for `forecast`. Only the refinement knobs, which nothing else in CommonOptions
// wants: a boundary is a located event, and how tightly to place one is a choice.
struct ForecastOptions {
  bool refine = false;
  double tolerance = 1.0e-6;
};

// Who leads, and when it changes. The same response table `rank` builds, read down the time
// axis instead of across it.
int runForecast(const CommonOptions& options, const ForecastOptions& forecast, const char* argv0) {
  std::string storePath;
  const NuclearData data = openStore(options, argv0, storePath);
  const Inventory inventory = loadInventory(options, data);
  const PackHolder packs(options, data, inventory);

  const std::vector<double> times = timesFrom(options);
  if (times.size() < 2) {
    throw InputError("time: a forecast needs a range of times -- try --times 1h:100y:log:60");
  }
  const DecayResult result = decay(data, inventory, times, decayOptionsFrom(options));

  ResponseSpec spec;
  spec.metric = metricFrom(options);
  spec.aggregate = aggregateFrom(options.aggregate);
  spec.unit = requireUnit(options.unit, spec.metric, Domain::Instant);
  spec.pack = packs.resolved();
  spec.geometry = geometryFrom(options);
  const ResponseTable table = buildResponse(data, result, spec);

  // Refinement re-solves inside each boundary's bracket instead of interpolating across it.
  // The evaluator outlives the call that uses it, which is all the lifetime this needs.
  EventTolerance tolerance;
  tolerance.relative = forecast.tolerance;
  std::optional<ResponseEvaluator> evaluator;
  if (forecast.refine) {
    evaluator.emplace(data, inventory, spec, decayOptionsFrom(options));
  }
  const std::vector<DominanceWindow> windows =
      evaluator ? dominanceWindows(table, *evaluator, /*minSamples=*/2, tolerance)
                : dominanceWindows(table);
  // topN doubles as the depth of "ever near the top"; 0 means no limit, which for a forecast
  // would print every nuclide in the chain, so it falls back to a readable default.
  const std::vector<std::int64_t> pins = pinsFrom(options, table);
  const std::vector<RankTrack> tracks =
      unionTopN(table, options.topN > 0 ? options.topN : 10, pins);

  ReportFormat format = ReportFormat::Text;
  parseReportFormat(options.format, format);
  ReportContext context = contextFor(data, storePath, inventory, table, packs);
  context.geometry = describeGeometry(options, spec.metric, spec.unit);

  OutputStream out(options.output);
  writeForecast(out.get(), windows, tracks, table, context, format);
  return 0;
}

int runDecay(const CommonOptions& options, const char* argv0) {
  std::string storePath;
  const NuclearData data = openStore(options, argv0, storePath);
  const Inventory inventory = loadInventory(options, data);
  const PackHolder packs(options, data, inventory);

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
  // The fourth gap, and the only one that costs atoms rather than photons: a spontaneous-fission
  // branch with no yield set removes the parent and produces nothing. Counted by the store
  // itself so the number here and the one staging warns about are the same number.
  const std::vector<Zai> sfWithoutYields = data.spontaneousFissionWithoutYields();

  std::printf("\nunstable nuclides:              %d\n", unstable);
  std::printf("  with discrete photon lines:   %d\n", withLines);
  std::printf("    of those, >5%% continuum:    %d\n", partialContinuum);
  std::printf("    of those, clamped lines:    %d\n", withClampedLines);
  std::printf("  emit photons, no spectrum:    %d\n", noSpectrumAtAll);
  std::printf("  SF branch, no yields staged:  %zu\n", sfWithoutYields.size());
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
  if (!sfWithoutYields.empty()) {
    std::printf(
        "\n%zu nuclides have a spontaneous-fission branch but no fission-yield set of any\n"
        "energy. The branch removes atoms from the parent and produces nothing in their\n"
        "place, so a chain passing through one of them does not conserve atoms. Negligible\n"
        "for a fission-product source; worth knowing for an actinide inventory. Staging the\n"
        "SFY sublibrary (--sfy-dir) closes the gap. Affected:\n   ",
        sfWithoutYields.size());
    constexpr std::size_t kMaxNamed = 8;
    const std::size_t named = std::min(sfWithoutYields.size(), kMaxNamed);
    for (std::size_t k = 0; k < named; ++k) {
      std::printf("%s%s", k == 0 ? " " : ", ", formatNuclideName(sfWithoutYields[k]).c_str());
    }
    if (sfWithoutYields.size() > named) {
      std::printf(", and %zu more", sfWithoutYields.size() - named);
    }
    std::printf("\n");
  }
  if (!data.hasPhotonLines()) {
    std::printf(
        "\nNote: this store carries no photon lines, so exposure and photon metrics are "
        "unavailable.\nStage from ENDF decay tapes to add them.\n");
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

// --- locating events on a curve ----------------------------------------------

// Options for `when`. Separate from CommonOptions because nothing else takes them, and folding
// a level into the options every subcommand shares would put `--level` in the help of six
// commands that ignore it.
struct WhenOptions {
  double level = 0.0;
  bool hasLevel = false;
  std::string of;
  std::string ratio;
  bool peaks = false;

  // The task form: the curve becomes what a job of this length accrues, plotted against when
  // the job STARTS, and --level becomes the budget it has to fit inside.
  std::string task;

  // Which side of the level the windows are wanted on. Empty takes the default for the curve
  // in force, and the two defaults differ because the questions do: a rate is asked about
  // above a level -- how long am I taking this -- and a task below a budget -- when may I
  // start.
  std::string windows;

  bool refine = false;
  double tolerance = 1.0e-6;
};

std::string trimmedText(std::string_view text) {
  std::size_t first = text.find_first_not_of(" \t");
  if (first == std::string_view::npos) {
    return {};
  }
  std::size_t last = text.find_last_not_of(" \t");
  return std::string(text.substr(first, last - first + 1));
}

// requirePin does the naming half -- accepting "Cs-137", "A=140" or "Cs" against whatever
// aggregate the table was built with -- and returns a key. A series needs the column.
int columnFor(const ResponseTable& table, const std::string& name) {
  const std::int64_t key = requirePin(table, name);
  for (int c = 0; c < table.contributorCount(); ++c) {
    if (table.contributors[static_cast<std::size_t>(c)].key == key) {
      return c;
    }
  }
  // requirePin refuses a name this table does not carry, so arriving here is a bug rather than
  // bad input, and the exit code should say so.
  throw NusiftError("trajectory: \"" + name + "\" resolved to a key the table does not hold");
}

int runWhen(const CommonOptions& options, const WhenOptions& when, const char* argv0) {
  std::string storePath;
  const NuclearData data = openStore(options, argv0, storePath);
  const Inventory inventory = loadInventory(options, data);
  const PackHolder packs(options, data, inventory);

  const std::vector<double> times = timesFrom(options);
  if (times.size() < 2) {
    throw InputError("time: an event has to sit between two samples -- try --times 1h:100y:log:60");
  }

  // A task accrues over its window, so its curve is in INTERVAL units -- roentgen rather than
  // roentgen per hour, decays rather than becquerel. Resolving the unit against the wrong
  // domain is how a budget silently becomes a rate.
  const bool isTask = !when.task.empty();
  const double taskSeconds = isTask ? parseDuration(when.task) : 0.0;
  const Domain domain = isTask ? Domain::Interval : Domain::Instant;

  ResponseSpec spec;
  spec.metric = metricFrom(options);
  spec.aggregate = aggregateFrom(options.aggregate);
  spec.unit = requireUnit(options.unit, spec.metric, domain);
  spec.pack = packs.resolved();
  spec.geometry = geometryFrom(options);

  const DecayOptions decayOptions = decayOptionsFrom(options);
  EventTolerance tolerance;
  tolerance.relative = when.tolerance;

  EventReport report;
  report.metric = metricName(spec.metric);
  report.gridStartSeconds = times.front();
  report.gridEndSeconds = times.back();
  report.gridPoints = static_cast<int>(times.size());
  report.hasLevel = when.hasLevel;
  report.level = when.level;

  // Windows default to the side the curve is usually asked about, and --windows overrides.
  bool windowsBelow = isTask;
  if (when.windows == "below") {
    windowsBelow = true;
  } else if (when.windows == "above") {
    windowsBelow = false;
  }

  // Both paths end with the same searches over the same kind of series. What differs is what
  // the series IS: a sampled rate, or a curve every point of which is an interval integral.
  EventSeries series;
  ResponseTable contextTable;

  // Declared out here because a series built with an evaluator holds a reference to it, and
  // the searches below run after this scope would have ended. Unused without --refine, which
  // is the one thing it must not silently be.
  std::optional<ResponseEvaluator> evaluator;

  if (isTask) {
    if (!when.of.empty() || !when.ratio.empty()) {
      throw InputError(
          "trajectory: --task follows what the whole inventory accrues over the window. A "
          "single contributor's share of a task is a ranking question -- try `integrate`");
    }
    series = taskSeries(data, inventory, spec, times, taskSeconds, decayOptions, when.refine);

    // The caveats a report carries -- which emitters have unmodelled photon continua -- belong
    // to the inventory and not to a particular window, so the first one speaks for all of them.
    std::vector<std::int64_t> keys;
    const std::vector<double> integral = intervalIntegral(
        data, inventory, times.front(), times.front() + taskSeconds, &keys, decayOptions);
    contextTable = buildIntervalResponse(data, keys, integral, times.front(),
                                         times.front() + taskSeconds, spec);

    report.curve = "a " + formatDuration(taskSeconds) + " task, by when it starts";
    report.unit = unitName(spec.unit);
  } else {
    const DecayResult result = decay(data, inventory, times, decayOptions);
    contextTable = buildResponse(data, result, spec);
    const ResponseTable& table = contextTable;

    if (when.refine) {
      evaluator.emplace(data, inventory, spec, decayOptions);
    }

    if (!when.ratio.empty()) {
      const std::size_t slash = when.ratio.find('/');
      if (slash == std::string::npos) {
        throw InputError(
            "trajectory: --ratio names two contributors as A/B, e.g. --ratio Zr-95/Nb-95");
      }
      const std::string top = trimmedText(std::string_view(when.ratio).substr(0, slash));
      const std::string bottom = trimmedText(std::string_view(when.ratio).substr(slash + 1));
      const int numerator = columnFor(table, top);
      const int denominator = columnFor(table, bottom);
      series = evaluator ? ratioSeries(table, numerator, denominator, *evaluator)
                         : ratioSeries(table, numerator, denominator);
      report.curve = top + " / " + bottom;
      // A ratio of two quantities in the same unit is dimensionless, and printing the unit
      // would say it is a number of becquerel when it is a number of times.
      report.unit.clear();
    } else if (!when.of.empty()) {
      const int column = columnFor(table, when.of);
      series = evaluator ? contributorSeries(table, column, *evaluator)
                         : contributorSeries(table, column);
      report.curve = table.labels[static_cast<std::size_t>(column)];
      report.unit = unitName(table.unit);
    } else {
      series = evaluator ? totalSeries(table, *evaluator) : totalSeries(table);
      report.curve = "the total";
      report.unit = unitName(table.unit);
    }
  }

  if (when.hasLevel) {
    report.events = crossings(series, when.level, tolerance);
    report.windows = windowsBelow ? nusift::windowsBelow(series, when.level, tolerance)
                                  : windowsAbove(series, when.level, tolerance);
    report.windowsBelowLevel = windowsBelow;
  }
  // With no level there is nothing to cross, so the turns are the whole answer rather than an
  // extra. With one they are an extra, and only if asked for.
  if (when.peaks || !when.hasLevel) {
    const std::vector<TrajectoryEvent> turns = extrema(series, tolerance);
    report.events.insert(report.events.end(), turns.begin(), turns.end());
    std::sort(report.events.begin(), report.events.end(),
              [](const TrajectoryEvent& a, const TrajectoryEvent& b) {
                return a.timeSeconds < b.timeSeconds;
              });
  }

  ReportFormat format = ReportFormat::Text;
  parseReportFormat(options.format, format);
  ReportContext context = contextFor(data, storePath, inventory, contextTable, packs);
  context.geometry = describeGeometry(options, spec.metric, spec.unit);

  OutputStream out(options.output);
  writeEvents(out.get(), report, context, format);
  return 0;
}

// --- how much is allowed ------------------------------------------------------

struct AllowableOptions {
  double limit = 0.0;
  std::string limitName;
  int limitingCount = 3;
};

int runAllowable(const CommonOptions& options, const AllowableOptions& allow, const char* argv0) {
  std::string storePath;
  const NuclearData data = openStore(options, argv0, storePath);
  const Inventory inventory = loadInventory(options, data);
  const PackHolder packs(options, data, inventory);

  const std::vector<double> times = timesFrom(options);
  const DecayResult result = decay(data, inventory, times, decayOptionsFrom(options));

  // One criterion, taking the metric, unit and geometry the command was given. Several criteria
  // over DIFFERENT metrics is where this capability earns its keep -- and is where the binding
  // one can change over time -- but expressing that needs a limits file this front end does not
  // have yet, so the library and the Python binding take a list and the CLI takes one.
  Criterion criterion;
  criterion.spec.metric = metricFrom(options);
  criterion.spec.aggregate = aggregateFrom(options.aggregate);
  criterion.spec.unit = requireUnit(options.unit, criterion.spec.metric, Domain::Instant);
  criterion.spec.pack = packs.resolved();
  criterion.spec.geometry = geometryFrom(options);
  criterion.limit = allow.limit;
  criterion.name = allow.limitName.empty()
                       ? std::string(metricName(criterion.spec.metric)) + " limit"
                       : allow.limitName;

  const std::vector<Criterion> criteria = {criterion};
  const std::vector<AllowableScale> scaled =
      allowableScale(data, result, criteria, allow.limitingCount);

  ReportFormat format = ReportFormat::Text;
  parseReportFormat(options.format, format);
  const ResponseTable table = buildResponse(data, result, criterion.spec);
  ReportContext context = contextFor(data, storePath, inventory, table, packs);
  context.geometry = describeGeometry(options, criterion.spec.metric, criterion.spec.unit);

  OutputStream out(options.output);
  writeAllowable(out.get(), scaled, criteria, context, format);

  // The single number a reader of this table came for, in the format a person reads. CSV and
  // JSON carry the curve it comes from, and a consumer wanting the instant locates it the same
  // way rather than being handed a derived field it cannot check.
  if (format == ReportFormat::Text && scaled.size() >= 2) {
    try {
      const std::optional<TrajectoryEvent> fits = firstCrossing(scaleSeries(scaled), 1.0);
      if (fits.has_value() && fits->kind == EventKind::Rising) {
        out.get() << "\n  The inventory as it stands first fits every limit at "
                  << formatDuration(fits->timeSeconds) << ", located within "
                  << formatDuration(fits->locatedToSeconds) << ".\n";
      }
    } catch (const InputError&) {
      // The scale is not one searchable curve over this grid: unbounded stretches, or too few
      // bounded samples. The table above already reports everything that was computed, and
      // inventing a date from a curve with a hole in it is exactly what scaleSeries refuses.
    }
  }
  return 0;
}

// --- counterfactual interventions ---------------------------------------------

struct IntervenOptions {
  std::vector<std::string> removals;  // "Cs", "Sr-90=0.5"
  std::string removeAt = "0";
  bool together = false;
};

// "Cs" or "Sr-90=0.5": a selector, optionally with the fraction removed. Without one the
// removal is complete, which is the upper bound on any real separation and the honest default
// for a question about what a separation could buy at most.
Removal parseRemoval(const std::string& text) {
  Removal removal;
  const std::size_t equals = text.rfind('=');
  // "Z=55" and "A=137" are selectors that contain their own '=', so a suffix only counts as a
  // fraction when it actually parses as one.
  if (equals != std::string::npos) {
    const std::string tail = text.substr(equals + 1);
    try {
      std::size_t used = 0;
      const double fraction = std::stod(tail, &used);
      if (used == tail.size()) {
        removal.selector = text.substr(0, equals);
        removal.fraction = fraction;
        return removal;
      }
    } catch (const std::exception&) {
      // Not a number, so the '=' belongs to the selector.
    }
  }
  removal.selector = text;
  removal.fraction = 1.0;
  return removal;
}

int runIntervene(const CommonOptions& options, const IntervenOptions& intervene,
                 const char* argv0) {
  std::string storePath;
  const NuclearData data = openStore(options, argv0, storePath);
  const Inventory inventory = loadInventory(options, data);
  const PackHolder packs(options, data, inventory);

  const std::vector<double> times = timesFrom(options);
  if (times.size() != 1) {
    throw InputError(
        "time: intervene reports one response time; pass a single --at (the benefit is measured "
        "against the response at that instant)");
  }
  const double removeAt = parseDuration(intervene.removeAt);

  ResponseSpec spec;
  spec.metric = metricFrom(options);
  spec.aggregate = aggregateFrom(options.aggregate);
  spec.unit = requireUnit(options.unit, spec.metric, Domain::Instant);
  spec.pack = packs.resolved();
  spec.geometry = geometryFrom(options);

  std::vector<Intervention> interventions;
  if (intervene.together) {
    // One intervention doing everything: the combined separation, rather than a comparison of
    // separations.
    Intervention combined;
    for (const std::string& text : intervene.removals) {
      const Removal removal = parseRemoval(text);
      combined.name += combined.name.empty() ? "" : " + ";
      combined.name += removal.selector;
      combined.removals.push_back(removal);
    }
    interventions.push_back(std::move(combined));
  } else {
    // The default is a COMPARISON: each removal is its own alternative, measured against the
    // same baseline. That is the question "which separation is worth doing", which is the one
    // people arrive with; --together answers "what does doing all of them buy".
    for (const std::string& text : intervene.removals) {
      const Removal removal = parseRemoval(text);
      Intervention one;
      one.name = removal.selector;
      one.removals.push_back(removal);
      interventions.push_back(std::move(one));
    }
  }

  const InterventionStudy study = compareInterventions(
      data, inventory, removeAt, times.front(), spec, interventions, decayOptionsFrom(options));

  ReportFormat format = ReportFormat::Text;
  parseReportFormat(options.format, format);
  ReportContext context;
  context.storePath = storePath;
  context.storeLibrary = data.provenance().library;
  context.storeCreatedUtc = data.provenance().createdUtc;
  context.storeNuclideCount = data.size();
  context.seedProvenance = study.seedProvenance;
  context.geometry = describeGeometry(options, spec.metric, spec.unit);

  OutputStream out(options.output);
  writeInterventions(out.get(), study, context, format);
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
  CLI::App* spectrumCmd = app.add_subcommand(
      "spectrum", "Top contributing photon lines (by line; exposure by default)");
  addCommonOptions(spectrumCmd, spectrumOptions, /*wantsTimes=*/true, /*wantsIntervals=*/false);

  CommonOptions attributeOptions;
  CLI::App* attributeCmd =
      app.add_subcommand("attribute", "Which SEEDED nuclides a response is riding on, at one time");
  addCommonOptions(attributeCmd, attributeOptions, /*wantsTimes=*/true, /*wantsIntervals=*/false);

  CommonOptions forecastOptions;
  ForecastOptions forecastExtra;
  CLI::App* forecastCmd =
      app.add_subcommand("forecast", "Who dominates, and over which time windows");
  addCommonOptions(forecastCmd, forecastOptions, /*wantsTimes=*/true, /*wantsIntervals=*/false);
  forecastCmd->add_flag("--refine", forecastExtra.refine,
                        "Place each boundary by re-solving inside its bracket, rather than "
                        "interpolating across it (costs a few solves per boundary)");
  forecastCmd
      ->add_option("--tolerance", forecastExtra.tolerance,
                   "Relative tolerance for --refine (default 1e-6)")
      ->check(CLI::PositiveNumber);

  CommonOptions whenOptions;
  WhenOptions whenExtra;
  CLI::App* whenCmd =
      app.add_subcommand("when", "When a curve crosses a level, peaks, or holds above it");
  addCommonOptions(whenCmd, whenOptions, /*wantsTimes=*/true, /*wantsIntervals=*/false);
  CLI::Option* levelOption =
      whenCmd->add_option("--level", whenExtra.level,
                          "Level to cross, in the metric's units; without one only turns are "
                          "reported");
  whenCmd->add_option("--of", whenExtra.of,
                      "Follow one contributor instead of the total, e.g. --of Cs-137");
  whenCmd->add_option("--ratio", whenExtra.ratio,
                      "Follow a ratio of two contributors, e.g. --ratio Zr-95/Nb-95");
  whenCmd->add_flag("--peaks", whenExtra.peaks, "Report turns as well as crossings");
  whenCmd->add_option("--task", whenExtra.task,
                      "Follow what a job of this length accrues, against when it starts "
                      "(e.g. 1h). --level is then the budget it has to fit inside");
  whenCmd
      ->add_option("--windows", whenExtra.windows,
                   "Which side of the level to report windows on; defaults to above for a "
                   "rate and below for a --task budget")
      ->check(CLI::IsMember({"above", "below"}));
  whenCmd->add_flag("--refine", whenExtra.refine,
                    "Place each event by re-solving inside its bracket, rather than "
                    "interpolating across it (costs a few solves per event)");
  whenCmd
      ->add_option("--tolerance", whenExtra.tolerance,
                   "Relative tolerance for --refine (default 1e-6)")
      ->check(CLI::PositiveNumber);

  CommonOptions allowableOptions;
  AllowableOptions allowableExtra;
  CLI::App* allowableCmd = app.add_subcommand(
      "allowable", "By what factor this inventory can be scaled before a limit binds");
  addCommonOptions(allowableCmd, allowableOptions, /*wantsTimes=*/true, /*wantsIntervals=*/false);
  allowableCmd->add_option("--limit", allowableExtra.limit, "The limit, in the metric's units")
      ->required()
      ->check(CLI::PositiveNumber);
  allowableCmd->add_option("--limit-name", allowableExtra.limitName,
                           "What the limit is, named in the report (e.g. \"A2 transport\")");
  allowableCmd->add_option("--limiting", allowableExtra.limitingCount,
                           "How many contributors to name as driving the binding criterion");

  CommonOptions intervenOptions;
  IntervenOptions intervenExtra;
  CLI::App* intervenCmd = app.add_subcommand(
      "intervene", "What removing something on a given date is worth to a later response");
  addCommonOptions(intervenCmd, intervenOptions, /*wantsTimes=*/true, /*wantsIntervals=*/false);
  intervenCmd
      ->add_option("--remove", intervenExtra.removals,
                   "Take this out: a nuclide, an element, Z=55 or A=137, optionally with the "
                   "fraction removed (Cs, Sr-90=0.5). Repeatable; each is a separate "
                   "alternative unless --together")
      ->required();
  intervenCmd->add_option("--remove-at", intervenExtra.removeAt,
                          "When the removal happens (default now)");
  intervenCmd->add_flag("--together", intervenExtra.together,
                        "Apply every --remove as one combined intervention rather than "
                        "comparing them");

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
    if (attributeCmd->parsed()) {
      return runAttribute(attributeOptions, argv0);
    }
    if (forecastCmd->parsed()) {
      return runForecast(forecastOptions, forecastExtra, argv0);
    }
    if (whenCmd->parsed()) {
      // Whether a level was GIVEN, not whether it is non-zero: zero is a perfectly good level
      // to ask about, and a default of zero would otherwise be indistinguishable from one.
      whenExtra.hasLevel = levelOption->count() > 0;
      return runWhen(whenOptions, whenExtra, argv0);
    }
    if (allowableCmd->parsed()) {
      return runAllowable(allowableOptions, allowableExtra, argv0);
    }
    if (intervenCmd->parsed()) {
      return runIntervene(intervenOptions, intervenExtra, argv0);
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
