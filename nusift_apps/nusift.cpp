// nusift -- the command-line driver.
//
// Exit codes:
//   0  success
//   1  runtime failure (a solve diverged, a store is unreadable)
//   2  bad input (an unparseable argument, a malformed inventory row)
//
#include <CLI/CLI.hpp>

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstdio>
#include <exception>
#include <fstream>
#include <iostream>
#include <memory>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "nusift/core/error.hpp"
#include "nusift/core/nuclide_name.hpp"
#include "nusift/engine/decay_engine.hpp"
#include "nusift/engine/inventory.hpp"
#include "nusift/engine/reconcile.hpp"
#include "nusift/exposure/air_coefficients.hpp"
#include "nusift/exposure/point_source.hpp"
#include "nusift/io/inventory_io.hpp"
#include "nusift/io/report.hpp"
#include "nusift/io/source_report.hpp"
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
#include "nusift/triage/spectrum.hpp"
#include "nusift/triage/task_plan.hpp"
#include "nusift/triage/triage_set.hpp"
#include "nusift/triage/uncertainty.hpp"
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
  // The area, volume or mass a concentration pack's inventory is spread through. A bare number:
  // which unit it is in is the pack's to say, and the report prints it back so a volume given
  // to a per-square-metre pack is visible rather than silently accepted.
  double packExtent = 0.0;

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

// Where the inventory comes from, and nothing else. Split out of addCommonOptions because
// `source` needs exactly this much of it: an inventory has to be named the same way on every
// command, but a command that computes an emission spectrum has no use for a distance, a
// metric or a ranking depth, and offering flags it then ignores is how a user comes to believe
// a number was computed with a geometry it never saw.
void addInventoryOptions(CLI::App* app, CommonOptions& options) {
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
}

// How the solve is run, as opposed to what is asked of it. Shared for the same reason: the
// cost and accuracy knobs mean the same thing on every command that decays anything.
void addSolverOptions(CLI::App* app, CommonOptions& options) {
  app->add_option("--cram-order", options.cramOrder, "CRAM order: 16 for screening, 48 default")
      ->check(CLI::IsMember({16, 48}));
  app->add_flag("--no-prune", options.noPrune,
                "Solve the whole chain instead of the seed's forward closure");
  app->add_option("--threads", options.threads,
                  "Worker threads for the per-time solves; 0 uses every core")
      ->check(CLI::NonNegativeNumber);
}

void addCommonOptions(CLI::App* app, CommonOptions& options, bool wantsTimes, bool wantsIntervals) {
  addInventoryOptions(app, options);

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
  app->add_option("--extent", options.packExtent,
                  "For a pack whose coefficients are per unit concentration: the volume, area "
                  "or mass the inventory is spread through, in the unit the pack declares")
      ->check(CLI::PositiveNumber);
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
  addSolverOptions(app, options);
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
    PackExtent extent;
    if (options.packExtent > 0.0) {
      extent.value = options.packExtent;
      // The unit is the pack's, so --extent takes a bare number and cannot disagree with it.
      // A pack that wants no extent refuses one, which is where a misplaced flag surfaces.
      extent.unit = pack_->provenance().per;
    }
    resolved_ = resolvePack(*pack_, data, seed, extent);
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
    // The extent is not a setting, it is half of what a concentration answer means: the same
    // inventory in ten times the volume gives a tenth the dose rate.
    if (resolved_ && resolved_->extent.value > 0.0) {
      char spread[64];
      std::snprintf(spread, sizeof(spread), "%.4g", resolved_->extent.value);
      text += ", spread through ";
      text += spread;
      text += " ";
      text += resolved_->extent.unit;
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

// Options for `source`. All of them are about the ENERGY axis, which is the axis that
// distinguishes this command from `spectrum`: one ranks the lines, the other resolves them onto
// a grid a transport code can read.
struct SourceOptions {
  int bins = 100;
  std::string scale = "linear";
  // keV on the way in, because that is the unit a decay photon is spoken in and the one the
  // line ranking already prints. Each deck converts on the way out to what its own code reads
  // -- MeV for MCNP, eV for OpenMC -- which is exactly the sort of conversion a user should not
  // be doing by hand at two in the morning.
  double minKeV = 0.0;
  double maxKeV = 0.0;
  std::string edgesKeV;
  std::string edgesFile;
  std::string format = "text";
};

// Explicit edges, from a comma-separated list or from a file of them. Both are in keV and both
// end up in the same place; the file exists because a real group structure has two hundred
// boundaries and does not belong on a command line.
std::vector<double> edgesFrom(const SourceOptions& source) {
  std::vector<double> edges;
  auto push = [&edges](const std::string& text) {
    if (text.empty()) {
      return;
    }
    try {
      std::size_t used = 0;
      const double value = std::stod(text, &used);
      if (used != text.size()) {
        throw std::invalid_argument("trailing");
      }
      edges.push_back(value * 1000.0);
    } catch (const std::exception&) {
      throw InputError("source: \"" + text + "\" is not an energy in keV");
    }
  };

  if (!source.edgesKeV.empty()) {
    std::string field;
    for (const char c : source.edgesKeV) {
      if (c == ',') {
        push(field);
        field.clear();
      } else if (!std::isspace(static_cast<unsigned char>(c))) {
        field += c;
      }
    }
    push(field);
  }
  if (!source.edgesFile.empty()) {
    std::ifstream file(source.edgesFile);
    if (!file) {
      throw InputError("source: cannot read edges from \"" + source.edgesFile + "\"");
    }
    std::string line;
    while (std::getline(file, line)) {
      const std::size_t hash = line.find('#');
      if (hash != std::string::npos) {
        line.erase(hash);
      }
      std::istringstream fields(line);
      std::string field;
      while (fields >> field) {
        push(field);
      }
    }
  }
  return edges;
}

BinningSpec binningFrom(const SourceOptions& source) {
  BinningSpec binning;
  binning.count = source.bins;
  binning.minEv = source.minKeV * 1000.0;
  binning.maxEv = source.maxKeV * 1000.0;
  binning.edgesEv = edgesFrom(source);
  if (!parseBinScale(source.scale, binning.scale)) {
    throw InputError("source: \"" + source.scale + "\" is not a bin scale (linear or log)");
  }
  // The refusal of edges-plus-a-range lives in resolveBinEdges, so this front end and the
  // binding cannot come to differ about which combination is meaningful.
  return binning;
}

// The binned photon emission of an inventory, as a report or as a transport code's source card.
//
// Deliberately NOT a mode of `spectrum`. That command ranks discrete lines, and its answer is
// read by a person deciding what to look at; this one hands a source term to another code,
// which will use whatever it is given without noticing what is missing from it. The two want
// different truncation, different formats and different warnings, and the one thing they must
// not share is the impression that they are the same answer at different resolutions.
int runSource(const CommonOptions& options, const SourceOptions& source, const char* argv0) {
  const bool hasInterval = !options.intervals.empty();
  if (hasInterval && !options.atTimes.empty()) {
    throw InputError(
        "time: a source is emitted either AT an instant or OVER an interval, and the two are "
        "different quantities -- photons per second against a count of photons. Give one");
  }
  if (!hasInterval && options.atTimes.size() != 1) {
    throw InputError("time: source takes exactly one --at, or one --interval T1,T2");
  }
  if (options.intervals.size() > 1) {
    throw InputError("time: source takes one --interval; a deck describes one source");
  }

  std::string storePath;
  const NuclearData data = openStore(options, argv0, storePath);
  const Inventory inventory = loadInventory(options, data);
  const BinningSpec binning = binningFrom(source);

  BinnedSpectrum spectrum;
  if (hasInterval) {
    const auto [t1, t2] = parseInterval(options.intervals.front());
    std::vector<std::int64_t> keys;
    const std::vector<double> integral =
        intervalIntegral(data, inventory, t1, t2, &keys, decayOptionsFrom(options));
    spectrum = binnedIntervalSpectrum(data, keys, integral, t1, t2, binning);
  } else {
    const std::vector<double> times{parseDuration(options.atTimes.front())};
    const DecayResult result = decay(data, inventory, times, decayOptionsFrom(options));
    spectrum = binnedSpectrum(data, result, 0, binning);
  }

  SourceFormat format = SourceFormat::Text;
  if (!parseSourceFormat(source.format, format)) {
    throw InputError("source: \"" + source.format +
                     "\" is not a format (text, csv, json, mcnp, openmc)");
  }

  // Built here rather than through contextFor(), which reads its continuum footnote off a
  // response table. There is no table: the spectrum carries its own flagged emitters, scanned
  // over the whole inventory for the same reason and reported the same way.
  ReportContext context;
  context.storePath = storePath;
  context.storeLibrary = data.provenance().library;
  context.storeCreatedUtc = data.provenance().createdUtc;
  context.storeNuclideCount = data.stagedCount();
  context.seedProvenance = inventory.provenance();

  OutputStream out(options.output);
  writeSourceSpectrum(out.get(), spectrum, context, format);
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

// --- what the assay's uncertainty does to the answer ---------------------------

// The one uncertainty question that needs no evaluated data. `attribute` prints the importance
// column already; this multiplies it by what the sheet said about each row and reports which
// measurement the error bar is actually resting on.
int runUncertainty(const CommonOptions& options, const char* argv0) {
  if (options.inventoryPath.empty()) {
    throw InputError(
        "uncertainty: give an inventory whose rows carry uncertainties with --inventory. A "
        "fission seed has no assay to propagate -- its uncertainty is in the yields, which is "
        "the evaluated-data half of the question and is not this");
  }

  std::string storePath;
  const NuclearData data = openStore(options, argv0, storePath);
  const PackHolder* packs = nullptr;

  InventoryReadOptions readOptions;
  readOptions.ignoreUnknown = options.ignoreUnknown;
  readOptions.warnings = &std::cerr;
  const DatedInventory dated = readInventoryDated(options.inventoryPath, data, readOptions);

  const bool anySigma =
      std::any_of(dated.groups.begin(), dated.groups.end(),
                  [](const AssayGroup& g) { return g.inventory.hasUncertainties(); });
  if (!anySigma) {
    throw InputError("uncertainty: no row in \"" + options.inventoryPath +
                     "\" states one. Add an `uncertainty` column -- absolute in the row's own "
                     "unit, or relative as 5% -- naming it in a header row if the file carries "
                     "no assay dates");
  }

  const double epoch = dated.dated ? latestAssayDate(dated.groups) : 0.0;
  const std::vector<double> times = timesFrom(options);
  if (times.size() != 1) {
    throw InputError("time: an error bar is on the response at ONE instant; give one --at");
  }

  // The pack is resolved against the MERGED seed, since a fold is a question about the whole
  // inventory rather than about one sheet.
  const Inventory merged = readInventory(options.inventoryPath, data, readOptions);
  const PackHolder holder(options, data, merged);
  packs = &holder;

  ResponseSpec spec;
  spec.metric = metricFrom(options);
  spec.aggregate = Aggregate::Nuclide;
  spec.unit = requireUnit(options.unit, spec.metric, Domain::Instant);
  spec.pack = holder.resolved();
  spec.geometry = geometryFrom(options);

  const ResponseUncertainty uncertainty = responseUncertainty(
      data, dated.groups, epoch, times.front(), spec, decayOptionsFrom(options));

  ReportFormat format = ReportFormat::Text;
  parseReportFormat(options.format, format);
  ReportContext context;
  context.storePath = storePath;
  context.storeLibrary = data.provenance().library;
  context.storeCreatedUtc = data.provenance().createdUtc;
  context.storeNuclideCount = data.stagedCount();
  context.seedProvenance = merged.provenance();
  context.geometry = describeGeometry(options, spec.metric, spec.unit);
  context.pack = packs->describe();

  OutputStream out(options.output);
  writeUncertainty(out.get(), uncertainty, context, format);
  return 0;
}

// --- a job with a shape --------------------------------------------------------

struct PlanOptions {
  std::vector<std::string> legs;
  double budget = 0.0;
};

// "name,duration,distance[,occupancy]", or "name,duration" for a break. Positional rather than
// keyed because a leg has exactly these parts and every one of them is required to say what the
// leg IS -- a keyword form would let a distance be forgotten, and a forgotten distance would
// silently take the default.
PlanLeg parseLeg(const std::string& text, const CommonOptions& options, bool isBreak) {
  std::vector<std::string> fields;
  std::string field;
  for (const char c : text) {
    if (c == ',') {
      fields.push_back(trimmedText(field));
      field.clear();
    } else {
      field += c;
    }
  }
  fields.push_back(trimmedText(field));

  PlanLeg leg;
  // A break needs no distance and is spelled as just a duration, since the whole point of it is
  // that nowhere is where you are.
  if (isBreak) {
    if (fields.size() == 1) {
      leg.name = "break";
      leg.durationSeconds = parseDuration(fields[0]);
    } else if (fields.size() == 2) {
      leg.name = fields[0];
      leg.durationSeconds = parseDuration(fields[1]);
    } else {
      throw InputError("plan: a break is \"duration\" or \"name,duration\", got \"" + text + "\"");
    }
    leg.occupancy = 0.0;
    return leg;
  }

  if (fields.size() < 3 || fields.size() > 4) {
    throw InputError(
        "plan: a leg is \"name,duration,distance\" with an optional occupancy "
        "(e.g. \"valve work,20m,0.8,0.6\"), got \"" +
        text + "\"");
  }
  leg.name = fields[0];
  if (leg.name.empty()) {
    throw InputError(
        "plan: a leg needs a name; \"leg 3 costs 60% of your dose\" is not an "
        "answer anybody can act on");
  }
  leg.durationSeconds = parseDuration(fields[1]);

  leg.geometry = geometryFrom(options);
  try {
    std::size_t used = 0;
    leg.geometry.distanceM = std::stod(fields[2], &used);
    if (used != fields[2].size()) {
      throw std::invalid_argument("trailing");
    }
  } catch (const std::exception&) {
    throw InputError("plan: \"" + fields[2] + "\" is not a distance in metres");
  }
  if (!(leg.geometry.distanceM > 0.0)) {
    throw InputError("plan: leg \"" + leg.name + "\" needs a positive distance");
  }
  if (fields.size() == 4) {
    try {
      std::size_t used = 0;
      leg.occupancy = std::stod(fields[3], &used);
      if (used != fields[3].size()) {
        throw std::invalid_argument("trailing");
      }
    } catch (const std::exception&) {
      throw InputError("plan: \"" + fields[3] + "\" is not an occupancy fraction");
    }
  }
  return leg;
}

// A job as it is actually done, rather than as one window at one distance. Which LEG costs the
// dose is the question a plan answers and a single interval cannot, because the dose is not
// distributed the way the durations are.
int runPlan(const CommonOptions& options, const PlanOptions& planOptions, const char* argv0) {
  if (planOptions.legs.empty()) {
    throw InputError("plan: give at least one --leg name,duration,distance");
  }

  std::string storePath;
  const NuclearData data = openStore(options, argv0, storePath);
  const Inventory inventory = loadInventory(options, data);
  const PackHolder packs(options, data, inventory);

  const std::vector<double> times = timesFrom(options);
  if (times.size() != 1) {
    throw InputError("time: a plan starts at one instant; give one --at");
  }

  // Legs are stored with their order marker stripped, in the order given: the sequence IS the
  // plan, so sorting or de-duplicating them would be rewriting the job.
  std::vector<PlanLeg> legs;
  for (const std::string& text : planOptions.legs) {
    const bool isBreak = text.rfind("break:", 0) == 0;
    legs.push_back(parseLeg(isBreak ? text.substr(6) : text, options, isBreak));
  }

  ResponseSpec spec;
  spec.metric = metricFrom(options);
  spec.aggregate = Aggregate::Nuclide;
  spec.unit = requireUnit(options.unit, spec.metric, Domain::Interval);
  spec.pack = packs.resolved();

  const TaskPlan plan = runTaskPlan(data, inventory, spec, times.front(), legs, planOptions.budget,
                                    decayOptionsFrom(options));

  ReportFormat format = ReportFormat::Text;
  parseReportFormat(options.format, format);
  ReportContext context;
  context.storePath = storePath;
  context.storeLibrary = data.provenance().library;
  context.storeCreatedUtc = data.provenance().createdUtc;
  context.storeNuclideCount = data.stagedCount();
  context.seedProvenance = inventory.provenance();
  // The distance is per leg and is in the table, so the model line names everything else the
  // geometry carries and deliberately not a distance that would be wrong for most rows.
  context.geometry = describeGeometry(options, spec.metric, spec.unit);
  const std::size_t at = context.geometry.find("point source at");
  if (at != std::string::npos) {
    context.geometry.replace(at, std::string("point source at").size(),
                             "point source, distance per leg;");
    const std::size_t metres = context.geometry.find(" m,", at);
    if (metres != std::string::npos) {
      context.geometry.erase(
          at + std::string("point source, distance per leg;").size(),
          metres + 3 - at - std::string("point source, distance per leg;").size());
    }
  }

  OutputStream out(options.output);
  writeTaskPlan(out.get(), plan, context, format);
  return 0;
}

// --- the smallest list that works everywhere -----------------------------------

struct ShortlistOptions {
  // Comma-separated, because several metrics at once is the point: a list chosen for activity
  // alone guarantees nothing about the exposure floor at a later time. Empty takes --metric.
  std::string metrics;
};

// `rank --coverage 0.95` is this question at one instant for one metric. This is the same
// question asked of every metric at every time simultaneously, which is a covering problem
// rather than a longer sort -- and the union of the per-time answers is not it.
int runShortlist(const CommonOptions& options, const ShortlistOptions& shortlist,
                 const char* argv0) {
  const double fraction = options.coverage > 0.0 ? options.coverage : 0.95;

  std::string storePath;
  const NuclearData data = openStore(options, argv0, storePath);
  const Inventory inventory = loadInventory(options, data);
  const PackHolder packs(options, data, inventory);

  const std::vector<double> times = timesFrom(options);
  const DecayResult result = decay(data, inventory, times, decayOptionsFrom(options));

  // Which metrics to be good for. A pack is one of them like any other, which is what makes
  // "hold 95% of the activity AND of the A2 index, always" a single question.
  std::vector<std::string> names;
  if (shortlist.metrics.empty()) {
    names.push_back(options.metric);
  } else {
    std::string field;
    for (const char c : shortlist.metrics) {
      if (c == ',') {
        names.push_back(trimmedText(field));
        field.clear();
      } else {
        field += c;
      }
    }
    names.push_back(trimmedText(field));
  }
  if (!options.packPath.empty() && shortlist.metrics.empty()) {
    names.clear();
  }

  // The tables outlive the requirements that point at them, which is why they are held here
  // rather than built inside the loop that fills the requirement list.
  std::vector<ResponseTable> tables;
  std::vector<CoverageRequirement> requirements;
  const Aggregate aggregate = aggregateFrom(options.aggregate);
  tables.reserve(names.size() + 1);

  for (const std::string& name : names) {
    ResponseSpec spec;
    spec.metric = metricFrom(name);
    spec.aggregate = aggregate;
    spec.unit = requireUnit("", spec.metric, Domain::Instant);
    spec.geometry = geometryFrom(options);
    tables.push_back(buildResponse(data, result, spec));
    CoverageRequirement requirement;
    requirement.fraction = fraction;
    requirement.label = std::string(metricName(spec.metric)) + " (" + unitName(spec.unit) + ")";
    requirements.push_back(requirement);
  }
  if (packs.resolved() != nullptr) {
    ResponseSpec spec;
    spec.metric = Metric::Pack;
    spec.aggregate = aggregate;
    spec.unit = Unit::PackDefined;
    spec.pack = packs.resolved();
    spec.geometry = geometryFrom(options);
    tables.push_back(buildResponse(data, result, spec));
    CoverageRequirement requirement;
    requirement.fraction = fraction;
    requirement.label = packs.resolved()->pack->provenance().name;
    requirements.push_back(requirement);
  }
  if (requirements.empty()) {
    throw InputError("shortlist: name at least one metric with --metrics or --pack");
  }
  // Filled after the vector has stopped growing: a pointer taken before a reallocation would
  // outlive the table it named.
  for (std::size_t i = 0; i < requirements.size(); ++i) {
    requirements[i].table = &tables[i];
  }

  const TriageSet set = robustTriageSet(requirements);

  ReportFormat format = ReportFormat::Text;
  parseReportFormat(options.format, format);
  ReportContext context = contextFor(data, storePath, inventory, tables.front(), packs);
  context.geometry = describeGeometry(
      options, requirements.size() == 1 ? metricFrom(names.front()) : Metric::Exposure,
      tables.front().unit);

  OutputStream out(options.output);
  writeTriageSet(out.get(), set, requirements, context, format);
  return 0;
}

// --- bringing several assays to one date --------------------------------------

struct ReconcileOptions {
  // An epoch later than every assay is allowed and is sometimes what is wanted -- "what will
  // these drums look like on the shipping date". An earlier one is refused by the library, for
  // the reason reconcile.hpp gives at length.
  std::string epoch;
  std::string writePath;
  std::string writeUnit = "atoms";
};

// Reading a dated inventory reconciles it whatever the command, because a date the tool ignored
// would be worse than one it refused. This command is the one that SHOWS the reconciliation --
// which sheet contributed what, how far each was carried -- and can write the merged inventory
// out as an ordinary undated one for the next run to seed from.
int runReconcile(const CommonOptions& options, const ReconcileOptions& reconcileOptions,
                 const char* argv0) {
  if (options.inventoryPath.empty()) {
    throw InputError("reconcile: give the assays to merge with --inventory");
  }

  std::string storePath;
  const NuclearData data = openStore(options, argv0, storePath);

  InventoryReadOptions readOptions;
  readOptions.ignoreUnknown = options.ignoreUnknown;
  readOptions.warnings = &std::cerr;
  const DatedInventory dated = readInventoryDated(options.inventoryPath, data, readOptions);
  if (!dated.dated) {
    throw InputError("reconcile: \"" + options.inventoryPath +
                     "\" carries no assay dates, so there is nothing to bring to a common one. "
                     "Add an `assayed` column of ISO dates (2024-03-15), or use `inventory "
                     "convert` to change its units");
  }

  const double epoch = reconcileOptions.epoch.empty() ? latestAssayDate(dated.groups)
                                                      : parseCalendarDate(reconcileOptions.epoch);
  const Reconciliation reconciled = reconcile(data, dated.groups, epoch, decayOptionsFrom(options));

  ReportFormat format = ReportFormat::Text;
  parseReportFormat(options.format, format);
  ReportContext context;
  context.storePath = storePath;
  context.storeLibrary = data.provenance().library;
  context.storeCreatedUtc = data.provenance().createdUtc;
  context.storeNuclideCount = data.stagedCount();
  context.seedProvenance = options.inventoryPath;

  OutputStream out(options.output);
  writeReconciliation(out.get(), reconciled, data, context, format);

  // The checkpoint. What comes out is an ORDINARY inventory with no dates in it, because at the
  // epoch there is only one date left and carrying it as a column would invite a second
  // reconciliation of an already-reconciled file.
  if (!reconcileOptions.writePath.empty()) {
    Quantity unit = Quantity::Atoms;
    if (!parseQuantity(reconcileOptions.writeUnit, unit)) {
      throw InputError("reconcile: \"" + reconcileOptions.writeUnit + "\" is not a unit");
    }
    std::ofstream file(reconcileOptions.writePath);
    if (!file) {
      throw InputError("reconcile: cannot write to \"" + reconcileOptions.writePath + "\"");
    }
    Inventory merged = reconciled.inventory;
    merged.setProvenance(options.inventoryPath + ", " +
                         std::to_string(reconciled.contributions.size()) +
                         " assays reconciled to " + formatCalendarDate(epoch));
    const std::string& path = reconcileOptions.writePath;
    if (path.size() >= 5 && path.compare(path.size() - 5, 5, ".json") == 0) {
      writeInventoryJson(file, merged, data, unit);
    } else {
      writeInventoryCsv(file, merged, data, unit);
    }
    std::cerr << "  wrote " << merged.size() << " nuclides to " << path << "\n";
  }
  return 0;
}

// --- how long may I stay ------------------------------------------------------

// Options for `stay`. The budget is spelled `--budget` rather than `--level` because it is not
// a level: `when --level` asks where a curve meets a value, and this asks how much of a
// quantity may be spent. Naming them the same would suggest the two searches are the same one.
struct StayOptions {
  double budget = 0.0;
  std::string maxStay = "24h";
  double tolerance = 1.0e-6;
};

// The converse of `when --task`: there the length is fixed and the start is the unknown, here
// the start is fixed and the LENGTH is. Two inversions of one exact integral, and neither is
// obtainable by sampling the other.
int runStay(const CommonOptions& options, const StayOptions& stayOptions, const char* argv0) {
  if (!(stayOptions.budget > 0.0)) {
    throw InputError("stay: give the budget a stay has to fit inside with --budget");
  }
  const double maxDuration = parseDuration(stayOptions.maxStay);

  std::string storePath;
  const NuclearData data = openStore(options, argv0, storePath);
  const Inventory inventory = loadInventory(options, data);
  const PackHolder packs(options, data, inventory);

  ResponseSpec spec;
  spec.metric = metricFrom(options);
  // Nuclide aggregate always: a stay time is a root on the TOTAL, and which nuclide is
  // responsible does not change how long the budget lasts. `rank` is where that question lives.
  spec.aggregate = Aggregate::Nuclide;
  // The interval domain is what makes this a budget rather than a rate. requireUnit refuses a
  // rate unit here with the same message it refuses one on `integrate`, which is the refusal a
  // user who typed Sv/h should see.
  spec.unit = requireUnit(options.unit, spec.metric, Domain::Interval);
  spec.pack = packs.resolved();
  spec.geometry = geometryFrom(options);

  EventTolerance tolerance;
  tolerance.relative = stayOptions.tolerance;

  StayReport report;
  report.metric = metricName(spec.metric);
  report.unit =
      spec.metric == Metric::Pack ? packs.resolved()->pack->provenance().unit : unitName(spec.unit);
  report.budget = stayOptions.budget;
  report.maxDurationSeconds = maxDuration;
  for (const double start : timesFrom(options)) {
    report.stays.push_back(stayTime(data, inventory, spec, start, stayOptions.budget, maxDuration,
                                    decayOptionsFrom(options), tolerance));
  }

  ReportFormat format = ReportFormat::Text;
  parseReportFormat(options.format, format);
  // A table is needed for the continuum footnote and nothing else, so it is built at the first
  // start time: which emitters carry an unmodelled continuum is a property of the inventory,
  // not of how long anyone stands next to it.
  const DecayResult result =
      decay(data, inventory, std::vector<double>{report.stays.front().startSeconds},
            decayOptionsFrom(options));
  ResponseSpec instant = spec;
  instant.unit = defaultUnit(spec.metric, Domain::Instant);
  ReportContext context =
      contextFor(data, storePath, inventory, buildResponse(data, result, instant), packs);
  context.geometry = describeGeometry(options, spec.metric, spec.unit);

  OutputStream out(options.output);
  writeStayTimes(out.get(), report, context, format);
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

  CommonOptions sourceOptions;
  SourceOptions sourceExtra;
  CLI::App* sourceCmd = app.add_subcommand(
      "source", "The binned photon emission spectrum, as a transport code's source definition");
  addInventoryOptions(sourceCmd, sourceOptions);
  sourceCmd->add_option("--at", sourceOptions.atTimes, "Cooling time (e.g. 30d, 1.5y)");
  sourceCmd->add_option("--interval", sourceOptions.intervals,
                        "Emit over a window T1,T2 instead: a count of photons rather than a "
                        "rate (e.g. 1h,30d)");
  sourceCmd->add_option("--bins", sourceExtra.bins, "Number of energy bins (default 100)")
      ->check(CLI::PositiveNumber);
  sourceCmd->add_option("--bin-scale", sourceExtra.scale, "Bin spacing: linear or log")
      ->check(CLI::IsMember({"linear", "lin", "log", "logarithmic"}));
  sourceCmd
      ->add_option("--min-kev", sourceExtra.minKeV,
                   "Lowest energy to bin, in keV; default is the softest line present")
      ->check(CLI::PositiveNumber);
  sourceCmd
      ->add_option("--max-kev", sourceExtra.maxKeV,
                   "Highest energy to bin, in keV; default is the hardest line present")
      ->check(CLI::PositiveNumber);
  sourceCmd->add_option("--edges-kev", sourceExtra.edgesKeV,
                        "Explicit bin boundaries in keV, comma separated");
  sourceCmd->add_option("--edges-file", sourceExtra.edgesFile,
                        "Explicit bin boundaries in keV, whitespace separated, # comments ok");
  sourceCmd->add_option("-o,--output", sourceOptions.output, "Write here instead of stdout");
  sourceCmd
      ->add_option("--format", sourceExtra.format,
                   "text, csv, json, mcnp (an SDEF card), or openmc (a Python snippet)")
      ->check(CLI::IsMember({"text", "csv", "json", "mcnp", "sdef", "openmc"}));
  addSolverOptions(sourceCmd, sourceOptions);

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

  CommonOptions uncertaintyOptions;
  CLI::App* uncertaintyCmd = app.add_subcommand(
      "uncertainty", "The error bar the assay's own uncertainties put on a response");
  addCommonOptions(uncertaintyCmd, uncertaintyOptions, /*wantsTimes=*/true,
                   /*wantsIntervals=*/false);

  CommonOptions planOptions;
  PlanOptions planExtra;
  CLI::App* planCmd = app.add_subcommand(
      "plan", "What a job costs leg by leg, each with its own distance and occupancy");
  addCommonOptions(planCmd, planOptions, /*wantsTimes=*/true, /*wantsIntervals=*/false);
  planCmd->add_option("--leg", planExtra.legs,
                      "A stretch of the job: name,duration,distance[,occupancy]. Repeatable, and "
                      "the ORDER is the plan (e.g. \"valve work,20m,0.8,0.6\"). A break is "
                      "\"break:15m\" -- spelled as a leg so it keeps its place in the sequence");
  planCmd
      ->add_option("--budget", planExtra.budget,
                   "Say whether the plan fits this, and where it runs out if not")
      ->check(CLI::PositiveNumber);

  CommonOptions shortlistOptions;
  ShortlistOptions shortlistExtra;
  CLI::App* shortlistCmd = app.add_subcommand(
      "shortlist", "The smallest set of contributors holding a coverage floor at every time");
  addCommonOptions(shortlistCmd, shortlistOptions, /*wantsTimes=*/true, /*wantsIntervals=*/false);
  shortlistCmd->add_option("--metrics", shortlistExtra.metrics,
                           "Be good for all of these at once, comma separated "
                           "(e.g. activity,exposure,photon). Defaults to --metric");

  CommonOptions reconcileOptions;
  ReconcileOptions reconcileExtra;
  CLI::App* reconcileCmd = app.add_subcommand(
      "reconcile", "Bring assays taken on different dates to one date, and merge them");
  addInventoryOptions(reconcileCmd, reconcileOptions);
  reconcileCmd->add_option("--epoch", reconcileExtra.epoch,
                           "Bring everything to this ISO date instead of the latest assay. "
                           "Must not be earlier than any assay");
  reconcileCmd->add_option("--write", reconcileExtra.writePath,
                           "Also write the merged inventory here, as a reusable inventory file");
  reconcileCmd->add_option("--write-units", reconcileExtra.writeUnit,
                           "Units for --write (atoms, g, Bq, Ci, ...)");
  reconcileCmd->add_option("-o,--output", reconcileOptions.output,
                           "Write the report here instead of stdout");
  reconcileCmd->add_option("--format", reconcileOptions.format, "text, csv, or json")
      ->check(CLI::IsMember({"text", "csv", "json"}));
  addSolverOptions(reconcileCmd, reconcileOptions);

  CommonOptions stayOptions;
  StayOptions stayExtra;
  CLI::App* stayCmd =
      app.add_subcommand("stay", "How long a stay beginning at a time can run on a dose budget");
  addCommonOptions(stayCmd, stayOptions, /*wantsTimes=*/true, /*wantsIntervals=*/false);
  stayCmd
      ->add_option("--budget", stayExtra.budget,
                   "The total a stay has to fit inside, in an INTERVAL unit (Sv, R, Gy, decays)")
      ->required()
      ->check(CLI::PositiveNumber);
  stayCmd->add_option("--max-stay", stayExtra.maxStay,
                      "The longest stay to consider (default 24h). A budget not spent inside it "
                      "is reported as not spent, rather than as a very long duration");
  stayCmd->add_option("--tolerance", stayExtra.tolerance, "Relative tolerance on the duration")
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
    if (uncertaintyCmd->parsed()) {
      return runUncertainty(uncertaintyOptions, argv0);
    }
    if (planCmd->parsed()) {
      return runPlan(planOptions, planExtra, argv0);
    }
    if (shortlistCmd->parsed()) {
      return runShortlist(shortlistOptions, shortlistExtra, argv0);
    }
    if (reconcileCmd->parsed()) {
      return runReconcile(reconcileOptions, reconcileExtra, argv0);
    }
    if (stayCmd->parsed()) {
      return runStay(stayOptions, stayExtra, argv0);
    }
    if (sourceCmd->parsed()) {
      return runSource(sourceOptions, sourceExtra, argv0);
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
