#include "nusift/triage/response.hpp"

#include <algorithm>
#include <cctype>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <string_view>

#include "nusift/core/element_symbols.hpp"
#include "nusift/core/error.hpp"
#include "nusift/core/nuclide_name.hpp"
#include "nusift/nucdata/nuclear_data.hpp"
#include "nusift/nucdata/photon_lines.hpp"
#include "nusift/units.hpp"

namespace nusift {
namespace {

constexpr const char* kModule = "response";
constexpr const char* kUnitsModule = "units";
constexpr const char* kPinModule = "pin";

// Above this fraction of a nuclide's photon energy sitting in an unmodeled continuum, the
// emitter is flagged. Named once because three places test it -- the two table assemblers and
// the caveats the adjoint path reports -- and a threshold that drifted between them would flag
// an emitter in one report and not in another for the same store.
constexpr double kUnmodeledContinuumFlag = 0.05;

// Every unit, in the order the help text lists them. The one place the set is enumerated, so
// parseUnit and the error message it raises cannot come to disagree about what exists.
constexpr Unit kAllUnits[] = {
    Unit::PackDefined,
    Unit::Becquerel,
    Unit::Curie,
    Unit::Decays,
    Unit::RoentgenPerHour,
    Unit::GrayPerHour,
    Unit::SievertPerHour,
    Unit::Roentgen,
    Unit::Gray,
    Unit::Sievert,
    Unit::PhotonsPerSecond,
    Unit::Photons,
    Unit::PhotonsPerSquareMeterPerSecond,
    Unit::PhotonsPerSquareMeter,
    Unit::Watt,
    Unit::Joule,
};

// Case-insensitive ASCII equality. Unit spellings are ASCII by construction -- they come from
// unitName() -- so there is no locale question here to get wrong.
bool equalsIgnoreCase(std::string_view a, std::string_view b) {
  if (a.size() != b.size()) {
    return false;
  }
  for (std::size_t i = 0; i < a.size(); ++i) {
    if (std::tolower(static_cast<unsigned char>(a[i])) !=
        std::tolower(static_cast<unsigned char>(b[i]))) {
      return false;
    }
  }
  return true;
}

// "Bq, Ci, decays" for activity; "R/h, Gy/h, Sv/h, R, Gy, Sv" for exposure. Built from the
// same predicate that enforces the pairing, so the message can never offer a unit the check
// would then refuse.
std::string spellingsFor(Metric metric) {
  std::string list;
  for (const Unit unit : kAllUnits) {
    if (!unitSuitsMetric(unit, metric)) {
      continue;
    }
    if (!list.empty()) {
      list += ", ";
    }
    list += unitName(unit);
  }
  return list;
}

// Scale from the metric's natural unit to the requested one. Each metric is computed in one
// base unit -- becquerel for activity, roentgen per hour for exposure -- and converted once
// at the very end, so the conversion cannot creep into the physics.
double unitScale(Unit unit) {
  switch (unit) {
    // A pack's coefficients are already in the pack's unit. There is nothing to convert to and
    // nothing this file could convert with.
    case Unit::PackDefined:
      return 1.0;
    case Unit::Curie:
      return 1.0 / units::kBqPerCi;
    case Unit::GrayPerHour:
    case Unit::Gray:
      return units::kGyPerR;
    // No conversion: a sievert is computed as a sievert. The weight below is already ICRP 116
    // effective dose per becquerel, not an air kerma waiting to be relabelled.
    case Unit::SievertPerHour:
    case Unit::Sievert:
      return 1.0;
    case Unit::Becquerel:
    case Unit::Decays:
    case Unit::RoentgenPerHour:
    case Unit::Roentgen:
    // Photon units are the metric's natural units: a weight in photons per decay times atoms
    // (or atom-seconds) is already the number printed. No conversion creeps in.
    case Unit::PhotonsPerSecond:
    case Unit::Photons:
    case Unit::PhotonsPerSquareMeterPerSecond:
    case Unit::PhotonsPerSquareMeter:
    // Heat is computed in watts per atom, which is joules per second: against atom-seconds the
    // same weight is joules, with no hour to take back out. Same arrangement as the photon
    // units and for the same reason -- the weight was never quoted per hour.
    case Unit::Watt:
    case Unit::Joule:
      return 1.0;
  }
  return 1.0;
}

// Reconcile the base unit's time with the domain's. Exposure is computed per HOUR, because
// that is the unit a rate is quoted in, while an interval weights atom-SECONDS -- so an
// integrated exposure carries an extra factor of an hour that has to come back out. Activity
// has no such mismatch: becquerel against atom-seconds is already a plain count of decays.
//
// Applied at the same point as the unit conversion, and nowhere else, for the same reason:
// one place where units are reconciled is one place where they can be wrong.
double domainScale(Metric metric, Domain domain) {
  if (metric == Metric::Exposure && domain == Domain::Interval) {
    return 1.0 / units::kSecondsPerHour;
  }
  return 1.0;
}

// The per-nuclide weight. This function IS the metric definition -- everything else in this
// file is bookkeeping over index spaces.
//
// Every metric is lambda times something: activity stops there, exposure carries on into the
// photon spectrum with its energy deposition factor, and photon stops at the photons
// themselves -- the same spectrum, with or without the geometry that says where they get to.
// The shared factor is not a coincidence: every metric NuSIFT reports is per-decay, so it is
// proportional to the decay rate, and the metric is what each decay is worth.
// A kernel pack's weight for one nuclide, and how much of it came from photon lines INSIDE the
// curve's tabulated range. The second number is what a coverage figure is built from: a
// spectrum sitting mostly off the end of a published curve is not an answer that curve can
// give, and clamping it quietly would hide exactly that.
struct KernelWeight {
  double total = 0.0;
  double inRange = 0.0;
};

KernelWeight kernelWeightFor(const ResponseSpec& spec, const NuclearData& data, int index) {
  const CoefficientPack& pack = *spec.pack->pack;
  const double lambda = data.decayConstant(index);
  KernelWeight weight;
  if (lambda <= 0.0) {
    return weight;
  }
  // The same construction the built-in photon metrics use, with the curve read from a file
  // instead of compiled in: each line's fluence at the point, times what the curve says a
  // photon of that energy is worth, summed with the geometry inside the sum.
  for (const GammaLine& line : data.lines(index)) {
    if (!(line.energyEv > 0.0 && line.intensity > 0.0)) {
      continue;
    }
    const double contribution = lambda * line.intensity *
                                exposure::pointFluenceCoeff(line.energyEv, spec.geometry) *
                                pack.kernelAt(line.energyEv);
    weight.total += contribution;
    if (!pack.kernelClamps(line.energyEv)) {
      weight.inRange += contribution;
    }
  }
  return weight;
}

double weightFor(const ResponseSpec& spec, const NuclearData& data, int index) {
  const double lambda = data.decayConstant(index);
  switch (spec.metric) {
    case Metric::Pack:
      // A kernel is evaluated here, against this nuclide's lines and this spec's geometry,
      // exactly as exposure is. A per-nuclide pack was resolved earlier -- basis and folds
      // included -- because which coefficient applies to a folded daughter is a question about
      // the seed, and the seed is not visible from here.
      if (spec.pack->pack->provenance().shape == PackShape::Kernel) {
        return kernelWeightFor(spec, data, index).total;
      }
      return spec.pack->weights[static_cast<std::size_t>(index)];
    case Metric::Activity:
      // Against atoms this is a rate in Bq; against atom-seconds it is a count of decays.
      // Same weight, different domain -- which is why Domain is a separate axis from Metric
      // rather than two metrics.
      return lambda;
    case Metric::Exposure:
      // Which QUANTITY this is depends on the unit, and that is the point rather than a
      // complication: roentgen and gray name air kerma, and a sievert names effective dose to
      // a person. They are different physics -- air's energy absorption against a phantom
      // calculation -- and a factor of five apart on a soft emitter, so one weight cannot
      // serve both by scaling.
      //
      // lambda * (per becquerel), in R/h or Sv/h per atom. Against atom-seconds it is the
      // accrued total, once domainScale() has taken the hour back out. Either per-becquerel
      // factor sums over the nuclide's photon lines WITH the geometry inside the sum, which is
      // why no single per-nuclide constant could stand in for it.
      if (isEffectiveDoseUnit(spec.unit)) {
        return lambda * exposure::effectiveDoseRatePerBecquerel(data.lines(index), spec.geometry);
      }
      return lambda * exposure::exposureRatePerBecquerel(data.lines(index), spec.geometry);
    case Metric::Photon:
      // lambda * (photons per decay): the source strength, in photons/s per atom, independent
      // of any geometry -- a property of the inventory alone. Against atom-seconds the same
      // weight is the photons emitted over the window, and there is no hour to take back out,
      // because the weight was never quoted per hour.
      if (isFluenceUnit(spec.unit)) {
        // The fluence rate at the point, in photons/(m^2*s) per atom. As for exposure the
        // attenuation sits inside the sum over lines, so the geometry is part of the weight
        // and part of what the answer means.
        return lambda * exposure::fluenceRatePerBecquerel(data.lines(index), spec.geometry);
      }
      return lambda * totalPhotonYield(data.lines(index));
    case Metric::Heat:
      // lambda * (recoverable energy per decay), in watts per atom. The three MT457 averages
      // summed: a decay's whole energy budget less the neutrinos the evaluation excludes.
      // No geometry and no spectrum, which makes this the cheapest metric in the file after
      // activity -- and, unlike exposure, a per-nuclide constant that a pack COULD have carried
      // if anyone published one. It is built in rather than packed because it is computed from
      // staged evaluated data rather than read from a published table.
      return lambda * data.decayEnergyEv(index) * units::kEvToJ;
  }
  return 0.0;
}

// Which aggregate bucket a nuclide falls into.
std::int64_t bucketKey(Aggregate aggregate, const Zai& zai) {
  switch (aggregate) {
    case Aggregate::Nuclide:
      return zai.key();
    case Aggregate::MassChain:
      return zai.a;
    case Aggregate::Element:
      return zai.z;
    case Aggregate::GammaLine:
      break;  // lines do not bucket by nuclide; see assembleLines
  }
  return zai.key();
}

// "Ba-137m 661.7 keV". keV rather than eV because that is how spectroscopy is spoken, and one
// decimal because that is the resolution at which lines are distinguished by name.
std::string lineLabel(const Zai& emitter, double energyEv) {
  char buffer[64];
  std::snprintf(buffer, sizeof(buffer), "%s %.1f keV", formatNuclideName(emitter).c_str(),
                energyEv / 1000.0);
  return buffer;
}

std::string bucketLabel(Aggregate aggregate, std::int64_t key, std::int64_t dominant) {
  switch (aggregate) {
    case Aggregate::Nuclide:
      return formatNuclideName(Zai::fromKey(key));
    case Aggregate::MassChain: {
      std::string label = "A=" + std::to_string(key);
      // Naming the dominant member is what makes the row actionable: "A=140" alone does not
      // tell anyone what to look at, while "A=140 (La-140)" does.
      if (dominant != 0) {
        label += " (" + formatNuclideName(Zai::fromKey(dominant)) + ")";
      }
      return label;
    }
    case Aggregate::Element: {
      std::string label = elementSymbol(static_cast<int>(key));
      if (dominant != 0) {
        label += " (" + formatNuclideName(Zai::fromKey(dominant)) + ")";
      }
      return label;
    }
    case Aggregate::GammaLine:
      break;  // labeled at construction, where the energy is in hand
  }
  return {};
}

// Accumulates columns for one aggregate bucket across every time.
struct Bucket {
  std::size_t column = 0;
  double peak = 0.0;              // largest single-nuclide contribution seen
  std::int64_t dominant = 0;      // the nuclide named in the label
  double dominantHalfLife = 0.0;  // its half-life, for breaking near-ties
  int flags = kFlagNone;
};

// Shared core: given per-time per-nuclide weighted values in the result's index space, fold
// them into aggregate columns. Both the instantaneous and interval builders funnel through
// here so the two can never disagree about how aggregation works.
ResponseTable assemble(const NuclearData& data, std::span<const std::int64_t> keys,
                       const std::vector<std::vector<double>>& weightedByTime,
                       const ResponseSpec& spec, Domain domain) {
  const std::size_t nNuc = keys.size();
  const std::size_t nT = weightedByTime.size();

  std::map<std::int64_t, Bucket> buckets;
  std::vector<std::size_t> columnOf(nNuc, 0);

  for (std::size_t i = 0; i < nNuc; ++i) {
    const Zai zai = Zai::fromKey(keys[i]);
    const std::int64_t bucket = bucketKey(spec.aggregate, zai);
    auto [it, inserted] = buckets.try_emplace(bucket);
    if (inserted) {
      it->second.column = buckets.size() - 1;
    }
    columnOf[i] = it->second.column;

    // Largest contribution this nuclide makes at any time, used to name the bucket. Taken
    // over the whole grid rather than at one time so the label does not change identity
    // partway down a column.
    double peak = 0.0;
    for (std::size_t k = 0; k < nT; ++k) {
      peak = std::max(peak, weightedByTime[k][i]);
    }

    const int dataIndex = data.indexOfKey(zai.key());
    const double halfLife = dataIndex >= 0 ? data.halfLifeSeconds(dataIndex) : 0.0;

    // Break a near-tie toward the longer-lived member. A chain in secular equilibrium has
    // every member at essentially the same activity, and whichever edges ahead numerically is
    // arbitrary -- but the answer is not: the long-lived parent is what controls the chain and
    // what anyone acting on the ranking would actually address. "A=90 (Sr-90)" is useful;
    // "A=90 (Y-90)" points at the 64-hour daughter that merely follows it.
    constexpr double kTieBand = 1.01;
    const bool clearlyLarger = peak > it->second.peak * kTieBand;
    const bool nearTieButLongerLived =
        peak * kTieBand >= it->second.peak && halfLife > it->second.dominantHalfLife;
    if (it->second.dominant == 0 || clearlyLarger || nearTieButLongerLived) {
      it->second.peak = std::max(it->second.peak, peak);
      it->second.dominant = zai.key();
      it->second.dominantHalfLife = halfLife;
    }

    // Photon metrics only, matching what unmodeledEnergyFraction below is gated on. The flag
    // says a photon spectrum is incomplete: that understates an exposure and a photon count
    // alike, and says nothing whatever about a count of decays, so an activity report carrying
    // it would end with a paragraph about a metric it never computed.
    if ((spec.metric == Metric::Exposure || spec.metric == Metric::Photon) && dataIndex >= 0 &&
        data.unmodeledPhotonFraction(dataIndex) > kUnmodeledContinuumFlag) {
      it->second.flags |= kFlagUnmodeledContinuum;
    }

    // What the pack has to say about this nuclide. A bucket carries the flag if ANY of its
    // members does, which is the honest direction for an aggregate: a mass chain half of which
    // the pack does not cover is a mass chain whose number is incomplete.
    if (spec.metric == Metric::Pack && dataIndex >= 0) {
      switch (spec.pack->coverage[static_cast<std::size_t>(dataIndex)]) {
        case PackCoverage::None:
          it->second.flags |= kFlagNotInPack;
          break;
        case PackCoverage::Folded:
          it->second.flags |= kFlagFoldedInPack;
          break;
        case PackCoverage::Own:
          break;
      }
    }
  }

  const std::size_t nC = buckets.size();

  ResponseTable table;
  table.metric = spec.metric;
  table.aggregate = spec.aggregate;
  table.domain = domain;
  table.unit = spec.unit;
  table.contributors.resize(nC);
  table.labels.resize(nC);
  table.flags.assign(nC, kFlagNone);
  table.values.assign(nT * nC, 0.0);
  table.totals.assign(nT, 0.0);

  for (const auto& [key, bucket] : buckets) {
    table.contributors[bucket.column] = ContributorId{key, bucket.dominant};
    table.labels[bucket.column] = bucketLabel(spec.aggregate, key, bucket.dominant);
    table.flags[bucket.column] = bucket.flags;
  }

  const double scale = unitScale(spec.unit) * domainScale(spec.metric, domain);
  for (std::size_t k = 0; k < nT; ++k) {
    const std::size_t base = k * nC;
    double total = 0.0;
    for (std::size_t i = 0; i < nNuc; ++i) {
      const double value = weightedByTime[k][i] * scale;
      table.values[base + columnOf[i]] += value;
      total += value;
    }
    table.totals[k] = total;
  }

  return table;
}

// One column per discrete photon line, rather than per nuclide.
//
// A line's weight is lambda_i * y_ij times the metric's kernel for that energy: k(E_j) for
// exposure, the fluence coefficient for a photon fluence unit, and one for photon strength.
// Multiplying by the emitter's atom count gives the response that one line contributes -- so
// the table is built from the same atoms as every other aggregate, only weighted more finely.
//
// A full evaluation carries 86000 lines, and a fission seed reaches thousands of emitters, so
// the columns are thresholded: a line contributing less than kLineFloor of its own emitter's
// exposure is dropped. Relative to the EMITTER rather than to the global total, deliberately --
// a global threshold would erase the entire spectrum of every minor nuclide, and the question
// "which line dominates THIS nuclide" is one people ask.
//
// Dropped lines still count toward the TOTAL. A total is over everything that contributes,
// and a line below the floor contributes -- it is only a column this table declines to carry.
// Keeping it in the total is what makes the total by line the total by nuclide exactly, rather
// than a number that falls short of it by however many trace lines the spectra happened to
// hold, and it is what lets the coverage a line ranking reports stay a true fraction.
ResponseTable assembleLines(const NuclearData& data, std::span<const std::int64_t> keys,
                            const std::vector<std::vector<double>>& atomsByTime,
                            const ResponseSpec& spec, Domain domain) {
  constexpr double kLineFloor = 1.0e-6;

  const std::size_t nNuc = keys.size();
  const std::size_t nT = atomsByTime.size();

  struct Column {
    std::int64_t emitterKey = 0;
    double energyEv = 0.0;
    double weight = 0.0;      // lambda * intensity * the metric's kernel for the energy
    std::size_t nuclide = 0;  // index into the result's space
    int flags = kFlagNone;
  };
  std::vector<Column> columns;
  // Per nuclide, lambda times what its sub-floor lines deliver per becquerel: the weight of
  // everything this table does not carry as a column.
  std::vector<double> droppedWeight(nNuc, 0.0);

  for (std::size_t i = 0; i < nNuc; ++i) {
    const Zai zai = Zai::fromKey(keys[i]);
    const int dataIndex = data.indexOfKey(zai.key());
    if (dataIndex < 0) {
      continue;
    }
    const double lambda = data.decayConstant(dataIndex);
    if (lambda <= 0.0) {
      continue;
    }
    const LineSpectrum lines = data.lines(dataIndex);
    if (lines.empty()) {
      continue;
    }

    // The emitter's total and the per-line coefficient in the SAME weight, so the relative
    // floor below drops a line by its share of its emitter whatever the metric.
    const bool fluence = isFluenceUnit(spec.unit);
    const double emitterTotal = spec.metric == Metric::Exposure
                                    ? exposure::exposureRatePerBecquerel(lines, spec.geometry)
                                : fluence ? exposure::fluenceRatePerBecquerel(lines, spec.geometry)
                                          : totalPhotonYield(lines);
    const double floor = emitterTotal * kLineFloor;
    const int flags = data.unmodeledPhotonFraction(dataIndex) > kUnmodeledContinuumFlag
                          ? kFlagUnmodeledContinuum
                          : kFlagNone;

    for (const GammaLine& line : lines) {
      const double perBecquerel =
          spec.metric == Metric::Exposure
              ? line.intensity * exposure::pointExposureCoeff(line.energyEv, spec.geometry)
          : fluence ? line.intensity * exposure::pointFluenceCoeff(line.energyEv, spec.geometry)
                    : line.intensity;
      if (perBecquerel <= 0.0) {
        continue;
      }
      if (perBecquerel < floor) {
        droppedWeight[i] += lambda * perBecquerel;
        continue;
      }
      columns.push_back(Column{zai.key(), line.energyEv, lambda * perBecquerel, i, flags});
    }
  }

  const std::size_t nC = columns.size();

  ResponseTable table;
  table.metric = spec.metric;
  table.aggregate = spec.aggregate;
  table.domain = domain;
  table.unit = spec.unit;
  table.contributors.resize(nC);
  table.labels.resize(nC);
  table.flags.assign(nC, kFlagNone);
  table.values.assign(nT * nC, 0.0);
  table.totals.assign(nT, 0.0);

  for (std::size_t c = 0; c < nC; ++c) {
    const Column& column = columns[c];
    const Zai emitter = Zai::fromKey(column.emitterKey);
    table.contributors[c] = ContributorId{column.emitterKey, column.emitterKey, column.energyEv};
    table.labels[c] = lineLabel(emitter, column.energyEv);
    table.flags[c] = column.flags;
  }

  const double scale = unitScale(spec.unit) * domainScale(spec.metric, domain);
  for (std::size_t k = 0; k < nT; ++k) {
    const std::size_t base = k * nC;
    double total = 0.0;
    for (std::size_t c = 0; c < nC; ++c) {
      const Column& column = columns[c];
      const double value = column.weight * atomsByTime[k][column.nuclide] * scale;
      table.values[base + c] = value;
      total += value;
    }
    for (std::size_t i = 0; i < nNuc; ++i) {
      total += droppedWeight[i] * atomsByTime[k][i] * scale;
    }
    table.totals[k] = total;
  }

  return table;
}

// Fraction of the emitted photon energy rate that lives in a spectrum NuSIFT does not model.
// Activity-weighted, so a nuclide with a large unmodeled fraction but negligible activity
// contributes negligibly -- which is the whole point of reporting a magnitude rather than a
// count of flagged nuclides.
std::vector<double> unmodeledEnergyFractions(const NuclearData& data,
                                             std::span<const std::int64_t> keys,
                                             const std::vector<std::vector<double>>& atomsByTime) {
  const std::size_t nNuc = keys.size();
  std::vector<double> modeled(nNuc, 0.0);
  std::vector<double> missing(nNuc, 0.0);

  for (std::size_t i = 0; i < nNuc; ++i) {
    const int index = data.indexOfKey(keys[i]);
    if (index < 0) {
      continue;
    }
    const double lambda = data.decayConstant(index);
    if (lambda <= 0.0) {
      continue;
    }
    modeled[i] = lambda * discretePhotonEnergyEv(data.lines(index));
    missing[i] = lambda * data.continuumPhotonEv(index);
  }

  std::vector<double> fractions(atomsByTime.size(), 0.0);
  for (std::size_t k = 0; k < atomsByTime.size(); ++k) {
    double modelledTotal = 0.0;
    double missingTotal = 0.0;
    for (std::size_t i = 0; i < nNuc; ++i) {
      const double atoms = atomsByTime[k][i];
      modelledTotal += modeled[i] * atoms;
      missingTotal += missing[i] * atoms;
    }
    const double total = modelledTotal + missingTotal;
    fractions[k] = total > 0.0 ? missingTotal / total : 0.0;
  }
  return fractions;
}

// The air path's optical depth at each time, averaged over the nuclides by the exposure each
// delivers. `weightedByTime` holds lambda times the per-becquerel exposure times atoms (or
// atom-seconds) -- each nuclide's share of the total -- so the unit scale cancels and the
// domain does not matter.
std::vector<double> meanOpticalDepths(const NuclearData& data, std::span<const std::int64_t> keys,
                                      const std::vector<std::vector<double>>& weightedByTime,
                                      const exposure::PointSourceGeometry& geometry) {
  const std::size_t nNuc = keys.size();
  std::vector<double> depth(nNuc, 0.0);
  for (std::size_t i = 0; i < nNuc; ++i) {
    const int index = data.indexOfKey(keys[i]);
    if (index >= 0) {
      depth[i] = exposure::meanOpticalDepth(data.lines(index), geometry);
    }
  }

  std::vector<double> means(weightedByTime.size(), 0.0);
  for (std::size_t k = 0; k < weightedByTime.size(); ++k) {
    double weighted = 0.0;
    double total = 0.0;
    for (std::size_t i = 0; i < nNuc; ++i) {
      const double share = weightedByTime[k][i];
      if (share > 0.0) {
        weighted += share * depth[i];
        total += share;
      }
    }
    means[k] = total > 0.0 ? weighted / total : 0.0;
  }
  return means;
}

// The caveats a table carries depend on which part of the photon model it actually used. The
// unmodeled continuum understates every photon metric -- an exposure and a photon count
// alike -- so both carry the energy fraction. The air-path depth is a property only of an
// answer that USED the geometry: a photon-strength number at no distance has no path to be
// thick, and footnoting one with a depth computed from a geometry it never used would report
// a limit of a model it did not run.
void fillPhotonCaveats(ResponseTable& table, const NuclearData& data,
                       std::span<const std::int64_t> keys,
                       const std::vector<std::vector<double>>& rawAtoms,
                       const std::vector<std::vector<double>>& weighted, const ResponseSpec& spec) {
  if (spec.metric != Metric::Exposure && spec.metric != Metric::Photon) {
    return;
  }
  table.unmodeledEnergyFraction = unmodeledEnergyFractions(data, keys, rawAtoms);
  if (spec.metric == Metric::Exposure || isFluenceUnit(spec.unit)) {
    table.meanOpticalDepth = meanOpticalDepths(data, keys, weighted, spec.geometry);
  }
}

// --- naming a contributor to pin ---------------------------------------------

// What fraction of the quantity the pack's basis measures is actually carried by the pack, at
// each time. Measured in the BASIS rather than in the response's own unit, because the response
// is zero for exactly the nuclides in question: weighting the missing ones by the coefficient
// they do not have would make coverage 100% by construction.
std::vector<double> packCoverageByTime(const NuclearData& data, std::span<const std::int64_t> keys,
                                       const std::vector<std::vector<double>>& atomsByTime,
                                       const ResponseSpec& spec) {
  const std::size_t nNuc = keys.size();
  const std::size_t nT = atomsByTime.size();

  // A kernel pack asks a different coverage question from a per-nuclide one. There is no list
  // of nuclides it does or does not carry -- it applies to every photon line there is -- so
  // what it can fail to speak for is an ENERGY: a line outside the curve's tabulated range
  // takes a clamped value. Coverage is therefore the share of the response that came from
  // lines inside the range, which is the same thing the built-in kernels report about their
  // own clamps.
  if (spec.pack->pack->provenance().shape == PackShape::Kernel) {
    std::vector<double> total(nNuc, 0.0);
    std::vector<double> inRange(nNuc, 0.0);
    for (std::size_t i = 0; i < nNuc; ++i) {
      const int index = data.indexOfKey(keys[i]);
      if (index < 0) {
        continue;
      }
      const KernelWeight weight = kernelWeightFor(spec, data, index);
      total[i] = weight.total;
      inRange[i] = weight.inRange;
    }
    std::vector<double> coverage(nT, 1.0);
    for (std::size_t k = 0; k < nT; ++k) {
      double all = 0.0;
      double carried = 0.0;
      for (std::size_t i = 0; i < nNuc; ++i) {
        all += total[i] * atomsByTime[k][i];
        carried += inRange[i] * atomsByTime[k][i];
      }
      coverage[k] = all > 0.0 ? carried / all : 1.0;
    }
    return coverage;
  }

  // The basis quantity per atom, which is what a coefficient would have multiplied.
  std::vector<double> perAtom(nNuc, 0.0);
  std::vector<char> counted(nNuc, 0);
  for (std::size_t i = 0; i < nNuc; ++i) {
    // The store indexes with int across its whole interface, so this is the one conversion the
    // loop cannot avoid -- and it is a different index space from `i`, which is why it is named.
    const int index = data.indexOfKey(keys[i]);
    if (index < 0) {
      continue;
    }
    switch (spec.pack->pack->provenance().basis) {
      case PackBasis::Activity:
        perAtom[i] = data.decayConstant(index);
        break;
      case PackBasis::Atoms:
        perAtom[i] = 1.0;
        break;
      case PackBasis::Mass:
        perAtom[i] = data.molarMassGPerMol(index) / units::kAvogadro;
        break;
      case PackBasis::Concentration:
        // Coverage asks which nuclides the pack can speak for, and the extent divides every
        // term of that ratio equally, so it cancels. Activity is the quantity being shared out.
        perAtom[i] = data.decayConstant(index);
        break;
    }
    // Folded counts as covered: the daughter's contribution is inside a parent's coefficient,
    // so it is accounted for even though it carries no weight of its own. Reporting it as
    // missing is how a coverage figure lies about the nuclides that usually dominate.
    counted[i] = static_cast<char>(spec.pack->coverage[static_cast<std::size_t>(index)] !=
                                   PackCoverage::None);
  }

  std::vector<double> coverage(nT, 0.0);
  for (std::size_t k = 0; k < nT; ++k) {
    double total = 0.0;
    double carried = 0.0;
    for (std::size_t i = 0; i < nNuc; ++i) {
      const double quantity = perAtom[i] * atomsByTime[k][i];
      total += quantity;
      if (counted[i] != 0) {
        carried += quantity;
      }
    }
    // An empty inventory covers nothing and misses nothing; reporting 0% would read as a
    // warning about a pack that has not been asked anything yet.
    coverage[k] = total > 0.0 ? carried / total : 1.0;
  }
  return coverage;
}

std::string_view trimmed(std::string_view text) {
  const auto space = [](char c) { return std::isspace(static_cast<unsigned char>(c)) != 0; };
  while (!text.empty() && space(text.front())) {
    text.remove_prefix(1);
  }
  while (!text.empty() && space(text.back())) {
    text.remove_suffix(1);
  }
  return text;
}

// Consume a leading "A=" or "Z=", case-insensitively. The qualified forms exist so a pin can
// say which number it means; the bare number is accepted too because in a table ranked by mass
// chain there is nothing else "140" could be.
bool stripQualifier(std::string_view& text, char letter) {
  if (text.size() >= 2 && std::tolower(static_cast<unsigned char>(text.front())) == letter &&
      text[1] == '=') {
    text.remove_prefix(2);
    return true;
  }
  return false;
}

std::optional<int> wholeNumber(std::string_view text) {
  if (text.empty()) {
    return std::nullopt;
  }
  int value = 0;
  for (const char c : text) {
    if (c < '0' || c > '9') {
      return std::nullopt;
    }
    value = value * 10 + (c - '0');
    if (value > 10'000'000) {
      return std::nullopt;  // far past anything either field can hold
    }
  }
  return value;
}

// The bucket key a pin spelling names, in the key space `aggregate` buckets into.
//
// Parsed here rather than through requireNuclideName so a failure can say what the spelling
// should have been for THIS aggregate: told "Cs-13" while ranking by mass chain, a message
// about nuclide names names only half the ways the pin could have been written.
std::int64_t pinKey(Aggregate aggregate, std::string_view text) {
  const std::string_view whole = trimmed(text);
  std::string_view rest = whole;
  if (rest.empty()) {
    throw InputError(tagged(kPinModule, "a pin needs a contributor to name"));
  }

  // The forms that would have worked, named only if none of them did.
  const char* expected = "a nuclide name";
  switch (aggregate) {
    case Aggregate::MassChain: {
      expected = "a mass number, A=140, or a nuclide name";
      const bool qualified = stripQualifier(rest, 'a');
      if (const std::optional<int> number = wholeNumber(rest);
          number.has_value() && *number >= 1 && *number <= kMaxMassNumber) {
        return *number;
      }
      // A nuclide name, or the raw-key form -- both of which know their own mass number.
      if (!qualified) {
        if (const std::optional<Zai> zai = parseNuclideName(rest); zai.has_value()) {
          return zai->a;
        }
      }
      break;
    }
    case Aggregate::Element: {
      expected = "an element symbol, Z=55, or a nuclide name";
      const bool qualified = stripQualifier(rest, 'z');
      if (const std::optional<int> number = wholeNumber(rest);
          number.has_value() && *number >= 1 && *number <= kMaxAtomicNumber) {
        return *number;
      }
      if (!qualified) {
        if (const int z = atomicNumber(rest); z > 0) {
          return z;
        }
        if (const std::optional<Zai> zai = parseNuclideName(rest); zai.has_value()) {
          return zai->z;
        }
      }
      break;
    }
    case Aggregate::Nuclide:
    case Aggregate::GammaLine:
      if (const std::optional<Zai> zai = parseNuclideName(rest); zai.has_value()) {
        return zai->key();
      }
      break;
  }

  throw InputError(tagged(kPinModule, "\"" + std::string(whole) + "\" is not " + expected +
                                          " (this ranking is by " + aggregateName(aggregate) +
                                          ")"));
}

void requireUsableSpec(const NuclearData& data, const ResponseSpec& spec, Domain domain) {
  if (spec.metric == Metric::Pack) {
    if (spec.pack == nullptr || spec.pack->pack == nullptr) {
      throw InputError(tagged(kModule, "this metric is a coefficient pack, and no pack was given"));
    }
    if (spec.pack->weights.size() != static_cast<std::size_t>(data.size())) {
      throw InputError(tagged(kModule,
                              "this pack was resolved against a different data store than the "
                              "one being solved, so its weights do not line up with the "
                              "nuclides they would weight"));
    }
    // The pack states which domains its quantity is meaningful in, and a rate-like coefficient
    // against an interval integral is a different quantity rather than the same one summed.
    const PackDomains domains = spec.pack->pack->provenance().domains;
    const bool instantOk = domains != PackDomains::IntervalOnly;
    const bool intervalOk = domains != PackDomains::InstantOnly;
    if ((domain == Domain::Instant && !instantOk) || (domain == Domain::Interval && !intervalOk)) {
      throw InputError(tagged(kModule, "pack \"" + spec.pack->pack->provenance().name +
                                           "\" declares itself " +
                                           (domains == PackDomains::InstantOnly ? "instantaneous"
                                                                                : "an "
                                                                                  "interval "
                                                                                  "quantity") +
                                           ", so it cannot answer the other"));
    }
  }
  if (!unitSuitsMetric(spec.unit, spec.metric)) {
    throw InputError(tagged(kModule, std::string("unit ") + unitName(spec.unit) +
                                         " does not measure " + metricName(spec.metric)));
  }
  if (!unitSuitsDomain(spec.unit, domain)) {
    throw InputError(
        tagged(kModule, std::string("unit ") + unitName(spec.unit) +
                            (domain == Domain::Interval
                                 ? " is a rate and cannot express a time-integrated total"
                                 : " is a total and cannot express an instantaneous value")));
  }
  // A store with no photon lines cannot answer a photon question at all. Returning zeros
  // would be indistinguishable from "nothing here emits photons", which is a different and
  // much more alarming statement.
  const bool isPhotonMetric = spec.metric == Metric::Exposure || spec.metric == Metric::Photon;
  if (spec.aggregate == Aggregate::GammaLine && !isPhotonMetric) {
    throw InputError(tagged(kModule,
                            "ranking by gamma line only makes sense for exposure or photon "
                            "output -- a photon line has no activity of its own, it is a way "
                            "its emitter's decays get out"));
  }
  if (isPhotonMetric && !data.hasPhotonLines()) {
    throw InputError(tagged(kModule, "this data store carries no photon lines, so " +
                                         std::string(metricName(spec.metric)) +
                                         " cannot be computed. Stage from ENDF decay tapes, "
                                         "which carry the discrete spectra, or rank by "
                                         "activity instead"));
  }
}

}  // namespace

bool isEffectiveDoseUnit(Unit unit) {
  return unit == Unit::SievertPerHour || unit == Unit::Sievert;
}

bool isFluenceUnit(Unit unit) {
  return unit == Unit::PhotonsPerSquareMeterPerSecond || unit == Unit::PhotonsPerSquareMeter;
}

bool unitSuitsDomain(Unit unit, Domain domain) {
  switch (unit) {
    // The pack decides, not the unit: `# domain:` in its header says whether its quantity is
    // instantaneous, integrated or both, and requireUsableSpec() enforces that against the
    // pack in hand. There is nothing for a unit with no fixed meaning to rule on here.
    case Unit::PackDefined:
      return true;
    case Unit::Becquerel:
    case Unit::Curie:
    case Unit::RoentgenPerHour:
    case Unit::GrayPerHour:
    case Unit::SievertPerHour:
    case Unit::PhotonsPerSecond:
    case Unit::PhotonsPerSquareMeterPerSecond:
    case Unit::Watt:
      return domain == Domain::Instant;
    case Unit::Decays:
    case Unit::Roentgen:
    case Unit::Gray:
    case Unit::Sievert:
    case Unit::Photons:
    case Unit::PhotonsPerSquareMeter:
    case Unit::Joule:
      return domain == Domain::Interval;
  }
  return false;
}

bool unitSuitsMetric(Unit unit, Metric metric) {
  switch (unit) {
    case Unit::PackDefined:
      return metric == Metric::Pack;
    case Unit::Becquerel:
    case Unit::Curie:
    case Unit::Decays:
      return metric == Metric::Activity;
    case Unit::RoentgenPerHour:
    case Unit::GrayPerHour:
    case Unit::SievertPerHour:
    case Unit::Roentgen:
    case Unit::Gray:
    case Unit::Sievert:
      return metric == Metric::Exposure;
    case Unit::PhotonsPerSecond:
    case Unit::Photons:
    case Unit::PhotonsPerSquareMeterPerSecond:
    case Unit::PhotonsPerSquareMeter:
      return metric == Metric::Photon;
    case Unit::Watt:
    case Unit::Joule:
      return metric == Metric::Heat;
  }
  return false;
}

const char* unitName(Unit unit) {
  switch (unit) {
    // A placeholder, and one a report should not print: ResponseTable::unitLabel carries the
    // spelling the pack declared, and every writer prefers it when it is set.
    case Unit::PackDefined:
      return "pack-defined";
    case Unit::Becquerel:
      return "Bq";
    case Unit::Curie:
      return "Ci";
    case Unit::Decays:
      return "decays";
    case Unit::RoentgenPerHour:
      return "R/h";
    case Unit::GrayPerHour:
      return "Gy/h";
    case Unit::SievertPerHour:
      return "Sv/h";
    case Unit::Roentgen:
      return "R";
    case Unit::Gray:
      return "Gy";
    case Unit::Sievert:
      return "Sv";
    case Unit::PhotonsPerSecond:
      return "photons/s";
    case Unit::Photons:
      return "photons";
    case Unit::PhotonsPerSquareMeterPerSecond:
      return "photons/m2/s";
    case Unit::PhotonsPerSquareMeter:
      return "photons/m2";
    case Unit::Watt:
      return "W";
    case Unit::Joule:
      return "J";
  }
  return "?";
}

bool parseUnit(std::string_view text, Unit& out) {
  // Matched against unitName() rather than against a table of its own: the accepted spellings
  // ARE the printed ones, so adding a unit cannot leave it unparseable.
  for (const Unit unit : kAllUnits) {
    if (equalsIgnoreCase(text, unitName(unit))) {
      out = unit;
      return true;
    }
  }
  return false;
}

Unit defaultUnit(Metric metric, Domain domain) {
  if (metric == Metric::Pack) {
    return Unit::PackDefined;
  }
  if (metric == Metric::Exposure) {
    return domain == Domain::Interval ? Unit::Roentgen : Unit::RoentgenPerHour;
  }
  if (metric == Metric::Photon) {
    // The source strength, not the fluence: the default is the geometry-free number, so a
    // photon answer that does not name a distance is not silently one at some distance.
    return domain == Domain::Interval ? Unit::Photons : Unit::PhotonsPerSecond;
  }
  if (metric == Metric::Heat) {
    return domain == Domain::Interval ? Unit::Joule : Unit::Watt;
  }
  return domain == Domain::Interval ? Unit::Decays : Unit::Becquerel;
}

Unit requireUnit(std::string_view text, Metric metric, Domain domain) {
  if (text.empty()) {
    return defaultUnit(metric, domain);
  }
  Unit unit = Unit::Becquerel;
  if (parseUnit(text, unit)) {
    return unit;
  }
  throw InputError(tagged(kUnitsModule, "\"" + std::string(text) + "\" is not a unit (activity: " +
                                            spellingsFor(Metric::Activity) +
                                            "; exposure: " + spellingsFor(Metric::Exposure) +
                                            "; photon: " + spellingsFor(Metric::Photon) +
                                            "; decay heat: " + spellingsFor(Metric::Heat) + ")"));
}

const char* metricName(Metric metric) {
  switch (metric) {
    case Metric::Activity:
      return "activity";
    case Metric::Exposure:
      return "exposure";
    case Metric::Photon:
      return "photon";
    case Metric::Heat:
      return "decay heat";
    // Deliberately generic. What this metric IS lives in the pack's own header -- its quantity,
    // its version, its scenario -- and a report prints those rather than this word.
    case Metric::Pack:
      return "pack";
  }
  return "?";
}

const char* aggregateName(Aggregate aggregate) {
  switch (aggregate) {
    case Aggregate::Nuclide:
      return "nuclide";
    case Aggregate::MassChain:
      return "mass chain";
    case Aggregate::Element:
      return "element";
    case Aggregate::GammaLine:
      return "gamma line";
  }
  return "?";
}

std::int64_t requirePin(const ResponseTable& table, std::string_view text) {
  const std::int64_t key = pinKey(table.aggregate, text);
  for (const ContributorId& id : table.contributors) {
    if (id.key == key) {
      return key;
    }
  }

  // Refused rather than pinned to a row of zeros. A pin that silently resolves to nothing is
  // the worst possible answer to "where does Cs-137 stand": it looks like the ranking was
  // asked and replied "nowhere", when in fact the question never reached the table.
  const std::string reading = table.aggregate == Aggregate::GammaLine
                                  ? formatNuclideName(Zai::fromKey(key))
                                  : bucketLabel(table.aggregate, key, 0);
  throw InputError(tagged(kPinModule, "\"" + std::string(trimmed(text)) + "\" names " + reading +
                                          ", which this table does not carry -- it is ranked by " +
                                          aggregateName(table.aggregate) +
                                          " and nothing in the inventory's chain reaches it"));
}

ResolvedPack resolvePack(const CoefficientPack& pack, const NuclearData& data,
                         const Inventory& seed, const PackExtent& extent) {
  ResolvedPack resolved;
  resolved.pack = &pack;
  resolved.extent = extent;
  resolved.weights = pack.weights(data, seed, extent);
  resolved.coverage = pack.covered(data, seed);
  return resolved;
}

std::vector<double> responseWeights(const NuclearData& data, const ResponseSpec& spec) {
  // Instant rather than the caller's domain: unitScale and domainScale are applied by the
  // table builders, and a weight is the same per-atom quantity either way. What this checks is
  // the part that does not depend on domain -- a store with no photon lines cannot answer an
  // exposure question however it is asked.
  requireUsableSpec(data, spec, Domain::Instant);
  // Scaled into spec.unit here, not left in the base unit, so a caller that dots this against
  // an inventory gets the number the report would print. The alternative -- returning becquerel
  // per atom and making every caller remember the conversion -- is the creep this file's
  // unitScale() comment exists to prevent.
  const double scale = unitScale(spec.unit) * domainScale(spec.metric, Domain::Instant);
  const int n = data.size();
  std::vector<double> weight(static_cast<std::size_t>(n), 0.0);
  for (int i = 0; i < n; ++i) {
    weight[static_cast<std::size_t>(i)] = weightFor(spec, data, i) * scale;
  }
  return weight;
}

std::vector<double> weightDecayDerivatives(const NuclearData& data, const ResponseSpec& spec) {
  const std::vector<double> weights = responseWeights(data, spec);
  std::vector<double> derivatives(weights.size(), 0.0);

  // Whether this spec's weight carries lambda at all. Every built-in does; a pack does only on
  // the two bases whose coefficient multiplies an activity.
  bool scalesWithLambda = true;
  if (spec.metric == Metric::Pack) {
    const PackBasis basis = spec.pack->pack->provenance().basis;
    scalesWithLambda = basis == PackBasis::Activity || basis == PackBasis::Concentration;
  }
  if (!scalesWithLambda) {
    return derivatives;
  }

  for (std::size_t i = 0; i < weights.size(); ++i) {
    const double lambda = data.decayConstant(static_cast<int>(i));
    if (lambda > 0.0) {
      derivatives[i] = weights[i] / lambda;
    }
  }
  return derivatives;
}

ExposureCaveats exposureCaveats(const NuclearData& data, std::span<const std::int64_t> keys,
                                std::span<const double> atoms, const ResponseSpec& spec) {
  ExposureCaveats out;
  if (spec.metric != Metric::Exposure && spec.metric != Metric::Photon) {
    return out;
  }
  if (keys.size() != atoms.size()) {
    throw NusiftError(tagged(kModule, "caveats: atoms do not match their index space"));
  }
  // Wrapped in a one-row outer vector because the two helpers below are written for a whole
  // time grid, and the alternative -- a second single-time copy of each -- is exactly the
  // duplication this function exists to remove.
  std::vector<std::vector<double>> raw{std::vector<double>(atoms.begin(), atoms.end())};
  std::vector<std::vector<double>> weighted{std::vector<double>(atoms.size(), 0.0)};

  std::set<std::string> emitters;
  for (std::size_t i = 0; i < keys.size(); ++i) {
    const int index = data.indexOfKey(keys[i]);
    if (index < 0) {
      continue;
    }
    weighted.front()[i] = weightFor(spec, data, index) * atoms[i];
    if (data.unmodeledPhotonFraction(index) > kUnmodeledContinuumFlag) {
      emitters.insert(formatNuclideName(Zai::fromKey(keys[i])));
    }
  }

  out.unmodeledEnergyFraction = unmodeledEnergyFractions(data, keys, raw).front();
  // Same rule as fillPhotonCaveats(): the depth only describes an answer that used the path.
  if (spec.metric == Metric::Exposure || isFluenceUnit(spec.unit)) {
    out.meanOpticalDepth = meanOpticalDepths(data, keys, weighted, spec.geometry).front();
  }
  out.unmodeledContinuum.assign(emitters.begin(), emitters.end());
  return out;
}

ResponseTable buildResponse(const NuclearData& data, const DecayResult& result,
                            const ResponseSpec& spec) {
  requireUsableSpec(data, spec, Domain::Instant);

  const std::size_t nNuc = result.nuclideKeys.size();
  const std::size_t nT = result.times.size();

  std::vector<double> weight(nNuc, 0.0);
  for (std::size_t i = 0; i < nNuc; ++i) {
    // indexOfKey answers in the STORE's index space, which is int-based and is not this loop's.
    const int index = data.indexOfKey(result.nuclideKeys[i]);
    weight[i] = index >= 0 ? weightFor(spec, data, index) : 0.0;
  }

  std::vector<std::vector<double>> weighted(nT);
  std::vector<std::vector<double>> rawAtoms(nT);
  for (std::size_t k = 0; k < nT; ++k) {
    const std::span<const double> atoms = result.atomsAt(static_cast<int>(k));
    std::vector<double> row(nNuc, 0.0);
    std::vector<double> raw(atoms.begin(), atoms.end());
    for (std::size_t i = 0; i < nNuc; ++i) {
      row[i] = weight[i] * atoms[i];
    }
    weighted[k] = std::move(row);
    rawAtoms[k] = std::move(raw);
  }

  // The line assembler applies its own per-line weights, so it takes atoms directly rather
  // than the per-nuclide weighted values every other aggregate folds together.
  ResponseTable table =
      spec.aggregate == Aggregate::GammaLine
          ? assembleLines(data, result.nuclideKeys, rawAtoms, spec, Domain::Instant)
          : assemble(data, result.nuclideKeys, weighted, spec, Domain::Instant);
  table.times = result.times;
  table.geometry = spec.geometry;
  fillPhotonCaveats(table, data, result.nuclideKeys, rawAtoms, weighted, spec);
  if (spec.metric == Metric::Pack) {
    table.unitLabel = spec.pack->pack->provenance().unit;
    table.packCoverage = packCoverageByTime(data, result.nuclideKeys, rawAtoms, spec);
  }
  return table;
}

ResponseTable buildIntervalResponse(const NuclearData& data, std::span<const std::int64_t> keys,
                                    std::span<const double> integral, double t1, double t2,
                                    const ResponseSpec& spec) {
  requireUsableSpec(data, spec, Domain::Interval);
  if (keys.size() != integral.size()) {
    throw NusiftError(tagged(kModule, "interval integral does not match its index space"));
  }

  const std::size_t nNuc = keys.size();
  const std::vector<double> integratedAtoms(integral.begin(), integral.end());

  // lambda against atom-seconds: a dimensionless count of decays over the window, or the
  // roentgen accrued once domainScale() has taken the hour back out.
  std::vector<std::vector<double>> weighted{std::vector<double>(nNuc)};
  for (std::size_t i = 0; i < nNuc; ++i) {
    const int index = data.indexOfKey(keys[i]);
    const double weight = index >= 0 ? weightFor(spec, data, index) : 0.0;
    weighted[0][i] = weight * integratedAtoms[i];
  }

  // Same split as the instantaneous builder, for the same reason: the line assembler applies
  // its own per-line weights, so it takes atom-seconds directly rather than the per-nuclide
  // weighted value every other aggregate folds together.
  ResponseTable table =
      spec.aggregate == Aggregate::GammaLine
          ? assembleLines(data, keys, std::vector<std::vector<double>>{integratedAtoms}, spec,
                          Domain::Interval)
          : assemble(data, keys, weighted, spec, Domain::Interval);

  table.times = {t1};
  table.timeEnds = {t2};
  table.geometry = spec.geometry;
  if (spec.metric == Metric::Pack) {
    table.unitLabel = spec.pack->pack->provenance().unit;
    // Over the integral rather than over an instant: the coverage of an accrued answer is the
    // coverage of what accrued it.
    table.packCoverage =
        packCoverageByTime(data, keys, std::vector<std::vector<double>>{integratedAtoms}, spec);
  }
  fillPhotonCaveats(table, data, keys, std::vector<std::vector<double>>{integratedAtoms}, weighted,
                    spec);
  return table;
}

}  // namespace nusift
