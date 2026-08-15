#include "nusift/triage/response.hpp"

#include <algorithm>
#include <cctype>
#include <map>
#include <optional>
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

// Every unit, in the order the help text lists them. The one place the set is enumerated, so
// parseUnit and the error message it raises cannot come to disagree about what exists.
constexpr Unit kAllUnits[] = {
    Unit::Becquerel,      Unit::Curie,    Unit::Decays, Unit::RoentgenPerHour, Unit::GrayPerHour,
    Unit::SievertPerHour, Unit::Roentgen, Unit::Gray,   Unit::Sievert,
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
    case Unit::Curie:
      return 1.0 / units::kBqPerCi;
    case Unit::GrayPerHour:
    case Unit::SievertPerHour:
    case Unit::Gray:
    case Unit::Sievert:
      return exposure::kGrayPerRoentgen;
    case Unit::Becquerel:
    case Unit::Decays:
    case Unit::RoentgenPerHour:
    case Unit::Roentgen:
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
// Both metrics are lambda times something: activity stops there, exposure carries on into the
// photon spectrum. That shared factor is not a coincidence -- every metric NuSIFT reports is
// per-decay, so it is proportional to the decay rate, and the metric is what each decay is
// worth.
double weightFor(const ResponseSpec& spec, const NuclearData& data, int index) {
  const double lambda = data.decayConstant(index);
  switch (spec.metric) {
    case Metric::Activity:
      // Against atoms this is a rate in Bq; against atom-seconds it is a count of decays.
      // Same weight, different domain -- which is why Domain is a separate axis from Metric
      // rather than two metrics.
      return lambda;
    case Metric::Exposure:
      // lambda * (exposure per becquerel), in R/h per atom. Against atom-seconds it is
      // roentgen accrued, once domainScale() has taken the hour back out. The per-becquerel
      // factor sums over the nuclide's photon lines WITH air attenuation inside the sum,
      // which is why it depends on the geometry and why no single per-nuclide constant could
      // stand in for it.
      return lambda * exposure::exposureRatePerBecquerel(data.lines(index), spec.geometry);
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
      break;  // labelled at construction, where the energy is in hand
  }
  return {};
}

// Accumulates columns for one aggregate bucket across every time.
struct Bucket {
  int column = 0;
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
  const int nNuc = static_cast<int>(keys.size());
  const int nT = static_cast<int>(weightedByTime.size());

  std::map<std::int64_t, Bucket> buckets;
  std::vector<int> columnOf(static_cast<std::size_t>(nNuc), -1);

  for (int i = 0; i < nNuc; ++i) {
    const Zai zai = Zai::fromKey(keys[static_cast<std::size_t>(i)]);
    const std::int64_t bucket = bucketKey(spec.aggregate, zai);
    auto [it, inserted] = buckets.try_emplace(bucket);
    if (inserted) {
      it->second.column = static_cast<int>(buckets.size()) - 1;
    }
    columnOf[static_cast<std::size_t>(i)] = it->second.column;

    // Largest contribution this nuclide makes at any time, used to name the bucket. Taken
    // over the whole grid rather than at one time so the label does not change identity
    // partway down a column.
    double peak = 0.0;
    for (int k = 0; k < nT; ++k) {
      peak =
          std::max(peak, weightedByTime[static_cast<std::size_t>(k)][static_cast<std::size_t>(i)]);
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

    // Exposure only, matching what unmodeledEnergyFraction below is gated on. The flag says a
    // photon spectrum is incomplete: that understates an exposure and says nothing whatever
    // about a count of decays, so an activity report carrying it would end with a paragraph
    // about a metric it never computed.
    if (spec.metric == Metric::Exposure && dataIndex >= 0 &&
        data.unmodeledPhotonFraction(dataIndex) > 0.05) {
      it->second.flags |= kFlagUnmodeledContinuum;
    }
  }

  const int nC = static_cast<int>(buckets.size());

  ResponseTable table;
  table.metric = spec.metric;
  table.aggregate = spec.aggregate;
  table.domain = domain;
  table.unit = spec.unit;
  table.contributors.resize(static_cast<std::size_t>(nC));
  table.labels.resize(static_cast<std::size_t>(nC));
  table.flags.assign(static_cast<std::size_t>(nC), kFlagNone);
  table.values.assign(static_cast<std::size_t>(nT) * static_cast<std::size_t>(nC), 0.0);
  table.totals.assign(static_cast<std::size_t>(nT), 0.0);

  for (const auto& [key, bucket] : buckets) {
    const std::size_t c = static_cast<std::size_t>(bucket.column);
    table.contributors[c] = ContributorId{key, bucket.dominant};
    table.labels[c] = bucketLabel(spec.aggregate, key, bucket.dominant);
    table.flags[c] = bucket.flags;
  }

  const double scale = unitScale(spec.unit) * domainScale(spec.metric, domain);
  for (int k = 0; k < nT; ++k) {
    const std::size_t base = static_cast<std::size_t>(k) * static_cast<std::size_t>(nC);
    double total = 0.0;
    for (int i = 0; i < nNuc; ++i) {
      const double value =
          weightedByTime[static_cast<std::size_t>(k)][static_cast<std::size_t>(i)] * scale;
      table.values[base + static_cast<std::size_t>(columnOf[static_cast<std::size_t>(i)])] += value;
      total += value;
    }
    table.totals[static_cast<std::size_t>(k)] = total;
  }

  return table;
}

// One column per discrete photon line, rather than per nuclide.
//
// A line's weight is lambda_i * y_ij * k(E_j): the emitter's decay rate, the photons per decay
// at that energy, and the geometry coefficient for that energy. Multiplying by the emitter's
// atom count gives the exposure that one line contributes -- so the table is built from the
// same atoms as every other aggregate, only weighted more finely.
//
// A full evaluation carries 86000 lines, and a fission seed reaches thousands of emitters, so
// the columns are thresholded: a line contributing less than kLineFloor of its own emitter's
// exposure is dropped. Relative to the EMITTER rather than to the global total, deliberately --
// a global threshold would erase the entire spectrum of every minor nuclide, and the question
// "which line dominates THIS nuclide" is one people ask.
ResponseTable assembleLines(const NuclearData& data, std::span<const std::int64_t> keys,
                            const std::vector<std::vector<double>>& atomsByTime,
                            const ResponseSpec& spec, Domain domain) {
  constexpr double kLineFloor = 1.0e-6;

  const int nNuc = static_cast<int>(keys.size());
  const int nT = static_cast<int>(atomsByTime.size());

  struct Column {
    std::int64_t emitterKey = 0;
    double energyEv = 0.0;
    double weight = 0.0;  // lambda * intensity * k(E)
    int nuclide = 0;      // index into the result's space
    int flags = kFlagNone;
  };
  std::vector<Column> columns;

  for (int i = 0; i < nNuc; ++i) {
    const Zai zai = Zai::fromKey(keys[static_cast<std::size_t>(i)]);
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

    const double emitterTotal = exposure::exposureRatePerBecquerel(lines, spec.geometry);
    const double floor = emitterTotal * kLineFloor;
    const int flags =
        data.unmodeledPhotonFraction(dataIndex) > 0.05 ? kFlagUnmodeledContinuum : kFlagNone;

    for (const GammaLine& line : lines) {
      const double perBecquerel =
          line.intensity * exposure::pointExposureCoeff(line.energyEv, spec.geometry);
      if (perBecquerel <= 0.0 || perBecquerel < floor) {
        continue;
      }
      columns.push_back(Column{zai.key(), line.energyEv, lambda * perBecquerel, i, flags});
    }
  }

  const int nC = static_cast<int>(columns.size());

  ResponseTable table;
  table.metric = spec.metric;
  table.aggregate = spec.aggregate;
  table.domain = domain;
  table.unit = spec.unit;
  table.contributors.resize(static_cast<std::size_t>(nC));
  table.labels.resize(static_cast<std::size_t>(nC));
  table.flags.assign(static_cast<std::size_t>(nC), kFlagNone);
  table.values.assign(static_cast<std::size_t>(nT) * static_cast<std::size_t>(nC), 0.0);
  table.totals.assign(static_cast<std::size_t>(nT), 0.0);

  for (int c = 0; c < nC; ++c) {
    const Column& column = columns[static_cast<std::size_t>(c)];
    const Zai emitter = Zai::fromKey(column.emitterKey);
    table.contributors[static_cast<std::size_t>(c)] =
        ContributorId{column.emitterKey, column.emitterKey, column.energyEv};
    table.labels[static_cast<std::size_t>(c)] = lineLabel(emitter, column.energyEv);
    table.flags[static_cast<std::size_t>(c)] = column.flags;
  }

  const double scale = unitScale(spec.unit) * domainScale(spec.metric, domain);
  for (int k = 0; k < nT; ++k) {
    const std::size_t base = static_cast<std::size_t>(k) * static_cast<std::size_t>(nC);
    double total = 0.0;
    for (int c = 0; c < nC; ++c) {
      const Column& column = columns[static_cast<std::size_t>(c)];
      const double value =
          column.weight *
          atomsByTime[static_cast<std::size_t>(k)][static_cast<std::size_t>(column.nuclide)] *
          scale;
      table.values[base + static_cast<std::size_t>(c)] = value;
      total += value;
    }
    table.totals[static_cast<std::size_t>(k)] = total;
  }

  return table;
}

// Fraction of the emitted photon energy rate that lives in a spectrum NuSIFT does not model.
// Activity-weighted, so a nuclide with a large unmodelled fraction but negligible activity
// contributes negligibly -- which is the whole point of reporting a magnitude rather than a
// count of flagged nuclides.
std::vector<double> unmodeledEnergyFractions(const NuclearData& data,
                                             std::span<const std::int64_t> keys,
                                             const std::vector<std::vector<double>>& atomsByTime) {
  const int nNuc = static_cast<int>(keys.size());
  std::vector<double> modelled(static_cast<std::size_t>(nNuc), 0.0);
  std::vector<double> missing(static_cast<std::size_t>(nNuc), 0.0);

  for (int i = 0; i < nNuc; ++i) {
    const int index = data.indexOfKey(keys[static_cast<std::size_t>(i)]);
    if (index < 0) {
      continue;
    }
    const double lambda = data.decayConstant(index);
    if (lambda <= 0.0) {
      continue;
    }
    modelled[static_cast<std::size_t>(i)] = lambda * discretePhotonEnergyEv(data.lines(index));
    missing[static_cast<std::size_t>(i)] = lambda * data.continuumPhotonEv(index);
  }

  std::vector<double> fractions(atomsByTime.size(), 0.0);
  for (std::size_t k = 0; k < atomsByTime.size(); ++k) {
    double modelledTotal = 0.0;
    double missingTotal = 0.0;
    for (int i = 0; i < nNuc; ++i) {
      const double atoms = atomsByTime[k][static_cast<std::size_t>(i)];
      modelledTotal += modelled[static_cast<std::size_t>(i)] * atoms;
      missingTotal += missing[static_cast<std::size_t>(i)] * atoms;
    }
    const double total = modelledTotal + missingTotal;
    fractions[k] = total > 0.0 ? missingTotal / total : 0.0;
  }
  return fractions;
}

// --- naming a contributor to pin ---------------------------------------------

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
  // A store with no photon lines cannot answer an exposure question at all. Returning zeros
  // would be indistinguishable from "nothing here emits photons", which is a different and
  // much more alarming statement.
  if (spec.aggregate == Aggregate::GammaLine && spec.metric != Metric::Exposure) {
    throw InputError(tagged(kModule,
                            "ranking by gamma line only makes sense for exposure -- a photon "
                            "line has no activity of its own, it is a way its emitter's decays "
                            "get out"));
  }
  if (spec.metric == Metric::Exposure && !data.hasPhotonLines()) {
    throw InputError(tagged(kModule,
                            "this data store carries no photon lines, so exposure cannot be "
                            "computed. Stage from ENDF decay tapes, which carry the discrete "
                            "spectra, or rank by activity instead"));
  }
}

}  // namespace

bool unitSuitsDomain(Unit unit, Domain domain) {
  switch (unit) {
    case Unit::Becquerel:
    case Unit::Curie:
    case Unit::RoentgenPerHour:
    case Unit::GrayPerHour:
    case Unit::SievertPerHour:
      return domain == Domain::Instant;
    case Unit::Decays:
    case Unit::Roentgen:
    case Unit::Gray:
    case Unit::Sievert:
      return domain == Domain::Interval;
  }
  return false;
}

bool unitSuitsMetric(Unit unit, Metric metric) {
  switch (unit) {
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
  }
  return false;
}

const char* unitName(Unit unit) {
  switch (unit) {
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
  if (metric == Metric::Exposure) {
    return domain == Domain::Interval ? Unit::Roentgen : Unit::RoentgenPerHour;
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
                                            "; exposure: " + spellingsFor(Metric::Exposure) + ")"));
}

const char* metricName(Metric metric) {
  switch (metric) {
    case Metric::Activity:
      return "activity";
    case Metric::Exposure:
      return "exposure";
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

ResponseTable buildResponse(const NuclearData& data, const DecayResult& result,
                            const ResponseSpec& spec) {
  requireUsableSpec(data, spec, Domain::Instant);

  const int nNuc = result.nuclideCount();
  const int nT = result.timeCount();

  std::vector<double> weight(static_cast<std::size_t>(nNuc), 0.0);
  for (int i = 0; i < nNuc; ++i) {
    const int index = data.indexOfKey(result.nuclideKeys[static_cast<std::size_t>(i)]);
    weight[static_cast<std::size_t>(i)] = index >= 0 ? weightFor(spec, data, index) : 0.0;
  }

  std::vector<std::vector<double>> weighted(static_cast<std::size_t>(nT));
  std::vector<std::vector<double>> rawAtoms(static_cast<std::size_t>(nT));
  for (int k = 0; k < nT; ++k) {
    const std::span<const double> atoms = result.atomsAt(k);
    std::vector<double> row(static_cast<std::size_t>(nNuc), 0.0);
    std::vector<double> raw(atoms.begin(), atoms.end());
    for (int i = 0; i < nNuc; ++i) {
      row[static_cast<std::size_t>(i)] =
          weight[static_cast<std::size_t>(i)] * atoms[static_cast<std::size_t>(i)];
    }
    weighted[static_cast<std::size_t>(k)] = std::move(row);
    rawAtoms[static_cast<std::size_t>(k)] = std::move(raw);
  }

  // The line assembler applies its own per-line weights, so it takes atoms directly rather
  // than the per-nuclide weighted values every other aggregate folds together.
  ResponseTable table =
      spec.aggregate == Aggregate::GammaLine
          ? assembleLines(data, result.nuclideKeys, rawAtoms, spec, Domain::Instant)
          : assemble(data, result.nuclideKeys, weighted, spec, Domain::Instant);
  table.times = result.times;
  if (spec.metric == Metric::Exposure) {
    table.unmodeledEnergyFraction = unmodeledEnergyFractions(data, result.nuclideKeys, rawAtoms);
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

  const int nNuc = static_cast<int>(keys.size());
  const std::vector<double> integratedAtoms(integral.begin(), integral.end());

  ResponseTable table;
  if (spec.aggregate == Aggregate::GammaLine) {
    // Same split as the instantaneous builder, for the same reason: the line assembler
    // applies its own per-line weights, so it takes atom-seconds directly rather than the
    // per-nuclide weighted value every other aggregate folds together.
    table = assembleLines(data, keys, std::vector<std::vector<double>>{integratedAtoms}, spec,
                          Domain::Interval);
  } else {
    std::vector<double> row(static_cast<std::size_t>(nNuc), 0.0);
    for (int i = 0; i < nNuc; ++i) {
      const int index = data.indexOfKey(keys[static_cast<std::size_t>(i)]);
      // lambda against atom-seconds: a dimensionless count of decays over the window.
      const double weight = index >= 0 ? weightFor(spec, data, index) : 0.0;
      row[static_cast<std::size_t>(i)] = weight * integratedAtoms[static_cast<std::size_t>(i)];
    }
    table = assemble(data, keys, std::vector<std::vector<double>>{std::move(row)}, spec,
                     Domain::Interval);
  }

  table.times = {t1};
  table.timeEnds = {t2};
  if (spec.metric == Metric::Exposure) {
    table.unmodeledEnergyFraction =
        unmodeledEnergyFractions(data, keys, std::vector<std::vector<double>>{integratedAtoms});
  }
  return table;
}

}  // namespace nusift
