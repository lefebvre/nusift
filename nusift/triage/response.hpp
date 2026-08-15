#pragma once
/**
 * @file
 * @brief Turns a decay result into a table of per-contributor response values.
 * @ingroup triage
 */
//
// This is where the engine's two matrices become the numbers a user asked for, and it is
// deliberately the ONLY place that knows how a metric is computed.
//
// Every response is
//
//     values(k, c) = w_c * X(k, parent(c))
//
// where X is `atoms` for an instantaneous metric and `integratedAtoms` for a time-integrated
// one, and w_c is a fixed per-contributor weight. Activity uses w = lambda; exposure will use
// w = lambda * (a photon sum). Aggregating by mass chain or element is a sum of nuclide
// columns, not a different calculation.
//
// The units fall out of the domain rather than being chosen: lambda*n is a rate in Bq, while
// lambda times an integral over time is a COUNT of decays. Conflating the two is the easiest
// mistake to make in a tool that reports both, so Domain determines which units are even
// offered.
//
#include <cstdint>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "nusift/core/nuclide.hpp"
#include "nusift/engine/decay_result.hpp"
#include "nusift/exposure/point_source.hpp"

namespace nusift {

class NuclearData;

enum class Metric {
  Activity,  // decays per second, or a count of decays over an interval
  Exposure,  // photon exposure rate, or exposure accrued over an interval
};

enum class Aggregate {
  Nuclide,
  MassChain,  // isobars: everything sharing a mass number A
  Element,    // everything sharing an atomic number Z
  // One column per discrete photon line. Exposure only -- a line has no activity of its own,
  // it is a way its emitter's decays get out. This is the aggregate no tool that collapses a
  // spectrum to a single per-nuclide constant can offer, and it is what a shielding or
  // detector question actually needs: not "which nuclide", but "which energy".
  GammaLine,
};

enum class Domain {
  Instant,   // evaluated at a point in time
  Interval,  // integrated over a window
};

enum class Unit {
  // Activity. A rate and a count, which is the same distinction Domain draws.
  Becquerel,
  Curie,
  Decays,  // interval only

  // Exposure. Likewise a rate and an accrued total. Gray and sievert are the same number here
  // -- air kerma with a photon radiation weighting factor of 1 -- but they name different
  // quantities, and a report that says "Sv" when it means absorbed dose in air invites exactly
  // the misreading the distinction exists to prevent.
  RoentgenPerHour,
  GrayPerHour,
  SievertPerHour,
  Roentgen,  // interval only
  Gray,      // interval only
  Sievert,   // interval only
};

// Which metric a unit can express. Reporting exposure in becquerel is not a rounding error,
// it is a category error, so it is refused at the boundary alongside the domain check.
bool unitSuitsMetric(Unit unit, Metric metric);

// Whether a unit can express a given domain. Bq is a rate and cannot describe an interval
// total; a count of decays cannot describe an instant. Checked at the boundary so the error
// names the mismatch instead of silently reporting a number in the wrong dimension.
bool unitSuitsDomain(Unit unit, Domain domain);
const char* unitName(Unit unit);

// Parse a unit spelling, case-insensitively against the names unitName() prints -- so "Bq",
// "bq", "gy/h" and "Sv" all resolve. This lives in the library rather than in each front end
// because the CLI and the Python binding both need it, and two tables drift: `--units bq`
// worked while `units="bq"` raised, for no reason a user could see.
bool parseUnit(std::string_view text, Unit& out);

// The unit a metric and domain fall back to when the user names none. It has to satisfy both,
// or the most ordinary invocation would fail on a unit nobody chose.
Unit defaultUnit(Metric metric, Domain domain);

// Resolve a user-supplied spelling, taking the default for empty text. Throws InputError
// naming the spellings that would have worked.
Unit requireUnit(std::string_view text, Metric metric, Domain domain);

const char* metricName(Metric metric);
const char* aggregateName(Aggregate aggregate);

// Identity of one column. Which fields are meaningful depends on the aggregate, and the rest
// stay zero rather than being overloaded with a second meaning.
struct ContributorId {
  std::int64_t key = 0;  // ZAI key, or a mass number A, or an atomic number Z
  // For an aggregate, the single nuclide contributing most within it at the time the table
  // was built. A mass chain named only by its number tells the user nothing actionable;
  // "A=140 (La-140)" tells them what to look at. For a gamma line, the emitting nuclide.
  std::int64_t dominantMemberKey = 0;
  // GammaLine only: the photon energy this column represents. Carried separately from the key
  // because a line is identified by its emitter AND its energy, and a consumer binning a
  // spectrum needs the number rather than the label.
  double lineEnergyEv = 0.0;
};

// Set when a contributor's value is known to understate the truth. Reported alongside the
// ranking rather than folded into it, because silently adjusting a number the data does not
// support would be worse than flagging it.
enum ContributorFlag : int {
  kFlagNone = 0,
  kFlagUnmodeledContinuum = 1 << 0,  // photon energy in a continuum NuSIFT does not model
};

struct ResponseTable {
  Metric metric = Metric::Activity;
  Aggregate aggregate = Aggregate::Nuclide;
  Domain domain = Domain::Instant;
  Unit unit = Unit::Becquerel;

  std::vector<double> times;     // [nT]
  std::vector<double> timeEnds;  // [nT], interval domain only

  std::vector<ContributorId> contributors;  // [nC]
  std::vector<std::string> labels;          // [nC], display names
  std::vector<int> flags;                   // [nC]

  std::vector<double> values;  // [nT * nC], row-major by time, already in `unit`
  // Sum over ALL contributors, never over a truncated prefix. Ranking reports coverage
  // against this, so a top-N view can never misrepresent itself as the whole.
  std::vector<double> totals;  // [nT]

  // Exposure only, empty otherwise: the fraction of emitted photon ENERGY that NuSIFT does not
  // model, at each time. [nT]
  //
  // This is what turns "356 nuclides are understated" into a number anyone can act on. A count
  // says nothing about magnitude -- three hundred negligible nuclides and three dominant ones
  // look identical -- whereas this says how much of the photon output is missing.
  //
  // It is an ENERGY fraction, not an exposure fraction, and deliberately so. Bounding the
  // missing exposure would require assuming an energy for photons whose energies are precisely
  // what is unknown. Emitted energy is exactly computable from what the evaluation does give,
  // and exposure tracks it closely enough over the range that matters to make it the right
  // indicator: a 10% energy shortfall means the exposure is understated by roughly that much.
  std::vector<double> unmodeledEnergyFraction;

  int timeCount() const { return static_cast<int>(times.size()); }
  int contributorCount() const { return static_cast<int>(contributors.size()); }

  std::span<const double> valuesAt(int timeIndex) const {
    const std::size_t n = contributors.size();
    return std::span<const double>(values.data() + static_cast<std::size_t>(timeIndex) * n, n);
  }
};

// Resolve a contributor named the way a user writes one, against the aggregate `table` was
// built with: "Cs-137" for a nuclide, "A=140" or "140" for a mass chain, "Cs" or "55" for an
// element. A nuclide name is accepted whatever the aggregate and resolves to the bucket that
// nuclide falls in, so `--pin Cs-137` means the same thing whether the table ranks nuclides or
// the mass chains they sit in -- someone who knows a nuclide name should not have to work out
// which isobar it belongs to in order to follow it.
//
// Returns the contributor key, which is what RankRequest::pinned holds. Throws InputError when
// the spelling names nothing of the right kind, and when it names something this table does not
// carry: an inventory whose chain never reaches Cs-137 cannot pin it, and saying so is better
// than a row of zeros that reads like an answer.
//
// A gamma-line table has one column per line and several per emitter, so a key there names an
// EMITTER, and pinning it pins every line of that emitter the table carries. Pinning one line
// out of a spectrum is not offered: what is understated or overlooked is the nuclide, and which
// of its lines you are looking at does not change that.
std::int64_t requirePin(const ResponseTable& table, std::string_view text);

struct ResponseSpec {
  Metric metric = Metric::Activity;
  Aggregate aggregate = Aggregate::Nuclide;
  Unit unit = Unit::Becquerel;

  // Used only by Metric::Exposure. Held here rather than passed separately because the
  // geometry is part of what the resulting numbers MEAN -- an exposure table without the
  // distance it was computed at is not interpretable, and keeping them together makes it hard
  // to report one without the other.
  exposure::PointSourceGeometry geometry;
};

// Build a table of instantaneous responses at each time in `result`.
ResponseTable buildResponse(const NuclearData& data, const DecayResult& result,
                            const ResponseSpec& spec);

// Build a single-row table for one integrated interval. `integral` is per-nuclide
// atom-seconds over [t1, t2] in the index space of `keys`, exactly as intervalIntegral
// produces them.
ResponseTable buildIntervalResponse(const NuclearData& data, std::span<const std::int64_t> keys,
                                    std::span<const double> integral, double t1, double t2,
                                    const ResponseSpec& spec);

}  // namespace nusift
