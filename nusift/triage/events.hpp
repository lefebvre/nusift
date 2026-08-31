#pragma once
/**
 * @file
 * @brief When a trajectory crosses a level, peaks, or holds inside a window.
 * @ingroup triage
 */
//
// Ranking answers "what dominates now". Forecasting answers "who leads, and when that
// changes". This answers the third question a trajectory carries and the one people arrive
// with a number in hand for: WHEN does a quantity reach a value someone cares about. When does
// the total fall below a release limit, when does an ingrowth-fed curve stop getting worse,
// how long does a field stay above a level, when does a measured ratio leave the band a
// scaling factor was calibrated in.
//
// Every answer here is a LOCATED EVENT, and a located event is the one quantity in NuSIFT that
// is not exact. An inventory and an interval integral are closed form within the decay model.
// A crossing is a root of a curve that was SAMPLED, and two consequences run through this
// whole interface:
//
//   The grid decides what is seen. An event is searched for only inside a bracket the samples
//   actually straddle. An excursion that rises and falls back between two consecutive samples
//   leaves no sign change behind and is not found -- not "unlikely to be found", not found --
//   and no amount of refinement inside the intervals that were observed will reveal it. Two
//   crossings inside one interval cancel the same way. The remedy is a denser grid, and
//   choosing it is the caller's decision, so nothing here silently subdivides on its own.
//
//   The bracket is the honest uncertainty. Every event carries the grid interval that observed
//   it and the width its location was narrowed to. Reporting a crossing as a bare instant
//   would claim a precision the sampling does not support.
//
// The same rule the forecast boundaries follow: refine inside an observed bracket, never
// invent one. dominanceWindows() in forecast.hpp locates a leader change this way already;
// this generalises it to a level, an extremum, and a window, and shares the refinement so the
// two can never disagree about where a crossing is.
//
#include <cstdint>
#include <functional>
#include <optional>
#include <span>
#include <vector>

#include "nusift/engine/decay_engine.hpp"
#include "nusift/engine/inventory.hpp"
#include "nusift/triage/response.hpp"

namespace nusift {

// The curve an event search runs over: what the grid saw, and -- when the caller can supply
// it -- a way to ask what the curve does between two samples.
struct EventSeries {
  std::vector<double> times;   // [nT], strictly increasing
  std::vector<double> values;  // [nT]

  // Evaluate the same quantity at an arbitrary time within the grid's span.
  //
  // Supplied, an event is refined by real evaluations. For a decay response each one is a
  // single-time solve, which costs what any other single-time answer costs and parallelises
  // the same way, so a crossing can be placed far more tightly than the grid that found it.
  //
  // Left empty, the event is INTERPOLATED inside its bracket instead and reported as known
  // only to that bracket's width. That is the difference between narrowing a location and
  // guessing within an interval, and `TrajectoryEvent::refined` says which happened rather
  // than leaving the caller to infer it from a tolerance.
  std::function<double(double)> evaluate;
};

// How tightly to place an event, and how hard to try.
//
// The tolerance is relative to the event time because response grids are log-spaced and span
// decades: a fixed second is absurdly tight at thirty years and uselessly loose at thirty
// seconds. The absolute floor exists for grids that start at zero.
struct EventTolerance {
  double relative = 1.0e-6;
  double absoluteSeconds = 0.0;
  int maxIterations = 100;
};

enum class EventKind {
  Rising,   // the series crossed the level going up
  Falling,  // going down
  Maximum,
  Minimum,
};

const char* eventKindName(EventKind kind);

struct TrajectoryEvent {
  EventKind kind = EventKind::Rising;
  double timeSeconds = 0.0;
  double value = 0.0;  // the series value at timeSeconds

  // The consecutive grid samples that observed the event. The search never leaves them, so
  // this pair is both the provenance of the event and the bound on how wrong it can be.
  double bracketStartSeconds = 0.0;
  double bracketEndSeconds = 0.0;

  // The width the location was narrowed to -- the honest error bar on `timeSeconds`. Equal to
  // the whole bracket when the series carried no evaluator, because an interpolated event is
  // known no better than the interval it sits in.
  double locatedToSeconds = 0.0;

  // True when real evaluations narrowed the bracket, false when the location is an
  // interpolation between two samples.
  bool refined = false;

  // False when refinement hit maxIterations before reaching the tolerance. The event is still
  // real and still inside its bracket -- only more loosely placed than was asked for.
  bool converged = true;
};

// A stretch over which the series holds on one side of a level.
struct LevelWindow {
  double startSeconds = 0.0;
  double endSeconds = 0.0;

  // False when the grid began already inside the window: the entry happened before the first
  // sample, or never happened at all. The start is then the grid's first time, which is a
  // bound and not a crossing, and saying so is the difference between "it was already above
  // when we started looking" and "it rose above at this instant".
  bool entryObserved = false;

  // False when the grid ended still inside. The window is open at that end: it has not been
  // observed to close, which is not the same as closing at the last sample.
  bool exitObserved = false;
};

// Every level crossing the grid observed, in time order.
//
// A crossing is recorded where consecutive samples fall on opposite sides of `level`, which
// makes "the series is below the level" a boolean that flips exactly once per reported event.
// A sample sitting exactly on the level counts as not-below, so a curve that touches and
// retreats reports the touch once rather than twice or not at all.
std::vector<TrajectoryEvent> crossings(const EventSeries& series, double level,
                                       const EventTolerance& tolerance = {});

// The first and last crossings, empty when the grid observed none. Thin wrappers over
// crossings() rather than cheaper searches of their own: the cost is the sampling, which has
// already happened, and one code path cannot disagree with itself about where a root is.
std::optional<TrajectoryEvent> firstCrossing(const EventSeries& series, double level,
                                             const EventTolerance& tolerance = {});
std::optional<TrajectoryEvent> lastCrossing(const EventSeries& series, double level,
                                            const EventTolerance& tolerance = {});

// Every interior maximum and minimum the grid resolved, in time order.
//
// An extremum needs three samples to be seen at all -- a turn is a sample higher (or lower)
// than both its neighbours -- so a peak in the first or last interval of the grid is outside
// what the sampling can distinguish from a monotone run, and is not reported. A curve that is
// flat across three samples is not turning and is not an extremum.
//
// This is the "is waiting actually helping" question: the sign change in the derivative of an
// ingrowth-fed curve is a minimum, and freshly separated Sr-90, sealed Ra-226 and Pu-241 all
// have one.
std::vector<TrajectoryEvent> extrema(const EventSeries& series,
                                     const EventTolerance& tolerance = {});

// The windows over which the series holds at or above `level`, with each boundary refined the
// way a crossing is. Windows open at the grid's edge are flagged rather than clipped silently.
//
// This is the stay-time and re-entry shape: "when is it above the level I care about, and for
// how long", answered as intervals rather than as a list of instants the caller has to pair up.
std::vector<LevelWindow> windowsAbove(const EventSeries& series, double level,
                                      const EventTolerance& tolerance = {});

// The same below the level -- the complement, and the form a "when is it safe" question takes.
std::vector<LevelWindow> windowsBelow(const EventSeries& series, double level,
                                      const EventTolerance& tolerance = {});

// --- series drawn from a response table --------------------------------------
//
// The table already holds the sampled curve, so these are the cheap path: they build a series
// with no evaluator, and every event they produce is interpolated inside its grid bracket. To
// refine instead, fill in EventSeries::evaluate with something that re-solves at a time.

// The response total against time -- the curve a limit, a release threshold or a re-entry
// criterion is compared with.
EventSeries totalSeries(const ResponseTable& table);

// One contributor's column.
EventSeries contributorSeries(const ResponseTable& table, int contributor);

// The ratio of two contributors -- the clock form. Zr-95/Nb-95 or Ba-140/La-140 for time since
// fission, Cs-137/Co-60 for a scaling factor's drift.
//
// Samples before the denominator becomes positive are dropped rather than reported as infinite
// or zero: a ratio clock does not start until the nuclide it divides by exists, and the first
// samples of an ingrowth are exactly that case. If the denominator returns to zero later the
// ratio is discontinuous inside the grid, no bracket across the gap would mean anything, and
// this throws InputError rather than searching over it.
EventSeries ratioSeries(const ResponseTable& table, int numerator, int denominator);

// --- refining an event by re-solving ------------------------------------------
//
// Everything above locates an event inside a bracket the grid observed. Whether that location
// is NARROWED or merely INTERPOLATED turns on one thing: whether EventSeries::evaluate can say
// what the curve does between two samples. For a response curve that evaluation is a
// single-time solve, costing what any other single-time answer costs, and this is the object
// that performs it.
//
// It is deliberately not automatic. Refinement multiplies solves by the number of events and
// by the iterations each one takes, and whether that is worth spending is the caller's
// judgement rather than this layer's. What is not the caller's judgement is being honest about
// which happened: `TrajectoryEvent::refined` says so either way, and an interpolated event is
// reported as known no better than its bracket.
//
// The evaluator has to be built from the same data, inventory and spec the table was, or it
// would narrow an event on one curve using values from another. The builders below check what
// can be checked -- metric, aggregate, unit, geometry -- and refuse a mismatch rather than
// silently mixing two curves.
//
// LIFETIME. A series built with an evaluator holds a reference to it, and the evaluator holds
// references to the data and the inventory. All three have to outlive every search that uses
// the series.
class ResponseEvaluator {
public:
  ResponseEvaluator(const NuclearData& data, const Inventory& inventory, const ResponseSpec& spec,
                    const DecayOptions& options = {});

  // The response total at `timeSeconds`, in the spec's unit.
  double total(double timeSeconds) const;

  // What one contributor holds at `timeSeconds`, named by the identity the table carries rather
  // than by a column index. Two tables built from the same seed agree on identities whatever
  // times they were built at; a column index means nothing away from the table it came from.
  double value(const ContributorId& id, double timeSeconds) const;

  const ResponseSpec& spec() const { return spec_; }

  // Solves spent so far. What refinement cost is not visible in the events it produced, so a
  // caller that wants to report the price has to ask for it.
  int solves() const { return solves_; }

private:
  const ResponseTable& tableAt(double timeSeconds) const;

  const NuclearData& data_;
  const Inventory& inventory_;
  ResponseSpec spec_;
  DecayOptions options_;

  // The last instant solved for. A ratio asks for two columns at the same time, and every
  // search asks about the same instant more than once; without this the same solve would be
  // paid for twice or more.
  mutable ResponseTable cached_;
  mutable double cachedTime_ = 0.0;
  mutable bool hasCache_ = false;
  mutable int solves_ = 0;
};

// The same three series, able to refine. Identical sampling to the versions above -- the values
// still come from the table, so nothing about what the grid SAW changes -- with `evaluate`
// filled in, so an event inside a bracket can be narrowed instead of interpolated.
EventSeries totalSeries(const ResponseTable& table, const ResponseEvaluator& evaluator);
EventSeries contributorSeries(const ResponseTable& table, int contributor,
                              const ResponseEvaluator& evaluator);
EventSeries ratioSeries(const ResponseTable& table, int numerator, int denominator,
                        const ResponseEvaluator& evaluator);

// --- the fixed-duration task --------------------------------------------------
//
// "What does a one-hour job cost, and when should it be done" is not a question about the rate
// curve. It is a question about the total ACCRUED over a window of fixed length, as a function
// of when that window starts, and no point of that curve is a point of the rate curve.
//
// Every sample is an exact interval integral, so the curve inherits the property that makes
// inverting it worth doing at all: there is no quadrature error inside a window to reason
// about, only the sampling of the start times. And once it is a series, every search above
// applies to it unchanged:
//
//   extrema()      -- the best and the worst time to start. For an ingrowth-fed mixture the
//                     best start is not the latest one, which is the whole reason to ask.
//   crossings()    -- when the job first fits a budget, and when it stops fitting.
//   windowsBelow() -- every stretch of start times the budget allows.
//
// The unit has to be an interval unit, because the values are accrued totals: roentgen rather
// than roentgen per hour, decays rather than becquerel. buildIntervalResponse() refuses the
// rest, and this leaves that refusal where it already lives.
//
// COST, stated because it is unlike everything else here: each sample is an interval integral,
// which is two to three solves rather than one, and no work is shared between samples. A
// sixty-point start grid is a couple of hundred solves -- seconds rather than milliseconds --
// and `refine` buys tighter events with more of them. Nothing about that is hidden by making
// the curve look like any other series.
EventSeries taskSeries(const NuclearData& data, const Inventory& inventory,
                       const ResponseSpec& spec, std::span<const double> startTimes,
                       double durationSeconds, const DecayOptions& options = {},
                       bool refine = false);

}  // namespace nusift
