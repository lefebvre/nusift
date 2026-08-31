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
#include <vector>

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

}  // namespace nusift
