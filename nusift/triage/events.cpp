#include "nusift/triage/events.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <string>
#include <vector>

#include "nusift/core/error.hpp"

namespace nusift {
namespace {

constexpr const char* kModule = "trajectory";

// A grid has to be usable as a set of brackets before anything can be located in it: two
// samples at minimum, paired values, and strictly increasing times. A repeated time would make
// a zero-width bracket that no refinement can narrow and no interpolation can divide by.
void requireUsableSeries(const EventSeries& series) {
  if (series.times.size() != series.values.size()) {
    throw InputError(tagged(kModule, "series times and values are different lengths"));
  }
  if (series.times.size() < 2) {
    throw InputError(tagged(kModule,
                            "an event needs at least two samples to sit between; a single "
                            "time can be evaluated but nothing can be located in it"));
  }
  for (std::size_t k = 1; k < series.times.size(); ++k) {
    if (!(series.times[k] > series.times[k - 1])) {
      throw InputError(tagged(kModule, "series times must be strictly increasing"));
    }
  }
}

// The width an event's location has to be narrowed to. Relative to the time itself because
// response grids span decades: one second is far tighter than a thirty-year answer can
// support and far looser than a thirty-second one deserves.
double toleranceAt(const EventTolerance& tolerance, double time) {
  return std::max(tolerance.absoluteSeconds, tolerance.relative * std::abs(time));
}

// Where the series meets `level` between two samples, without evaluating anything.
//
// Over one grid interval a decay response is close to exponential, so log(value) is close to
// linear in time and the crossing of a positive level has a closed form there. This is the
// same reasoning -- and the same accuracy -- as the interpolated boundaries in
// dominanceWindows(). Falls back to a straight line when a log is not available, which is the
// case for a level or a sample at or below zero, and for a ratio series that legitimately
// passes through zero.
double interpolateCrossing(double t0, double t1, double v0, double v1, double level) {
  if (level > 0.0 && v0 > 0.0 && v1 > 0.0) {
    const double span = std::log(v0 / v1);
    if (span != 0.0) {
      const double fraction = std::log(v0 / level) / span;
      return t0 + (t1 - t0) * std::clamp(fraction, 0.0, 1.0);
    }
  }
  const double span = v1 - v0;
  if (span == 0.0) {
    return 0.5 * (t0 + t1);
  }
  const double fraction = (level - v0) / span;
  return t0 + (t1 - t0) * std::clamp(fraction, 0.0, 1.0);
}

// Narrow a bracketed crossing by evaluating the curve inside it.
//
// Illinois -- regula falsi with the stale endpoint's residual halved whenever the same side is
// kept twice. The bracket is maintained at every step, so the root never escapes the interval
// the grid observed, and halving cures the one-sided stalling that makes plain false position
// converge arbitrarily slowly on a convex curve. A decay curve is convex in exactly that way,
// so the fix is not academic here.
//
// `located` receives the width of the final bracket: the honest error bar, not the tolerance
// that was asked for.
double refineCrossing(const std::function<double(double)>& evaluate, double t0, double t1,
                      double v0, double v1, double level, const EventTolerance& tolerance,
                      double& located, bool& converged, int& evaluations) {
  double a = t0;
  double b = t1;
  double fa = v0 - level;
  double fb = v1 - level;
  double best = interpolateCrossing(t0, t1, v0, v1, level);
  int side = 0;
  converged = false;
  evaluations = 0;

  for (int iteration = 0; iteration < tolerance.maxIterations; ++iteration) {
    if (b - a <= toleranceAt(tolerance, 0.5 * (a + b))) {
      converged = true;
      break;
    }

    double c = 0.0;
    const double denominator = fb - fa;
    if (denominator == 0.0 || !std::isfinite(denominator)) {
      c = 0.5 * (a + b);
    } else {
      c = (a * fb - b * fa) / denominator;
    }
    // Keep every step strictly inside the bracket. A false-position step cannot leave it in
    // exact arithmetic, but rounding at the ends of a decade-wide interval can put it on a
    // boundary, and a step that lands on an endpoint stops making progress.
    const double margin = 0.25 * (b - a);
    c = std::clamp(c, a + 1.0e-12 * margin, b - 1.0e-12 * margin);
    if (!std::isfinite(c) || !(c > a && c < b)) {
      c = 0.5 * (a + b);
    }

    const double fc = evaluate(c) - level;
    ++evaluations;
    best = c;
    if (fc == 0.0) {
      a = c;
      b = c;
      converged = true;
      break;
    }

    if ((fc < 0.0) == (fa < 0.0)) {
      a = c;
      fa = fc;
      if (side == -1) {
        fb *= 0.5;
      }
      side = -1;
    } else {
      b = c;
      fb = fc;
      if (side == 1) {
        fa *= 0.5;
      }
      side = 1;
    }
  }

  if (b - a <= toleranceAt(tolerance, 0.5 * (a + b))) {
    converged = true;
  }
  located = b - a;
  return std::clamp(best, a, b);
}

// Build one crossing event from the interval [k, k+1], refining it when the series can be
// evaluated and interpolating it when it cannot.
TrajectoryEvent locateCrossing(const EventSeries& series, int k, double level, bool rising,
                               const EventTolerance& tolerance) {
  const std::size_t i = static_cast<std::size_t>(k);
  const double t0 = series.times[i];
  const double t1 = series.times[i + 1];
  const double v0 = series.values[i];
  const double v1 = series.values[i + 1];

  TrajectoryEvent event;
  event.kind = rising ? EventKind::Rising : EventKind::Falling;
  event.bracketStartSeconds = t0;
  event.bracketEndSeconds = t1;
  event.value = level;

  if (series.evaluate) {
    double located = t1 - t0;
    bool converged = true;
    int evaluations = 0;
    event.timeSeconds = refineCrossing(series.evaluate, t0, t1, v0, v1, level, tolerance, located,
                                       converged, evaluations);
    event.locatedToSeconds = located;
    // A bracket already inside the tolerance is returned without ever evaluating, and that
    // location is an interpolation like any other. Saying otherwise would make the flag mean
    // "an evaluator was available" rather than "an evaluation placed this".
    event.refined = evaluations > 0;
    event.converged = converged;
  } else {
    event.timeSeconds = interpolateCrossing(t0, t1, v0, v1, level);
    // An interpolated event is known no better than the interval it sits in, whatever the
    // interpolation suggests to three decimal places.
    event.locatedToSeconds = t1 - t0;
    event.refined = false;
    event.converged = true;
  }
  return event;
}

// The abscissa a three-point parabola is fitted in.
//
// A parabola through three points is best conditioned when they are EVENLY spaced, and which
// abscissa achieves that is a property of the grid: response grids are log-spaced, where
// log-time evens them out, while a short or hand-built grid is usually arithmetic, where it
// does the opposite. So the choice follows the samples rather than being fixed -- each spacing
// is scored by how uneven it leaves the three points, as a dimensionless fraction of the span,
// and the more even one wins. A grid reaching zero has no log and is linear by default.
bool useLogAbscissa(double t0, double t1, double t2) {
  if (!(t0 > 0.0 && t1 > 0.0 && t2 > 0.0)) {
    return false;
  }
  const double span = t2 - t0;
  const double logSpan = std::log(t2 / t0);
  if (!(span > 0.0) || !(logSpan > 0.0)) {
    return false;
  }
  const double linearSkew = std::abs((t2 - t1) - (t1 - t0)) / span;
  const double logSkew = std::abs(std::log(t2 / t1) - std::log(t1 / t0)) / logSpan;
  return logSkew < linearSkew;
}

// The turning point of the parabola through three samples, clamped into the bracket. Used when
// there is no evaluator: a peak located by interpolation alone.
double interpolateExtremum(double t0, double t1, double t2, double y0, double y1, double y2,
                           double& valueAtVertex) {
  const bool useLog = useLogAbscissa(t0, t1, t2);
  const double x0 = useLog ? std::log(t0) : t0;
  const double x1 = useLog ? std::log(t1) : t1;
  const double x2 = useLog ? std::log(t2) : t2;

  const double d1 = x1 - x0;
  const double d2 = x1 - x2;
  const double num = d1 * d1 * (y1 - y2) - d2 * d2 * (y1 - y0);
  const double den = d1 * (y1 - y2) - d2 * (y1 - y0);

  double xv = x1;
  if (den != 0.0 && std::isfinite(num / den)) {
    xv = std::clamp(x1 - 0.5 * num / den, std::min(x0, x2), std::max(x0, x2));
  }

  // Lagrange value at the vertex, so the reported peak is the parabola's and not the middle
  // sample's -- the whole point of interpolating is that the true peak is not at a sample.
  const double l0 = ((xv - x1) * (xv - x2)) / ((x0 - x1) * (x0 - x2));
  const double l1 = ((xv - x0) * (xv - x2)) / ((x1 - x0) * (x1 - x2));
  const double l2 = ((xv - x0) * (xv - x1)) / ((x2 - x0) * (x2 - x1));
  valueAtVertex = l0 * y0 + l1 * y1 + l2 * y2;

  return useLog ? std::exp(xv) : xv;
}

// Narrow a bracketed turn by evaluating inside it.
//
// Golden section rather than a derivative method: the series is a callable with no gradient,
// and the section search needs only that the bracket hold one turn -- which is exactly what
// three samples straddling a peak establish. It reduces the bracket by a fixed factor every
// step and cannot leave it.
double refineExtremum(const std::function<double(double)>& evaluate, double lo, double hi,
                      bool maximum, const EventTolerance& tolerance, double& valueAtBest,
                      double& located, bool& converged) {
  const double invPhi = 1.0 / std::numbers::phi;
  const auto score = [&](double t) {
    const double value = evaluate(t);
    return maximum ? value : -value;
  };

  double x1 = hi - invPhi * (hi - lo);
  double x2 = lo + invPhi * (hi - lo);
  double f1 = score(x1);
  double f2 = score(x2);
  converged = false;

  for (int iteration = 0; iteration < tolerance.maxIterations; ++iteration) {
    if (hi - lo <= toleranceAt(tolerance, 0.5 * (lo + hi))) {
      converged = true;
      break;
    }
    if (f1 < f2) {
      lo = x1;
      x1 = x2;
      f1 = f2;
      x2 = lo + invPhi * (hi - lo);
      f2 = score(x2);
    } else {
      hi = x2;
      x2 = x1;
      f2 = f1;
      x1 = hi - invPhi * (hi - lo);
      f1 = score(x1);
    }
  }
  if (hi - lo <= toleranceAt(tolerance, 0.5 * (lo + hi))) {
    converged = true;
  }

  const double best = f1 >= f2 ? x1 : x2;
  const double bestScore = std::max(f1, f2);
  valueAtBest = maximum ? bestScore : -bestScore;
  located = hi - lo;
  return best;
}

TrajectoryEvent locateExtremum(const EventSeries& series, int k, bool maximum,
                               const EventTolerance& tolerance) {
  const std::size_t i = static_cast<std::size_t>(k);
  const double t0 = series.times[i - 1];
  const double t1 = series.times[i];
  const double t2 = series.times[i + 1];

  TrajectoryEvent event;
  event.kind = maximum ? EventKind::Maximum : EventKind::Minimum;
  event.bracketStartSeconds = t0;
  event.bracketEndSeconds = t2;

  if (series.evaluate) {
    double value = series.values[i];
    double located = t2 - t0;
    bool converged = true;
    const double best =
        refineExtremum(series.evaluate, t0, t2, maximum, tolerance, value, located, converged);
    // Never report a turn worse than a sample the grid already holds: the middle sample is a
    // real evaluation, and a section search that lands beside it is not an improvement on it.
    const double sampled = series.values[i];
    const bool sampleWins = maximum ? sampled > value : sampled < value;
    event.timeSeconds = sampleWins ? t1 : best;
    event.value = sampleWins ? sampled : value;
    event.locatedToSeconds = located;
    event.refined = true;
    event.converged = converged;
  } else {
    double value = series.values[i];
    event.timeSeconds = interpolateExtremum(t0, t1, t2, series.values[i - 1], series.values[i],
                                            series.values[i + 1], value);
    event.value = value;
    event.locatedToSeconds = t2 - t0;
    event.refined = false;
    event.converged = true;
  }
  return event;
}

// Windows on one side of a level, sharing the crossing refinement with crossings() so a window
// edge and a crossing can never be reported at different instants.
std::vector<LevelWindow> windowsOn(const EventSeries& series, double level,
                                   const EventTolerance& tolerance, bool above) {
  requireUsableSeries(series);
  const int nT = static_cast<int>(series.times.size());

  // "Inside" is at-or-above for windowsAbove and strictly below for windowsBelow, so the two
  // partition the grid: every sample is inside exactly one of them and a window edge belongs
  // to the window that contains the level.
  const auto inside = [&](int k) {
    const double value = series.values[static_cast<std::size_t>(k)];
    return above ? !(value < level) : value < level;
  };

  std::vector<LevelWindow> windows;
  LevelWindow current;
  bool open = false;

  if (inside(0)) {
    current = LevelWindow{};
    current.startSeconds = series.times.front();
    current.entryObserved = false;
    open = true;
  }

  for (int k = 0; k + 1 < nT; ++k) {
    if (inside(k) == inside(k + 1)) {
      continue;
    }
    const bool entering = inside(k + 1);
    // Rising or falling is a property of the CURVE, not of which side the window is on: going
    // below a level enters a windowsBelow window and is still a falling crossing. Taken from
    // the samples exactly as crossings() takes it, so the two agree on a shared edge.
    const bool rising = series.values[static_cast<std::size_t>(k)] < level;
    const TrajectoryEvent edge = locateCrossing(series, k, level, rising, tolerance);
    if (entering) {
      current = LevelWindow{};
      current.startSeconds = edge.timeSeconds;
      current.entryObserved = true;
      open = true;
    } else if (open) {
      current.endSeconds = edge.timeSeconds;
      current.exitObserved = true;
      windows.push_back(current);
      open = false;
    }
  }

  if (open) {
    current.endSeconds = series.times.back();
    current.exitObserved = false;
    windows.push_back(current);
  }
  return windows;
}

}  // namespace

const char* eventKindName(EventKind kind) {
  switch (kind) {
    case EventKind::Rising:
      return "rising";
    case EventKind::Falling:
      return "falling";
    case EventKind::Maximum:
      return "maximum";
    case EventKind::Minimum:
      return "minimum";
  }
  return "?";
}

std::vector<TrajectoryEvent> crossings(const EventSeries& series, double level,
                                       const EventTolerance& tolerance) {
  requireUsableSeries(series);
  const int nT = static_cast<int>(series.times.size());

  std::vector<TrajectoryEvent> events;
  for (int k = 0; k + 1 < nT; ++k) {
    // "Below the level" as a boolean, so a sample sitting exactly on the level is not-below
    // and the flag flips exactly once per crossing. Comparing signs of a difference instead
    // would report a touch twice, or not at all.
    const bool below0 = series.values[static_cast<std::size_t>(k)] < level;
    const bool below1 = series.values[static_cast<std::size_t>(k) + 1] < level;
    if (below0 == below1) {
      continue;
    }
    events.push_back(locateCrossing(series, k, level, /*rising=*/below0, tolerance));
  }
  return events;
}

std::optional<TrajectoryEvent> firstCrossing(const EventSeries& series, double level,
                                             const EventTolerance& tolerance) {
  const std::vector<TrajectoryEvent> events = crossings(series, level, tolerance);
  if (events.empty()) {
    return std::nullopt;
  }
  return events.front();
}

std::optional<TrajectoryEvent> lastCrossing(const EventSeries& series, double level,
                                            const EventTolerance& tolerance) {
  const std::vector<TrajectoryEvent> events = crossings(series, level, tolerance);
  if (events.empty()) {
    return std::nullopt;
  }
  return events.back();
}

std::vector<TrajectoryEvent> extrema(const EventSeries& series, const EventTolerance& tolerance) {
  requireUsableSeries(series);
  const int nT = static_cast<int>(series.times.size());
  std::vector<TrajectoryEvent> events;
  if (nT < 3) {
    return events;
  }

  for (int k = 1; k + 1 < nT; ++k) {
    const double before = series.values[static_cast<std::size_t>(k) - 1];
    const double here = series.values[static_cast<std::size_t>(k)];
    const double after = series.values[static_cast<std::size_t>(k) + 1];

    // A turn has to be strict on at least one side. Requiring it on both would miss a peak
    // the grid happened to sample twice at the same height; requiring it on neither would
    // report every flat stretch as an extremum.
    const bool maximum = here >= before && here >= after && (here > before || here > after);
    const bool minimum = here <= before && here <= after && (here < before || here < after);
    if (!maximum && !minimum) {
      continue;
    }
    events.push_back(locateExtremum(series, k, maximum, tolerance));
  }
  return events;
}

std::vector<LevelWindow> windowsAbove(const EventSeries& series, double level,
                                      const EventTolerance& tolerance) {
  return windowsOn(series, level, tolerance, /*above=*/true);
}

std::vector<LevelWindow> windowsBelow(const EventSeries& series, double level,
                                      const EventTolerance& tolerance) {
  return windowsOn(series, level, tolerance, /*above=*/false);
}

// --- series drawn from a response table --------------------------------------

EventSeries totalSeries(const ResponseTable& table) {
  EventSeries series;
  series.times = table.times;
  series.values = table.totals;
  return series;
}

EventSeries contributorSeries(const ResponseTable& table, int contributor) {
  if (contributor < 0 || contributor >= table.contributorCount()) {
    throw InputError(tagged(kModule, "contributor index " + std::to_string(contributor) +
                                         " is outside this table's " +
                                         std::to_string(table.contributorCount()) + " columns"));
  }
  EventSeries series;
  series.times = table.times;
  series.values.reserve(table.times.size());
  for (int k = 0; k < table.timeCount(); ++k) {
    series.values.push_back(table.valuesAt(k)[static_cast<std::size_t>(contributor)]);
  }
  return series;
}

EventSeries ratioSeries(const ResponseTable& table, int numerator, int denominator) {
  const EventSeries top = contributorSeries(table, numerator);
  const EventSeries bottom = contributorSeries(table, denominator);

  // Skip the leading samples where the denominator does not yet exist. A ratio clock does not
  // start until the nuclide it divides by has grown in, and reporting infinity for the first
  // sample of an ingrowth would put a spurious crossing in front of every real one.
  std::size_t start = 0;
  while (start < bottom.values.size() && !(bottom.values[start] > 0.0)) {
    ++start;
  }

  EventSeries series;
  for (std::size_t k = start; k < bottom.values.size(); ++k) {
    if (!(bottom.values[k] > 0.0)) {
      // Positive, then not, then positive again: the ratio has a pole inside the grid. No
      // bracket spanning it means anything, so this refuses rather than searching across it.
      throw InputError(tagged(
          kModule, "the denominator returns to zero at t=" + std::to_string(bottom.times[k]) +
                       " s, so the ratio is discontinuous inside this grid. Restrict "
                       "the time range to where the denominator is present"));
    }
    series.times.push_back(top.times[k]);
    series.values.push_back(top.values[k] / bottom.values[k]);
  }

  if (series.times.size() < 2) {
    throw InputError(tagged(kModule,
                            "the denominator is present for fewer than two samples, which is "
                            "not enough grid for a ratio to be located in"));
  }
  return series;
}

// --- refining an event by re-solving -----------------------------------------

namespace {

// An evaluator narrows an event by producing values the grid never sampled, so it has to be
// producing values of the SAME curve. Metric, aggregate and unit are the identity of what is
// plotted; the geometry is part of what an exposure or a fluence MEANS, and two distances are
// two different curves however alike they look.
void requireSameCurve(const ResponseTable& table, const ResponseEvaluator& evaluator) {
  const ResponseSpec& spec = evaluator.spec();
  const exposure::PointSourceGeometry& sampled = table.geometry;
  const exposure::PointSourceGeometry& asked = spec.geometry;
  const bool sameGeometry =
      sampled.distanceM == asked.distanceM && sampled.airDensityKgM3 == asked.airDensityKgM3 &&
      sampled.airAttenuation == asked.airAttenuation && sampled.buildup == asked.buildup &&
      sampled.irradiation == asked.irradiation;
  if (table.metric != spec.metric || table.aggregate != spec.aggregate || table.unit != spec.unit ||
      !sameGeometry) {
    throw InputError(tagged(kModule,
                            "this evaluator was built for a different response than the table it "
                            "would refine: it would narrow an event on one curve using values "
                            "taken from another"));
  }
}

}  // namespace

ResponseEvaluator::ResponseEvaluator(const NuclearData& data, const Inventory& inventory,
                                     const ResponseSpec& spec, const DecayOptions& options)
    : data_(data), inventory_(inventory), spec_(spec), options_(options) {}

const ResponseTable& ResponseEvaluator::tableAt(double timeSeconds) const {
  if (hasCache_ && cachedTime_ == timeSeconds) {
    return cached_;
  }
  const double times[] = {timeSeconds};
  const DecayResult result = decay(data_, inventory_, std::span<const double>(times, 1), options_);
  cached_ = buildResponse(data_, result, spec_);
  cachedTime_ = timeSeconds;
  hasCache_ = true;
  ++solves_;
  return cached_;
}

double ResponseEvaluator::total(double timeSeconds) const {
  return tableAt(timeSeconds).totals.front();
}

double ResponseEvaluator::value(const ContributorId& id, double timeSeconds) const {
  const ResponseTable& table = tableAt(timeSeconds);
  for (int c = 0; c < table.contributorCount(); ++c) {
    const ContributorId& other = table.contributors[static_cast<std::size_t>(c)];
    if (other.key != id.key) {
      continue;
    }
    // A line table carries several columns per emitter, so the key alone does not name one.
    if (table.aggregate == Aggregate::GammaLine && other.lineEnergyEv != id.lineEnergyEv) {
      continue;
    }
    return table.valuesAt(0)[static_cast<std::size_t>(c)];
  }
  // The columns of a response are the buckets of the seed's forward closure, which does not
  // depend on time, and a line's share of its own emitter is a fixed fraction -- so the same
  // contributors exist at every instant. Arriving here means the evaluator and the table were
  // not built from the same seed, which is a programming error rather than bad input.
  throw NusiftError(tagged(kModule,
                           "the response at this time carries no column for a contributor the "
                           "table holds, so the two were not built from the same inventory"));
}

EventSeries totalSeries(const ResponseTable& table, const ResponseEvaluator& evaluator) {
  requireSameCurve(table, evaluator);
  EventSeries series = totalSeries(table);
  series.evaluate = [&evaluator](double t) { return evaluator.total(t); };
  return series;
}

EventSeries contributorSeries(const ResponseTable& table, int contributor,
                              const ResponseEvaluator& evaluator) {
  requireSameCurve(table, evaluator);
  // The sampling half validates the index, so the identity below is only read once it is known
  // to exist.
  EventSeries series = contributorSeries(table, contributor);
  const ContributorId id = table.contributors[static_cast<std::size_t>(contributor)];
  series.evaluate = [&evaluator, id](double t) { return evaluator.value(id, t); };
  return series;
}

EventSeries ratioSeries(const ResponseTable& table, int numerator, int denominator,
                        const ResponseEvaluator& evaluator) {
  requireSameCurve(table, evaluator);
  EventSeries series = ratioSeries(table, numerator, denominator);
  const ContributorId top = table.contributors[static_cast<std::size_t>(numerator)];
  const ContributorId bottom = table.contributors[static_cast<std::size_t>(denominator)];
  series.evaluate = [&evaluator, top, bottom](double t) {
    const double below = evaluator.value(bottom, t);
    // The sampled series already refused a denominator that vanishes at any sample, so a zero
    // here sits between two samples that both saw it present. Returning an infinity would drag
    // the root to whichever end of the bracket the search last tried; saying what happened is
    // the only honest answer.
    if (!(below > 0.0)) {
      throw InputError(
          tagged(kModule, "the ratio's denominator vanishes at t=" + std::to_string(t) +
                              " s, between two samples that both found it present"));
    }
    return evaluator.value(top, t) / below;
  };
  return series;
}

EventSeries taskSeries(const NuclearData& data, const Inventory& inventory,
                       const ResponseSpec& spec, std::span<const double> startTimes,
                       double durationSeconds, const DecayOptions& options, bool refine) {
  if (!(durationSeconds > 0.0)) {
    throw InputError(tagged(kModule,
                            "a task needs a positive duration: a window of no length accrues "
                            "nothing, whenever it starts"));
  }
  if (startTimes.size() < 2) {
    throw InputError(tagged(kModule,
                            "a task curve needs at least two start times for an event to sit "
                            "between; one start time is a single answer, not a curve"));
  }
  for (std::size_t k = 1; k < startTimes.size(); ++k) {
    if (!(startTimes[k] > startTimes[k - 1])) {
      throw InputError(tagged(kModule, "task start times must be strictly increasing"));
    }
  }

  // Copies rather than captured references: this callable outlives the call when `refine` puts
  // it in the series, and a spec or a set of solver options that died with the argument list
  // would be read during the search. The data and the inventory are the caller's to keep alive,
  // which the header says.
  const ResponseSpec heldSpec = spec;
  const DecayOptions heldOptions = options;
  const auto accrued = [&data, &inventory, heldSpec, heldOptions, durationSeconds](double start) {
    std::vector<std::int64_t> keys;
    const std::vector<double> integral =
        intervalIntegral(data, inventory, start, start + durationSeconds, &keys, heldOptions);
    return buildIntervalResponse(data, keys, integral, start, start + durationSeconds, heldSpec)
        .totals.front();
  };

  EventSeries series;
  series.times.assign(startTimes.begin(), startTimes.end());
  series.values.reserve(startTimes.size());
  for (const double start : startTimes) {
    series.values.push_back(accrued(start));
  }
  if (refine) {
    series.evaluate = accrued;
  }
  return series;
}

StayTime stayTime(const NuclearData& data, const Inventory& inventory, const ResponseSpec& spec,
                  double startSeconds, double budget, double maxDurationSeconds,
                  const DecayOptions& options, const EventTolerance& tolerance) {
  if (!(startSeconds >= 0.0)) {
    throw InputError(tagged(kModule, "a stay cannot begin before the inventory exists"));
  }
  if (!(budget > 0.0)) {
    throw InputError(tagged(kModule,
                            "a stay needs a positive budget: a budget of zero is spent before "
                            "the door opens, which is not an answer about a duration"));
  }
  if (!(maxDurationSeconds > 0.0)) {
    throw InputError(tagged(kModule,
                            "a stay needs a positive ceiling to search inside; how long a stay "
                            "would be considered at all is part of the question"));
  }

  const auto accrued = [&data, &inventory, &spec, &options, startSeconds](double duration) {
    const double end = startSeconds + duration;
    std::vector<std::int64_t> keys;
    const std::vector<double> integral =
        intervalIntegral(data, inventory, startSeconds, end, &keys, options);
    return buildIntervalResponse(data, keys, integral, startSeconds, end, spec).totals.front();
  };

  StayTime stay;
  stay.startSeconds = startSeconds;
  stay.budget = budget;
  stay.maxDurationSeconds = maxDurationSeconds;
  stay.accruedAtMax = accrued(maxDurationSeconds);
  stay.samples = 1;

  // Not reached inside the window asked about. Reported as such rather than as a very large
  // duration, for the reason allowableScale() reports an unbounded criterion rather than a huge
  // scale: "the budget is not spent in a day" and "you may stay 8.6e17 seconds" are different
  // statements, and only one of them is something this tool observed.
  if (stay.accruedAtMax < budget) {
    return stay;
  }
  // Spent exactly at the ceiling. Located to nothing because no search happened, which is the
  // truthful width rather than a tolerance nobody used.
  if (stay.accruedAtMax == budget) {
    stay.bounded = true;
    stay.durationSeconds = maxDurationSeconds;
    return stay;
  }

  // A(0) = 0 exactly and needs no solve, so the bracket is free at one end. Illinois takes it
  // from there -- the same refinement every located event in this file is narrowed by, on a
  // curve that happens to guarantee the root inside it is the only one.
  int evaluations = 0;
  stay.durationSeconds =
      refineCrossing(accrued, 0.0, maxDurationSeconds, 0.0, stay.accruedAtMax, budget, tolerance,
                     stay.locatedToSeconds, stay.converged, evaluations);
  stay.samples += evaluations;
  stay.bounded = true;
  return stay;
}

}  // namespace nusift
