#include "nusift/triage/forecast.hpp"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <limits>
#include <map>
#include <optional>
#include <span>
#include <vector>

namespace nusift {
namespace {

// The contributor holding the largest value at `timeIndex`, or -1 if nothing does. Ties break
// on the contributor key for the same reason ranking does: without it the leader can appear to
// change at a crossover purely from sort order.
int leaderAt(const ResponseTable& table, int timeIndex) {
  const std::span<const double> values = table.valuesAt(timeIndex);
  int best = -1;
  for (int c = 0; c < table.contributorCount(); ++c) {
    const double value = values[static_cast<std::size_t>(c)];
    if (value <= 0.0) {
      continue;
    }
    if (best < 0 || value > values[static_cast<std::size_t>(best)]) {
      best = c;
      continue;
    }
    // Key then energy, so two lines of one emitter tying exactly resolve the way a ranking does.
    const ContributorId& mine = table.contributors[static_cast<std::size_t>(c)];
    const ContributorId& held = table.contributors[static_cast<std::size_t>(best)];
    if (value == values[static_cast<std::size_t>(best)] &&
        (mine.key < held.key || (mine.key == held.key && mine.lineEnergyEv < held.lineEnergyEv))) {
      best = c;
    }
  }
  return best;
}

// Where two contributors' values cross between consecutive samples.
//
// A leader change is a crossing of the ratio a/b through one, so it is handed to the event
// engine rather than located here: one code path places a crossing, and a boundary in a
// forecast and a crossing from `when` can never disagree about the same instant. The engine's
// log-linear interpolation of a positive level reduces, for a ratio against one, to exactly
// the closed form this function used to carry -- both contenders are close to exponential over
// one interval -- so the merge does not move an unrefined boundary.
//
// With an evaluator the same bracket is narrowed by real solves instead. Empty when the ratio
// does not cross one inside [t_k, t_k+1], or cannot be formed there.
std::optional<double> crossingIn(const ResponseTable& table, int leaving, int arriving, int k,
                                 const ResponseEvaluator* evaluator,
                                 const EventTolerance& tolerance) {
  const double t0 = table.times[static_cast<std::size_t>(k)];
  const double t1 = table.times[static_cast<std::size_t>(k) + 1];

  const std::span<const double> before = table.valuesAt(k);
  const std::span<const double> after = table.valuesAt(k + 1);
  const double a0 = before[static_cast<std::size_t>(leaving)];
  const double b0 = before[static_cast<std::size_t>(arriving)];
  const double a1 = after[static_cast<std::size_t>(leaving)];
  const double b1 = after[static_cast<std::size_t>(arriving)];
  if (!(a0 > 0.0 && b0 > 0.0 && a1 > 0.0 && b1 > 0.0)) {
    return std::nullopt;
  }

  EventSeries pair;
  pair.times = {t0, t1};
  pair.values = {a0 / b0, a1 / b1};
  if (evaluator != nullptr) {
    const ContributorId top = table.contributors[static_cast<std::size_t>(leaving)];
    const ContributorId bottom = table.contributors[static_cast<std::size_t>(arriving)];
    pair.evaluate = [evaluator, top, bottom](double t) {
      const double below = evaluator->value(bottom, t);
      // Both contenders are positive at both ends of this bracket, so a zero in between is a
      // curve that dipped out and back. The ratio is genuinely unbounded there and the leader
      // is genuinely the numerator, which is what an infinity says; the search reads it as a
      // residual with a sign and bisects, rather than stepping to a nonsense time.
      if (!(below > 0.0)) {
        return std::numeric_limits<double>::infinity();
      }
      return evaluator->value(top, t) / below;
    };
  }

  const std::vector<TrajectoryEvent> found = crossings(pair, 1.0, tolerance);
  if (found.empty()) {
    return std::nullopt;
  }
  return found.front().timeSeconds;
}

// The boundary between two windows, searched over every interval from the last sample the
// outgoing contributor actually led to the first the incoming one leads. The outgoing one leads
// at one end and the incoming one at the other, so their ratio crosses one somewhere inside.
//
// Searching only the interval after the outgoing WINDOW ended -- which is what this did -- is
// wrong whenever a short run in between was absorbed: for a leader sequence A, B, C with B
// lasting one sample, the A/C crossing lies before B's sample and the interval after it holds
// none, so the boundary fell back to a midpoint in the wrong interval, refined or not.
double boundaryTime(const ResponseTable& table, int leaving, int arriving, int lastLed,
                    int firstLed, const ResponseEvaluator* evaluator,
                    const EventTolerance& tolerance) {
  for (int k = lastLed; k < firstLed; ++k) {
    if (const std::optional<double> at =
            crossingIn(table, leaving, arriving, k, evaluator, tolerance)) {
      return *at;
    }
  }
  // Reached only when one of the two is zero across the span, so no ratio exists to cross.
  return 0.5 * (table.times[static_cast<std::size_t>(lastLed)] +
                table.times[static_cast<std::size_t>(firstLed)]);
}

RankTrack trackOf(const ResponseTable& table, int contributor) {
  RankTrack track;
  track.id = table.contributors[static_cast<std::size_t>(contributor)];
  track.label = table.labels[static_cast<std::size_t>(contributor)];
  track.rank.assign(static_cast<std::size_t>(table.timeCount()), 0);
  track.fraction.assign(static_cast<std::size_t>(table.timeCount()), 0.0);

  for (int k = 0; k < table.timeCount(); ++k) {
    const std::span<const double> values = table.valuesAt(k);
    const double mine = values[static_cast<std::size_t>(contributor)];
    const double total = table.totals[static_cast<std::size_t>(k)];
    const double fraction = total > 0.0 ? mine / total : 0.0;
    track.fraction[static_cast<std::size_t>(k)] = fraction;

    int rank = 1;
    for (int c = 0; c < table.contributorCount(); ++c) {
      if (values[static_cast<std::size_t>(c)] > mine) {
        ++rank;
      }
    }
    track.rank[static_cast<std::size_t>(k)] = mine > 0.0 ? rank : 0;

    if (fraction > track.peakFraction) {
      track.peakFraction = fraction;
      track.peakTimeIndex = k;
    }
  }
  return track;
}

// The body both entry points share. The evaluator is a pointer rather than two copies of two
// hundred lines: everything about which runs exist and how they are absorbed is identical, and
// only how tightly the boundaries between them are placed differs.
std::vector<DominanceWindow> dominanceWindowsImpl(const ResponseTable& table, int minSamples,
                                                  const ResponseEvaluator* evaluator,
                                                  const EventTolerance& tolerance) {
  const int nT = table.timeCount();
  if (nT == 0 || table.contributorCount() == 0) {
    return {};
  }

  std::vector<int> leaders(static_cast<std::size_t>(nT));
  for (int k = 0; k < nT; ++k) {
    leaders[static_cast<std::size_t>(k)] = leaderAt(table, k);
  }

  // Coalesce consecutive samples with the same leader into runs.
  // `lastLed` is the last sample this contributor actually led. It equals lastIndex until a
  // short run after it is absorbed, which extends the window but not the samples it leads.
  struct Run {
    int contributor;
    int firstIndex;
    int lastIndex;
    int lastLed;
  };
  std::vector<Run> runs;
  for (int k = 0; k < nT; ++k) {
    const int leader = leaders[static_cast<std::size_t>(k)];
    if (leader < 0) {
      continue;
    }
    if (!runs.empty() && runs.back().contributor == leader) {
      runs.back().lastIndex = k;
      runs.back().lastLed = k;
    } else {
      runs.push_back(Run{leader, k, k, k});
    }
  }
  if (runs.empty()) {
    return {};
  }

  // Absorb runs the grid barely resolved. Near a crossover two contenders can trade places
  // sample to sample within numerical noise, and reporting six one-sample windows is less
  // truthful than reporting one boundary.
  if (minSamples > 1 && runs.size() > 1) {
    std::vector<Run> kept;
    for (const Run& run : runs) {
      const int samples = run.lastIndex - run.firstIndex + 1;
      if (samples < minSamples && !kept.empty()) {
        kept.back().lastIndex = run.lastIndex;
        continue;
      }
      kept.push_back(run);
    }
    // Merge neighbours that ended up the same contributor after absorbing.
    runs.clear();
    for (const Run& run : kept) {
      if (!runs.empty() && runs.back().contributor == run.contributor) {
        runs.back().lastIndex = run.lastIndex;
        runs.back().lastLed = run.lastLed;
      } else {
        runs.push_back(run);
      }
    }
  }

  // Each interior boundary once, so a window's end and the next one's start are the same number.
  std::vector<double> boundaries;
  for (std::size_t r = 0; r + 1 < runs.size(); ++r) {
    boundaries.push_back(boundaryTime(table, runs[r].contributor, runs[r + 1].contributor,
                                      runs[r].lastLed, runs[r + 1].firstIndex, evaluator,
                                      tolerance));
  }

  std::vector<DominanceWindow> windows;
  windows.reserve(runs.size());
  for (std::size_t r = 0; r < runs.size(); ++r) {
    const Run& run = runs[r];
    DominanceWindow window;
    window.id = table.contributors[static_cast<std::size_t>(run.contributor)];
    window.label = table.labels[static_cast<std::size_t>(run.contributor)];

    // The first window starts at the grid's start and the last ends at its end; interior
    // boundaries are the interpolated crossings.
    window.startSeconds = r == 0 ? table.times.front() : boundaries[r - 1];
    window.endSeconds = r + 1 == runs.size() ? table.times.back() : boundaries[r];

    for (int k = run.firstIndex; k <= run.lastIndex; ++k) {
      const double total = table.totals[static_cast<std::size_t>(k)];
      if (total > 0.0) {
        window.peakFraction =
            std::max(window.peakFraction,
                     table.valuesAt(k)[static_cast<std::size_t>(run.contributor)] / total);
      }
    }
    windows.push_back(std::move(window));
  }
  return windows;
}

}  // namespace

std::vector<DominanceWindow> dominanceWindows(const ResponseTable& table, int minSamples) {
  return dominanceWindowsImpl(table, minSamples, /*evaluator=*/nullptr, EventTolerance{});
}

std::vector<DominanceWindow> dominanceWindows(const ResponseTable& table,
                                              const ResponseEvaluator& evaluator, int minSamples,
                                              const EventTolerance& tolerance) {
  return dominanceWindowsImpl(table, minSamples, &evaluator, tolerance);
}

std::vector<RankTrack> unionTopN(const ResponseTable& table, int n) {
  if (n <= 0 || table.timeCount() == 0) {
    return {};
  }

  std::vector<char> everInTop(static_cast<std::size_t>(table.contributorCount()), 0);
  for (int k = 0; k < table.timeCount(); ++k) {
    const std::span<const double> values = table.valuesAt(k);
    std::vector<int> order(static_cast<std::size_t>(table.contributorCount()));
    for (int c = 0; c < table.contributorCount(); ++c) {
      order[static_cast<std::size_t>(c)] = c;
    }
    const int depth = std::min(n, table.contributorCount());
    std::partial_sort(order.begin(), order.begin() + depth, order.end(), [&](int a, int b) {
      const double va = values[static_cast<std::size_t>(a)];
      const double vb = values[static_cast<std::size_t>(b)];
      if (va != vb) {
        return va > vb;
      }
      // Key then energy, as the ranking breaks ties, so the same line enters the top n here.
      const ContributorId& ia = table.contributors[static_cast<std::size_t>(a)];
      const ContributorId& ib = table.contributors[static_cast<std::size_t>(b)];
      if (ia.key != ib.key) {
        return ia.key < ib.key;
      }
      return ia.lineEnergyEv < ib.lineEnergyEv;
    });
    for (int i = 0; i < depth; ++i) {
      if (values[static_cast<std::size_t>(order[static_cast<std::size_t>(i)])] > 0.0) {
        everInTop[static_cast<std::size_t>(order[static_cast<std::size_t>(i)])] = 1;
      }
    }
  }

  std::vector<RankTrack> tracks;
  for (int c = 0; c < table.contributorCount(); ++c) {
    if (everInTop[static_cast<std::size_t>(c)] != 0) {
      tracks.push_back(trackOf(table, c));
    }
  }
  // By peak share, so the row that matters most somewhere comes first -- ordering by value at
  // any single time would bury exactly the contributor a forecast exists to surface.
  std::sort(tracks.begin(), tracks.end(), [](const RankTrack& a, const RankTrack& b) {
    if (a.peakFraction != b.peakFraction) {
      return a.peakFraction > b.peakFraction;
    }
    if (a.id.key != b.id.key) {
      return a.id.key < b.id.key;
    }
    return a.id.lineEnergyEv < b.id.lineEnergyEv;
  });
  return tracks;
}

std::vector<RankTrack> unionTopN(const ResponseTable& table, int n,
                                 std::span<const std::int64_t> pinned) {
  std::vector<RankTrack> tracks = unionTopN(table, n);
  if (pinned.empty() || table.timeCount() == 0) {
    return tracks;
  }

  std::vector<RankTrack> extras;
  for (int c = 0; c < table.contributorCount(); ++c) {
    const ContributorId& id = table.contributors[static_cast<std::size_t>(c)];
    if (std::find(pinned.begin(), pinned.end(), id.key) == pinned.end()) {
      continue;
    }
    // A gamma-line table has one column per line, so identity there is the emitter AND the
    // energy; comparing keys alone would drop every line of an emitter after its first.
    const bool already = std::any_of(tracks.begin(), tracks.end(), [&](const RankTrack& track) {
      return track.id.key == id.key && track.id.lineEnergyEv == id.lineEnergyEv;
    });
    if (already) {
      continue;
    }
    RankTrack track = trackOf(table, c);
    track.pinned = true;
    extras.push_back(std::move(track));
  }

  // Ordered among themselves the way the tracks above them are, and appended rather than merged:
  // a pinned track is in the list because it was asked for, and sorting it in with the ones the
  // forecast found would erase that distinction at exactly the moment it matters.
  std::sort(extras.begin(), extras.end(), [](const RankTrack& a, const RankTrack& b) {
    if (a.peakFraction != b.peakFraction) {
      return a.peakFraction > b.peakFraction;
    }
    if (a.id.key != b.id.key) {
      return a.id.key < b.id.key;
    }
    return a.id.lineEnergyEv < b.id.lineEnergyEv;
  });
  tracks.insert(tracks.end(), std::make_move_iterator(extras.begin()),
                std::make_move_iterator(extras.end()));
  return tracks;
}

std::vector<RankTrack> persistentTopN(const ResponseTable& table, int n) {
  std::vector<RankTrack> tracks = unionTopN(table, n);
  tracks.erase(std::remove_if(tracks.begin(), tracks.end(),
                              [n](const RankTrack& track) {
                                return std::any_of(track.rank.begin(), track.rank.end(),
                                                   [n](int rank) { return rank == 0 || rank > n; });
                              }),
               tracks.end());
  return tracks;
}

}  // namespace nusift
