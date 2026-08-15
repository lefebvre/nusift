#include "nusift/triage/forecast.hpp"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <map>
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
    if (best < 0 || value > values[static_cast<std::size_t>(best)] ||
        (value == values[static_cast<std::size_t>(best)] &&
         table.contributors[static_cast<std::size_t>(c)].key <
             table.contributors[static_cast<std::size_t>(best)].key)) {
      best = c;
    }
  }
  return best;
}

// Where two contributors' values cross between consecutive samples.
//
// Both are close to exponential over a grid interval, so log(a/b) is close to linear in time
// and its zero is the crossing. Falls back to the interval midpoint when the ratio does not
// actually change sign, which happens only if the caller asked about the wrong interval.
double crossingTime(const ResponseTable& table, int leaving, int arriving, int k) {
  const double t0 = table.times[static_cast<std::size_t>(k)];
  const double t1 = table.times[static_cast<std::size_t>(k) + 1];

  const std::span<const double> before = table.valuesAt(k);
  const std::span<const double> after = table.valuesAt(k + 1);
  const double a0 = before[static_cast<std::size_t>(leaving)];
  const double b0 = before[static_cast<std::size_t>(arriving)];
  const double a1 = after[static_cast<std::size_t>(leaving)];
  const double b1 = after[static_cast<std::size_t>(arriving)];
  if (!(a0 > 0.0 && b0 > 0.0 && a1 > 0.0 && b1 > 0.0)) {
    return 0.5 * (t0 + t1);
  }

  const double d0 = std::log(a0 / b0);
  const double d1 = std::log(a1 / b1);
  if (d0 == d1 || (d0 > 0.0) == (d1 > 0.0)) {
    return 0.5 * (t0 + t1);
  }
  return t0 + (t1 - t0) * d0 / (d0 - d1);
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

}  // namespace

std::vector<DominanceWindow> dominanceWindows(const ResponseTable& table, int minSamples) {
  const int nT = table.timeCount();
  if (nT == 0 || table.contributorCount() == 0) {
    return {};
  }

  std::vector<int> leaders(static_cast<std::size_t>(nT));
  for (int k = 0; k < nT; ++k) {
    leaders[static_cast<std::size_t>(k)] = leaderAt(table, k);
  }

  // Coalesce consecutive samples with the same leader into runs.
  struct Run {
    int contributor;
    int firstIndex;
    int lastIndex;
  };
  std::vector<Run> runs;
  for (int k = 0; k < nT; ++k) {
    const int leader = leaders[static_cast<std::size_t>(k)];
    if (leader < 0) {
      continue;
    }
    if (!runs.empty() && runs.back().contributor == leader) {
      runs.back().lastIndex = k;
    } else {
      runs.push_back(Run{leader, k, k});
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
      } else {
        runs.push_back(run);
      }
    }
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
    window.startSeconds = r == 0 ? table.times.front()
                                 : crossingTime(table, runs[r - 1].contributor, run.contributor,
                                                runs[r - 1].lastIndex);
    window.endSeconds = r + 1 == runs.size() ? table.times.back()
                                             : crossingTime(table, run.contributor,
                                                            runs[r + 1].contributor, run.lastIndex);

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
      return table.contributors[static_cast<std::size_t>(a)].key <
             table.contributors[static_cast<std::size_t>(b)].key;
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
    return a.id.key < b.id.key;
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
