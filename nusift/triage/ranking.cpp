#include "nusift/triage/ranking.hpp"

#include <algorithm>
#include <cstdint>
#include <string>
#include <vector>

#include "nusift/core/error.hpp"

namespace nusift {
namespace {

constexpr const char* kModule = "ranking";

}  // namespace

Ranking rank(const ResponseTable& table, int timeIndex, const RankRequest& request) {
  if (timeIndex < 0 || timeIndex >= table.timeCount()) {
    throw NusiftError(tagged(kModule, "time index " + std::to_string(timeIndex) +
                                          " out of range [0, " + std::to_string(table.timeCount()) +
                                          ")"));
  }

  Ranking ranking;
  ranking.metric = table.metric;
  ranking.aggregate = table.aggregate;
  ranking.domain = table.domain;
  ranking.unit = table.unit;
  ranking.time = table.times[static_cast<std::size_t>(timeIndex)];
  if (!table.timeEnds.empty()) {
    ranking.timeEnd = table.timeEnds[static_cast<std::size_t>(timeIndex)];
  }
  ranking.total = table.totals[static_cast<std::size_t>(timeIndex)];
  if (!table.unmodeledEnergyFraction.empty()) {
    ranking.unmodeledEnergyFraction =
        table.unmodeledEnergyFraction[static_cast<std::size_t>(timeIndex)];
  }

  const std::span<const double> values = table.valuesAt(timeIndex);
  const int nC = table.contributorCount();

  std::vector<int> order(static_cast<std::size_t>(nC));
  for (int c = 0; c < nC; ++c) {
    order[static_cast<std::size_t>(c)] = c;
  }
  std::sort(order.begin(), order.end(), [&](int a, int b) {
    const double va = values[static_cast<std::size_t>(a)];
    const double vb = values[static_cast<std::size_t>(b)];
    if (va != vb) {
      return va > vb;
    }
    // Deterministic tiebreak; see the header.
    return table.contributors[static_cast<std::size_t>(a)].key <
           table.contributors[static_cast<std::size_t>(b)].key;
  });

  // A zero or negative total makes every fraction meaningless. That happens legitimately --
  // an inventory of nothing but stable nuclides has no activity at all -- so it is reported
  // as an empty ranking rather than treated as an error or divided by.
  const bool haveTotal = ranking.total > 0.0;

  // Walk the whole ordering once, before deciding what to return: where each contributor
  // stands, what it is worth, and what everything down to it covers.
  //
  // The selection below only ever needs the prefix of this, but a pinned row needs all three
  // for a place it was cut from -- and reporting a pinned contributor's rank and the coverage
  // above it is what makes it an answer ("Cs-137 is 37th, and the rows over it are 99.8% of
  // the total") rather than a number floating outside the ranking that excluded it.
  std::vector<int> trueRank(static_cast<std::size_t>(nC), 0);
  std::vector<double> fractionOf(static_cast<std::size_t>(nC), 0.0);
  std::vector<double> cumulativeThrough(static_cast<std::size_t>(nC), 0.0);
  int contributing = 0;
  double running = 0.0;
  for (const int c : order) {
    const double value = values[static_cast<std::size_t>(c)];
    if (value <= 0.0) {
      break;  // sorted descending, so nothing after this contributes
    }
    const double fraction = haveTotal ? value / ranking.total : 0.0;
    running += fraction;
    fractionOf[static_cast<std::size_t>(c)] = fraction;
    trueRank[static_cast<std::size_t>(c)] = ++contributing;
    cumulativeThrough[static_cast<std::size_t>(c)] = running;
  }

  std::vector<char> shown(static_cast<std::size_t>(nC), 0);
  const auto append = [&](int c, bool pinned) {
    Contributor contributor;
    contributor.id = table.contributors[static_cast<std::size_t>(c)];
    contributor.label = table.labels[static_cast<std::size_t>(c)];
    contributor.value = values[static_cast<std::size_t>(c)];
    contributor.fraction = fractionOf[static_cast<std::size_t>(c)];
    contributor.cumulativeFraction = cumulativeThrough[static_cast<std::size_t>(c)];
    contributor.rank = trueRank[static_cast<std::size_t>(c)];
    contributor.flags = table.flags[static_cast<std::size_t>(c)];
    contributor.pinned = pinned;
    ranking.contributors.push_back(std::move(contributor));
    shown[static_cast<std::size_t>(c)] = 1;
  };

  double covered = 0.0;
  int kept = 0;
  for (const int c : order) {
    if (trueRank[static_cast<std::size_t>(c)] == 0) {
      break;  // nothing from here on contributes
    }
    if (fractionOf[static_cast<std::size_t>(c)] < request.minFraction) {
      break;
    }
    if (request.topN > 0 && kept >= request.topN) {
      break;
    }

    append(c, /*pinned=*/false);
    ++kept;
    covered = cumulativeThrough[static_cast<std::size_t>(c)];

    // Coverage is checked AFTER appending, so the returned prefix is the smallest one that
    // reaches the requested fraction rather than the largest one that stays below it.
    if (request.coverage > 0.0 && covered >= request.coverage) {
      break;
    }
  }

  // Pinned contributors the prefix did not already reach, appended below it rather than mixed
  // into it. Ordered by where they stand, so a pinned tail reads the same way the ranking above
  // it does.
  if (!request.pinned.empty()) {
    std::vector<int> extras;
    for (int c = 0; c < nC; ++c) {
      const std::int64_t key = table.contributors[static_cast<std::size_t>(c)].key;
      if (shown[static_cast<std::size_t>(c)] == 0 &&
          std::find(request.pinned.begin(), request.pinned.end(), key) != request.pinned.end()) {
        extras.push_back(c);
      }
    }
    std::sort(extras.begin(), extras.end(), [&](int a, int b) {
      // Contributors holding no place at all go last: they are the answer "nothing, at this
      // time", and putting them among ranked rows would suggest otherwise. Ties are possible
      // only between the lines of one emitter in a gamma-line table, where energy separates
      // them -- and a tie there is exact, since the columns share a nuclide's atom count.
      const int ra = trueRank[static_cast<std::size_t>(a)];
      const int rb = trueRank[static_cast<std::size_t>(b)];
      if (ra != rb) {
        return (ra == 0 ? nC + 1 : ra) < (rb == 0 ? nC + 1 : rb);
      }
      return table.contributors[static_cast<std::size_t>(a)].lineEnergyEv <
             table.contributors[static_cast<std::size_t>(b)].lineEnergyEv;
    });
    for (const int c : extras) {
      append(c, /*pinned=*/true);
      covered += fractionOf[static_cast<std::size_t>(c)];
    }
  }

  ranking.coveredFraction = covered;

  // Count only contributors that were left out AND actually contribute. A chain always
  // carries stable terminators with exactly zero activity; reporting "3 further contributors
  // omitted" alongside "shown rows cover 100%" is a contradiction, and it invites the reader
  // to go looking for something that is not there. A pinned row that contributes nothing was
  // never in that count either, so it cannot come out of it.
  int omitted = contributing;
  for (const Contributor& contributor : ranking.contributors) {
    if (contributor.rank > 0) {
      --omitted;
    }
  }
  ranking.omittedCount = omitted;
  return ranking;
}

std::vector<Ranking> rankAll(const ResponseTable& table, const RankRequest& request) {
  std::vector<Ranking> rankings;
  rankings.reserve(static_cast<std::size_t>(table.timeCount()));
  for (int k = 0; k < table.timeCount(); ++k) {
    rankings.push_back(rank(table, k, request));
  }
  return rankings;
}

}  // namespace nusift
