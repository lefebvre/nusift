#include "nusift/triage/ranking.hpp"

#include <algorithm>
#include <string>

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

  double cumulative = 0.0;
  int kept = 0;
  for (const int c : order) {
    const double value = values[static_cast<std::size_t>(c)];
    if (value <= 0.0) {
      break;  // sorted descending, so nothing after this contributes
    }
    const double fraction = haveTotal ? value / ranking.total : 0.0;
    if (fraction < request.minFraction) {
      break;
    }
    if (request.topN > 0 && kept >= request.topN) {
      break;
    }

    cumulative += fraction;
    Contributor contributor;
    contributor.id = table.contributors[static_cast<std::size_t>(c)];
    contributor.label = table.labels[static_cast<std::size_t>(c)];
    contributor.value = value;
    contributor.fraction = fraction;
    contributor.cumulativeFraction = cumulative;
    contributor.rank = ++kept;
    contributor.flags = table.flags[static_cast<std::size_t>(c)];
    ranking.contributors.push_back(std::move(contributor));

    // Coverage is checked AFTER appending, so the returned prefix is the smallest one that
    // reaches the requested fraction rather than the largest one that stays below it.
    if (request.coverage > 0.0 && cumulative >= request.coverage) {
      break;
    }
  }

  ranking.coveredFraction = cumulative;

  // Count only contributors that were left out AND actually contribute. A chain always
  // carries stable terminators with exactly zero activity; reporting "3 further contributors
  // omitted" alongside "shown rows cover 100%" is a contradiction, and it invites the reader
  // to go looking for something that is not there.
  int omitted = 0;
  for (int c = 0; c < nC; ++c) {
    if (values[static_cast<std::size_t>(c)] > 0.0) {
      ++omitted;
    }
  }
  ranking.omittedCount = omitted - static_cast<int>(ranking.contributors.size());
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
