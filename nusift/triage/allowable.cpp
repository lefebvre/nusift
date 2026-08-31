#include "nusift/triage/allowable.hpp"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "nusift/core/error.hpp"
#include "nusift/nucdata/nuclear_data.hpp"
#include "nusift/triage/ranking.hpp"

namespace nusift {
namespace {

constexpr const char* kModule = "allowable scale";

// Every way a criterion can fail to mean anything, checked before a single solve is weighted so
// the error names the criterion rather than surfacing from inside the response builder with no
// idea which of several it came from.
void requireUsableCriteria(std::span<const Criterion> criteria) {
  if (criteria.empty()) {
    throw InputError(tagged(kModule,
                            "no criteria given, so there is nothing for a scale to be limited "
                            "by -- an unconstrained inventory has no maximum"));
  }
  for (std::size_t q = 0; q < criteria.size(); ++q) {
    const Criterion& criterion = criteria[q];
    const std::string where = "criterion " + std::to_string(q);
    if (criterion.name.empty()) {
      throw InputError(tagged(kModule, where +
                                           " has no name. Naming the binding criterion is most "
                                           "of the answer, and an unnamed one cannot be named"));
    }
    if (!(criterion.limit > 0.0)) {
      throw InputError(tagged(kModule, "\"" + criterion.name + "\" has a limit of " +
                                           std::to_string(criterion.limit) +
                                           ". A limit of zero forbids the material outright "
                                           "rather than scaling it, and a negative one is not a "
                                           "limit"));
    }
    if (!unitSuitsDomain(criterion.spec.unit, Domain::Instant)) {
      throw InputError(
          tagged(kModule, "\"" + criterion.name + "\" is limited in " +
                              std::string(unitName(criterion.spec.unit)) +
                              ", which is a total accrued over a window. A possession or "
                              "transport limit constrains a quantity held at an instant, so "
                              "the criterion needs a rate unit"));
    }
  }
}

// The contributors deciding one criterion's total, most first. Taken from the criterion's own
// table so the aggregate the caller asked for is the aggregate they are named in: a criterion
// ranked by mass chain names chains, one ranked by gamma line names lines.
std::vector<LimitingContributor> limitingFor(const ResponseTable& table, int timeIndex,
                                             int limitingCount) {
  std::vector<LimitingContributor> limiting;
  if (limitingCount <= 0) {
    return limiting;
  }
  RankRequest request;
  request.topN = limitingCount;
  const Ranking ranking = rank(table, timeIndex, request);

  limiting.reserve(ranking.contributors.size());
  for (const Contributor& contributor : ranking.contributors) {
    limiting.push_back(
        LimitingContributor{contributor.id, contributor.label, contributor.fraction});
  }
  return limiting;
}

}  // namespace

std::vector<AllowableScale> allowableScale(const NuclearData& data, const DecayResult& result,
                                           std::span<const Criterion> criteria, int limitingCount) {
  requireUsableCriteria(criteria);

  // One table per criterion, built through the same path a report uses. Dotting weights against
  // atoms here instead would be cheaper by a constant and would let the number checked against
  // the limit drift from the number a ranking of the same spec would print -- which is the one
  // divergence this whole layer cannot afford, since the two are read side by side.
  std::vector<ResponseTable> tables;
  tables.reserve(criteria.size());
  for (const Criterion& criterion : criteria) {
    tables.push_back(buildResponse(data, result, criterion.spec));
  }

  const int nT = result.timeCount();
  std::vector<AllowableScale> out;
  out.reserve(static_cast<std::size_t>(nT));

  for (int k = 0; k < nT; ++k) {
    AllowableScale at;
    at.timeSeconds = result.times[static_cast<std::size_t>(k)];
    at.criteria.reserve(criteria.size());

    int binding = -1;
    double best = 0.0;
    for (std::size_t q = 0; q < criteria.size(); ++q) {
      const double response = tables[q].totals[static_cast<std::size_t>(k)];

      CriterionHeadroom headroom;
      headroom.name = criteria[q].name;
      headroom.limit = criteria[q].limit;
      headroom.response = response;

      if (response > 0.0) {
        headroom.fraction = response / criteria[q].limit;
        headroom.scale = criteria[q].limit / response;
        if (binding < 0 || headroom.scale < best) {
          best = headroom.scale;
          binding = static_cast<int>(q);
        }
      } else {
        // Nothing of this quantity is present, so no multiple of the inventory can exceed the
        // limit on it. Not a very large scale -- no scale at all.
        headroom.unbounded = true;
      }
      at.criteria.push_back(std::move(headroom));
    }

    if (binding >= 0) {
      at.bounded = true;
      at.scale = best;
      at.bindingIndex = binding;
      at.criteria[static_cast<std::size_t>(binding)].binding = true;
      at.limiting = limitingFor(tables[static_cast<std::size_t>(binding)], k, limitingCount);
    }
    out.push_back(std::move(at));
  }
  return out;
}

EventSeries scaleSeries(std::span<const AllowableScale> scaled) {
  std::size_t start = 0;
  while (start < scaled.size() && !scaled[start].bounded) {
    ++start;
  }

  EventSeries series;
  for (std::size_t k = start; k < scaled.size(); ++k) {
    if (!scaled[k].bounded) {
      // Bounded, then not, then bounded again. Interpolating across the gap would invent a
      // scale for a stretch where nothing constrained the inventory, so this refuses instead.
      for (std::size_t rest = k; rest < scaled.size(); ++rest) {
        if (scaled[rest].bounded) {
          throw InputError(
              tagged(kModule, "nothing constrains the inventory between t=" +
                                  std::to_string(scaled[k].timeSeconds) +
                                  " s and t=" + std::to_string(scaled[rest].timeSeconds) +
                                  " s, so the scale is not one curve over this grid. "
                                  "Restrict the time range to where a criterion binds"));
        }
      }
      break;
    }
    series.times.push_back(scaled[k].timeSeconds);
    series.values.push_back(scaled[k].scale);
  }

  if (series.times.size() < 2) {
    throw InputError(tagged(kModule,
                            "fewer than two times have a bounded scale, which is not enough "
                            "grid for a curve to be located in"));
  }
  return series;
}

}  // namespace nusift
