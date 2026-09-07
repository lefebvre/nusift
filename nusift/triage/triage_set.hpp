#pragma once
/**
 * @file
 * @brief The smallest set of contributors that stays above a coverage floor everywhere at once.
 * @ingroup triage
 */
//
// `rank --coverage 0.95` answers "which nuclides are 95% of this, HERE" -- one metric, one
// instant. This answers the question a monitoring list actually poses: which nuclides are 95% of
// it EVERYWHERE, across every metric that matters and every time on the grid, as one list.
//
// The two are not the same question and the second is not obtainable by repeating the first. The
// union of each time's top-N is a set nobody chose: it is as large as the worst time makes it,
// it carries whatever the other times happened to rank highly, and nothing about it is minimal.
// Worse, a union of per-metric lists guarantees nothing about the metrics jointly -- a nuclide
// dropped from the activity list because it ranked eleventh there may be the one holding up the
// exposure floor at a later time.
//
// So it is a covering problem, over the values matrices that already exist:
//
//     choose S minimising |S|  subject to  sum over S of value(r, t, c) >= f_r * total(r, t)
//                                          for every requirement r and every time t
//
// Set cover is NP-hard and this is a deterministic greedy, which is the honest thing to run at
// this size and the honest thing to say about it: the result is a SMALL set that meets the
// floor, not a proof that none smaller exists. Greedy is within a ln(n) factor of optimal and in
// practice much closer, and the alternative -- an exact solver over ten thousand gamma lines --
// would buy a nuclide or two at a cost nobody asked for.
//
// What the result is NOT is a truncation. Every other knob in ranking.hpp cuts a list down; this
// one builds a list up until a stated property holds, and the property is checked against the
// full totals rather than against the rows it selected.
//
#include <cstdint>
#include <span>
#include <string>
#include <vector>

#include "nusift/triage/response.hpp"

namespace nusift {

// One thing the set has to be good for: a sampled response, and the share of its total the
// chosen contributors have to account for at EVERY time in it.
//
// Several of these is the point. One requirement per metric asks for a list that serves them
// jointly; the same metric against two evaluated libraries asks for one that survives the
// disagreement between them.
struct CoverageRequirement {
  const ResponseTable* table = nullptr;
  // In (0, 1]. The floor, applied at every time independently -- not to an average over them,
  // which would let a list fail badly somewhere and pass on the strength of everywhere else.
  double fraction = 0.95;
  // What this requirement is, for the report: "exposure at 2 m", "ENDF/B-VIII.1 activity".
  // Required, because a shortfall reported against requirement 2 of 4 names nothing.
  std::string label;
};

// A contributor in the set, and why it is there.
struct SetMember {
  ContributorId id;
  std::string label;

  // The order greedy chose it in, from 1. Not a ranking: the second member is the one that most
  // improved the worst-covered constraints GIVEN the first, which is usually not the second
  // largest contributor to anything. Reported because the order is the argument for the set.
  int order = 0;

  // What this member closed when it was chosen, summed over every constraint still unmet at
  // that moment. The marginal value that earned it its place.
  //
  // Each term is a FRACTION of that constraint's own required coverage, not an amount in the
  // requirement's unit, so the sum is dimensionless and a run in Ci reports what the same run
  // in Bq does. It is therefore comparable between members of one set, and not between sets
  // built from different numbers of requirements -- an unmet constraint can contribute at most
  // 1 to it.
  double closedShortfall = 0.0;

  // The largest share of any single requirement's total this contributor holds, over all times.
  // What a reader checks the list against: a member whose peak share is 0.1% is there to hold up
  // some particular instant, not because it is ever large.
  double peakFraction = 0.0;
};

// How a finished set fares against one requirement at one time.
struct CoveragePoint {
  std::string requirement;
  double timeSeconds = 0.0;
  double achieved = 0.0;  // fraction of that requirement's total the set accounts for
  double required = 0.0;
};

struct TriageSet {
  std::vector<SetMember> members;

  // Where the set is closest to failing: the smallest achieved-minus-required margin over every
  // requirement and time. This is the number that says whether the answer is comfortable or
  // exactly on the edge, and it is the one a reader should look at first.
  CoveragePoint binding;

  // Constraints the set does not meet, empty when it meets them all. Not a failure of the
  // search: a gamma-line table's total includes lines below the column floor, so selecting
  // EVERY column can still fall short of 100%, and a floor set above what the table can express
  // is unreachable by any set. Reported rather than thrown, because the set is still the best
  // one available and the reader needs to see how far short it lands.
  std::vector<CoveragePoint> shortfalls;

  // Contributors considered, and constraints checked. The size of the problem the answer came
  // from, which is what makes "23 out of 3800" a different statement from "23 out of 25".
  int candidateCount = 0;
  int constraintCount = 0;
};

// The smallest set greedy finds. Members come back in the order chosen.
//
// Every requirement's table must share one aggregate: a nuclide column and a gamma-line column
// are not the same kind of thing, and a set mixing them would satisfy neither question. Throws
// InputError for that, for an empty requirement list, for a null table, for a fraction outside
// (0, 1], and for an unlabelled requirement.
//
// Ties are broken by contributor key so the result is reproducible, on the same terms the
// ranking sort is: an answer that depended on iteration order would be a different monitoring
// list on a different machine.
TriageSet robustTriageSet(std::span<const CoverageRequirement> requirements);

}  // namespace nusift
