#pragma once
/**
 * @file
 * @brief Which contributors lead, and over which time windows.
 * @ingroup triage
 */
//
// Ranking answers "what dominates at 30 days". Forecasting answers the question people
// actually arrive with: "what dominates, and when does that change". A run over sixty time
// points already contains that answer; it is spread across sixty separate tables and nobody
// reads it that way.
//
// The output is a handful of windows -- Ba-140 leads from three days to six weeks, then Cs-137
// from two years onward -- which is both what a reader remembers and what decides when a
// measurement or a shielding assumption stops being valid.
//
#include <string>
#include <vector>

#include "nusift/triage/response.hpp"

namespace nusift {

// One stretch of time over which a single contributor leads.
struct DominanceWindow {
  ContributorId id;
  std::string label;
  double startSeconds = 0.0;
  double endSeconds = 0.0;
  // Largest share of the total this contributor reaches while it leads. A leader holding 80%
  // is a different situation from one holding 21% of a flat field, and the window alone does
  // not say which.
  double peakFraction = 0.0;
};

// One contributor followed across the whole grid.
struct RankTrack {
  ContributorId id;
  std::string label;
  std::vector<int> rank;         // [nT], 1-based; 0 where the contributor is absent
  std::vector<double> fraction;  // [nT]
  double peakFraction = 0.0;
  int peakTimeIndex = 0;
};

// The windows over which each successive leader holds first place.
//
// Boundaries are interpolated rather than snapped to grid points. Between two samples the
// contenders' values are close to exponential, so the log of their ratio is close to linear in
// time and crosses zero at a well-defined instant. Reporting the grid point instead would put
// the crossover wherever the grid happened to fall, which for a log-spaced grid at late times
// can be years away from the truth.
//
// Runs shorter than `minSamples` grid points are absorbed into their neighbours: near a
// crossover the leader can flicker between two contenders within numerical noise, and six
// one-sample windows are less true than one boundary.
//
// The threshold counts SAMPLES rather than elapsed time, which matters on a log-spaced grid.
// Measured as a fraction of the linear span, every window before the last decade is
// vanishingly short -- a 33-minute nuclide's genuine reign over the first hours is under
// 0.001% of a hundred-year axis -- so a duration test absorbs the entire early history into
// whichever contributor happened to lead first. Sample count is grid-agnostic and says what
// was meant: a run the grid barely resolved.
std::vector<DominanceWindow> dominanceWindows(const ResponseTable& table, int minSamples = 2);

// Every contributor that reaches the top `n` at any time, ordered by peak share. This is what
// a forecast prints rows for: a nuclide that only matters at thirty years still belongs in the
// table, and a nuclide that is never in the top `n` anywhere does not.
std::vector<RankTrack> unionTopN(const ResponseTable& table, int n);

// Contributors in the top `n` at EVERY time. The complement of the above, and the answer to
// "what do I always have to account for" as opposed to "what could ever matter".
std::vector<RankTrack> persistentTopN(const ResponseTable& table, int n);

}  // namespace nusift
