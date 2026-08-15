#pragma once
/**
 * @file
 * @brief Ranks the contributors in a response table and reports how much of the total they
 *        actually account for.
 * @ingroup triage
 */
//
// The question NuSIFT exists to answer -- which isotopes or mass chains dominate -- is a
// ranking question, and a ranking that omits how much it left out is misleading. So every
// Ranking carries the total over ALL contributors and the fraction the returned rows cover,
// regardless of how it was truncated. A top-10 that accounts for 40% of the activity and a
// top-10 that accounts for 99% look identical without that number, and they mean entirely
// different things.
//
#include <string>
#include <vector>

#include "nusift/triage/response.hpp"

namespace nusift {

struct RankRequest {
  // Largest N to return; 0 means no limit.
  int topN = 10;

  // If > 0, return the SMALLEST prefix whose cumulative fraction reaches this -- the literal
  // "which nuclides are 95% of the activity". Applied after topN, so setting both gives
  // whichever is shorter.
  double coverage = 0.0;

  // Drop contributors below this fraction of the total. Trims the long tail of nuclides that
  // are present but irrelevant.
  double minFraction = 0.0;

  // Contributor keys that appear whatever they rank, as resolved by requirePin(). Every other
  // knob here is a way of truncating the ranking; this is the one that reaches past the cut.
  //
  // A pin is not a filter and does not reorder anything: the ranking is still the top of the
  // table, with the pinned contributors appended below it carrying the rank they actually
  // hold. That is the only way "and where does Cs-137 stand?" can be asked without either
  // printing the whole chain or guessing a --top large enough to reach it.
  std::vector<std::int64_t> pinned;
};

struct Contributor {
  ContributorId id;
  std::string label;
  double value = 0.0;               // in the table's unit
  double fraction = 0.0;            // of the total over all contributors
  double cumulativeFraction = 0.0;  // including every higher-ranked contributor
  // 1-based, over ALL contributors rather than over the returned ones -- so a pinned row says
  // where it really stands. Zero means it contributes nothing at this time and so holds no
  // place in the ordering at all, which only a pinned row can be.
  int rank = 0;
  int flags = kFlagNone;
  // Present because it was pinned rather than because it ranked. Says nothing about the row's
  // numbers, which are the true ones either way -- it says the row is not part of the prefix
  // above it, and a report that ran the two together would misstate what the ranking covers.
  bool pinned = false;
};

struct Ranking {
  Metric metric = Metric::Activity;
  Aggregate aggregate = Aggregate::Nuclide;
  Domain domain = Domain::Instant;
  Unit unit = Unit::Becquerel;

  double time = 0.0;
  double timeEnd = 0.0;  // interval domain only

  // Over every contributor, not just the returned ones.
  double total = 0.0;
  // Sum of the returned contributors' fractions. The honesty number: a truncated ranking
  // reports this so it can never present itself as complete. Pinned rows count toward it like
  // any other returned row, which is why it is a sum over what was returned rather than the
  // cumulative of the last one -- with a pin below the cut those two are no longer the same
  // number, and only the sum describes what the reader can actually see.
  double coveredFraction = 0.0;
  int omittedCount = 0;

  // Exposure only: the fraction of emitted photon energy NuSIFT does not model at this time.
  // Carried on the ranking so a report can state the magnitude of what is missing without
  // needing the table it came from.
  double unmodeledEnergyFraction = 0.0;

  std::vector<Contributor> contributors;
};

// Rank one time slice.
//
// Ordering is by value descending, then by contributor key ascending. The tiebreak is
// specified rather than incidental: exactly-tied contributors are common with symmetric
// synthetic inputs, and without a deterministic secondary key their order varies with the
// sort implementation and golden tests become flaky across platforms.
//
// Any pinned contributor the ranking did not already reach follows the ranked rows, in the
// order it stands in -- so the returned list is a ranked prefix followed by a pinned tail, and
// never an interleaving of the two.
Ranking rank(const ResponseTable& table, int timeIndex, const RankRequest& request);

// One ranking per time in the table.
std::vector<Ranking> rankAll(const ResponseTable& table, const RankRequest& request);

}  // namespace nusift
