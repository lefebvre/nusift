#pragma once
/**
 * @file
 * @brief Attributes a response to the SEEDED nuclides it came from, rather than the ones
 *        emitting it now.
 * @ingroup triage
 */
//
// TWO ATTRIBUTIONS OF THE SAME NUMBER.
//
// triage/ranking.hpp answers "what is producing this response right now" -- it ranks the
// nuclides present at the response time, which is why a Cs-137 source's exposure lands on the
// Ba-137m daughter that actually emits the photons. That is the correct answer to the question
// it asks.
//
// This asks a different one: which of the nuclides I SEEDED is that response riding on. For an
// inventory read from a file the two questions are close. For a fission seed they are almost
// disjoint -- the emitters are what survives at the response time, while the seeds are the
// short-lived precursors whose yields determined how much of it there would be. Xe-140 is gone
// long before 30 days and is a top seed contributor at 30 days; Ba-140 is a top emitter and
// was barely seeded at all. Neither ranking can be derived from the other.
//
// It is an exact decomposition, not an estimate. Decay is linear, so the shares below sum to
// the same total the ordinary ranking reports for the same metric, time, and geometry -- which
// is why this carries a coveredFraction like any other ranking, and why it needs no
// uncertainty data to be meaningful. See engine/adjoint_engine.hpp for why one adjoint solve
// yields all of them at once.
//
// NOT MODELED HERE: interval domains. The shares are of an instantaneous response; attributing
// a time-integrated total to its seed needs the integrated adjoint, which is a different solve.
//
#include <cstdint>
#include <string>
#include <vector>

#include "nusift/engine/decay_engine.hpp"
#include "nusift/engine/inventory.hpp"
#include "nusift/triage/ranking.hpp"
#include "nusift/triage/response.hpp"

namespace nusift {

class NuclearData;

// One seeded nuclide's share of the response.
struct SeedShare {
  std::int64_t key = 0;  // ZAI key of the SEEDED nuclide
  std::string label;

  double seedAtoms = 0.0;  // n0_i, atoms placed in the inventory
  // dR/dn0_i, in the ranking's unit per seed atom. The intensive half of the share: a large
  // importance with a small seed is a nuclide worth seeding more accurately, which is not the
  // same statement as a large share.
  double importance = 0.0;
  double value = 0.0;  // seedAtoms * importance -- this seed's exact share of the total

  double fraction = 0.0;            // of the total over all seeds
  double cumulativeFraction = 0.0;  // including every higher-ranked seed
  int rank = 0;                     // 1-based over ALL seeds, as Contributor::rank is
  bool pinned = false;
};

struct SeedAttribution {
  Metric metric = Metric::Activity;
  Unit unit = Unit::Becquerel;
  double time = 0.0;

  // Over every seeded nuclide. Equal, to rounding, to the total a Ranking reports for the same
  // spec and time -- the two are partitions of one number.
  double total = 0.0;
  double coveredFraction = 0.0;
  // Over the seeds that CONTRIBUTE. A seeded nuclide worth nothing to this metric at this time
  // -- a stable seed, or a pure beta emitter under an exposure metric -- is not a contributor
  // that happened to be cut, so counting it here would report an omission the reader can do
  // nothing about and that would not change the total if it were shown.
  int omittedCount = 0;

  // The same two caveats a Ranking carries for the same spec and time, and the same numbers:
  // an exposure attributed to its seed is the identical figure the ranking reports, so it is
  // understated by the identical amount. Zero for Metric::Activity. See ExposureCaveats.
  double unmodeledEnergyFraction = 0.0;
  double meanOpticalDepth = 0.0;
  double buildup = 1.0;
  std::vector<std::string> unmodeledContinuum;

  std::string seedProvenance;
  std::vector<SeedShare> shares;
};

// Attribute R at `time` to the nuclides in `inventory`.
//
// `request` is honoured exactly as it is for a Ranking: topN, coverage, and minFraction
// truncate, and `pinned` reaches past the cut carrying the rank it really holds. A pin here
// names a SEEDED nuclide, so pinning one the inventory does not carry is refused rather than
// answered with a row of zeros.
//
// A seed worth nothing to this metric at this time holds no place in the ordering and is not
// among the shares -- it would occupy a topN slot ahead of a real contributor and inflate the
// omitted count without adding anything a reader could act on. Pinning one still returns it,
// rankless, because "Cs-133 contributes nothing here" is a real answer to a question that was
// actually asked.
//
// For Metric::Exposure this also runs the forward solve the caveats need, so an attributed
// exposure carries the same warnings the ranking of the same number does. That is a second
// CRAM solve on top of the adjoint one; activity attribution does only the adjoint.
//
// Throws InputError for an interval unit (the shares are instantaneous -- see the file note),
// for a gamma-line aggregate (a photon line has no seed), and on the same terms buildResponse
// does for a spec the store cannot answer.
SeedAttribution attributeToSeed(const NuclearData& data, const Inventory& inventory, double time,
                                const ResponseSpec& spec, const RankRequest& request,
                                const DecayOptions& options = {});

// Resolve a nuclide spelling against the seeds `inventory` carries, for RankRequest::pinned.
// Throws InputError naming what went wrong: an unparseable name, or one the inventory does not
// seed.
std::int64_t requireSeedPin(const NuclearData& data, const Inventory& inventory,
                            std::string_view text);

}  // namespace nusift
