#pragma once
/**
 * @file
 * @brief How much of this inventory is allowed, against limits the caller supplies.
 * @ingroup triage
 */
//
// Ranking answers what dominates, forecasting when that changes, and the event engine when a
// curve reaches a value. This answers the question a shipper, a holder or a waste generator
// asks instead: HOW MUCH of this is allowed. By what factor could the inventory be multiplied
// before the first limit binds -- and, run over a grid, how that factor grows as the material
// decays.
//
//   s_max(t) = min_q  L_q / R_q(t)
//
// over the criteria q the caller supplies, each a quantity and a limit on it. The criterion
// achieving the minimum is the one that binds, and the contributors driving its response are
// the nuclides that decide the answer.
//
// This is EXACT in the sense the methodology docs use, and for the same reason the attribution
// shares are: every response is a linear functional of the inventory, so multiplying the
// inventory by s multiplies every R_q by exactly s. There is no search, no iteration, and no
// convergence criterion -- the scale at which a criterion binds is a division. What is not
// exact is anything the grid had to observe, which is why s_max(t) is reported per sample and
// located events over it are the event engine's business rather than this file's.
//
// What this deliberately does NOT know is whether a criterion is legitimate, whether it is
// genuinely linear, or how a regulation composes several of them. A sum-of-fractions index,
// a package-type condition, a per-table threshold, a fissile exception: those live in a rule
// layer above these weights, versioned with them. Here a criterion is a quantity and a number,
// and the caller is the one asserting that dividing by it means something.
//
#include <cstdint>
#include <span>
#include <string>
#include <vector>

#include "nusift/engine/decay_result.hpp"
#include "nusift/triage/events.hpp"
#include "nusift/triage/response.hpp"

namespace nusift {

class NuclearData;

// One thing the inventory is limited by: a quantity, and how much of it is allowed.
struct Criterion {
  // Named rather than numbered because naming the binding one is most of the answer. Required:
  // a report that says the binding criterion is "" has not said anything.
  std::string name;

  // The quantity limited, exactly as a ranking would compute it -- so the number checked
  // against the limit is the number a report would print for the same spec. Instantaneous
  // only: an interval unit describes an accrued total over a window, and a window is not what
  // a possession or transport limit constrains.
  ResponseSpec spec;

  // In spec.unit. Must be positive: a limit of zero forbids the material outright rather than
  // scaling it, and a negative one is not a limit.
  double limit = 0.0;
};

// What one criterion permits at one time.
struct CriterionHeadroom {
  std::string name;
  double response = 0.0;  // R_q(t), in the criterion's unit
  double limit = 0.0;

  // R/L, the sum-of-fractions idiom: the share of its limit this criterion already uses.
  double fraction = 0.0;
  // L/R, the same number inverted: the factor this criterion ALONE would permit.
  double scale = 0.0;

  bool binding = false;

  // The response is zero, so this criterion permits any scale and constrains nothing. Reported
  // rather than folded away as a very large number: "nothing here is limited by transport
  // index" and "the transport index allows 1e18 times this" are different statements, and only
  // one of them is true.
  bool unbounded = false;
};

// A contributor driving the binding criterion's response.
struct LimitingContributor {
  ContributorId id;
  std::string label;
  double fraction = 0.0;  // of the binding criterion's total
};

// The whole answer at one time.
struct AllowableScale {
  double timeSeconds = 0.0;

  // The largest factor the inventory can be multiplied by with every criterion still satisfied.
  // Meaningful only when `bounded`; left at zero otherwise.
  double scale = 0.0;

  // False when every criterion is unbounded, i.e. nothing constrains the inventory at this
  // time at all. A caller that treats an unset scale as zero would invert the meaning, so the
  // flag is checked rather than the number.
  bool bounded = false;

  // Index into the criteria the caller passed, or -1 when nothing binds.
  int bindingIndex = -1;

  std::vector<CriterionHeadroom> criteria;

  // The contributors deciding the binding criterion, most first. Empty when nothing binds.
  std::vector<LimitingContributor> limiting;
};

// s_max at every time on the decay grid, with the binding criterion and its limiting
// contributors named at each one.
//
// `limitingCount` caps how many contributors are named per time; 0 names none, which is the
// cheap path when only the curve is wanted.
//
// Throws InputError for an empty criteria list, a criterion with no name, a non-positive
// limit, or a unit that cannot express an instantaneous value.
std::vector<AllowableScale> allowableScale(const NuclearData& data, const DecayResult& result,
                                           std::span<const Criterion> criteria,
                                           int limitingCount = 3);

// s_max(t) as a curve the event engine can search: when the scale first reaches 1 is the date
// the inventory as it stands becomes shippable, holdable or releasable, and that is a crossing
// like any other.
//
// Unbounded samples carry no scale, so they are not part of the curve. The bounded samples are
// taken as one contiguous run: a grid where nothing constrains the inventory, then something
// does again, has no single curve to search, and that is refused rather than stitched across.
//
// The series carries no evaluator, so events located on it are interpolated inside their grid
// brackets. Refining one means re-solving the decay at a candidate time and re-forming the
// minimum there, which is a caller's decision to pay for -- attach EventSeries::evaluate and
// every event over the curve is refined instead.
EventSeries scaleSeries(std::span<const AllowableScale> scaled);

}  // namespace nusift
