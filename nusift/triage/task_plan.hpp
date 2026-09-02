#pragma once
/**
 * @file
 * @brief A job as a sequence of legs, each with its own distance and occupancy.
 * @ingroup triage
 */
//
// `taskSeries()` treats a job as one window at one distance, which is what a screening question
// needs and not what a job is. Real work has a shape: walk in, spend twenty minutes at arm's
// length from the thing that matters, step back to wait on a tool, walk out. Those legs are at
// different distances and the dose is not distributed the way their durations are, so a plan
// answers a question the single window cannot: WHICH LEG costs the dose, and would moving it or
// shortening it help.
//
// It is the interval integral applied leg by leg and summed:
//
//     accrued = sum over legs k of  occupancy_k * w(geometry_k) . integral over [t_k, t_k+D_k]
//
// Exact, on the same terms as everything else built on that integral: the atom-seconds over each
// leg's window are closed form, and the geometry is a fixed weight post-multiplied on top -- so a
// leg at 0.8 m and a leg at 3 m cost one solve each rather than one model each. Since item 1 the
// weight can be an effective-dose kernel, so this is a worker-DOSE plan rather than an exposure
// plan the day it is built.
//
// WHAT IT DOES NOT MODEL, and both matter more here than in a screening answer:
//
//   * No shielding. Every leg is an unshielded point source in air. A plan whose middle leg is
//     behind a wall cannot say so, and the number for that leg is an overestimate of unknown
//     size. This is the shielding item and not this one.
//   * Occupancy is a fraction, and multiplying by it assumes the presence is spread EVENLY over
//     the leg. Over a leg short against the decay time constant that is exact to floating point;
//     over a long one it is not, and the truth lies between "all of it at the start" and "all of
//     it at the end". The exact alternative needs no new machinery -- split the leg into the
//     stretches actually spent there, which a plan already expresses.
//
#include <span>
#include <string>
#include <vector>

#include "nusift/engine/decay_engine.hpp"
#include "nusift/engine/inventory.hpp"
#include "nusift/exposure/point_source.hpp"
#include "nusift/triage/response.hpp"

namespace nusift {

class NuclearData;

// One stretch of the job.
struct PlanLeg {
  // Required. A plan whose legs are numbered rather than named answers "leg 3 costs 60% of your
  // dose", which is not an answer anybody can act on.
  std::string name;
  double durationSeconds = 0.0;

  // The fraction of the leg actually spent in the field. Zero makes the leg a BREAK: the clock
  // runs, nothing accrues, and no solve is spent on it -- which is not an optimisation but the
  // definition, since a window nobody is standing in contributes nothing whatever the source is
  // doing.
  double occupancy = 1.0;

  // This leg's own geometry. Held per leg rather than taken from the spec because the distance
  // is the thing a plan varies; ResponseSpec::geometry is ignored entirely by this path, and
  // saying so here is cheaper than a reader wondering which of the two won.
  exposure::PointSourceGeometry geometry;
};

struct LegResult {
  std::string name;
  double startSeconds = 0.0;
  double endSeconds = 0.0;
  double occupancy = 1.0;
  double distanceM = 0.0;

  // In the plan's unit, already scaled by occupancy.
  double accrued = 0.0;
  double fraction = 0.0;    // of the plan total
  double cumulative = 0.0;  // including every earlier leg
  double cumulativeFraction = 0.0;

  // Accrued divided by the time actually spent there. The column that separates "this leg is
  // expensive because it is long" from "this leg is expensive because it is close", which is the
  // whole reason to break a job into legs rather than average it.
  double meanRate = 0.0;

  bool isBreak = false;
};

struct TaskPlan {
  Unit unit = Unit::Roentgen;
  double startSeconds = 0.0;
  double endSeconds = 0.0;

  double total = 0.0;
  double elapsedSeconds = 0.0;  // wall clock, breaks included
  double exposedSeconds = 0.0;  // sum of duration * occupancy, which is what earns the dose

  std::vector<LegResult> legs;

  // Optional budget accounting. `budget` is what the caller set; zero means none was asked
  // about and the two fields below are meaningless.
  double budget = 0.0;
  bool budgetSpent = false;
  // Where it runs out: the index of the leg the budget is exhausted in, and how far into that
  // leg. Located by the same root-find `stayTime()` uses, on that leg's own window -- so a plan
  // says not merely "this does not fit" but "you have to be out 6.2 minutes into the valve
  // work".
  int spentInLeg = -1;
  double spentAtSeconds = 0.0;
};

// Run the plan. `startSeconds` is when the first leg begins; each leg follows the one before it
// with no gap, so a wait between legs is a break leg rather than an implied silence.
//
// `spec` supplies the metric, the unit and any pack; its geometry is NOT used. The unit has to be
// an interval unit, because a leg accrues a total rather than holding a rate, and
// buildIntervalResponse() refuses the rest where that refusal already lives.
//
// `budget` is optional: pass a positive value to have the plan say where it runs out.
//
// COST: one interval integral per working leg -- two to three solves each -- and none at all for
// a break. A six-leg plan is a dozen or so solves, which is the same order as one stay time and
// far less than a task curve.
//
// Throws InputError for an empty plan, an unnamed leg, a non-positive duration, an occupancy
// outside [0, 1], or a negative start.
TaskPlan runTaskPlan(const NuclearData& data, const Inventory& inventory, const ResponseSpec& spec,
                     double startSeconds, std::span<const PlanLeg> legs, double budget = 0.0,
                     const DecayOptions& options = {});

}  // namespace nusift
