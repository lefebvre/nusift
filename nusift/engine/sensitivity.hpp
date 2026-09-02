#pragma once
/**
 * @file
 * @brief How much a response depends on the EVALUATED data, rather than on the inventory.
 * @ingroup engine
 */
//
// adjoint_engine.hpp answers dR/dn0: how the response moves with the seed, which is a partition
// and needs no evaluated data to mean anything. This answers the other one -- dR/dlambda, how it
// moves with the half-lives the evaluation supplied -- and it is a genuinely different
// calculation rather than the same solve read differently.
//
// THREE TERMS, AND ONLY ONE OF THEM COMES FROM THE ADJOINT.
//
//   dR/dlambda_i  =  implicit          the response moving because the MATRIX moved
//                  + explicit          ... because the WEIGHT moved, since w = lambda * k
//                  + basis             ... because the SEED moved, for a row given as activity
//
// cram defines R = <w, n(T)> with w held fixed, so its quadrature returns the implicit term
// alone. NuSIFT's weight IS lambda times something -- activity is lambda, exposure is lambda
// times a photon sum -- so the explicit term (dw_i/dlambda_i) * n_i(T) is real and is the same
// order as the first. And an inventory row written in becquerel fixes A0 rather than n0, so
// n0 = A0/lambda moves too and contributes dR/dn0_i * (-n0_i/lambda_i). A tool that reported one
// or two of the three would be answering a question nobody asked.
//
// REPORT THE ELASTICITY, NOT THE DERIVATIVE. At secular equilibrium the implicit and explicit
// terms cancel across seven orders of magnitude -- Ba-137m's activity is pinned by its parent's
// feed rate, so moving its own lambda barely moves R -- and a relative error against that
// difference is meaningless. The elasticity lambda_i * (dR/dlambda_i) / R stays interpretable
// through the cancellation, and is what an error budget wants anyway: a fractional response to a
// fractional perturbation, directly multiplicable by a relative sigma.
//
// THE QUADRATURE IS A FOOTGUN AT ITS DEFAULTS, and this is where that is handled rather than
// left to the caller. The integrand is a product of two trajectories, each carrying a boundary
// layer at the interval ends, and cram's default refinement suits "a thermal pin with a few
// intervals of days to months". A decay chain spans microsecond isomers to primordials: a single
// thirty-day interval at the default is 37% wrong for the Cs-137/Ba-137m pair, silently. cram's
// own documented rule fixes it -- refine until the smallest piece is below the shortest removal
// time in the problem -- and refinement costs two more pieces each, not twice as many, so
// satisfying it is cheap. chooseEndRefinements() applies it from the pruned set.
//
#include <cstdint>
#include <span>
#include <string>
#include <vector>

#include "nusift/engine/decay_engine.hpp"
#include "nusift/engine/inventory.hpp"
#include "nusift/triage/response.hpp"

namespace nusift {

class NuclearData;

// One nuclide's decay constant, and what the response does when it moves.
struct DecaySensitivity {
  std::int64_t key = 0;
  std::string label;

  double decayConstant = 0.0;
  double halfLifeSeconds = 0.0;

  // The three terms, kept apart because they are separately meaningful and because their SUM is
  // the thing that cancels. A reader who sees only the total cannot tell a genuinely
  // insensitive nuclide from two large terms that happened to annihilate.
  double implicit = 0.0;
  double explicitWeight = 0.0;
  double basis = 0.0;
  double total = 0.0;  // in the response's unit per (1/s)

  // lambda * (dR/dlambda) / R. Dimensionless, sign-carrying, and the number to read: a
  // fractional change in the response per fractional change in the half-life, so multiplying it
  // by a relative sigma gives a relative contribution to the error bar directly.
  double elasticity = 0.0;

  // The evaluated 1-sigma on this nuclide's half-life, as a FRACTION, or zero when the store
  // stages none. ENDF/B-VIII.1 states one for 85% of staged nuclides.
  double relativeUncertainty = 0.0;
  // |elasticity| * relativeUncertainty: this nuclide's relative contribution to sigma_R. Zero
  // where no uncertainty is evaluated, which is not the same as a zero contribution and is why
  // the report counts how much of the answer that covers.
  double sigmaContribution = 0.0;
};

// One decay mode's branching fraction, and what the response does when it moves.
struct BranchingSensitivity {
  std::int64_t parentKey = 0;
  std::string parent;
  std::string daughter;
  std::string mode;  // the ENDF RTYP, spelled

  double branching = 0.0;
  double sigma = 0.0;  // absolute 1-sigma on the branching, or 0 when none is evaluated

  double total = 0.0;       // dR/db, in the response's unit
  double elasticity = 0.0;  // b * (dR/db) / R
};

// What the branchings of ONE nuclide contribute to the variance, both ways.
//
// The constraint is per nuclide -- a nuclide's modes sum to one, and nothing couples two
// nuclides' branchings -- so the covariance is block diagonal with one small block each, and
// the contraction is exact rather than approximated.
struct BranchingBlock {
  std::int64_t parentKey = 0;
  std::string parent;
  int modes = 0;

  // Sum over modes of (dR/db_i * sigma_i)^2: what a DIAGONAL treatment would report. Kept so the
  // report can show what the constraint changes rather than asserting that it matters.
  double diagonalVariance = 0.0;
  // e^T C e with C the constrained covariance. This is the honest figure.
  double constrainedVariance = 0.0;
};

struct DecaySensitivities {
  double time = 0.0;
  double response = 0.0;  // R, from the same forward march the quadrature used

  // Ranked by |sigmaContribution|, or by |elasticity| when nothing carries an uncertainty.
  std::vector<DecaySensitivity> nuclides;

  // Root-sum-square of the per-nuclide contributions, as a FRACTION of R. This is a sensitivity
  // norm and NOT an error budget: it takes Sigma diagonal, and evaluated half-lives are not
  // independent of the branchings and yields fitted alongside them. See attribution.md section 6.
  double relativeNorm = 0.0;
  // Share of |dR/dlambda| carried by nuclides whose half-life uncertainty IS evaluated. A norm
  // over 40% of the sensitivity is not a norm over the answer.
  double coveredFraction = 0.0;
  int withUncertainty = 0;
  int withoutUncertainty = 0;

  // What the quadrature actually did, because its cost and its correctness are the same knob.
  // --- branchings, and the one correlation that is derivable -----------------
  //
  // ENDF carries no covariance for decay data at all: the tapes hold MF1 and MF8 and nothing
  // else. But a nuclide's branching fractions SUM TO ONE, and that constraint is a fact about
  // the data model rather than an evaluated quantity -- so the correlation it induces can be
  // derived rather than imported. Imposing it on the stated sigmas gives
  //
  //     C = D - (D u u^T D) / (u^T D u),      D = diag(sigma^2),  u = (1, 1, ... 1)
  //
  // the nearest covariance consistent with both, and C u = 0 exactly, so a perturbation can
  // never take the branchings off the simplex they live on.
  //
  // Two consequences are worth expecting. A nuclide with ONE mode has b = 1 by construction and
  // contributes exactly nothing, whatever sigma the evaluation states for it -- a diagonal
  // treatment credits it with a variance it cannot have. And a nuclide with TWO modes is forced
  // entirely: db_1 = -db_2, so the two sigmas must be equal and the correlation is exactly -1,
  // which the projection recovers rather than assumes.
  std::vector<BranchingSensitivity> branchings;
  std::vector<BranchingBlock> branchingBlocks;
  // sqrt(sum of constrained variances) / R, and the same with the constraint ignored. Both are
  // reported because the difference between them IS the result: it says what treating the
  // branchings as independent would have cost.
  double branchingNorm = 0.0;
  double branchingNormDiagonal = 0.0;

  int endRefinements = 0;
  int solves = 0;
  // True when the refinement rule asked for more than `maxEndRefinements` allowed. The answer is
  // then UNDER-REFINED by cram's own criterion and may be wrong by tens of percent -- reported
  // rather than silently returned, since nothing about the numbers would show it.
  bool refinementCapped = false;
  double shortestRemovalSeconds = 0.0;
};

struct SensitivityOptions {
  // The schedule the quadrature integrates over. More intervals resolve interior variation;
  // they do NOT substitute for refinement, which the spike measured directly -- an eight-point
  // log grid at the default refinement is worse than one interval properly refined, for four
  // times the work.
  int scheduleIntervals = 1;
  // A ceiling on the automatic refinement, defaulting to cram's own hard limit of 30. A fission
  // source reaches it: the chain carries microsecond isomers, so resolving every boundary layer
  // in it would need more pieces than cram will take. That is reported rather than worked
  // around -- see DecaySensitivities::refinementCapped.
  int maxEndRefinements = 30;
  int gaussPoints = 5;
  int subIntervals = 4;
};

// The refinement cram's own rule asks for: enough that the smallest quadrature piece falls below
// `shortestRemovalSeconds`, given `intervalSeconds` split into `subIntervals` pieces. Exposed
// because it is the one number that decides whether the answer is right, and a caller that wants
// to check it against a finer rule has to be able to ask what was used.
int chooseEndRefinements(double intervalSeconds, int subIntervals, double shortestRemovalSeconds,
                         int cap);

// dR/dlambda for every nuclide the seed can reach, at `time`.
//
// COST, and it is unlike anything else in the library. This is quadrature over the forward and
// adjoint trajectories rather than a matrix exponential: roughly
// (subIntervals + 2*endRefinements) * gaussPoints * 2 solves per schedule interval, which for a
// decay problem needing twelve refinements is some three hundred solves where a ranking costs
// one. Seconds, not milliseconds, and `solves` reports the count.
//
// Throws InputError for a non-positive time (the derivative at t = 0 is identically zero and the
// schedule would be empty), an interval unit, a gamma-line aggregate, and on the same terms
// buildResponse does for a spec the store cannot answer.
DecaySensitivities decaySensitivities(const NuclearData& data, const Inventory& inventory,
                                      double time, const ResponseSpec& spec,
                                      const SensitivityOptions& sensitivity = {},
                                      const DecayOptions& options = {});

}  // namespace nusift
