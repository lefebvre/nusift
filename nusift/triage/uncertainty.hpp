#pragma once
/**
 * @file
 * @brief The error bar an assay's uncertainty puts on a response, and which row dominates it.
 * @ingroup triage
 */
//
// The one uncertainty question NuSIFT can answer without any evaluated data at all.
//
// `attribute` already computes the importance vector g_i = dR/dn0_i, and R is LINEAR in n0, so
// given a covariance on the seed the propagation is exact rather than first-order:
//
//     sigma_R^2 = g^T Sigma g
//
// with no series expansion, no assumption of small errors, and no derivative that had to be
// estimated. Nothing else in the tool's error story is this clean, and the reason is the reason
// everything else here is exact: decay is linear and the metric is a fixed weight on top of it.
//
// WHAT THE ROW UNCERTAINTIES ARE, AND WHAT THEY ARE NOT. A per-row sigma supplies only the
// DIAGONAL of Sigma, which asserts that the assay errors are independent. Two aliquots counted
// on the same detector against the same standard are not independent, and a sheet whose rows
// were fitted to a total is not either. The report says so; it cannot detect it.
//
// UNITS CARRY THROUGH FOR FREE, which was not expected. toAtoms() is value * k in every branch
// -- 1 for atoms, N_A for moles, N_A*g/M for a mass, 1/lambda for an activity -- so it is
// strictly linear with no offset, and the same call converts a sigma as correctly as it converts
// a quantity. The measurement basis therefore does NOT have to be retained for this question.
// It does have to be retained for the OTHER half of the uncertainty item, dR/dlambda, where
// n0 = A0/lambda moves with the parameter being differentiated -- see attribution.md section 6.
// Those are different questions and only one of them needs the basis.
//
// ASSAYS TAKEN ON DIFFERENT DATES NEED NO NEW MACHINERY, which was also not expected. A sigma
// cannot ride on a reconciled inventory -- a diagonal Sigma at assay becomes D Sigma D^T at the
// epoch and D is not diagonal -- so the propagation has to reach back to assay time. But
//
//     R = <w, exp(A T) sum_a exp(A tau_a) n_a> = sum_a <exp(A^T (T + tau_a)) w, n_a>
//
// because the two exponentials share one matrix and commute. So the importance of an assay
// carried forward by tau is just the ordinary adjoint run for a longer time, T + tau. One
// existing solve per assay, and R = sum_a <g_a, n_a> is an identity a test can check.
//
#include <cstdint>
#include <span>
#include <string>
#include <vector>

#include "nusift/engine/decay_engine.hpp"
#include "nusift/engine/reconcile.hpp"
#include "nusift/triage/response.hpp"

namespace nusift {

class NuclearData;

// One seeded row's contribution to the answer and to its error bar.
struct SeedUncertainty {
  std::int64_t key = 0;
  std::string label;

  // Which assay this row came from, and how far it was carried to reach the epoch. Both are
  // reported because an old assay's uncertainty is carried the whole way with it: a row
  // measured to 5% two years ago is still 5% of whatever it has become.
  std::string assay;
  double carriedSeconds = 0.0;

  double seedAtoms = 0.0;
  double sigmaAtoms = 0.0;  // 1-sigma on the seed, in atoms; zero when the row stated none

  double importance = 0.0;  // dR/dn0 for THIS row, at T + carried
  double share = 0.0;       // importance * seedAtoms -- this row's exact share of R

  // |importance| * sigmaAtoms: the standard deviation this row alone puts on R. Not additive --
  // the variance fractions below are what sum to one.
  double sigmaContribution = 0.0;
  // (sigmaContribution / sigma_R)^2. THE assay-planning number: it says which measurement to
  // improve, and it is a share of the variance rather than of the answer, so a row with a small
  // share of R can dominate it.
  double varianceFraction = 0.0;
};

struct ResponseUncertainty {
  Metric metric = Metric::Activity;
  Unit unit = Unit::Becquerel;
  double time = 0.0;

  // R itself, formed as sum over assays of <g_a, n_a>. Equal to what a ranking reports for the
  // same spec and time, which is the identity that says the per-assay adjoints were run at the
  // right times.
  double response = 0.0;

  double sigma = 0.0;     // sqrt(sum of (g_i sigma_i)^2)
  double relative = 0.0;  // sigma / response, when response > 0

  // How much of the answer rests on rows that stated no uncertainty at all. The figure that
  // decides whether the error bar means anything: a 3% error bar over rows carrying 40% of the
  // response is not a 3% error bar on the response.
  double coveredFraction = 0.0;
  int rowsWithSigma = 0;
  int rowsWithoutSigma = 0;

  // Ranked by varianceFraction, largest first.
  std::vector<SeedUncertainty> seeds;
};

// Propagate the assay uncertainties in `assays` to the response at `time` AFTER `epochSeconds`.
//
// Costs one adjoint solve per assay -- the same solve `attribute` runs, at a later time for an
// assay carried further. An undated inventory is one assay carried nothing and costs exactly
// what attribution costs.
//
// Throws InputError for an empty assay list, an epoch earlier than an assay (the same refusal
// reconcile() makes and for the same reason), a negative time, an interval unit (the shares are
// instantaneous, as attributeToSeed's are), and a gamma-line aggregate.
ResponseUncertainty responseUncertainty(const NuclearData& data, std::span<const AssayGroup> assays,
                                        double epochSeconds, double time, const ResponseSpec& spec,
                                        const DecayOptions& options = {});

}  // namespace nusift
