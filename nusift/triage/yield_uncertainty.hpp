#pragma once
/**
 * @file
 * @brief The error bar the evaluated fission yields put on a response, correlation included.
 * @ingroup triage
 */
//
// The other half of the error-bar question. uncertainty.hpp propagates what the ASSAY said, at
// fixed nuclear data; this propagates what the EVALUATION said about the yields, at a fixed
// number of fissions.
//
// It is the same propagation and deliberately so. A fission seed is n0_i = N_f * Y_i with N_f
// exact by assumption, so the seed's covariance is N_f^2 Sigma_Y and
//
//     sigma_R^2 = g^T (N_f^2 Sigma_Y) g,        Sigma_Y[i][j] = rho_ij sigma_i sigma_j
//
// with g the importance vector `attribute` already computes. Exact, not first-order, for the
// reason everything else here is exact: R is linear in the seed. The only new ingredient is
// rho, which no evaluation publishes and which yield_covariance.hpp imports.
//
// THREE FIGURES, NOT TWO, and the reason is coverage. The correlated figure can only be formed
// over products the imported matrix carries. So:
//
//   sigmaDiagonal         over every seeded product stating a sigma -- what the shipped diagonal
//                         treatment reports, and what the correlated figure is an improvement on
//   sigmaDiagonalMatched  the same sum restricted to products the matrix carries
//   sigmaCorrelated       that same restricted set, with the off-diagonal
//
// The pair to compare is the last two: they run over one set of products, so their whole
// difference is the correlation and nothing else. The first says what the matrix left out.
//
// EXPECT THE CORRELATION TO WIDEN THE ERROR BAR, which is the opposite of what the branching
// constraint does and is not a contradiction. A branching constraint is LOCAL -- a handful of
// modes forced to sum to one, anti-correlated, variance falls. The yield constraint is GLOBAL:
// one sum over a thousand products, so its anti-correlation is spread thin, while what is
// concentrated is the strong POSITIVE correlation inside a mass chain. A response that sums a
// chain sums those correlations with it. On U-235 thermal at 30 days the measured effect is a
// factor of 5.2 in variance, 2.06% against 4.72%.
//
// THE QUADRATIC FORM'S SIGN IS CHECKED RATHER THAN ASSUMED. The published matrices are not
// positive semi-definite -- see yield_covariance.hpp -- so g^T Sigma g can come out negative.
// When it does, no standard deviation is reported for it: a negative variance has no square
// root, and printing one would be inventing a number to fill a column.
//
#include <cstdint>
#include <string>
#include <vector>

#include "nusift/engine/decay_engine.hpp"
#include "nusift/nucdata/yield_covariance.hpp"
#include "nusift/seed/seed_fission.hpp"
#include "nusift/triage/response.hpp"

namespace nusift {

class NuclearData;

// One fission product's yield, and what that yield's uncertainty does to the answer.
struct YieldContribution {
  std::int64_t key = 0;
  std::string label;

  double yield = 0.0;       // atoms per fission, as the evaluation tabulated it
  double sigmaYield = 0.0;  // absolute 1-sigma on that yield; 0 when the evaluation states none
  double seedAtoms = 0.0;   // N_f * yield

  double importance = 0.0;  // dR/dn0 for this product, at `time`
  double share = 0.0;       // importance * seedAtoms -- this product's exact share of R

  // |importance| * N_f * sigmaYield: the standard deviation this product alone puts on R,
  // ignoring every correlation. Not additive; the variance fractions below are what sum to one.
  double sigmaContribution = 0.0;
  double varianceFraction = 0.0;  // of the DIAGONAL variance, which is the additive one

  // Whether the imported matrix carries this product. A product it does not carry keeps its
  // diagonal term and is absent from the correlated figure, which is what the coverage
  // fractions below exist to quantify.
  bool correlated = false;
};

struct YieldUncertainty {
  Metric metric = Metric::Activity;
  Unit unit = Unit::Becquerel;
  double time = 0.0;

  double fissions = 0.0;
  double response = 0.0;  // R = <g, n0>, equal to what a ranking reports for the same spec

  // The three figures, and their relatives against R. See the header note for which pair to
  // compare with which.
  double sigmaDiagonal = 0.0;
  double relativeDiagonal = 0.0;
  double sigmaDiagonalMatched = 0.0;
  double relativeDiagonalMatched = 0.0;
  double sigmaCorrelated = 0.0;
  double relativeCorrelated = 0.0;

  // The raw quadratic form g^T (N_f^2 Sigma_Y) g, SIGNED. Kept because the published matrices
  // are indefinite and a negative value is a real outcome that a square root would erase.
  double varianceCorrelated = 0.0;
  bool varianceNegative = false;
  // Ratio of the correlated variance to the matched diagonal one: what the correlation costs,
  // as one number. Zero when the form came out negative.
  double varianceRatio = 0.0;

  // The most negative eigenvalue of the correlation block actually contracted. Negative means
  // the published matrix is indefinite on this seed's products, which is expected and is
  // reported rather than repaired.
  double smallestEigenvalue = 0.0;

  // Share of R carried by products the matrix covers, and by products stating a sigma. A 5%
  // error bar over products carrying 60% of the response is not a 5% error bar on the response,
  // and neither figure can be inferred from the other.
  double coveredFraction = 0.0;
  double sigmaCoveredFraction = 0.0;
  int productsMatched = 0;
  int productsUnmatched = 0;
  int productsWithSigma = 0;
  int productsWithoutSigma = 0;
  // Products the matrix carries that this seed does not produce at all. Slack rather than a
  // problem, and the only visible sign of the edition gap when it is small.
  int productsUnusedInMatrix = 0;

  // The declared pairing, both halves. These differ by construction -- the matrices are built
  // for ENDF/B-VIII.0 and the store is VIII.1 -- and the report says so rather than letting the
  // agreement be assumed.
  YieldCovarianceProvenance provenance;
  std::string storeLibrary;
  double incidentEnergyEv = 0.0;

  // Ranked by varianceFraction, largest first.
  std::vector<YieldContribution> products;
};

// Propagate the evaluated yield uncertainties of `seed` to the response at `time`, using
// `correlation` for the off-diagonal.
//
// Costs exactly one adjoint solve -- the same solve `attribute` runs -- plus one dense
// contraction and one symmetric eigensolve over the matched products.
//
// Throws InputError when the store carries no yields for the seed (naming what it does carry),
// when no staged yield states an uncertainty at all, when the matrix and the seed share no
// product, for an interval unit or a gamma-line aggregate (as responseUncertainty refuses them,
// and for the same reasons), and for a negative time.
YieldUncertainty yieldUncertainty(const NuclearData& data, const seed::FissionSeed& seed,
                                  const YieldCorrelation& correlation, double time,
                                  const ResponseSpec& spec, const DecayOptions& options = {});

}  // namespace nusift
