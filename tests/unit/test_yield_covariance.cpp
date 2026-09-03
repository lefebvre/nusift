// The imported yield correlation, and the error bar it puts on a fission source.
//
// Three things are worth checking beyond the arithmetic, and each of them is something the
// published product does rather than something the code might get wrong on its own.
//
// The KEY TRANSLATION, because FYCoM's Z*10000 + M*1000 + A and NuSIFT's Z*10000 + A*10 + I
// collide on legal keys rather than on malformed ones. 290068 is Cu-68 in one convention and
// Zn-6m8 in nothing; a translation that silently did not happen would contract the matrix
// against the wrong products and report a plausible number.
//
// The CONTRACTION against a hand-computable case, since g^T Sigma g for two products with a
// known correlation is a closed form and there is no excuse for a tolerance.
//
// And the INDEFINITENESS, because the published matrices are not positive semi-definite and the
// whole design rests on saying so rather than repairing it. A matrix that makes the variance
// come out negative is constructed here deliberately: it is the case where a square root would
// be an invented number, and the test pins that no square root is offered.
//
#include <gtest/gtest.h>

#include <cmath>
#include <fstream>
#include <string>

#include "nusift/core/error.hpp"
#include "nusift/engine/adjoint_engine.hpp"
#include "nusift/nucdata/nuclear_data.hpp"
#include "nusift/nucdata/yield_covariance.hpp"
#include "nusift/seed/seed_fission.hpp"
#include "nusift/triage/yield_uncertainty.hpp"
#include "synthetic_chain.hpp"

namespace nusift {
namespace {

const Zai kU235{92, 235, 0};
const Zai kSn{50, 100, 0};  // first product of the synthetic yield set
const Zai kSb{51, 100, 0};  // second

// A chain whose one fissionable parent splits into the two products of the synthetic chain,
// each with a stated yield uncertainty. Small enough that every figure below has a closed form.
NuclearData chainWithYields(double sigmaSn = 0.12, double sigmaSb = 0.08) {
  StoreArrays arrays = synth::linearChain({1.0e-3, 5.0e-4});
  arrays.nfyParentKey = {kU235.key()};
  arrays.nfyEnergyEv = {0.0253};
  arrays.nfySetOffset = {0, 2};
  arrays.nfyProductKey = {kSn.key(), kSb.key()};
  arrays.nfyProductYield = {1.2, 0.8};  // sums to 2.0, as independent yields must
  arrays.nfyProductYieldUncertainty = {sigmaSn, sigmaSb};
  return NuclearData::fromArrays(std::move(arrays));
}

// Write a FYCoM-shaped correlation CSV over the two products, with `rho` off the diagonal.
// Keys are written in FYCoM's convention, which is the point of writing a file at all.
std::string writeCorrelation(const std::string& path, double rho) {
  std::ofstream file(path);
  file << ", 500100, 510100\n";
  file << "500100, 1.0, " << rho << "\n";
  file << "510100, " << rho << ", 1.0\n";
  return path;
}

std::string tempPath(const char* name) {
  return std::string(::testing::TempDir()) + name;
}

seed::FissionSeed seed(double fissions = 1.0e20) {
  seed::FissionSeed s;
  s.fissile = kU235;
  s.incidentEnergyEv = 0.0253;
  s.fissions = fissions;
  return s;
}

// --- the key translation -----------------------------------------------------

// FYCoM puts the isomeric state in the thousands digit and NuSIFT puts it in the units digit,
// so the two conventions produce different keys for the same nuclide and the SAME key for
// different ones. Cu-68m is 291068 there and 290681 here; read either as the other's and the
// matrix contracts against a nuclide nobody asked about.
TEST(YieldCovariance, TranslatesFycomKeys) {
  EXPECT_EQ(yieldKeyFromFycom(290068), (Zai{29, 68, 0}.key()));   // Cu-68
  EXPECT_EQ(yieldKeyFromFycom(291068), (Zai{29, 68, 1}.key()));   // Cu-68m
  EXPECT_EQ(yieldKeyFromFycom(550137), (Zai{55, 137, 0}.key()));  // Cs-137
  EXPECT_EQ(yieldKeyFromFycom(561137), (Zai{56, 137, 1}.key()));  // Ba-137m

  // The collision that makes this worth a function rather than an inline expression: 290068 is
  // a legal key in both conventions and means different things in each.
  EXPECT_NE(yieldKeyFromFycom(290068), 290068);
}

// --- reading -----------------------------------------------------------------

TEST(YieldCovariance, ReadsASymmetricCorrelation) {
  const std::string path = writeCorrelation(tempPath("nusift_corr_ok.csv"), 0.5);
  const YieldCorrelation matrix = YieldCorrelation::read(path, "ENDF/B-VIII.0");

  EXPECT_EQ(matrix.size(), 2u);
  EXPECT_EQ(matrix.indexOf(kSn.key()), 0);
  EXPECT_EQ(matrix.indexOf(kSb.key()), 1);
  EXPECT_EQ(matrix.indexOf(Zai{55, 137, 0}.key()), -1);  // absent, not zero
  EXPECT_DOUBLE_EQ(matrix.at(0, 1), 0.5);
  EXPECT_DOUBLE_EQ(matrix.at(0, 0), 1.0);

  // The provenance is recovered from the path and the caller, and both halves reach the answer.
  EXPECT_EQ(matrix.provenance().library, "ENDF/B-VIII.0");
  EXPECT_FALSE(matrix.provenance().citation.empty());
}

// A covariance file has variances on its diagonal, not ones. Handing one over is the easiest
// mistake to make -- the files sit in the same directory and differ by one word in the name --
// and it would silently scale every answer by the yields' variances twice.
TEST(YieldCovariance, RefusesACovarianceFile) {
  const std::string path = tempPath("nusift_corr_cov.csv");
  {
    std::ofstream file(path);
    file << ", 500100, 510100\n";
    file << "500100, 0.0144, 0.0048\n";
    file << "510100, 0.0048, 0.0064\n";
  }
  EXPECT_THROW(YieldCorrelation::read(path, "ENDF/B-VIII.0"), InputError);
}

TEST(YieldCovariance, RefusesARaggedOrAsymmetricMatrix) {
  const std::string ragged = tempPath("nusift_corr_ragged.csv");
  {
    std::ofstream file(ragged);
    file << ", 500100, 510100\n";
    file << "500100, 1.0\n";
  }
  EXPECT_THROW(YieldCorrelation::read(ragged, "x"), InputError);

  const std::string asymmetric = tempPath("nusift_corr_asym.csv");
  {
    std::ofstream file(asymmetric);
    file << ", 500100, 510100\n";
    file << "500100, 1.0, 0.5\n";
    file << "510100, 0.4, 1.0\n";
  }
  EXPECT_THROW(YieldCorrelation::read(asymmetric, "x"), InputError);

  // Rows and columns in different orders would contract against the wrong products while
  // looking perfectly well formed.
  const std::string reordered = tempPath("nusift_corr_reordered.csv");
  {
    std::ofstream file(reordered);
    file << ", 500100, 510100\n";
    file << "510100, 1.0, 0.5\n";
    file << "500100, 0.5, 1.0\n";
  }
  EXPECT_THROW(YieldCorrelation::read(reordered, "x"), InputError);
}

// --- the contraction ---------------------------------------------------------

// Two products, so the quadratic form is a closed form:
//
//   var = (g1 N s1)^2 + (g2 N s2)^2 + 2 rho (g1 N s1)(g2 N s2)
//
// Checked at rho = 0 against the diagonal figure, and at rho = +/-0.5 against that expression
// with the importances the tool itself reports -- which is not circular, because the assertion
// is on the CONTRACTION and the importances are an input to it either way.
TEST(YieldUncertainty, ContractsAgainstTheClosedForm) {
  const NuclearData data = chainWithYields();
  const double time = 500.0;
  ResponseSpec spec;  // activity, becquerel

  const std::string zero = writeCorrelation(tempPath("nusift_corr_rho0.csv"), 0.0);
  const YieldUncertainty uncorrelated =
      yieldUncertainty(data, seed(), YieldCorrelation::read(zero, "test"), time, spec);

  // With no off-diagonal the correlated figure IS the diagonal one, to every digit. This is the
  // identity that says the two are the same sum and differ only by what is imported.
  EXPECT_NEAR(uncorrelated.sigmaCorrelated, uncorrelated.sigmaDiagonalMatched,
              uncorrelated.sigmaDiagonalMatched * 1.0e-12);
  EXPECT_NEAR(uncorrelated.varianceRatio, 1.0, 1.0e-12);
  EXPECT_EQ(uncorrelated.productsMatched, 2);
  EXPECT_EQ(uncorrelated.productsUnmatched, 0);

  // Both terms, from the report's own per-product figures.
  ASSERT_EQ(uncorrelated.products.size(), 2u);
  const double a = uncorrelated.products[0].sigmaContribution;
  const double b = uncorrelated.products[1].sigmaContribution;
  EXPECT_NEAR(uncorrelated.sigmaDiagonalMatched, std::hypot(a, b), std::hypot(a, b) * 1.0e-12);

  for (const double rho : {0.5, -0.5}) {
    const std::string path =
        writeCorrelation(tempPath(rho > 0 ? "nusift_corr_pos.csv" : "nusift_corr_neg.csv"), rho);
    const YieldUncertainty u =
        yieldUncertainty(data, seed(), YieldCorrelation::read(path, "test"), time, spec);
    // The signs of the two importances are both positive for an activity response, so a
    // positive correlation widens and a negative one narrows -- the mechanism the whole entry
    // turns on, in its two-product form.
    const double expected = a * a + b * b + 2.0 * rho * a * b;
    EXPECT_NEAR(u.varianceCorrelated, expected, std::abs(expected) * 1.0e-10) << "rho=" << rho;
    if (rho > 0) {
      EXPECT_GT(u.sigmaCorrelated, u.sigmaDiagonalMatched);
    } else {
      EXPECT_LT(u.sigmaCorrelated, u.sigmaDiagonalMatched);
    }
  }
}

// The response the propagation forms has to be the one a ranking reports, which is the identity
// that says the adjoint was run against the right seed at the right time.
TEST(YieldUncertainty, ResponseMatchesTheSeededTotal) {
  const NuclearData data = chainWithYields();
  const std::string path = writeCorrelation(tempPath("nusift_corr_resp.csv"), 0.3);
  const YieldCorrelation matrix = YieldCorrelation::read(path, "test");
  ResponseSpec spec;

  const Inventory inventory = seed::seedFromFission(data, seed());
  const std::vector<double> weight = responseWeights(data, spec);
  for (const double time : {0.0, 100.0, 5000.0}) {
    const YieldUncertainty u = yieldUncertainty(data, seed(), matrix, time, spec);
    const SeedImportance importance = seedImportance(data, inventory, weight, time);
    EXPECT_NEAR(u.response, importance.response, importance.response * 1.0e-12) << "at t=" << time;
  }
}

// A relative sigma on every yield gives that same relative sigma on the answer when the
// products are perfectly correlated, whatever the response weights are: rho = 1 makes the
// perturbation a common scaling of the whole seed, and R is linear in the seed. Exact, and a
// sharper check on the contraction than any tolerance-based one.
TEST(YieldUncertainty, PerfectCorrelationIsACommonScaling) {
  const double relative = 0.05;
  // sigma_i = 0.05 * Y_i for both products.
  const NuclearData data = chainWithYields(0.05 * 1.2, 0.05 * 0.8);
  const std::string path = writeCorrelation(tempPath("nusift_corr_one.csv"), 1.0);
  const YieldCorrelation matrix = YieldCorrelation::read(path, "test");
  ResponseSpec spec;

  for (const double time : {0.0, 1000.0, 20000.0}) {
    const YieldUncertainty u = yieldUncertainty(data, seed(), matrix, time, spec);
    EXPECT_NEAR(u.relativeCorrelated, relative, 1.0e-12) << "at t=" << time;
  }
}

// --- indefiniteness ----------------------------------------------------------

// The published matrices are not positive semi-definite, so a contraction can come out
// negative. The design decision is to report that rather than repair it, and this pins both
// halves: the diagnostic is negative, and no standard deviation is invented for a variance
// that has no square root.
TEST(YieldUncertainty, ReportsAnIndefiniteMatrixRatherThanRepairingIt) {
  const NuclearData data = chainWithYields();
  // |rho| > 1 is not a legal correlation and is exactly what a stochastic estimate can produce
  // once it has been post-processed. Its eigenvalues are 1 +/- |rho|, so this one is indefinite.
  // It has to be indefinite ENOUGH: the two products' terms are in a ratio of about three here,
  // and 1 + r^2 + 2 rho r stays positive until rho passes -(1 + r^2)/2r.
  const std::string path = writeCorrelation(tempPath("nusift_corr_indef.csv"), -2.0);
  const YieldCorrelation matrix = YieldCorrelation::read(path, "test");
  ResponseSpec spec;

  const YieldUncertainty u = yieldUncertainty(data, seed(), matrix, 500.0, spec);
  EXPECT_LT(u.smallestEigenvalue, 0.0);
  EXPECT_NEAR(u.smallestEigenvalue, -1.0, 1.0e-12);  // 1 - |rho|
  EXPECT_TRUE(u.varianceNegative);
  EXPECT_LT(u.varianceCorrelated, 0.0);
  // No square root offered, and no ratio either. The diagonal figure is untouched by any of it.
  EXPECT_DOUBLE_EQ(u.sigmaCorrelated, 0.0);
  EXPECT_DOUBLE_EQ(u.varianceRatio, 0.0);
  EXPECT_GT(u.sigmaDiagonalMatched, 0.0);
}

// --- coverage ----------------------------------------------------------------

// A matrix that carries only one of the seed's products still answers, over the products it
// covers, and says how much of the response that was. The diagonal figure stays over
// everything, so the two are not silently the same sum.
TEST(YieldUncertainty, ReportsPartialCoverage) {
  const NuclearData data = chainWithYields();
  const std::string path = tempPath("nusift_corr_partial.csv");
  {
    std::ofstream file(path);
    file << ", 500100\n";
    file << "500100, 1.0\n";
  }
  const YieldCorrelation matrix = YieldCorrelation::read(path, "test");
  ResponseSpec spec;

  const YieldUncertainty u = yieldUncertainty(data, seed(), matrix, 500.0, spec);
  EXPECT_EQ(u.productsMatched, 1);
  EXPECT_EQ(u.productsUnmatched, 1);
  EXPECT_GT(u.coveredFraction, 0.0);
  EXPECT_LT(u.coveredFraction, 1.0);
  // The restricted diagonal is strictly smaller than the full one, which is what makes the
  // three figures worth reporting separately.
  EXPECT_LT(u.sigmaDiagonalMatched, u.sigmaDiagonal);
}

// A matrix for the wrong fissioning system shares no product with the seed, and that is a
// refusal rather than a zero: contracting over an empty set would report an error bar of zero
// and look like a confident answer.
TEST(YieldUncertainty, RefusesAMatrixForAnotherSystem) {
  const NuclearData data = chainWithYields();
  const std::string path = tempPath("nusift_corr_other.csv");
  {
    std::ofstream file(path);
    file << ", 551370\n";
    file << "551370, 1.0\n";
  }
  const YieldCorrelation matrix = YieldCorrelation::read(path, "test");
  ResponseSpec spec;
  EXPECT_THROW(yieldUncertainty(data, seed(), matrix, 500.0, spec), InputError);
}

// A store staged before the uncertainty columns were carried has no sigma at all, and the
// honest answer is to say so rather than to report a zero error bar.
TEST(YieldUncertainty, RefusesAStoreWithNoStagedSigma) {
  const NuclearData data = chainWithYields(0.0, 0.0);
  const std::string path = writeCorrelation(tempPath("nusift_corr_nosigma.csv"), 0.5);
  const YieldCorrelation matrix = YieldCorrelation::read(path, "test");
  ResponseSpec spec;
  EXPECT_THROW(yieldUncertainty(data, seed(), matrix, 500.0, spec), InputError);
}

}  // namespace
}  // namespace nusift
