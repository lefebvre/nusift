#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "nusift/core/error.hpp"
#include "nusift/exposure/air_coefficients.hpp"
#include "nusift/exposure/point_source.hpp"
#include "nusift/nucdata/photon_lines.hpp"

namespace nusift::exposure {
namespace {

// Published specific gamma-ray constants are quoted in R*cm^2/(h*mCi); NuSIFT computes
// R*m^2/(h*Bq). 1 m^2 = 1e4 cm^2 and 1 mCi = 3.7e7 Bq.
constexpr double kToPublishedUnits = 1.0e4 * 3.7e7;

LineSpectrum spanOf(const std::vector<GammaLine>& lines) {
  return LineSpectrum(lines.data(), lines.size());
}

// --- the air tables --------------------------------------------------------

// Exact tabulated NIST values at their own grid points. If interpolation ever shifts these,
// every exposure number moves with them.
TEST(AirCoefficients, ReproducesNistTableNodesExactly) {
  EXPECT_NEAR(airMassEnergyAbsorption(1.0e6), 0.002789, 1e-9);
  EXPECT_NEAR(airMassEnergyAbsorption(6.0e5), 0.002953, 1e-9);
  EXPECT_NEAR(airMassAttenuation(1.0e6), 0.006358, 1e-9);
  EXPECT_NEAR(airMassAttenuation(1.0e4), 0.5120, 1e-9);
}

// Attenuation falls monotonically across the whole tabulated range; energy absorption does
// not -- it dips near 100 keV where Compton scattering takes over from the photoelectric
// effect, then rises again. Asserting the dip exists guards against a table entered in the
// wrong column, which monotonicity alone would not catch.
TEST(AirCoefficients, EnergyAbsorptionHasItsComptonMinimum) {
  const double at30keV = airMassEnergyAbsorption(3.0e4);
  const double at100keV = airMassEnergyAbsorption(1.0e5);
  const double at500keV = airMassEnergyAbsorption(5.0e5);
  EXPECT_LT(at100keV, at30keV) << "should fall from 30 keV toward the minimum";
  EXPECT_LT(at100keV, at500keV) << "should rise again above the minimum";
}

TEST(AirCoefficients, AttenuationFallsMonotonically) {
  double previous = airMassAttenuation(1.0e4);
  for (const double energy : {1.5e4, 3.0e4, 1.0e5, 5.0e5, 1.0e6, 5.0e6, 1.0e7}) {
    const double current = airMassAttenuation(energy);
    EXPECT_LT(current, previous) << "at " << energy << " eV";
    previous = current;
  }
}

// Outside the table the coefficients clamp rather than extrapolate. Below 10 keV a power-law
// extrapolation of the photoelectric rise is badly wrong, and clamping is the honest failure.
TEST(AirCoefficients, ClampsOutsideTheTabulatedRange) {
  EXPECT_DOUBLE_EQ(airMassEnergyAbsorption(1.0e3), airMassEnergyAbsorption(kMinTabulatedEv));
  EXPECT_DOUBLE_EQ(airMassAttenuation(1.0e8), airMassAttenuation(kMaxTabulatedEv));
  EXPECT_TRUE(isOutsideTabulatedRange(1.0e3));
  EXPECT_TRUE(isOutsideTabulatedRange(1.0e8));
  EXPECT_FALSE(isOutsideTabulatedRange(6.6e5));
}

// --- the physics gate ------------------------------------------------------
//
// Two nuclides whose gamma constants are published everywhere. These are the tests that say
// the exposure model is right rather than merely self-consistent, and they are the reason the
// NIST tables are entered at NIST's own grid.

// Co-60: 1.173 and 1.333 MeV, both at essentially unit intensity. Published Gamma is
// 13.2 R*cm^2/(h*mCi).
TEST(PointSource, Cobalt60GammaConstantMatchesPublishedValue) {
  const std::vector<GammaLine> lines = {
      {1173228.0, 0.9985, SpectrumType::Gamma},
      {1332492.0, 0.9998, SpectrumType::Gamma},
  };
  const double published = gammaConstant(spanOf(lines)) * kToPublishedUnits;

  // Comes out near 12.9 against a published 13.2. The residual is the model's, not an error:
  // published constants fold in a small scatter contribution that an uncollided calculation
  // omits, and differ slightly by which mu_en/rho evaluation they used.
  EXPECT_NEAR(published, 13.2, 13.2 * 0.05) << "computed " << published;
}

// Ba-137m carries the 661.657 keV line that everyone attributes to Cs-137. Per Ba-137m decay
// its intensity is 0.898, giving a gamma constant near 3.38.
TEST(PointSource, Barium137mGammaConstantMatchesPublishedValue) {
  const std::vector<GammaLine> lines = {{661657.0, 0.8994, SpectrumType::Gamma}};
  const double published = gammaConstant(spanOf(lines)) * kToPublishedUnits;
  EXPECT_NEAR(published, 3.38, 3.38 * 0.06) << "computed " << published;
}

// The reason attaching lines to the emitting nuclide is right rather than merely tidy.
//
// Published tables quote Gamma = 3.3 for "Cs-137", but Cs-137 emits almost no photons -- the
// 662 keV line comes from its Ba-137m daughter, populated in 94.7% of decays. NuSIFT puts the
// line on Ba-137m, where it physically belongs, and the published value falls out of the
// equilibrium activity ratio without anyone having to fold the branching into a constant.
TEST(PointSource, Caesium137SystemReproducesItsPublishedConstantThroughEquilibrium) {
  const std::vector<GammaLine> barium = {{661657.0, 0.8994, SpectrumType::Gamma}};
  const double gammaBa137m = gammaConstant(spanOf(barium)) * kToPublishedUnits;

  // In secular equilibrium A(Ba-137m) = 0.947 * A(Cs-137), so exposure per unit Cs-137
  // activity is the branching ratio times Ba-137m's constant.
  constexpr double kBranchingToIsomer = 0.947;
  const double perCaesiumActivity = kBranchingToIsomer * gammaBa137m;

  EXPECT_NEAR(perCaesiumActivity, 3.3, 3.3 * 0.06)
      << "computed " << perCaesiumActivity << " from Gamma(Ba-137m) = " << gammaBa137m;
}

// --- geometry --------------------------------------------------------------

// With attenuation off the model must be exactly inverse-square. Any error in the geometric
// term shows up here with nothing else to hide behind.
TEST(PointSource, IsExactlyInverseSquareInVacuum) {
  const std::vector<GammaLine> lines = {{661657.0, 0.9, SpectrumType::Gamma}};
  PointSourceGeometry geometry;
  geometry.airAttenuation = false;

  geometry.distanceM = 1.0;
  const double atOne = exposureRate(spanOf(lines), 1.0e9, geometry);
  geometry.distanceM = 2.0;
  const double atTwo = exposureRate(spanOf(lines), 1.0e9, geometry);
  geometry.distanceM = 10.0;
  const double atTen = exposureRate(spanOf(lines), 1.0e9, geometry);

  EXPECT_NEAR(atTwo, atOne / 4.0, atOne * 1e-14);
  EXPECT_NEAR(atTen, atOne / 100.0, atOne * 1e-14);
}

// The vacuum gamma constant is what the exposure rate reduces to at 1 m, which is what makes
// it comparable against a published table at all.
TEST(PointSource, GammaConstantIsTheVacuumRateAtOneMetre) {
  const std::vector<GammaLine> lines = {{1173228.0, 0.9985, SpectrumType::Gamma},
                                        {1332492.0, 0.9998, SpectrumType::Gamma}};
  PointSourceGeometry geometry;
  geometry.airAttenuation = false;
  geometry.distanceM = 1.0;

  EXPECT_NEAR(exposureRate(spanOf(lines), 1.0, geometry), gammaConstant(spanOf(lines)),
              gammaConstant(spanOf(lines)) * 1e-12);
}

// Air attenuation must reduce the rate, and more so with distance and density -- but only
// slightly at a metre, which is why published constants ignore it.
TEST(PointSource, AirAttenuationReducesTheRateAndScalesWithPathAndDensity) {
  const std::vector<GammaLine> lines = {{661657.0, 0.9, SpectrumType::Gamma}};
  PointSourceGeometry vacuum;
  vacuum.airAttenuation = false;
  PointSourceGeometry air;

  EXPECT_LT(exposureRate(spanOf(lines), 1.0e9, air), exposureRate(spanOf(lines), 1.0e9, vacuum));

  // At 1 m of air the correction is well under a percent, so a vacuum constant is a good
  // approximation there -- and at 100 m it is not, which is the whole point of modelling it.
  const double ratioAtOneMetre =
      exposureRate(spanOf(lines), 1.0e9, air) / exposureRate(spanOf(lines), 1.0e9, vacuum);
  EXPECT_GT(ratioAtOneMetre, 0.99);

  PointSourceGeometry farAir;
  farAir.distanceM = 100.0;
  PointSourceGeometry farVacuum = farAir;
  farVacuum.airAttenuation = false;
  const double ratioAtHundredMetres =
      exposureRate(spanOf(lines), 1.0e9, farAir) / exposureRate(spanOf(lines), 1.0e9, farVacuum);
  EXPECT_LT(ratioAtHundredMetres, 0.9) << "100 m of air should attenuate appreciably";

  PointSourceGeometry dense = air;
  dense.airDensityKgM3 = 2.0 * air.airDensityKgM3;
  EXPECT_LT(exposureRate(spanOf(lines), 1.0e9, dense), exposureRate(spanOf(lines), 1.0e9, air));
}

// THE structural claim: because mu_air depends on energy, attenuation cannot be factored out
// of the sum over lines. If it could, a per-nuclide constant would suffice and the store would
// not need to persist spectra at all.
//
// Two spectra with the same total photon energy but different hardness must diverge with
// distance -- identical in vacuum, different through air.
TEST(PointSource, SpectralHardnessChangesTheDistanceDependence) {
  // One 2 MeV photon versus four 500 keV photons: equal energy, very different penetration.
  const std::vector<GammaLine> hard = {{2.0e6, 1.0, SpectrumType::Gamma}};
  const std::vector<GammaLine> soft = {{5.0e5, 4.0, SpectrumType::Gamma}};

  PointSourceGeometry near;
  near.distanceM = 1.0;
  PointSourceGeometry far;
  far.distanceM = 200.0;

  const double hardNear = exposureRate(spanOf(hard), 1.0e9, near);
  const double softNear = exposureRate(spanOf(soft), 1.0e9, near);
  const double hardFar = exposureRate(spanOf(hard), 1.0e9, far);
  const double softFar = exposureRate(spanOf(soft), 1.0e9, far);

  // The harder spectrum loses proportionally less over the longer path, so the ratio between
  // the two shifts. A single distance-independent constant per nuclide cannot express this.
  const double nearRatio = hardNear / softNear;
  const double farRatio = hardFar / softFar;
  EXPECT_GT(farRatio, nearRatio * 1.05) << "near " << nearRatio << ", far " << farRatio
                                        << " -- attenuation must sit inside the sum over lines";
}

// --- linearity and composition ---------------------------------------------

TEST(PointSource, IsLinearInActivity) {
  const std::vector<GammaLine> lines = {{661657.0, 0.9, SpectrumType::Gamma}};
  const PointSourceGeometry geometry;
  const double single = exposureRate(spanOf(lines), 1.0e9, geometry);
  EXPECT_NEAR(exposureRate(spanOf(lines), 3.0e9, geometry), 3.0 * single, single * 1e-12);
}

TEST(PointSource, IsAdditiveOverLines) {
  const std::vector<GammaLine> first = {{661657.0, 0.9, SpectrumType::Gamma}};
  const std::vector<GammaLine> second = {{1173228.0, 0.5, SpectrumType::Gamma}};
  const std::vector<GammaLine> both = {{661657.0, 0.9, SpectrumType::Gamma},
                                       {1173228.0, 0.5, SpectrumType::Gamma}};
  const PointSourceGeometry geometry;
  EXPECT_NEAR(
      exposureRate(spanOf(both), 1.0e9, geometry),
      exposureRate(spanOf(first), 1.0e9, geometry) + exposureRate(spanOf(second), 1.0e9, geometry),
      exposureRate(spanOf(both), 1.0e9, geometry) * 1e-12);
}

TEST(PointSource, BuildupIsACleanMultiplicativeScale) {
  const std::vector<GammaLine> lines = {{661657.0, 0.9, SpectrumType::Gamma}};
  PointSourceGeometry plain;
  PointSourceGeometry scaled;
  scaled.buildup = 2.5;
  EXPECT_NEAR(exposureRate(spanOf(lines), 1.0e9, scaled),
              2.5 * exposureRate(spanOf(lines), 1.0e9, plain),
              exposureRate(spanOf(lines), 1.0e9, plain) * 1e-12);
}

TEST(PointSource, ANuclideWithNoLinesHasNoExposure) {
  const std::vector<GammaLine> none;
  const PointSourceGeometry geometry;
  EXPECT_DOUBLE_EQ(exposureRate(spanOf(none), 1.0e15, geometry), 0.0);
  EXPECT_DOUBLE_EQ(gammaConstant(spanOf(none)), 0.0);
}

// --- units and guards ------------------------------------------------------

TEST(PointSource, RoentgenConversions) {
  EXPECT_DOUBLE_EQ(roentgenToGray(1.0), 0.00876);
  EXPECT_DOUBLE_EQ(roentgenToSievert(1.0), 0.00876);
  // A photon radiation weighting factor of 1 makes the two numerically identical; they are
  // named separately because absorbed dose and equivalent dose are different quantities.
  EXPECT_DOUBLE_EQ(roentgenToGray(2.5), roentgenToSievert(2.5));
}

// A point source has no exposure rate at zero distance. Reporting an infinity would put a
// meaningless number into a report; refusing says what went wrong.
TEST(PointSource, RejectsAnImpossibleGeometry) {
  const std::vector<GammaLine> lines = {{661657.0, 0.9, SpectrumType::Gamma}};
  PointSourceGeometry zero;
  zero.distanceM = 0.0;
  EXPECT_THROW(exposureRate(spanOf(lines), 1.0e9, zero), InputError);

  PointSourceGeometry negative;
  negative.distanceM = -1.0;
  EXPECT_THROW(exposureRate(spanOf(lines), 1.0e9, negative), InputError);

  PointSourceGeometry badBuildup;
  badBuildup.buildup = 0.0;
  EXPECT_THROW(exposureRate(spanOf(lines), 1.0e9, badBuildup), InputError);
}

}  // namespace
}  // namespace nusift::exposure
