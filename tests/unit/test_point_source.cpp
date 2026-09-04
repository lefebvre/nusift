#include <gtest/gtest.h>

#include <cmath>
#include <numbers>
#include <vector>

#include "nusift/core/error.hpp"
#include "nusift/exposure/air_coefficients.hpp"
#include "nusift/exposure/dose_coefficients.hpp"
#include "nusift/exposure/point_source.hpp"
#include "nusift/nucdata/photon_lines.hpp"
#include "nusift/units.hpp"

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
//
// WHICH published value, though. Tabulated constants disagree with each other by more than
// any of them disagrees with NuSIFT, and the spread is convention rather than physics:
//
//   * The roentgen depends on W/e, revised from 33.7 to 33.85 J/C (ICRU 1979) and then to
//     33.97. Identical physics tabulated before that revision reads about 0.8% higher.
//   * An air KERMA rate constant uses the mass energy-TRANSFER coefficient, i.e. total kerma.
//     Exposure needs mass energy-ABSORPTION, i.e. collision kerma. Worth 0.3% at 1.25 MeV.
//   * Tabulations differ on which photons they include and on whose mu/rho evaluation.
//
// The reference used below is Ninkovic & Adrovic's recalculation, which exists precisely
// because "published data are in strong disagreement". Their 309.0 uGy*m^2/(GBq*h) for Co-60
// is 13.05 R*cm^2/(h*mCi) in the modern roentgen and 13.15 in the pre-1979 one -- the latter
// being the classic 13.2 that older tables quote. Same physics, different decade.
//
// None of these are scatter. Published constants are vacuum quantities by definition: a point
// source in a vacuum, no self-attenuation, no air scatter.

// Co-60: 1.173 and 1.333 MeV, both at essentially unit intensity.
//
// Nothing else in the spectrum matters. Gamma is 194 * sum(y E mu_en) to three digits, and it
// is pinned there: forcing both intensities to exactly 1.0 moves it +0.07%, collapsing both
// lines onto the 1.25 MeV NIST grid point (removing interpolation entirely) moves it -0.00%,
// and interpolating linearly instead of log-log moves it +0.12%. No plausible change to the
// decay data or to the interpolation reaches a quarter of a percent, so any disagreement with
// a published value lives in the conversion convention above rather than in this file.
TEST(PointSource, Cobalt60GammaConstantMatchesPublishedValue) {
  const std::vector<GammaLine> lines = {
      {1173228.0, 0.9985, SpectrumType::Gamma},
      {1332492.0, 0.9998, SpectrumType::Gamma},
  };
  const double published = gammaConstant(spanOf(lines)) * kToPublishedUnits;

  // 12.91 against 13.05. Of that gap, 0.3% is the transfer/absorption difference and the rest
  // is the air-coefficient evaluation -- and Ba-137m carries the same 0.9% offset, so it is
  // one uniform bias in the air table rather than anything specific to Co-60.
  constexpr double kModernRecalculation = 13.05;
  EXPECT_NEAR(published, kModernRecalculation, kModernRecalculation * 0.03)
      << "computed " << published;
}

// Ba-137m carries the 661.657 keV line that everyone attributes to Cs-137. Per Ba-137m decay
// its intensity is 0.899, giving 3.38 from the gamma alone.
//
// The shipped store also carries the Ba K X-rays near 32 keV, which lift the whole-spectrum
// constant to 3.47. That is why published Cs-137 values fall into a ~3.2 family (gamma only)
// and a ~3.3 family (X-rays included). This test deliberately isolates the gamma so that the
// number means one thing and does not move when the X-ray intensities are re-evaluated.
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
//
// This is an identity rather than an independent check, and deliberately so: a tabulated
// "Cs-137" constant IS the Ba-137m constant times the branch. Unger & Trubey say as much --
// their Cs-137 entry is the product of 94.6% and their computed Ba-137m value, added as a
// convenience, with a warning not to double-count it against a data set that also carries
// Ba-137m activity. Attributing each line to its emitter makes that error unrepresentable.
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
// slightly at a meter, which is why published constants ignore it.
TEST(PointSource, AirAttenuationReducesTheRateAndScalesWithPathAndDensity) {
  const std::vector<GammaLine> lines = {{661657.0, 0.9, SpectrumType::Gamma}};
  PointSourceGeometry vacuum;
  vacuum.airAttenuation = false;
  PointSourceGeometry air;

  EXPECT_LT(exposureRate(spanOf(lines), 1.0e9, air), exposureRate(spanOf(lines), 1.0e9, vacuum));

  // At 1 m of air the correction is well under a percent, so a vacuum constant is a good
  // approximation there -- and at 100 m it is not, which is the whole point of modeling it.
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

// --- the fluence kernel ------------------------------------------------------
//
// Fluence shares its transport with exposure: the same inverse-square spreading, the same air
// attenuation, the same buildup. The only difference is that the kerma-to-exposure conversion
// is left out. Every test in this section checks that identity directly, so a change to one
// kernel cannot move the other's numbers without being seen.

// Exposure is fluence times the kerma conversion, so the ratio of the two per-energy
// coefficients must be exactly that conversion at every energy. If the kernels ever
// diverged, every fluence number would be silently wrong in exactly the places where the
// exposure is still right.
TEST(PointSource, FluenceCoefficientIsTheExposureKernelWithoutKerma) {
  const PointSourceGeometry geometry;
  for (const double energy : {1.0e4, 1.0e5, 6.61657e5, 2.0e6, 1.0e7}) {
    const double conversion = energy * units::kEvToJ * airMassEnergyAbsorption(energy) *
                              units::kSecondsPerHour / units::kGyPerR;
    const double exposure = pointExposureCoeff(energy, geometry);
    EXPECT_NEAR(pointFluenceCoeff(energy, geometry) * conversion, exposure, exposure * 1e-12)
        << "at " << energy << " eV";
  }
}

// What a person can check on paper: the per-becquerel fluence rate is the sum over the
// spectrum of (1 / 4 pi r^2) times the exponential attenuation times the buildup, weighted
// by intensity -- and the rate scales linearly with the activity.
TEST(PointSource, FluenceRateMatchesSpreadingAttenuationAndBuildup) {
  const std::vector<GammaLine> lines = {{661657.0, 0.9, SpectrumType::Gamma},
                                        {1.0e6, 0.1, SpectrumType::Gamma}};
  PointSourceGeometry geometry;
  geometry.distanceM = 3.0;
  geometry.buildup = 2.5;

  const double d = geometry.distanceM;
  double expected = 0.0;
  for (const auto& line : lines) {
    expected += line.intensity / (4.0 * std::numbers::pi * d * d) *
                std::exp(-airMassAttenuation(line.energyEv) * geometry.airDensityKgM3 * d) *
                geometry.buildup;
  }

  EXPECT_DOUBLE_EQ(fluenceRatePerBecquerel(spanOf(lines), geometry), expected);
  // Two identical emitters double the fluence: the coefficient is linear in the spectrum.
  std::vector<GammaLine> doubled = lines;
  doubled.insert(doubled.end(), lines.begin(), lines.end());
  EXPECT_DOUBLE_EQ(fluenceRatePerBecquerel(spanOf(doubled), geometry), 2.0 * expected);
}

// With attenuation off, a unit-intensity line at 1 m is exactly 1 / (4 pi) photons per
// square meter per second per becquerel, and the inverse-square scaling is exact: there is
// nowhere for an error in the geometric term to hide.
TEST(PointSource, FluenceIsExactlyInverseSquareInVacuum) {
  const std::vector<GammaLine> lines = {{661657.0, 0.9, SpectrumType::Gamma}};
  const std::vector<GammaLine> unity = {{661657.0, 1.0, SpectrumType::Gamma}};
  PointSourceGeometry geometry;
  geometry.airAttenuation = false;

  geometry.distanceM = 1.0;
  const double atOne = fluenceRatePerBecquerel(spanOf(lines), geometry);
  geometry.distanceM = 2.0;
  const double atTwo = fluenceRatePerBecquerel(spanOf(lines), geometry);
  geometry.distanceM = 10.0;
  const double atTen = fluenceRatePerBecquerel(spanOf(lines), geometry);

  EXPECT_NEAR(atTwo, atOne / 4.0, atOne * 1e-14);
  EXPECT_NEAR(atTen, atOne / 100.0, atOne * 1e-14);
  geometry.distanceM = 1.0;
  EXPECT_DOUBLE_EQ(fluenceRatePerBecquerel(spanOf(unity), geometry),
                   1.0 / (4.0 * std::numbers::pi));
}

// The guards the exposure path has must hold on the fluence path as well: a zero distance
// would put an infinity into a report, and a zero buildup a number that looks meaningful.
TEST(PointSource, FluenceRefusesAnImpossibleGeometry) {
  const std::vector<GammaLine> lines = {{661657.0, 0.9, SpectrumType::Gamma}};
  PointSourceGeometry atSource;
  atSource.distanceM = 0.0;
  EXPECT_THROW(pointFluenceCoeff(661657.0, atSource), InputError);
  EXPECT_THROW(fluenceRatePerBecquerel(spanOf(lines), atSource), InputError);

  PointSourceGeometry noBuildup;
  noBuildup.buildup = 0.0;
  EXPECT_THROW(pointFluenceCoeff(661657.0, noBuildup), InputError);
  EXPECT_THROW(fluenceRatePerBecquerel(spanOf(lines), noBuildup), InputError);
}

TEST(PointSource, ANuclideWithNoLinesHasNoExposure) {
  const std::vector<GammaLine> none;
  const PointSourceGeometry geometry;
  EXPECT_DOUBLE_EQ(exposureRate(spanOf(none), 1.0e15, geometry), 0.0);
  EXPECT_DOUBLE_EQ(gammaConstant(spanOf(none)), 0.0);
  EXPECT_DOUBLE_EQ(meanOpticalDepth(spanOf(none), geometry), 0.0);
  EXPECT_DOUBLE_EQ(fluenceRatePerBecquerel(spanOf(none), geometry), 0.0);
}

// --- the air path ----------------------------------------------------------

// A single line's optical depth is mu(E) rho d, the quantity buildup is tabulated against; a
// spectrum's is the mean over its lines weighted by the exposure each delivers, so it sits
// between the lines' own depths and nearer the one that carries the exposure.
TEST(PointSource, MeanOpticalDepthIsExposureWeighted) {
  const std::vector<GammaLine> soft = {{1.0e5, 1.0, SpectrumType::Gamma}};
  const std::vector<GammaLine> hard = {{2.0e6, 1.0, SpectrumType::Gamma}};
  const std::vector<GammaLine> both = {{1.0e5, 1.0, SpectrumType::Gamma},
                                       {2.0e6, 1.0, SpectrumType::Gamma}};
  PointSourceGeometry geometry;
  geometry.distanceM = 50.0;

  const double softDepth = meanOpticalDepth(spanOf(soft), geometry);
  const double hardDepth = meanOpticalDepth(spanOf(hard), geometry);
  const double expectedSoft = airMassAttenuation(1.0e5) * geometry.airDensityKgM3 * 50.0;
  EXPECT_NEAR(softDepth, expectedSoft, expectedSoft * 1e-12);
  EXPECT_GT(softDepth, hardDepth) << "air is thicker to a 100 keV photon than to a 2 MeV one";

  // At 50 m the 2 MeV photon delivers many times the exposure of the 100 keV one, so the mean
  // sits much nearer its depth than the soft line's.
  const double mixed = meanOpticalDepth(spanOf(both), geometry);
  EXPECT_GT(mixed, hardDepth);
  EXPECT_LT(mixed, softDepth);
  EXPECT_LT(mixed - hardDepth, softDepth - mixed);

  // No attenuation means no path to be thick, whatever the distance.
  PointSourceGeometry vacuum = geometry;
  vacuum.airAttenuation = false;
  EXPECT_DOUBLE_EQ(meanOpticalDepth(spanOf(both), vacuum), 0.0);
}

// --- units and guards ------------------------------------------------------

TEST(PointSource, RoentgenConversions) {
  EXPECT_DOUBLE_EQ(roentgenToGray(1.0), 0.00876);
  EXPECT_DOUBLE_EQ(roentgenToGray(2.5), 2.5 * 0.00876);
  // There is no roentgenToSievert to test beside it. A sievert is effective dose and comes
  // from the ICRP 116 kernel below; it is not air kerma multiplied by anything.
}

// --- effective dose --------------------------------------------------------

// The cross-check that matters, and it is against ICRP's OWN second table rather than against
// a number this code produced. Publication 116 gives effective dose per fluence (Table A.1,
// which is what NuSIFT ships) and effective dose per air kerma free-in-air (Table A.2, which
// it does not). Dividing our dose kernel by our air-kerma kernel has to reproduce A.2 -- and
// the two paths share nothing but the geometry: one goes through ICRP's phantom coefficients,
// the other through NIST's mass energy-absorption coefficients for air.
//
// Agreement to a fraction of a percent across the range decay photons occupy therefore checks
// the transcription, the interpolation, the unit conversion and the kernel assembly at once.
TEST(PointSource, EffectiveDoseReproducesIcrpsOwnRatioToAirKerma) {
  struct Point {
    double mev;
    double svPerGy;  // ICRP 116 Table A.2, AP column
  };
  constexpr Point kPoints[] = {{0.02, 0.130}, {0.03, 0.423}, {0.05, 1.13}, {0.06, 1.33},
                               {0.1, 1.39},   {0.2, 1.17},   {0.5, 1.04},  {0.662, 1.02},
                               {1.0, 1.00},   {1.33, 0.996}, {2.0, 0.990}};

  PointSourceGeometry vacuum;
  vacuum.airAttenuation = false;

  for (const Point& point : kPoints) {
    const double energyEv = point.mev * 1.0e6;
    const double sievert = pointEffectiveDoseCoeff(energyEv, vacuum);
    const double gray = pointExposureCoeff(energyEv, vacuum) * units::kGyPerR;
    EXPECT_NEAR(sievert / gray, point.svPerGy, point.svPerGy * 0.01)
        << "at " << point.mev << " MeV, against ICRP 116 Table A.2";
  }
}

// Above a few MeV the same comparison drifts, and the drift is physics rather than a defect:
// ICRP's denominator is air KERMA -- energy transferred to charged particles, including what
// those particles later radiate away as bremsstrahlung -- while NIST's mu_en/rho counts only
// what is absorbed. The gap between them is air's radiative yield, about 1.8% at 5 MeV and
// 3.4% at 10 MeV, and our ratio sits high by very nearly that. It is asserted rather than
// merely tolerated, because a departure that stopped tracking the radiative yield would mean
// something else was wrong.
TEST(PointSource, AboveAFewMeVTheRatioDepartsByAirsRadiativeYield) {
  PointSourceGeometry vacuum;
  vacuum.airAttenuation = false;

  const auto ratioAt = [&vacuum](double mev) {
    const double energyEv = mev * 1.0e6;
    return (pointEffectiveDoseCoeff(energyEv, vacuum) /
            (pointExposureCoeff(energyEv, vacuum) * units::kGyPerR));
  };

  EXPECT_NEAR(ratioAt(5.0) / 0.943, 1.0 + 0.018, 0.01) << "air's radiative yield at 5 MeV";
  EXPECT_NEAR(ratioAt(10.0) / 0.848, 1.0 + 0.034, 0.01) << "and at 10 MeV";
}

// A tabulated energy comes back as tabulated rather than as an interpolation that lands
// nearby. ICRP put the lines this tool cares about on their own grid -- 0.511, 0.662, 1.117
// and 1.33 MeV are all grid points -- so the commonest answers involve no interpolation at
// all. Checked to round-off rather than to the bit, because a log-log interpolator returns
// exp(log(v)) at a grid point and that is not bitwise v.
TEST(PointSource, ATabulatedEnergyIsReturnedAsTabulated) {
  EXPECT_NEAR(effectiveDosePerFluence(0.662e6, Irradiation::AP), 3.17e-16, 3.17e-16 * 1.0e-12);
  EXPECT_NEAR(effectiveDosePerFluence(1.33e6, Irradiation::AP), 5.59e-16, 5.59e-16 * 1.0e-12);
  EXPECT_NEAR(effectiveDosePerFluence(0.511e6, Irradiation::ROT), 1.96e-16, 1.96e-16 * 1.0e-12);
  EXPECT_NEAR(effectiveDosePerFluence(1.0e4, Irradiation::ISO), 0.0288e-16, 0.0288e-16 * 1.0e-12);
}

// Below the table the coefficient is clamped rather than extrapolated. The curve is falling
// steeply there -- a factor of two per five keV -- so an extrapolation would be invention, and
// clamping is conservative in the direction that matters: it overstates a soft photon's dose
// rather than understating it. Am-241 is the nuclide this actually reaches.
TEST(PointSource, BelowTheTableTheCoefficientIsClampedNotExtrapolated) {
  const double atFloor = effectiveDosePerFluence(kMinTabulatedDoseEv, Irradiation::AP);
  EXPECT_DOUBLE_EQ(effectiveDosePerFluence(5.0e3, Irradiation::AP), atFloor);
  EXPECT_DOUBLE_EQ(effectiveDosePerFluence(1.0, Irradiation::AP), atFloor);
  EXPECT_TRUE(isOutsideTabulatedDoseRange(5.0e3));
  EXPECT_FALSE(isOutsideTabulatedDoseRange(kMinTabulatedDoseEv));

  const double atCeiling = effectiveDosePerFluence(kMaxTabulatedDoseEv, Irradiation::AP);
  EXPECT_DOUBLE_EQ(effectiveDosePerFluence(1.0e15, Irradiation::AP), atCeiling);
  EXPECT_TRUE(isOutsideTabulatedDoseRange(1.0e15));
}

// The irradiation geometry is not a label on the answer, it IS part of the answer. At 30 keV
// the same fluence delivers 3.3 times the effective dose to someone facing the source as to
// someone facing away, because at that energy the organs carrying the largest tissue weights
// are shielded by the whole body in PA and by nothing in AP.
TEST(PointSource, IrradiationGeometryChangesTheDoseSubstantially) {
  const double ap = effectiveDosePerFluence(3.0e4, Irradiation::AP);
  const double pa = effectiveDosePerFluence(3.0e4, Irradiation::PA);
  const double iso = effectiveDosePerFluence(3.0e4, Irradiation::ISO);
  EXPECT_NEAR(ap / pa, 0.313 / 0.0940, 1.0e-9);
  EXPECT_LT(iso, ap);
  EXPECT_GT(iso, pa);

  // And it converges at high energy, where a body is no longer able to shield itself.
  const double apHard = effectiveDosePerFluence(5.0e6, Irradiation::AP);
  const double paHard = effectiveDosePerFluence(5.0e6, Irradiation::PA);
  EXPECT_NEAR(apHard / paHard, 13.4 / 13.1, 1.0e-9);
}

TEST(PointSource, IrradiationNamesRoundTrip) {
  for (const Irradiation irradiation : {Irradiation::AP, Irradiation::PA, Irradiation::LLAT,
                                        Irradiation::RLAT, Irradiation::ROT, Irradiation::ISO}) {
    Irradiation parsed = Irradiation::PA;
    ASSERT_TRUE(parseIrradiation(irradiationName(irradiation), parsed));
    EXPECT_EQ(parsed, irradiation);
  }

  Irradiation lower = Irradiation::PA;
  EXPECT_TRUE(parseIrradiation("ap", lower));
  EXPECT_EQ(lower, Irradiation::AP);
  EXPECT_TRUE(parseIrradiation("RoT", lower));
  EXPECT_EQ(lower, Irradiation::ROT);
  EXPECT_FALSE(parseIrradiation("sideways", lower));
  EXPECT_FALSE(parseIrradiation("", lower));
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
