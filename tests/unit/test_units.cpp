#include <gtest/gtest.h>

#include <cmath>

#include "nusift/units.hpp"

namespace nusift::units {
namespace {

// Half-life -> decay constant. Cs-137 at 30.08 y is the reference case that appears in
// every downstream activity number, so it is pinned against a hand-computed value.
TEST(Units, DecayConstantFromHalfLife) {
  const double cs137HalfLife = 30.08 * kSecondsPerYear;
  const double lambda = decayConstant(cs137HalfLife);
  EXPECT_NEAR(lambda, 7.3e-10, 0.1e-10);
  // The defining relation, to machine precision.
  EXPECT_NEAR(std::exp(-lambda * cs137HalfLife), 0.5, 1e-15);
}

// Stable nuclides are encoded as a non-positive half-life throughout the store and the
// chain. They must yield lambda = 0 -- a terminator with no removal rate -- rather than an
// infinity that would poison the decay matrix.
TEST(Units, StableNuclideHasZeroDecayConstant) {
  EXPECT_EQ(decayConstant(0.0), 0.0);
  EXPECT_EQ(decayConstant(-1.0), 0.0);
  EXPECT_TRUE(std::isfinite(decayConstant(0.0)));
}

// ENDF gives the atomic weight ratio, not a molar mass. Getting this conversion wrong is
// the "use A as the molar mass" error the staged AWR field exists to prevent, so the test
// asserts the discrepancy is real and in the direction expected.
TEST(Units, MolarMassFromAwr) {
  // Cs-137: AWR = 135.8351 -> 137.02 g/mol.
  const double molar = molarMassFromAwr(135.8351);
  EXPECT_NEAR(molar, 137.02, 0.05);
  // Materially different from naively using A = 137, which is the whole point.
  EXPECT_GT(std::abs(molar - 137.0), 0.005);
}

TEST(Units, TimeConversions) {
  EXPECT_EQ(kSecondsPerMinute, 60.0);
  EXPECT_EQ(kSecondsPerHour, 3600.0);
  EXPECT_EQ(kSecondsPerDay, 86400.0);
  // The Julian year, 365.25 d -- stated in --help because the choice is arbitrary but its
  // consequences over a 100 y decay are not.
  EXPECT_NEAR(kSecondsPerYear, 31557600.0, 1e-9);
  EXPECT_NEAR(kSecondsPerYear / kSecondsPerDay, 365.25, 1e-12);
}

TEST(Units, ActivityConversion) {
  EXPECT_EQ(kBqPerCi, 3.7e10);  // exact by definition
}

// Exposure -> absorbed dose in air. NuSIFT reports one physical quantity and converts with
// a photon radiation weighting factor of 1, so Gy/R and Sv/R are numerically the same
// constant. Pinned so a future "improvement" to an ICRP-74 H*(10) coefficient has to be a
// deliberate, visible change.
TEST(Units, ExposureToDoseConversion) {
  EXPECT_EQ(kGyPerR, 0.00876);
  EXPECT_EQ(kRadiationWeightingPhoton, 1.0);
  const double exposureRPerHour = 1.0;
  EXPECT_NEAR(exposureRPerHour * kGyPerR * kRadiationWeightingPhoton, 0.00876, 1e-15);
}

TEST(Units, CrossSectionConversion) {
  EXPECT_EQ(kBarnToCm2, 1.0e-24);
  // A 37.18 barn Co-59 (n,gamma) at 1e14 n/cm^2/s gives a reaction rate near 3.7e-9 /s.
  const double rate = 37.18 * kBarnToCm2 * 1.0e14;
  EXPECT_NEAR(rate, 3.718e-9, 1e-12);
}

TEST(Units, FundamentalConstants) {
  EXPECT_EQ(kAvogadro, 6.02214076e23);  // exact, SI 2019
  EXPECT_EQ(kEvToJ, 1.602176634e-19);   // exact, SI 2019
  EXPECT_NEAR(kLn2, 0.6931471805599453, 1e-15);
}

}  // namespace
}  // namespace nusift::units
