#include <gtest/gtest.h>

#include <cmath>
#include <string>

#include "nusift/core/nuclide_name.hpp"
#include "nusift/exposure/point_source.hpp"
#include "nusift/nucdata/nuclear_data.hpp"
#include "nusift/units.hpp"
#include "validation_store.hpp"

namespace nusift::validation {
namespace {

// Published specific gamma-ray constants are quoted in R*cm^2/(h*mCi); NuSIFT computes
// R*m^2/(h*Bq). 1 m^2 = 1e4 cm^2 and 1 mCi = 3.7e7 Bq.
constexpr double kToPublishedUnits = 1.0e4 * 3.7e7;

int indexOfName(const NuclearData& data, const char* name) {
  return data.indexOf(requireNuclideName(name));
}

double gammaConstantOf(const NuclearData& data, const char* name) {
  const int i = indexOfName(data, name);
  return i >= 0 ? exposure::gammaConstant(data.lines(i)) * kToPublishedUnits : 0.0;
}

// --- published gamma constants, through the shipped store -------------------
//
// test_point_source.cpp already checks these against literal hand-written spectra, which is
// what isolates the exposure model from the data. This file asks the other half of the
// question: do the spectra NuSIFT actually ships reproduce the same published values. A
// staging bug that dropped a line, halved an intensity, or attached a spectrum to the wrong
// nuclide passes there and fails here.
//
// The reference is Ninkovic & Adrovic's recalculation -- 309.0 uGy*m^2/(GBq*h) for Co-60,
// which is 13.05 R*cm^2/(h*mCi) in the modern roentgen. See docs/exposure.md section 4 for why
// published tabulations disagree with each other by more than any of them disagrees with
// NuSIFT, and for the residual's decomposition. The band here is the same 3% the unit suite
// uses, for the same reason: it is wide enough to contain the convention spread and far too
// narrow to contain a data error.

TEST(PublishedConstants, Cobalt60FromTheStoreMatchesItsPublishedConstant) {
  const double computed = gammaConstantOf(committedStore(), "Co-60");
  EXPECT_NEAR(computed, 13.05, 13.05 * 0.03) << "computed " << computed;
}

// Ba-137m carries the 661.657 keV line that everyone attributes to Cs-137, and the store also
// carries the Ba K X-rays near 32 keV. The whole-spectrum constant is therefore about 3.47
// against the ~3.38 the gamma alone gives, which is exactly the difference between the two
// families of published Cs-137 values (~3.2 gamma-only, ~3.3 with X-rays). The band is 6%
// because the reference itself is quoted to two significant figures.
TEST(PublishedConstants, Barium137mFromTheStoreMatchesItsPublishedConstant) {
  const double computed = gammaConstantOf(committedStore(), "Ba-137m");
  EXPECT_NEAR(computed, 3.47, 3.47 * 0.06) << "computed " << computed;
}

// The reason attributing photons to their actual emitter is right rather than merely tidy.
// Cs-137 itself emits almost nothing -- the store gives it a single 283.5 keV line at 5.8e-6
// intensity -- so its own constant is four orders of magnitude below the tabulated "Cs-137"
// value, and that is not an error. The published number is the Ba-137m constant times the
// 94.7% branch, and it falls out of the equilibrium ratio without anyone folding a branching
// into a table.
TEST(PublishedConstants, Caesium137IsAlmostAllDaughter) {
  const double parent = gammaConstantOf(committedStore(), "Cs-137");
  const double daughter = gammaConstantOf(committedStore(), "Ba-137m");
  EXPECT_LT(parent, daughter * 1.0e-3) << "Cs-137 alone computes " << parent;

  constexpr double kBranchingToIsomer = 0.947;
  const double system = kBranchingToIsomer * daughter;
  EXPECT_NEAR(system, 3.3, 3.3 * 0.06)
      << "computed " << system << " from Gamma(Ba-137m) = " << daughter;
}

// --- staged decay data against evaluated values -----------------------------
//
// Half-lives are the single most consequential number in the store: every decay constant, and
// therefore every activity and every exposure, is ln(2) over one of these. They are checked
// against the evaluated values quoted in the standard compilations at 0.5%, which the staged
// data clears by three orders of magnitude -- ENDF/B-VIII.1 adopts the same evaluations. The
// point is not the margin, it is that a unit slip or a tape misread cannot hide.

struct HalfLifeCase {
  const char* name;
  double expected;
  double unitSeconds;
  const char* unit;
};

TEST(PublishedConstants, StagedHalfLivesMatchTheEvaluatedValues) {
  constexpr double kMinute = 60.0;
  constexpr double kHour = 3600.0;
  constexpr double kDay = 86400.0;
  const double year = units::kSecondsPerYear;

  const HalfLifeCase cases[] = {
      {"Cs-137", 30.08, year, "y"},   {"Co-60", 5.2711, year, "y"},
      {"Sr-90", 28.79, year, "y"},    {"H-3", 12.32, year, "y"},
      {"Am-241", 432.6, year, "y"},   {"I-131", 8.0252, kDay, "d"},
      {"Ir-192", 73.827, kDay, "d"},  {"Mo-99", 65.976, kHour, "h"},
      {"Tc-99m", 6.0067, kHour, "h"}, {"Ba-137m", 2.552, kMinute, "min"},
  };

  for (const HalfLifeCase& c : cases) {
    const int i = indexOfName(committedStore(), c.name);
    ASSERT_GE(i, 0) << c.name << " is not in the store";
    const double staged = committedStore().halfLifeSeconds(i) / c.unitSeconds;
    EXPECT_NEAR(staged, c.expected, c.expected * 0.005)
        << c.name << ": staged " << staged << " " << c.unit << ", evaluated " << c.expected;
  }
}

// Specific activity is ln(2) N_A / (T_half M), and it is what a user reads off a report, so it
// is checked as the quotient here even though the validation sweeps check both of its inputs
// separately against ENSDF and AME2020.
//
// Published tables are themselves this derivation, which is why the band is 2% rather than
// tighter: a table computed from a superseded half-life differs from this one by exactly that
// revision, and the tables disagree with each other by more. Sr-90 is the visible case,
// tabulated at 136 Ci/g from the older 29.1 y against the 28.79 y staged here.
TEST(PublishedConstants, SpecificActivitiesMatchPublishedTables) {
  struct Case {
    const char* name;
    double curiePerGram;
  };
  const Case cases[] = {
      {"H-3", 9650.0},  {"Co-60", 1131.0}, {"Cs-137", 87.0},
      {"Sr-90", 136.0}, {"Am-241", 3.43},  {"Ra-226", 0.9886},
  };

  for (const Case& c : cases) {
    const int i = indexOfName(committedStore(), c.name);
    ASSERT_GE(i, 0) << c.name << " is not in the store";
    const double molarMass = committedStore().molarMassGPerMol(i);
    ASSERT_GT(molarMass, 0.0) << c.name << " has no staged atomic weight";
    const double perGram =
        committedStore().decayConstant(i) * units::kAvogadro / molarMass / units::kBqPerCi;
    EXPECT_NEAR(perGram, c.curiePerGram, c.curiePerGram * 0.02)
        << c.name << ": computed " << perGram << " Ci/g, tabulated " << c.curiePerGram;
  }
}

}  // namespace
}  // namespace nusift::validation
