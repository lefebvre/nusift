#include <gtest/gtest.h>

#include <cmath>
#include <string>

#include "nusift/core/error.hpp"
#include "nusift/nucdata/nuclear_data.hpp"
#include "nusift/seed/seed_fission.hpp"
#include "synthetic_chain.hpp"

namespace nusift {
namespace {

// A chain with two fissionable parents and a handful of products, so lookup by parent and by
// energy can both be checked without a real evaluation.
NuclearData chainWithYields() {
  StoreArrays arrays = synth::linearChain({1.0e-3, 5.0e-4});
  const std::int64_t sn = Zai{50, 100, 0}.key();
  const std::int64_t sb = Zai{51, 100, 0}.key();

  // U-235-like parent, two tabulated energies; and a second parent at one energy.
  arrays.nfyParentKey = {Zai{92, 235, 0}.key(), Zai{92, 235, 0}.key(), Zai{94, 239, 0}.key()};
  arrays.nfyEnergyEv = {0.0253, 5.0e5, 0.0253};
  arrays.nfySetOffset = {0, 2, 4, 6};
  arrays.nfyProductKey = {sn, sb, sn, sb, sn, sb};
  // Thermal splits 1.2/0.8; fast splits 0.5/1.5; the second parent 1.0/1.0. Each sums to 2.0,
  // as independent yields must.
  arrays.nfyProductYield = {1.2, 0.8, 0.5, 1.5, 1.0, 1.0};
  return NuclearData::fromArrays(std::move(arrays));
}

const Zai kU235{92, 235, 0};
const Zai kPu239{94, 239, 0};

// --- energy conversions ------------------------------------------------------

// Glasstone's canonical figure. A kiloton is 1e12 calories by definition, and at 180 MeV per
// fission that is 1.45e23 fissions -- the number every weapons-effects text quotes, and the
// reason 180 rather than 200 is the default for an explosive yield.
TEST(FissionEnergy, OneKilotonIsGlasstonesFissionCount) {
  EXPECT_NEAR(seed::fissionsFromKt(1.0), 1.45e23, 1.45e23 * 0.01);
}

// The reactor convention counts the delayed beta and gamma energy that an explosive yield does
// not, so the same kiloton takes fewer fissions. An 11% difference in every downstream number,
// which is why the choice is explicit rather than assumed.
TEST(FissionEnergy, RecoverableConventionGivesFewerFissionsPerKiloton) {
  const double explosive = seed::fissionsFromKt(1.0);
  const double recoverable = seed::fissionsFromKt(1.0, seed::kMeVPerFissionRecoverable);
  EXPECT_NEAR(recoverable, 1.31e23, 1.31e23 * 0.01);
  EXPECT_LT(recoverable, explosive);
  EXPECT_NEAR(explosive / recoverable, 200.0 / 180.0, 1e-9);
}

TEST(FissionEnergy, KilotonsAndFissionsRoundTrip) {
  for (const double kt : {0.001, 1.0, 20.0, 1000.0}) {
    EXPECT_NEAR(seed::ktFromFissions(seed::fissionsFromKt(kt)), kt, kt * 1e-12);
  }
}

TEST(FissionEnergy, JoulesAgreeWithKilotons) {
  EXPECT_NEAR(seed::fissionsFromEnergyJ(seed::kJoulesPerKt), seed::fissionsFromKt(1.0),
              seed::fissionsFromKt(1.0) * 1e-12);
}

TEST(FissionEnergy, ParsesTheNamedConventions) {
  double meV = 0.0;
  EXPECT_TRUE(seed::parseMeVPerFission("explosive", meV));
  EXPECT_DOUBLE_EQ(meV, 180.0);
  EXPECT_TRUE(seed::parseMeVPerFission("recoverable", meV));
  EXPECT_DOUBLE_EQ(meV, 200.0);
  EXPECT_TRUE(seed::parseMeVPerFission("195", meV));
  EXPECT_DOUBLE_EQ(meV, 195.0);
  EXPECT_FALSE(seed::parseMeVPerFission("lots", meV));
}

TEST(IncidentEnergy, ParsesNamedEnergiesAndBareValues) {
  double energy = -1.0;
  EXPECT_TRUE(parseIncidentEnergy("thermal", energy));
  EXPECT_DOUBLE_EQ(energy, kThermalEv);
  EXPECT_TRUE(parseIncidentEnergy("fast", energy));
  EXPECT_DOUBLE_EQ(energy, kFastEv);
  EXPECT_TRUE(parseIncidentEnergy("14mev", energy));
  EXPECT_DOUBLE_EQ(energy, kFusionEv);
  EXPECT_TRUE(parseIncidentEnergy("spontaneous", energy));
  EXPECT_DOUBLE_EQ(energy, 0.0);
  EXPECT_TRUE(parseIncidentEnergy("2.53e-2", energy));
  EXPECT_NEAR(energy, 0.0253, 1e-12);
  EXPECT_FALSE(parseIncidentEnergy("warm", energy));
}

// --- seeding -----------------------------------------------------------------

// Every fission puts its yield's worth of atoms into the inventory, so the total is the
// fission count times the summed yield. Checking the total rather than individual products is
// what catches a normalisation error, which would scale everything equally.
TEST(SeedFission, TotalAtomsAreFissionsTimesTotalYield) {
  const NuclearData data = chainWithYields();
  seed::FissionSeed fissionSeed;
  fissionSeed.fissile = kU235;
  fissionSeed.incidentEnergyEv = kThermalEv;
  fissionSeed.fissions = 1.0e20;

  const Inventory inventory = seed::seedFromFission(data, fissionSeed);
  EXPECT_NEAR(inventory.totalAtoms(), 1.0e20 * 2.0, 1.0e20 * 2.0 * 1e-12);
  EXPECT_NEAR(inventory.atomsOf(Zai{50, 100, 0}), 1.2e20, 1.2e20 * 1e-12);
  EXPECT_NEAR(inventory.atomsOf(Zai{51, 100, 0}), 0.8e20, 0.8e20 * 1e-12);
}

// Yields differ by incident energy, and the lookup must pick the tabulated set nearest what
// was asked rather than defaulting to the first one.
TEST(SeedFission, SelectsTheYieldSetNearestTheRequestedEnergy) {
  const NuclearData data = chainWithYields();
  seed::FissionSeed fissionSeed;
  fissionSeed.fissile = kU235;
  fissionSeed.fissions = 1.0e20;

  fissionSeed.incidentEnergyEv = kThermalEv;
  const Inventory thermal = seed::seedFromFission(data, fissionSeed);
  EXPECT_NEAR(thermal.atomsOf(Zai{50, 100, 0}), 1.2e20, 1e8);

  fissionSeed.incidentEnergyEv = kFastEv;
  const Inventory fast = seed::seedFromFission(data, fissionSeed);
  EXPECT_NEAR(fast.atomsOf(Zai{50, 100, 0}), 0.5e20, 1e8);

  // An energy between the two snaps to whichever is closer, rather than interpolating -- the
  // evaluation has no data in between and inventing some would be worse than rounding.
  fissionSeed.incidentEnergyEv = 4.0e5;
  const Inventory between = seed::seedFromFission(data, fissionSeed);
  EXPECT_NEAR(between.atomsOf(Zai{50, 100, 0}), 0.5e20, 1e8) << "should snap to the fast set";
}

TEST(SeedFission, SelectsByParentAsWellAsEnergy) {
  const NuclearData data = chainWithYields();
  seed::FissionSeed fissionSeed;
  fissionSeed.incidentEnergyEv = kThermalEv;
  fissionSeed.fissions = 1.0e20;

  fissionSeed.fissile = kU235;
  const Inventory fromU = seed::seedFromFission(data, fissionSeed);
  fissionSeed.fissile = kPu239;
  const Inventory fromPu = seed::seedFromFission(data, fissionSeed);

  EXPECT_NEAR(fromU.atomsOf(Zai{50, 100, 0}), 1.2e20, 1e8);
  EXPECT_NEAR(fromPu.atomsOf(Zai{50, 100, 0}), 1.0e20, 1e8);
}

// The set of fissionable nuclides in an evaluation is small, so an unknown parent is answered
// with the list rather than leaving the user to guess a spelling.
TEST(SeedFission, UnknownParentNamesTheAvailableOnes) {
  const NuclearData data = chainWithYields();
  seed::FissionSeed fissionSeed;
  fissionSeed.fissile = Zai{55, 137, 0};  // Cs-137 does not fission
  fissionSeed.fissions = 1.0e20;

  try {
    seed::seedFromFission(data, fissionSeed);
    FAIL() << "expected InputError";
  } catch (const InputError& e) {
    const std::string what = e.what();
    EXPECT_NE(what.find("Cs-137"), std::string::npos) << what;
    EXPECT_NE(what.find("U-235"), std::string::npos)
        << "the message should list what the store does carry: " << what;
    EXPECT_NE(what.find("Pu-239"), std::string::npos) << what;
  }
}

TEST(SeedFission, RejectsANonPositiveFissionCount) {
  const NuclearData data = chainWithYields();
  seed::FissionSeed fissionSeed;
  fissionSeed.fissile = kU235;
  fissionSeed.fissions = 0.0;
  EXPECT_THROW(seed::seedFromFission(data, fissionSeed), InputError);
  fissionSeed.fissions = -1.0;
  EXPECT_THROW(seed::seedFromFission(data, fissionSeed), InputError);
}

// The provenance is what a report header shows, and it has to carry enough to reconstruct the
// run: which nuclide, at what energy, how many fissions, and the yield-sum integrity check.
TEST(SeedFission, ProvenanceRecordsTheSeedAndItsYieldSum) {
  const NuclearData data = chainWithYields();
  seed::FissionSeed fissionSeed;
  fissionSeed.fissile = kU235;
  fissionSeed.incidentEnergyEv = kThermalEv;
  fissionSeed.fissions = 1.0e20;

  const std::string provenance = seed::seedFromFission(data, fissionSeed).provenance();
  EXPECT_NE(provenance.find("U-235"), std::string::npos) << provenance;
  EXPECT_NE(provenance.find("thermal"), std::string::npos) << provenance;
  EXPECT_NE(provenance.find("1e+20"), std::string::npos) << provenance;
  EXPECT_NE(provenance.find("sum Y_indep = 2"), std::string::npos) << provenance;
  EXPECT_EQ(provenance.find("WARNING"), std::string::npos)
      << "a set summing to 2.0 is correct and should not be flagged: " << provenance;
}

// Independent yields sum to about 2.0 because fission makes two fragments. Cumulative yields
// sum to far more, so a wrong-MT staging shows up here -- and every downstream number would be
// wrong by that factor with nothing else to reveal it.
TEST(SeedFission, WarnsWhenTheYieldSumIsNotAboutTwo) {
  StoreArrays arrays = synth::linearChain({1.0e-3, 5.0e-4});
  arrays.nfyParentKey = {Zai{92, 235, 0}.key()};
  arrays.nfyEnergyEv = {0.0253};
  arrays.nfySetOffset = {0, 2};
  arrays.nfyProductKey = {Zai{50, 100, 0}.key(), Zai{51, 100, 0}.key()};
  arrays.nfyProductYield = {3.0, 2.5};  // sums to 5.5: not independent yields
  const NuclearData data = NuclearData::fromArrays(std::move(arrays));

  seed::FissionSeed fissionSeed;
  fissionSeed.fissile = kU235;
  fissionSeed.fissions = 1.0e20;
  const std::string provenance = seed::seedFromFission(data, fissionSeed).provenance();
  EXPECT_NE(provenance.find("WARNING"), std::string::npos) << provenance;
}

// A store staged without yields cannot seed from fission at all, and says so.
TEST(SeedFission, AStoreWithNoYieldsSaysSo) {
  const NuclearData data = NuclearData::fromArrays(synth::linearChain({1.0e-3}));
  ASSERT_TRUE(data.fissionYields().empty());

  seed::FissionSeed fissionSeed;
  fissionSeed.fissile = kU235;
  fissionSeed.fissions = 1.0e20;
  try {
    seed::seedFromFission(data, fissionSeed);
    FAIL() << "expected InputError";
  } catch (const InputError& e) {
    EXPECT_NE(std::string(e.what()).find("none at all"), std::string::npos) << e.what();
  }
}

}  // namespace
}  // namespace nusift
