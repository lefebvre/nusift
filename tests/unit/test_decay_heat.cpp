// Decay heat: the power an inventory releases as it decays.
//
// The metric is lambda times the sum of ENDF MT457's three average decay energies, which makes
// it the simplest weight in the file after activity -- no geometry, no spectrum, no kernel. What
// is worth testing is therefore not the arithmetic but the things the arithmetic is easy to get
// wrong around it:
//
//   that the three energies are SUMMED rather than one of them being taken for the total, which
//   a single-component test nuclide could never catch;
//
//   that watts against atoms and joules against atom-seconds are the SAME weight, since decay
//   heat is the one metric whose rate and total differ by no hour -- exposure's per-hour scaling
//   applied here would be wrong by 3600;
//
//   and that a nuclide the evaluation gave no energies is worth zero watts rather than being
//   quietly skipped, because "no staged energy" and "no heat" are different claims and only the
//   first is ever true of an unstable nuclide.
//
#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <vector>

#include "nusift/core/error.hpp"
#include "nusift/engine/decay_engine.hpp"
#include "nusift/nucdata/nuclear_data.hpp"
#include "nusift/triage/ranking.hpp"
#include "nusift/triage/response.hpp"
#include "nusift/units.hpp"
#include "synthetic_chain.hpp"

namespace nusift {
namespace {

constexpr double kLambda = 1.0e-6;
const Zai kParent{50, 100, 0};

// A two-member chain whose parent carries all three decay energies, deliberately unequal so that
// summing them is distinguishable from taking any one.
NuclearData chainWithEnergies(double em = 300000.0, double lp = 200000.0, double hp = 500000.0) {
  StoreArrays arrays = synth::linearChain({kLambda, 5.0e-7});
  arrays.emEnergyEv.assign(static_cast<std::size_t>(arrays.nuclideCount()), 0.0);
  arrays.lpEnergyEv.assign(static_cast<std::size_t>(arrays.nuclideCount()), 0.0);
  arrays.hpEnergyEv.assign(static_cast<std::size_t>(arrays.nuclideCount()), 0.0);
  arrays.emEnergyEv[0] = em;
  arrays.lpEnergyEv[0] = lp;
  arrays.hpEnergyEv[0] = hp;
  return NuclearData::fromArrays(std::move(arrays));
}

Inventory one(const Zai& zai, double atoms) {
  Inventory inv;
  inv.add(zai, atoms);
  return inv;
}

// --- the weight --------------------------------------------------------------

// P = lambda N E, with E the SUM of the three. The three are unequal here, so a weight that took
// any single one -- or any two -- differs from this by more than a tolerance could hide.
TEST(DecayHeat, IsLambdaTimesTheSummedDecayEnergy) {
  const NuclearData data = chainWithEnergies();
  const double atoms = 1.0e20;
  const double totalEv = 300000.0 + 200000.0 + 500000.0;
  const double expected = kLambda * atoms * totalEv * units::kEvToJ;

  ResponseSpec spec;
  spec.metric = Metric::Heat;
  spec.unit = Unit::Watt;

  const std::vector<double> when = {0.0};
  const DecayResult result = decay(data, one(kParent, atoms), when);
  const ResponseTable table = buildResponse(data, result, spec);
  EXPECT_NEAR(table.totals.front(), expected, expected * 1.0e-12);

  // And each component alone, so the summing is pinned rather than inferred from one total that
  // several wrong weights could also produce.
  const std::vector<std::array<double, 3>> single_component = {
      {1.0e6, 0.0, 0.0}, {0.0, 1.0e6, 0.0}, {0.0, 0.0, 1.0e6}};
  for (const std::array<double, 3>& e : single_component) {
    const double em = e[0], lp = e[1], hp = e[2];
    const NuclearData single = chainWithEnergies(em, lp, hp);
    const std::vector<double> at = {0.0};
    const DecayResult r = decay(single, one(kParent, atoms), at);
    const ResponseTable t = buildResponse(single, r, spec);
    const double want = kLambda * atoms * 1.0e6 * units::kEvToJ;
    EXPECT_NEAR(t.totals.front(), want, want * 1.0e-12);
  }
}

// The rate and the total are one weight with no hour between them, unlike exposure. Over a
// window the energy released is the integral of the power, and for a single nuclide decaying
// from N0 that integral has a closed form: E = N0 (1 - e^{-lambda T}) * E_decay -- the energy
// per decay times the number of decays, which is what makes this checkable to every digit.
TEST(DecayHeat, JoulesAreTheIntegralOfWatts) {
  const NuclearData data = chainWithEnergies();
  const double atoms = 1.0e20;
  const double totalEv = 1.0e6;
  const double window = 2.0e6;

  ResponseSpec spec;
  spec.metric = Metric::Heat;
  spec.unit = Unit::Joule;

  std::vector<std::int64_t> keys;
  const std::vector<double> integral =
      intervalIntegral(data, one(kParent, atoms), 0.0, window, &keys);
  const ResponseTable table = buildIntervalResponse(data, keys, integral, 0.0, window, spec);

  const double decays = atoms * (1.0 - std::exp(-kLambda * window));
  const double expected = decays * totalEv * units::kEvToJ;
  EXPECT_NEAR(table.totals.front(), expected, expected * 1.0e-10);
}

// --- what the units mean -----------------------------------------------------

// A watt is a rate and a joule is a total, and each refuses the other's domain -- the same line
// every metric here draws, and the reason Domain is an axis rather than two metrics.
TEST(DecayHeat, UnitsAreGatedByDomainAndMetric) {
  EXPECT_TRUE(unitSuitsDomain(Unit::Watt, Domain::Instant));
  EXPECT_FALSE(unitSuitsDomain(Unit::Watt, Domain::Interval));
  EXPECT_TRUE(unitSuitsDomain(Unit::Joule, Domain::Interval));
  EXPECT_FALSE(unitSuitsDomain(Unit::Joule, Domain::Instant));

  EXPECT_TRUE(unitSuitsMetric(Unit::Watt, Metric::Heat));
  EXPECT_TRUE(unitSuitsMetric(Unit::Joule, Metric::Heat));
  // Reporting decay heat in becquerel is a category error rather than a rounding one, and is
  // refused at the boundary like every other.
  EXPECT_FALSE(unitSuitsMetric(Unit::Becquerel, Metric::Heat));
  EXPECT_FALSE(unitSuitsMetric(Unit::Watt, Metric::Activity));
  EXPECT_FALSE(unitSuitsMetric(Unit::Watt, Metric::Exposure));

  Unit parsed = Unit::Becquerel;
  EXPECT_TRUE(parseUnit("W", parsed));
  EXPECT_EQ(parsed, Unit::Watt);
  EXPECT_TRUE(parseUnit("J", parsed));
  EXPECT_EQ(parsed, Unit::Joule);

  EXPECT_EQ(defaultUnit(Metric::Heat, Domain::Instant), Unit::Watt);
  EXPECT_EQ(defaultUnit(Metric::Heat, Domain::Interval), Unit::Joule);
}

// Decay heat needs no photon lines, which is what separates it from exposure: a store staged
// without spectra can still answer it, because MT457's averages are staged whether or not the
// discrete lines were.
TEST(DecayHeat, NeedsNoPhotonLines) {
  const NuclearData data = chainWithEnergies();
  ASSERT_FALSE(data.hasPhotonLines());

  ResponseSpec spec;
  spec.metric = Metric::Heat;
  spec.unit = Unit::Watt;
  const std::vector<double> when = {0.0};
  const DecayResult result = decay(data, one(kParent, 1.0e20), when);
  EXPECT_GT(buildResponse(data, result, spec).totals.front(), 0.0);

  // Where exposure on the same store is refused outright rather than answered with zeros.
  ResponseSpec exposure;
  exposure.metric = Metric::Exposure;
  exposure.unit = Unit::RoentgenPerHour;
  EXPECT_THROW(buildResponse(data, result, exposure), InputError);
}

// --- what a missing energy means ---------------------------------------------

// A nuclide the evaluation gave no average energies contributes zero watts. It is not skipped
// and it is not an error: 18 of the store's 3562 unstable nuclides are in this state, and the
// honest report of their heat is zero with the coverage said elsewhere.
TEST(DecayHeat, ANuclideWithNoStagedEnergyIsWorthZero) {
  const NuclearData data = chainWithEnergies(0.0, 0.0, 0.0);
  ResponseSpec spec;
  spec.metric = Metric::Heat;
  spec.unit = Unit::Watt;

  const std::vector<double> when = {0.0};
  const DecayResult result = decay(data, one(kParent, 1.0e20), when);
  const ResponseTable table = buildResponse(data, result, spec);
  EXPECT_DOUBLE_EQ(table.totals.front(), 0.0);
}

// A store staged before the two new columns existed carries neither, and every nuclide's heat is
// then whatever the electromagnetic term alone says -- zero here. The point is that it LOADS: an
// older store is readable rather than rejected, which is what the columns being optional means.
TEST(DecayHeat, AStoreWithoutTheNewColumnsStillLoads) {
  StoreArrays arrays = synth::linearChain({kLambda, 5.0e-7});
  arrays.emEnergyEv.assign(static_cast<std::size_t>(arrays.nuclideCount()), 1.0e6);
  // lpEnergyEv and hpEnergyEv left empty, as a store staged before them would have.
  const NuclearData data = NuclearData::fromArrays(std::move(arrays));

  EXPECT_DOUBLE_EQ(data.lpEnergyEv(0), 0.0);
  EXPECT_DOUBLE_EQ(data.hpEnergyEv(0), 0.0);
  EXPECT_DOUBLE_EQ(data.decayEnergyEv(0), 1.0e6);

  ResponseSpec spec;
  spec.metric = Metric::Heat;
  spec.unit = Unit::Watt;
  const std::vector<double> when = {0.0};
  const DecayResult result = decay(data, one(kParent, 1.0e20), when);
  const double expected = kLambda * 1.0e20 * 1.0e6 * units::kEvToJ;
  EXPECT_NEAR(buildResponse(data, result, spec).totals.front(), expected, expected * 1.0e-12);
}

}  // namespace
}  // namespace nusift
