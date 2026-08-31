#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "nusift/core/error.hpp"
#include "nusift/engine/decay_engine.hpp"
#include "nusift/io/time_spec.hpp"
#include "nusift/nucdata/nuclear_data.hpp"
#include "nusift/triage/allowable.hpp"
#include "nusift/triage/events.hpp"
#include "synthetic_chain.hpp"

namespace nusift {
namespace {

// One nuclide decaying to a stable daughter. Activity is lambda * N0 * e^{-lambda t}, so every
// scale below has a closed form.
NuclearData oneEmitter(double lambda) {
  StoreArrays a;
  a.provenance.version = 1;
  a.nuclideKey = {Zai{50, 100, 0}.key(), Zai{51, 100, 0}.key()};
  a.halfLife = {synth::halfLifeFor(lambda), 0.0};
  a.modeOffset = {0, 1, 1};
  a.modeRtyp = {synth::kBetaMinus};
  a.modeBranching = {1.0};
  a.modeFinalState = {0};
  a.modeIsFission = {0};
  return NuclearData::fromArrays(std::move(a));
}

// A pure beta parent feeding a photon-emitting daughter, so activity and exposure are carried
// by different members and two criteria over them can genuinely trade places.
NuclearData betaParentPhotonDaughter(double parent, double daughter) {
  StoreArrays a = synth::linearChain({parent, daughter});
  synth::addLines(a, 1, {661657.0}, {0.9});
  return NuclearData::fromArrays(std::move(a));
}

Criterion activityLimit(const std::string& name, double becquerel) {
  Criterion criterion;
  criterion.name = name;
  criterion.limit = becquerel;
  return criterion;  // ResponseSpec defaults to instantaneous activity in Bq
}

DecayResult decayOf(const NuclearData& data, double atoms, const std::vector<double>& times,
                    int z = 50) {
  Inventory inv;
  inv.add(Zai{z, 100, 0}, atoms);
  return decay(data, inv, times);
}

// --- the scale itself ---------------------------------------------------------

TEST(AllowableScale, ScaleIsTheLimitOverTheResponse) {
  const double lambda = 1.0e-4;
  const double atoms = 1.0e20;
  const double limit = 1.0e14;
  const NuclearData data = oneEmitter(lambda);

  const std::vector<double> times = logspace(1.0, 1.0e5, 12);
  const DecayResult result = decayOf(data, atoms, times);

  const std::vector<Criterion> criteria = {activityLimit("possession", limit)};
  const std::vector<AllowableScale> scaled = allowableScale(data, result, criteria);

  ASSERT_EQ(scaled.size(), times.size());
  for (std::size_t k = 0; k < times.size(); ++k) {
    const double expected = limit / (lambda * atoms * std::exp(-lambda * times[k]));
    EXPECT_TRUE(scaled[k].bounded);
    EXPECT_NEAR(scaled[k].scale, expected, expected * 1.0e-9);
    EXPECT_EQ(scaled[k].bindingIndex, 0);
    ASSERT_EQ(scaled[k].criteria.size(), 1u);
    EXPECT_TRUE(scaled[k].criteria[0].binding);
    EXPECT_NEAR(scaled[k].criteria[0].fraction, 1.0 / expected, 1.0e-9)
        << "the sum-of-fractions idiom is the same number inverted";
  }
}

// The exactness claim, checked the only way that means anything: multiply the inventory BY the
// answer, decay it again, and the binding criterion's response has to land on its limit. This
// is what "exact for the same reason the shares are" buys -- no search, no tolerance.
TEST(AllowableScale, ScalingTheInventoryByTheAnswerMeetsTheLimitExactly) {
  const double lambda = 1.0e-4;
  const double atoms = 1.0e20;
  const double limit = 3.7e13;
  const NuclearData data = oneEmitter(lambda);

  const std::vector<double> times = {0.0, 1.0e3, 1.0e4};
  const std::vector<Criterion> criteria = {activityLimit("possession", limit)};
  const std::vector<AllowableScale> scaled =
      allowableScale(data, decayOf(data, atoms, times), criteria);

  for (std::size_t k = 0; k < times.size(); ++k) {
    ASSERT_TRUE(scaled[k].bounded);
    // Decay the SCALED inventory and read the response the report would print for it.
    const DecayResult again = decayOf(data, atoms * scaled[k].scale, times);
    const ResponseTable table = buildResponse(data, again, criteria[0].spec);
    EXPECT_NEAR(table.totals[k], limit, limit * 1.0e-9)
        << "at t=" << times[k] << " the scaled inventory must sit exactly on the limit";
  }
}

TEST(AllowableScale, TheTightestCriterionBindsAndIsTheOnlyOneFlagged) {
  const NuclearData data = oneEmitter(1.0e-4);
  const DecayResult result = decayOf(data, 1.0e20, {0.0});

  const std::vector<Criterion> criteria = {
      activityLimit("loose", 1.0e18),
      activityLimit("tight", 1.0e14),
      activityLimit("middling", 1.0e16),
  };
  const std::vector<AllowableScale> scaled = allowableScale(data, result, criteria);

  ASSERT_EQ(scaled.size(), 1u);
  EXPECT_EQ(scaled[0].bindingIndex, 1);
  EXPECT_EQ(scaled[0].criteria[1].name, "tight");
  EXPECT_TRUE(scaled[0].criteria[1].binding);
  EXPECT_FALSE(scaled[0].criteria[0].binding);
  EXPECT_FALSE(scaled[0].criteria[2].binding);

  // The reported scale is the minimum, and every criterion still reports what it alone allows.
  EXPECT_DOUBLE_EQ(scaled[0].scale, scaled[0].criteria[1].scale);
  EXPECT_LT(scaled[0].criteria[1].scale, scaled[0].criteria[2].scale);
  EXPECT_LT(scaled[0].criteria[2].scale, scaled[0].criteria[0].scale);
}

// Two criteria over quantities carried by different members of the chain. At t=0 the daughter
// does not exist, so the exposure limit constrains nothing; once it grows in, exposure becomes
// the binding criterion. Which limit binds is not a property of the inventory alone.
TEST(AllowableScale, TheBindingCriterionChangesAsTheChainEvolves) {
  const NuclearData data = betaParentPhotonDaughter(1.0e-3, 1.0e-5);
  // t=0 included deliberately: it is the one sample where the daughter genuinely does not
  // exist, so the exposure criterion has nothing to constrain and says so.
  std::vector<double> times = {0.0};
  for (const double t : logspace(1.0, 1.0e6, 40)) {
    times.push_back(t);
  }
  const DecayResult result = decayOf(data, 1.0e20, times);

  Criterion exposure;
  exposure.name = "dose rate at 1 m";
  exposure.spec.metric = Metric::Exposure;
  exposure.spec.unit = Unit::RoentgenPerHour;
  exposure.spec.geometry.distanceM = 1.0;
  exposure.limit = 1.0e-2;

  const std::vector<Criterion> criteria = {activityLimit("possession", 1.0e15), exposure};
  const std::vector<AllowableScale> scaled = allowableScale(data, result, criteria);

  EXPECT_EQ(scaled.front().bindingIndex, 0) << "no daughter yet, so only activity constrains";
  EXPECT_TRUE(scaled.front().criteria[1].unbounded);

  const bool everExposure = std::any_of(
      scaled.begin(), scaled.end(), [](const AllowableScale& at) { return at.bindingIndex == 1; });
  EXPECT_TRUE(everExposure) << "once the emitter grows in, the exposure limit takes over";

  // Whichever binds, the reported scale is always the smallest any criterion permits.
  for (const AllowableScale& at : scaled) {
    if (!at.bounded) {
      continue;
    }
    for (const CriterionHeadroom& headroom : at.criteria) {
      if (!headroom.unbounded) {
        EXPECT_LE(at.scale, headroom.scale * (1.0 + 1.0e-12));
      }
    }
  }
}

// --- what "unbounded" means ---------------------------------------------------

// A criterion with nothing to constrain permits any scale. Reported as unbounded rather than
// as a very large number, because "nothing here is limited by this" and "this allows 1e18
// times the inventory" are different statements and only one of them is true.
TEST(AllowableScale, ACriterionWithNoResponseIsUnboundedRatherThanEnormous) {
  const NuclearData data = betaParentPhotonDaughter(1.0e-3, 1.0e-5);
  const DecayResult result = decayOf(data, 1.0e20, {0.0});

  Criterion exposure;
  exposure.name = "dose rate";
  exposure.spec.metric = Metric::Exposure;
  exposure.spec.unit = Unit::RoentgenPerHour;
  exposure.spec.geometry.distanceM = 1.0;
  exposure.limit = 1.0e-3;

  const std::vector<AllowableScale> scaled =
      allowableScale(data, result, std::vector<Criterion>{exposure});

  ASSERT_EQ(scaled.size(), 1u);
  EXPECT_TRUE(scaled[0].criteria[0].unbounded);
  EXPECT_DOUBLE_EQ(scaled[0].criteria[0].scale, 0.0) << "no scale, rather than a huge one";
  EXPECT_FALSE(scaled[0].bounded);
  EXPECT_EQ(scaled[0].bindingIndex, -1);
  EXPECT_TRUE(scaled[0].limiting.empty()) << "nothing binds, so nothing is limiting";
}

TEST(AllowableScale, AStableInventoryIsConstrainedByNoActivityLimitAtAll) {
  const NuclearData data = oneEmitter(1.0e-4);
  // Seed the stable daughter directly: it has no activity at any time.
  const DecayResult result = decayOf(data, 1.0e20, {0.0, 1.0e4}, /*z=*/51);

  const std::vector<Criterion> criteria = {activityLimit("possession", 1.0e10)};
  const std::vector<AllowableScale> scaled = allowableScale(data, result, criteria);

  for (const AllowableScale& at : scaled) {
    EXPECT_FALSE(at.bounded);
    EXPECT_TRUE(at.criteria[0].unbounded);
  }
}

// --- naming what drives the binding criterion ---------------------------------

TEST(AllowableScale, TheLimitingContributorsAreTheBindingCriterionsTopRows) {
  const NuclearData data = betaParentPhotonDaughter(1.0e-3, 1.0e-5);
  const DecayResult result = decayOf(data, 1.0e20, {1.0e3});

  const std::vector<Criterion> criteria = {activityLimit("possession", 1.0e15)};
  const std::vector<AllowableScale> scaled = allowableScale(data, result, criteria, 2);

  ASSERT_EQ(scaled.size(), 1u);
  ASSERT_FALSE(scaled[0].limiting.empty());
  EXPECT_LE(scaled[0].limiting.size(), 2u) << "limitingCount caps how many are named";

  // Ordered by share, and the shares are of the binding criterion's own total.
  for (std::size_t i = 1; i < scaled[0].limiting.size(); ++i) {
    EXPECT_GE(scaled[0].limiting[i - 1].fraction, scaled[0].limiting[i].fraction);
  }
  EXPECT_FALSE(scaled[0].limiting.front().label.empty());
  EXPECT_GT(scaled[0].limiting.front().fraction, 0.0);
}

TEST(AllowableScale, NamingCanBeTurnedOffForTheCurveAlone) {
  const NuclearData data = oneEmitter(1.0e-4);
  const DecayResult result = decayOf(data, 1.0e20, {0.0});
  const std::vector<Criterion> criteria = {activityLimit("possession", 1.0e15)};

  const std::vector<AllowableScale> scaled = allowableScale(data, result, criteria, 0);
  ASSERT_EQ(scaled.size(), 1u);
  EXPECT_TRUE(scaled[0].bounded) << "the scale is still computed";
  EXPECT_TRUE(scaled[0].limiting.empty());
}

// --- criteria that do not mean anything ---------------------------------------

TEST(AllowableScale, ACriterionThatCannotMeanAnythingIsRefused) {
  const NuclearData data = oneEmitter(1.0e-4);
  const DecayResult result = decayOf(data, 1.0e20, {0.0});

  EXPECT_THROW(allowableScale(data, result, std::vector<Criterion>{}), InputError)
      << "an unconstrained inventory has no maximum";

  Criterion unnamed = activityLimit("", 1.0e15);
  EXPECT_THROW(allowableScale(data, result, std::vector<Criterion>{unnamed}), InputError);

  EXPECT_THROW(allowableScale(data, result, std::vector<Criterion>{activityLimit("zero", 0.0)}),
               InputError);
  EXPECT_THROW(allowableScale(data, result, std::vector<Criterion>{activityLimit("neg", -1.0)}),
               InputError);

  // An interval unit is an accrued total over a window, which is not what a possession limit
  // constrains.
  Criterion accrued = activityLimit("decays", 1.0e15);
  accrued.spec.unit = Unit::Decays;
  EXPECT_THROW(allowableScale(data, result, std::vector<Criterion>{accrued}), InputError);
}

TEST(AllowableScale, TheRefusalNamesTheOffendingCriterion) {
  const NuclearData data = oneEmitter(1.0e-4);
  const DecayResult result = decayOf(data, 1.0e20, {0.0});
  const std::vector<Criterion> criteria = {activityLimit("fine", 1.0e15),
                                           activityLimit("the bad one", -3.0)};
  try {
    allowableScale(data, result, criteria);
    FAIL() << "expected a refusal";
  } catch (const InputError& error) {
    EXPECT_NE(std::string(error.what()).find("the bad one"), std::string::npos) << error.what();
  }
}

// --- the curve, and the event engine over it ----------------------------------

// s_max(t) = (L / lambda N0) e^{lambda t}, so it reaches 1 -- the instant the inventory as it
// stands becomes shippable -- at t = ln(lambda N0 / L) / lambda. The composition with the event
// engine is the point: a scale curve is a curve like any other.
TEST(AllowableScale, TheScaleCurveLocatesWhenTheInventoryBecomesShippable) {
  const double lambda = 1.0e-4;
  const double atoms = 1.0e20;
  const double limit = 1.0e13;
  const double expected = std::log(lambda * atoms / limit) / lambda;
  ASSERT_GT(expected, 0.0);

  const NuclearData data = oneEmitter(lambda);
  const DecayResult result = decayOf(data, atoms, logspace(1.0, 1.0e6, 60));
  const std::vector<Criterion> criteria = {activityLimit("transport", limit)};

  const EventSeries series = scaleSeries(allowableScale(data, result, criteria, 0));
  const std::optional<TrajectoryEvent> shippable = firstCrossing(series, 1.0);

  ASSERT_TRUE(shippable.has_value());
  EXPECT_EQ(shippable->kind, EventKind::Rising) << "headroom grows as the material decays";
  EXPECT_NEAR(shippable->timeSeconds, expected, expected * 0.01);
  EXPECT_GE(shippable->timeSeconds, shippable->bracketStartSeconds);
  EXPECT_LE(shippable->timeSeconds, shippable->bracketEndSeconds);
}

TEST(AllowableScale, TheScaleCurveRefusesAGridWithNothingToSearch) {
  std::vector<AllowableScale> none(3);  // all default-constructed: unbounded everywhere
  EXPECT_THROW(scaleSeries(none), InputError);

  std::vector<AllowableScale> one(3);
  one[1].timeSeconds = 1.0;
  one[1].bounded = true;
  one[1].scale = 2.0;
  EXPECT_THROW(scaleSeries(one), InputError) << "one bounded sample is not a curve";
}

// Bounded, then not, then bounded again has no single curve to search, and stitching across the
// gap would invent a scale for a stretch where nothing constrained the inventory.
TEST(AllowableScale, TheScaleCurveRefusesAGapInTheMiddle) {
  std::vector<AllowableScale> gapped(4);
  for (std::size_t k = 0; k < gapped.size(); ++k) {
    gapped[k].timeSeconds = static_cast<double>(k);
    gapped[k].bounded = k != 2;
    gapped[k].scale = 1.0 + static_cast<double>(k);
  }
  EXPECT_THROW(scaleSeries(gapped), InputError);
}

}  // namespace
}  // namespace nusift
