// dR/dlambda: how much the answer rests on the evaluated half-lives.
//
// Checked against closed forms rather than against a finer run of itself. For a single decaying
// nuclide with a stable daughter the activity is R = lambda N0 e^{-lambda T}, so
//
//     dR/dlambda = N0 e^{-lambda T} (1 - lambda T)      elasticity = 1 - lambda T
//
// and the two terms are separately known: the implicit one is -lambda N0 T e^{-lambda T} and the
// explicit one is N0 e^{-lambda T}. That pins each half of the derivative, not merely their sum,
// which matters because at equilibrium the sum is what cancels.
//
// The basis term has its own closed form and differs from the first by EXACTLY one. A row given
// in becquerel fixes A0, so R = A0 e^{-lambda T} with no lambda in front and the elasticity is
// -lambda T. Two expressions a unit apart is as sharp a test of that term as exists.
//
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <vector>

#include "nusift/core/error.hpp"
#include "nusift/engine/sensitivity.hpp"
#include "nusift/nucdata/nuclear_data.hpp"
#include "synthetic_chain.hpp"

namespace nusift {
namespace {

constexpr double kLambda = 1.0e-6;

// One unstable nuclide with a stable terminator: it holds all the activity itself, so every
// closed form above applies exactly.
NuclearData lone() {
  StoreArrays arrays = synth::linearChain({kLambda});
  arrays.awr = {99.1, 99.2};
  return NuclearData::fromArrays(std::move(arrays));
}

const Zai kParent{50, 100, 0};

const DecaySensitivity& parentOf(const DecaySensitivities& s) {
  const auto it = std::find_if(s.nuclides.begin(), s.nuclides.end(),
                               [](const DecaySensitivity& one) { return one.label == "Sn-100"; });
  EXPECT_NE(it, s.nuclides.end());
  return *it;
}

TEST(DecaySensitivity, MatchesTheClosedFormTermByTerm) {
  const NuclearData data = lone();
  const double atoms = 1.0e20;
  Inventory inv;
  inv.add(kParent, atoms);

  ResponseSpec spec;  // activity in becquerel
  for (const double time : {1.0e5, 1.0e6, 3.0e6}) {
    const DecaySensitivities s = decaySensitivities(data, inv, time, spec);
    const DecaySensitivity& one = parentOf(s);

    const double decayed = atoms * std::exp(-kLambda * time);
    // Both halves, not merely the sum: at equilibrium the sum is what cancels, so a test that
    // only checked the total could pass with both terms wrong by the same amount.
    EXPECT_NEAR(one.implicit, -kLambda * decayed * time, std::abs(one.implicit) * 1.0e-4)
        << "at t=" << time;
    EXPECT_NEAR(one.explicitWeight, decayed, decayed * 1.0e-9) << "at t=" << time;
    EXPECT_DOUBLE_EQ(one.basis, 0.0) << "an atoms-specified row does not move with lambda";
    // Judged against the scale of the TERMS, not of their difference. At lambda*T = 1 the two
    // are equal and opposite and the total is exactly zero, so a relative tolerance on the total
    // would be a tolerance on nothing -- which is the cancellation this file's header warns
    // about, in miniature and at a time this loop deliberately visits.
    EXPECT_NEAR(one.total, decayed * (1.0 - kLambda * time), decayed * 1.0e-4) << "at t=" << time;

    // The number the report leads with, and the one that survives the cancellation.
    EXPECT_NEAR(one.elasticity, 1.0 - kLambda * time, 1.0e-4) << "at t=" << time;
    EXPECT_NEAR(s.response, kLambda * decayed, kLambda * decayed * 1.0e-9);
  }
}

// The cancellation itself, at the instant it is total. lambda*T = 1 makes the implicit and
// explicit terms exactly equal and opposite, and what survives is fifteen orders of magnitude
// below either -- which is why the elasticity is what this reports and why a relative error
// against the total would be an error against nothing.
TEST(DecaySensitivity, ShowsTheTermsCancellingWhereTheyDo) {
  const NuclearData data = lone();
  Inventory inv;
  inv.add(kParent, 1.0e20);

  ResponseSpec spec;
  const DecaySensitivities s = decaySensitivities(data, inv, 1.0 / kLambda, spec);
  const DecaySensitivity& one = parentOf(s);

  EXPECT_NEAR(one.implicit, -one.explicitWeight, std::abs(one.explicitWeight) * 1.0e-6);
  EXPECT_LT(std::abs(one.total), std::abs(one.explicitWeight) * 1.0e-10)
      << "the sum is what cancels; both terms are large";
  EXPECT_NEAR(one.elasticity, 0.0, 1.0e-6)
      << "and the elasticity says so without dividing by a vanishing quantity";
}

// The basis term, isolated. The same nuclide and the same instant, specified two ways: as an
// atom count and as the activity that count represents. The response is identical and the
// elasticities differ by exactly one, which is the whole content of the third term.
TEST(DecaySensitivity, TheMeasurementBasisMovesTheDerivativeByExactlyOne) {
  const NuclearData data = lone();
  const double atoms = 1.0e20;
  const double time = 2.0e6;

  Inventory byAtoms;
  byAtoms.add(kParent, atoms);
  Inventory byActivity;
  // Same atoms, but declared to have come from an activity measurement -- which fixes A0 rather
  // than n0, so the seed itself moves when lambda does.
  byActivity.add(kParent, atoms, /*sigmaAtoms=*/0.0, /*atomsFromActivity=*/atoms);

  ResponseSpec spec;
  const DecaySensitivities a = decaySensitivities(data, byAtoms, time, spec);
  const DecaySensitivities b = decaySensitivities(data, byActivity, time, spec);

  EXPECT_NEAR(a.response, b.response, a.response * 1.0e-12) << "the same material either way";
  EXPECT_NEAR(parentOf(a).elasticity, 1.0 - kLambda * time, 1.0e-4);
  EXPECT_NEAR(parentOf(b).elasticity, -kLambda * time, 1.0e-4);
  EXPECT_NEAR(parentOf(a).elasticity - parentOf(b).elasticity, 1.0, 1.0e-6);
  EXPECT_LT(parentOf(b).basis, 0.0) << "a shorter half-life means fewer atoms for the same Bq";
}

// A pack whose coefficient multiplies ATOMS carries no lambda, so its weight does not move when
// a decay constant does. Assuming w/lambda there would invent a term that is not present, which
// is why the derivative is asked for rather than derived from the weight.
TEST(DecaySensitivity, KnowsWhichMetricsWeightsMoveWithLambda) {
  const NuclearData data = lone();
  ResponseSpec activity;
  const std::vector<double> byActivity = weightDecayDerivatives(data, activity);

  const int index = data.indexOf(kParent);
  ASSERT_GE(index, 0);
  // Activity's weight IS lambda, so its derivative is exactly one.
  EXPECT_NEAR(byActivity[static_cast<std::size_t>(index)], 1.0, 1.0e-12);

  // A stable nuclide has no lambda to differentiate with respect to.
  const int stable = data.indexOf(Zai{51, 100, 0});
  ASSERT_GE(stable, 0);
  EXPECT_DOUBLE_EQ(byActivity[static_cast<std::size_t>(stable)], 0.0);
}

// cram's own rule: refine until the smallest quadrature piece falls below the shortest removal
// time. Its default of 6 is for "a thermal pin with a few intervals of days to months", and a
// decay chain spanning minutes to millennia needs far more -- a single 30-day interval against
// Ba-137m's 221 s removal time needs twelve.
TEST(DecaySensitivity, ChoosesRefinementFromTheShortestRemovalTime) {
  // 30 days in four pieces is 6.48e5 s; reaching 221 s needs log2(2932) = 11.5 -> 12.
  EXPECT_EQ(chooseEndRefinements(30.0 * 86400.0, 4, 221.0, 40), 12);
  // A piece already shorter than the removal time needs no refinement at all.
  EXPECT_EQ(chooseEndRefinements(100.0, 4, 1000.0, 40), 0);
  // And the cap binds rather than running away on a problem spanning twenty decades.
  EXPECT_EQ(chooseEndRefinements(1.0e18, 4, 1.0e-6, 8), 8);
  EXPECT_EQ(chooseEndRefinements(0.0, 4, 1.0, 40), 0);
}

TEST(DecaySensitivity, RefusesAQuestionItCannotAnswer) {
  const NuclearData data = lone();
  Inventory inv;
  inv.add(kParent, 1.0e20);
  ResponseSpec spec;

  // At t = 0 the response is the seed and no half-life has acted on it, so the derivative is
  // identically zero and the schedule would be empty.
  EXPECT_THROW(decaySensitivities(data, inv, 0.0, spec), InputError);
  EXPECT_THROW(decaySensitivities(data, inv, -1.0, spec), InputError);

  ResponseSpec interval;
  interval.unit = Unit::Decays;
  EXPECT_THROW(decaySensitivities(data, inv, 1.0e6, interval), InputError);

  SensitivityOptions bad;
  bad.gaussPoints = 9;
  EXPECT_THROW(decaySensitivities(data, inv, 1.0e6, spec, bad), InputError);
}

}  // namespace
}  // namespace nusift
