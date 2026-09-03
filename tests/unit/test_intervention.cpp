#include <gtest/gtest.h>

#include <cmath>
#include <string>
#include <vector>

#include "nusift/core/error.hpp"
#include "nusift/engine/decay_engine.hpp"
#include "nusift/nucdata/nuclear_data.hpp"
#include "nusift/triage/intervention.hpp"
#include "nusift/triage/response.hpp"
#include "synthetic_chain.hpp"

namespace nusift {
namespace {

// Sn-100 -> Sb-100 -> Te-100 (stable). Three DIFFERENT elements on one mass chain, which is
// what makes "remove the cesium and the barium stays" testable at all.
NuclearData twoStepChain(double parent, double daughter) {
  StoreArrays arrays = synth::linearChain({parent, daughter});
  return NuclearData::fromArrays(std::move(arrays));
}

Inventory seededWith(std::int64_t key, double atoms) {
  Inventory inv;
  inv.addKey(key, atoms);
  return inv;
}

const std::int64_t kParent = Zai{50, 100, 0}.key();
const std::int64_t kDaughter = Zai{51, 100, 0}.key();

// R(T) the ordinary way: decay this inventory and read the total a ranking would print.
double responseAfter(const NuclearData& data, const Inventory& inv, double window,
                     const ResponseSpec& spec) {
  const DecayResult result = decay(data, inv, std::vector<double>{window});
  return buildResponse(data, result, spec).totals[0];
}

// The inventory at t0 with a fraction of one nuclide taken out, built by hand so the adjoint's
// answer can be checked against a full forward solve of the counterfactual rather than against
// itself.
Inventory withRemoval(const NuclearData& data, const Inventory& inv, double t0, std::int64_t key,
                      double fraction) {
  const DecayResult result = decay(data, inv, std::vector<double>{t0});
  Inventory at;
  for (int i = 0; i < result.nuclideCount(); ++i) {
    const std::int64_t here = result.nuclideKeys[static_cast<std::size_t>(i)];
    const double atoms = result.atomsAt(0)[static_cast<std::size_t>(i)];
    if (atoms > 0.0) {
      at.addKey(here, here == key ? atoms * (1.0 - fraction) : atoms);
    }
  }
  return at;
}

Intervention removing(const std::string& name, const std::string& selector, double fraction) {
  Intervention intervention;
  intervention.name = name;
  intervention.removals.push_back(Removal{selector, fraction});
  return intervention;
}

// --- the benefit is exact ------------------------------------------------------

// The whole claim: R is linear in n(t0), so the adjoint's answer for a counterfactual has to
// equal a full forward solve of that counterfactual. Not close to it -- equal, because neither
// side involves a search, a perturbation, or a tolerance.
TEST(Intervention, TheBenefitMatchesAFullForwardSolveOfTheCounterfactual) {
  const NuclearData data = twoStepChain(1.0e-4, 1.0e-5);
  const Inventory inv = seededWith(kParent, 1.0e20);
  const double t0 = 1.0e4;
  const double target = 1.0e6;

  for (const double fraction : {0.25, 0.5, 1.0}) {
    const std::vector<Intervention> plan = {removing("strip the parent", "Sn-100", fraction)};
    const InterventionStudy study =
        compareInterventions(data, inv, t0, target, ResponseSpec{}, plan);

    ASSERT_EQ(study.effects.size(), 1u);
    const double independently = responseAfter(data, withRemoval(data, inv, t0, kParent, fraction),
                                               target - t0, ResponseSpec{});

    EXPECT_NEAR(study.effects[0].response, independently, std::abs(independently) * 1.0e-9)
        << "fraction " << fraction;
    EXPECT_NEAR(study.effects[0].removed, study.baseline - independently,
                std::abs(study.baseline) * 1.0e-9);
  }
}

// The adjoint's baseline and the forward path's total are the same number by the duality
// identity. Asserted rather than trusted, because a units or index-space slip would show up
// here first.
TEST(Intervention, TheBaselineIsTheOrdinaryRankingTotal) {
  const NuclearData data = twoStepChain(1.0e-4, 1.0e-5);
  const Inventory inv = seededWith(kParent, 1.0e20);

  const std::vector<Intervention> plan = {removing("none really", "Sn-100", 0.0)};
  const InterventionStudy study =
      compareInterventions(data, inv, 1.0e4, 1.0e6, ResponseSpec{}, plan);

  const double forward = responseAfter(data, inv, 1.0e6, ResponseSpec{});
  EXPECT_NEAR(study.baseline, forward, forward * 1.0e-9);

  // Removing nothing is worth nothing, and leaves the response where it was.
  EXPECT_DOUBLE_EQ(study.effects[0].removed, 0.0);
  EXPECT_NEAR(study.effects[0].response, study.baseline, study.baseline * 1.0e-12);
}

// --- what removing a parent does, and does not do -------------------------------

// Removing a parent takes with it everything it would have gone on to produce. The importance
// vector already carries the whole forward evolution, so this needs no special handling -- but
// it is the non-obvious half of the answer and has to be true.
TEST(Intervention, RemovingAParentAlsoRemovesTheDaughtersItWouldHaveFed) {
  // A short-lived parent feeding a long-lived daughter: by the target time essentially all of
  // the response is daughter, and essentially all of the daughter came from the parent.
  const NuclearData data = twoStepChain(1.0e-3, 1.0e-8);
  const Inventory inv = seededWith(kParent, 1.0e20);

  const std::vector<Intervention> plan = {removing("strip the parent", "Sn-100", 1.0)};
  const InterventionStudy study = compareInterventions(data, inv, 0.0, 1.0e6, ResponseSpec{}, plan);

  ASSERT_EQ(study.effects.size(), 1u);
  // Everything was seeded as parent, so removing all of it at t=0 leaves nothing at all.
  EXPECT_NEAR(study.effects[0].removedFraction, 1.0, 1.0e-9);
  EXPECT_NEAR(study.effects[0].response, 0.0, study.baseline * 1.0e-9);

  // And the benefit is booked against the parent, which is what was removed -- not against the
  // daughter that would have carried the response.
  ASSERT_EQ(study.effects[0].contributors.size(), 1u);
  EXPECT_EQ(study.effects[0].contributors[0].label, "Sn-100");
}

// The other half, and the one that makes a separation a separation: the daughter standing there
// at t0 is a different element and is not removed with its parent.
TEST(Intervention, RemovingAnElementLeavesTheDaughterAlreadyPresent) {
  const NuclearData data = twoStepChain(1.0e-3, 1.0e-6);
  const Inventory inv = seededWith(kParent, 1.0e20);
  // Late enough that most of the parent has already become daughter.
  const double t0 = 1.0e5;

  const std::vector<Intervention> plan = {removing("separate Sn", "Sn", 1.0)};
  const InterventionStudy study =
      compareInterventions(data, inv, t0, t0 + 1.0, ResponseSpec{}, plan);

  ASSERT_EQ(study.effects.size(), 1u);
  // The response one second later is nearly all daughter, and the daughter is untouched, so
  // taking out every atom of tin buys almost nothing.
  EXPECT_LT(study.effects[0].removedFraction, 0.01)
      << "the Sb-100 already present is not tin and does not leave with it";
  EXPECT_GT(study.effects[0].response, 0.0);
}

// An element selector takes every isotope of that element, which is what a chemical separation
// does and what a nuclide selector deliberately does not.
TEST(Intervention, AnElementSelectorTakesEveryIsotopeAndANuclideOneTakesOne) {
  // Two tin isotopes, each with its own chain, so "Sn" and "Sn-100" differ.
  StoreArrays arrays;
  arrays.provenance.version = 1;
  arrays.nuclideKey = {Zai{50, 100, 0}.key(), Zai{50, 102, 0}.key(), Zai{51, 100, 0}.key(),
                       Zai{51, 102, 0}.key()};
  arrays.halfLife = {synth::halfLifeFor(1.0e-4), synth::halfLifeFor(2.0e-4), 0.0, 0.0};
  arrays.modeOffset = {0, 1, 2, 2, 2};
  arrays.modeRtyp = {synth::kBetaMinus, synth::kBetaMinus};
  arrays.modeBranching = {1.0, 1.0};
  arrays.modeFinalState = {0, 0};
  arrays.modeIsFission = {0, 0};
  const NuclearData data = NuclearData::fromArrays(std::move(arrays));

  Inventory inv;
  inv.addKey(Zai{50, 100, 0}.key(), 1.0e20);
  inv.addKey(Zai{50, 102, 0}.key(), 1.0e20);

  const std::vector<Intervention> plan = {removing("all tin", "Sn", 1.0),
                                          removing("one isotope", "Sn-100", 1.0)};
  const InterventionStudy study = compareInterventions(data, inv, 0.0, 1.0e3, ResponseSpec{}, plan);

  ASSERT_EQ(study.effects.size(), 2u);
  EXPECT_NEAR(study.effects[0].removedFraction, 1.0, 1.0e-9) << "every isotope of tin";
  EXPECT_EQ(study.effects[0].contributors.size(), 2u);

  EXPECT_LT(study.effects[1].removedFraction, 1.0) << "only one of the two";
  EXPECT_EQ(study.effects[1].contributors.size(), 1u);
  EXPECT_EQ(study.effects[1].contributors[0].label, "Sn-100");
}

TEST(Intervention, TheSameNuclideIsReachableByEveryEquivalentSpelling) {
  const NuclearData data = twoStepChain(1.0e-4, 1.0e-5);
  const Inventory inv = seededWith(kParent, 1.0e20);

  const std::vector<Intervention> plan = {removing("by symbol", "Sn", 1.0),
                                          removing("by number", "Z=50", 1.0),
                                          removing("by chain", "A=100", 1.0)};
  const InterventionStudy study =
      compareInterventions(data, inv, 1.0e4, 1.0e5, ResponseSpec{}, plan);

  ASSERT_EQ(study.effects.size(), 3u);
  EXPECT_DOUBLE_EQ(study.effects[0].removed, study.effects[1].removed)
      << "Sn and Z=50 name the same set";
  // A=100 is the whole chain, so it takes the daughter too and is worth strictly more.
  EXPECT_GT(study.effects[2].removed, study.effects[0].removed);
}

// --- shares --------------------------------------------------------------------

TEST(Intervention, ContributorsAreOrderedByWhatTheirRemovalWasWorth) {
  const NuclearData data = twoStepChain(1.0e-3, 1.0e-6);
  const Inventory inv = seededWith(kParent, 1.0e20);

  const std::vector<Intervention> plan = {removing("the whole chain", "A=100", 1.0)};
  const InterventionStudy study =
      compareInterventions(data, inv, 1.0e4, 1.0e6, ResponseSpec{}, plan);

  ASSERT_EQ(study.effects.size(), 1u);
  const std::vector<RemovedContributor>& rows = study.effects[0].contributors;
  ASSERT_GE(rows.size(), 2u);

  double sum = 0.0;
  for (std::size_t i = 0; i < rows.size(); ++i) {
    if (i > 0) {
      EXPECT_GE(rows[i - 1].value, rows[i].value);
    }
    EXPECT_GT(rows[i].atomsRemoved, 0.0);
    sum += rows[i].fraction;
  }
  EXPECT_NEAR(sum, 1.0, 1.0e-9) << "the shares partition what the intervention removed";
}

// A selector that resolves but finds nothing left is a real answer, not an error: "there is no
// tin by then" is what someone asking about a late separation needs to hear.
TEST(Intervention, ASelectorWithNothingLeftToRemoveIsAnAnswerNotAFailure) {
  const NuclearData data = twoStepChain(1.0, 1.0e-6);  // parent gone within seconds
  const Inventory inv = seededWith(kParent, 1.0e20);

  const std::vector<Intervention> plan = {removing("late separation", "Sn", 1.0)};
  const InterventionStudy study =
      compareInterventions(data, inv, 1.0e4, 1.0e5, ResponseSpec{}, plan);

  ASSERT_EQ(study.effects.size(), 1u);
  EXPECT_NEAR(study.effects[0].removed, 0.0, study.baseline * 1.0e-9);
  EXPECT_NEAR(study.effects[0].response, study.baseline, study.baseline * 1.0e-12);
}

// --- what is refused -------------------------------------------------------------

TEST(Intervention, RemovingSomethingAfterTheResponseIsRefused) {
  const NuclearData data = twoStepChain(1.0e-4, 1.0e-5);
  const Inventory inv = seededWith(kParent, 1.0e20);
  const std::vector<Intervention> plan = {removing("too late", "Sn-100", 1.0)};

  EXPECT_THROW(compareInterventions(data, inv, 1.0e6, 1.0e4, ResponseSpec{}, plan), InputError);
  EXPECT_THROW(compareInterventions(data, inv, -1.0, 1.0e4, ResponseSpec{}, plan), InputError);
}

TEST(Intervention, AnInterventionThatCannotMeanAnythingIsRefused) {
  const NuclearData data = twoStepChain(1.0e-4, 1.0e-5);
  const Inventory inv = seededWith(kParent, 1.0e20);

  EXPECT_THROW(
      compareInterventions(data, inv, 0.0, 1.0e4, ResponseSpec{}, std::vector<Intervention>{}),
      InputError);

  EXPECT_THROW(compareInterventions(data, inv, 0.0, 1.0e4, ResponseSpec{},
                                    std::vector<Intervention>{removing("", "Sn-100", 1.0)}),
               InputError)
      << "naming the intervention that wins is the answer";

  EXPECT_THROW(compareInterventions(data, inv, 0.0, 1.0e4, ResponseSpec{},
                                    std::vector<Intervention>{removing("over", "Sn-100", 1.5)}),
               InputError);
  EXPECT_THROW(compareInterventions(data, inv, 0.0, 1.0e4, ResponseSpec{},
                                    std::vector<Intervention>{removing("under", "Sn-100", -0.1)}),
               InputError);

  Intervention empty;
  empty.name = "nothing at all";
  EXPECT_THROW(
      compareInterventions(data, inv, 0.0, 1.0e4, ResponseSpec{}, std::vector<Intervention>{empty}),
      InputError);
}

// A removal matching nothing reads as "taking this out is worth nothing" when in fact the
// question never reached the chain -- the same trap a pin naming nothing sets.
TEST(Intervention, ASelectorNamingNothingTheChainReachesIsRefused) {
  const NuclearData data = twoStepChain(1.0e-4, 1.0e-5);
  const Inventory inv = seededWith(kParent, 1.0e20);

  EXPECT_THROW(compareInterventions(data, inv, 0.0, 1.0e4, ResponseSpec{},
                                    std::vector<Intervention>{removing("absent", "Cs-137", 1.0)}),
               InputError);
  EXPECT_THROW(compareInterventions(data, inv, 0.0, 1.0e4, ResponseSpec{},
                                    std::vector<Intervention>{removing("gibberish", "wat", 1.0)}),
               InputError);
}

// Half of something taken out twice is not a stated quantity, and guessing between 75% and 100%
// would be worse than refusing.
TEST(Intervention, RemovingTheSameNuclideTwiceInOneInterventionIsRefused) {
  const NuclearData data = twoStepChain(1.0e-4, 1.0e-5);
  const Inventory inv = seededWith(kParent, 1.0e20);

  Intervention overlapping;
  overlapping.name = "double-dipped";
  overlapping.removals.push_back(Removal{"Sn", 0.5});
  overlapping.removals.push_back(Removal{"Sn-100", 0.5});

  EXPECT_THROW(compareInterventions(data, inv, 0.0, 1.0e4, ResponseSpec{},
                                    std::vector<Intervention>{overlapping}),
               InputError);
}

TEST(Intervention, AGammaLineCannotBeRemovedFromAnInventory) {
  const NuclearData data = twoStepChain(1.0e-4, 1.0e-5);
  const Inventory inv = seededWith(kParent, 1.0e20);

  ResponseSpec spec;
  spec.aggregate = Aggregate::GammaLine;
  EXPECT_THROW(compareInterventions(data, inv, 0.0, 1.0e4, spec,
                                    std::vector<Intervention>{removing("x", "Sn-100", 1.0)}),
               InputError);
}

}  // namespace
}  // namespace nusift
