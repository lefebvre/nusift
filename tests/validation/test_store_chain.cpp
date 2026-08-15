#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

#include "nusift/core/nuclide_name.hpp"
#include "nusift/engine/decay_engine.hpp"
#include "nusift/engine/inventory.hpp"
#include "nusift/nucdata/nuclear_data.hpp"
#include "nusift/units.hpp"
#include "synthetic_chain.hpp"
#include "validation_store.hpp"

namespace nusift::validation {
namespace {

constexpr double kSeedAtoms = 1.0e20;

Zai zaiOf(const char* name) {
  return requireNuclideName(name);
}

double lambdaOf(const NuclearData& data, const char* name) {
  const int i = data.indexOf(zaiOf(name));
  return i >= 0 ? data.decayConstant(i) : 0.0;
}

// The result's index space is the PRUNED chain, not the store's, so a lookup has to go through
// the key rather than reusing a NuclearData index.
int resultIndexOf(const DecayResult& result, const char* name) {
  const std::int64_t key = zaiOf(name).key();
  const auto it = std::find(result.nuclideKeys.begin(), result.nuclideKeys.end(), key);
  return it == result.nuclideKeys.end()
             ? -1
             : static_cast<int>(std::distance(result.nuclideKeys.begin(), it));
}

DecayResult decaySingle(const char* parent, const std::vector<double>& times) {
  Inventory inventory;
  inventory.add(zaiOf(parent), kSeedAtoms);
  return decay(committedStore(), inventory, times);
}

// --- the engine on real chains against their closed forms -------------------
//
// The unit suite proves the solver against Bateman on synthetic chains, where every decay
// constant is chosen and every branch is 1.0. That leaves one thing unproven: whether the
// chain BUILT FROM THE SHIPPED STORE is the chain the physics describes. A branching fraction
// staged against the wrong daughter, a decay mode parsed as the wrong RTYP, or a half-life
// attached to the wrong isomer all produce a perfectly well-conditioned solve of the wrong
// matrix, and only a comparison on a real, named chain catches it.
//
// The closed forms come from synthetic_chain.hpp and are evaluated with the store's own decay
// constants, so what is being tested is the chain topology and the solver together, against
// arithmetic that knows nothing about either.

// Sr-90 -> Y-90 -> Zr-90 is the cleanest two-step chain in the fission-product inventory: both
// steps are beta-minus at unit branching, so the unmodified two-member Bateman solution applies
// with no branching factor at all.
TEST(StoreChain, Strontium90AndYttrium90FollowBateman) {
  const double lambdaSr = lambdaOf(committedStore(), "Sr-90");
  const double lambdaY = lambdaOf(committedStore(), "Y-90");
  ASSERT_GT(lambdaSr, 0.0);
  ASSERT_GT(lambdaY, 0.0);

  const std::vector<double> times = {86400.0, 30.0 * 86400.0, units::kSecondsPerYear,
                                     10.0 * units::kSecondsPerYear};
  const DecayResult result = decaySingle("Sr-90", times);
  const int iSr = resultIndexOf(result, "Sr-90");
  const int iY = resultIndexOf(result, "Y-90");
  ASSERT_GE(iSr, 0);
  ASSERT_GE(iY, 0);

  for (int k = 0; k < result.timeCount(); ++k) {
    const double t = times[static_cast<std::size_t>(k)];
    const double parent = synth::batemanN0(kSeedAtoms, lambdaSr, t);
    const double daughter = synth::batemanN1(kSeedAtoms, lambdaSr, lambdaY, t);
    EXPECT_NEAR(result.atomsAt(k)[iSr], parent, parent * 1e-9) << "Sr-90 at t = " << t;
    EXPECT_NEAR(result.atomsAt(k)[iY], daughter, daughter * 1e-9) << "Y-90 at t = " << t;
  }
}

// Secular equilibrium: with a parent 4000 times longer-lived than its daughter, the daughter's
// activity climbs to the parent's and stays there. This is the textbook statement, and it is
// the behaviour a Sr-90/Y-90 source actually has.
TEST(StoreChain, Yttrium90ReachesSecularEquilibriumWithItsParent) {
  const std::vector<double> times = {30.0 * 86400.0};
  const DecayResult result = decaySingle("Sr-90", times);
  const int iSr = resultIndexOf(result, "Sr-90");
  const int iY = resultIndexOf(result, "Y-90");
  ASSERT_GE(iSr, 0);
  ASSERT_GE(iY, 0);

  const double activitySr = lambdaOf(committedStore(), "Sr-90") * result.atomsAt(0)[iSr];
  const double activityY = lambdaOf(committedStore(), "Y-90") * result.atomsAt(0)[iY];
  EXPECT_NEAR(activityY / activitySr, 1.0, 0.005) << "ratio " << activityY / activitySr;
}

// Mo-99 -> Tc-99m is TRANSIENT equilibrium, and it is the one that exposes a branching error.
// Only 87.6% of Mo-99 decays populate the metastable state; the rest go straight to Tc-99. If
// the staged branch were 1.0 the Bateman curve below would be too high by a seventh, and no
// half-life check anywhere would notice.
TEST(StoreChain, Molybdenum99FeedsTechnetium99mThroughItsStagedBranch) {
  const double lambdaMo = lambdaOf(committedStore(), "Mo-99");
  const double lambdaTc = lambdaOf(committedStore(), "Tc-99m");
  ASSERT_GT(lambdaMo, 0.0);
  ASSERT_GT(lambdaTc, 0.0);

  const std::vector<double> times = {6.0 * 3600.0, 24.0 * 3600.0, 3.0 * 86400.0};
  const DecayResult result = decaySingle("Mo-99", times);
  const int iTc = resultIndexOf(result, "Tc-99m");
  ASSERT_GE(iTc, 0);

  // Recover the branch the store staged rather than asserting a literal, then check it against
  // the evaluated value. Doing it this way says which of the two is wrong when it fails.
  const double unbranched = synth::batemanN1(kSeedAtoms, lambdaMo, lambdaTc, times[0]);
  const double branch = result.atomsAt(0)[iTc] / unbranched;
  EXPECT_NEAR(branch, 0.876, 0.876 * 0.02) << "staged branch to the isomer is " << branch;

  for (int k = 0; k < result.timeCount(); ++k) {
    const double expected = branch * synth::batemanN1(kSeedAtoms, lambdaMo, lambdaTc,
                                                      times[static_cast<std::size_t>(k)]);
    EXPECT_NEAR(result.atomsAt(k)[iTc], expected, expected * 1e-9)
        << "Tc-99m at t = " << times[static_cast<std::size_t>(k)];
  }
}

// The equilibrium ratio the published Cs-137 gamma constant falls out of. Once Ba-137m's
// 2.55-minute transient has passed, its activity is the branch times its parent's, and that
// number -- not a folded-in table constant -- is what makes a Cs-137 source's exposure come
// out right. docs/exposure.md section 5 explains why the alternative is a live bug.
TEST(StoreChain, Barium137mHoldsTheStagedBranchOfItsParentsActivity) {
  const std::vector<double> times = {86400.0, units::kSecondsPerYear};
  const DecayResult result = decaySingle("Cs-137", times);
  const int iCs = resultIndexOf(result, "Cs-137");
  const int iBa = resultIndexOf(result, "Ba-137m");
  ASSERT_GE(iCs, 0);
  ASSERT_GE(iBa, 0);

  for (int k = 0; k < result.timeCount(); ++k) {
    const double activityCs = lambdaOf(committedStore(), "Cs-137") * result.atomsAt(k)[iCs];
    const double activityBa = lambdaOf(committedStore(), "Ba-137m") * result.atomsAt(k)[iBa];
    EXPECT_NEAR(activityBa / activityCs, 0.947, 0.947 * 0.01)
        << "ratio " << activityBa / activityCs << " at t = " << times[static_cast<std::size_t>(k)];
  }
}

// Beta decay moves a nuclide along a mass chain without changing A, so the total atom count on
// a pure beta chain is conserved, forever. Nothing about this depends on the decay constants,
// which is what makes it a check on the chain's TOPOLOGY -- production leaking into the wrong
// mass number, or a daughter registered but never produced into, breaks it while every
// individual Bateman curve above still passes.
//
// The tolerances differ per chain, and the reason is worth stating because it is a property of
// the evaluated data rather than of the solver. A chain conserves atoms exactly only if its
// branching fractions sum to exactly 1. Sr-90 and Y-90 each have a single mode at unit
// branching, so that chain conserves to round-off. Cs-137's two modes are staged as 0.05300549
// and 0.9469945, which sum to 0.99999999 -- ENDF quotes them to eight decimals and they do not
// quite close. NuSIFT stages what the evaluation says, so exactly 1e-8 of every Cs-137 decay
// goes nowhere, and by 300 years the chain is short by that same 1e-8.
//
// This was measured, not assumed: the deficit is identical at CRAM order 16 and 48, which rules
// out the approximation, and it tracks the decayed fraction exactly. Across the shipped store
// 226 nuclides have branchings that miss unity by more than 1e-12, and none miss by more than
// 1e-6, so 2e-8 here bounds the effect for these two chains while staying far below anything a
// real topology error would produce.
TEST(StoreChain, PureBetaMassChainsConserveTheirAtoms) {
  struct Case {
    const char* parent;
    double tolerance;
  };
  const Case cases[] = {
      {"Sr-90", 1e-10},  // unit branchings throughout
      {"Cs-137", 2e-8},  // staged branchings sum to 1 - 1e-8
  };

  for (const Case& c : cases) {
    const std::vector<double> times = {0.0, 86400.0, 100.0 * units::kSecondsPerYear};
    const DecayResult result = decaySingle(c.parent, times);
    const int massNumber = zaiOf(c.parent).a;

    for (int k = 0; k < result.timeCount(); ++k) {
      double total = 0.0;
      for (int i = 0; i < result.nuclideCount(); ++i) {
        const Zai zai = Zai::fromKey(result.nuclideKeys[static_cast<std::size_t>(i)]);
        if (zai.a == massNumber) {
          total += result.atomsAt(k)[i];
        }
      }
      EXPECT_NEAR(total, kSeedAtoms, kSeedAtoms * c.tolerance)
          << c.parent << " chain at t = " << times[static_cast<std::size_t>(k)];
    }
  }
}

}  // namespace
}  // namespace nusift::validation
