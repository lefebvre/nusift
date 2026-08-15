#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "nusift/core/error.hpp"
#include "nusift/engine/decay_engine.hpp"
#include "nusift/engine/inventory.hpp"
#include "nusift/nucdata/nuclear_data.hpp"
#include "synthetic_chain.hpp"

namespace nusift {
namespace {

using synth::batemanIntegralN0;
using synth::batemanIntegralN1;
using synth::batemanN0;
using synth::batemanN1;
using synth::batemanN2;
using synth::linearChain;

constexpr double kSeedAtoms = 1.0e20;

// Index of a nuclide within a DecayResult, which lives in the PRUNED index space and so
// cannot be indexed by NuclearData's ordering.
int resultIndex(const DecayResult& result, const Zai& zai) {
  for (int i = 0; i < result.nuclideCount(); ++i) {
    if (result.nuclideKeys[static_cast<std::size_t>(i)] == zai.key()) {
      return i;
    }
  }
  return -1;
}

Inventory seedOf(const Zai& zai, double atoms = kSeedAtoms) {
  Inventory inv;
  inv.add(zai, atoms);
  return inv;
}

// n(t) = n0 e^{-lambda t}, to near machine precision. The floor under everything else: if a
// single exponential is wrong, no chain result means anything.
TEST(DecayEngine, SingleNuclideIsExponential) {
  const double lambda = 1.0e-3;
  const NuclearData data = NuclearData::fromArrays(linearChain({lambda}));
  const Zai parent{50, 100, 0};

  const std::vector<double> times = {0.0, 100.0, 1000.0, 5000.0};
  const DecayResult result = decay(data, seedOf(parent), times);
  const int i = resultIndex(result, parent);
  ASSERT_GE(i, 0);

  for (int k = 0; k < result.timeCount(); ++k) {
    const double expected = batemanN0(kSeedAtoms, lambda, times[static_cast<std::size_t>(k)]);
    EXPECT_NEAR(result.atomsAt(k)[i], expected, expected * 1e-12)
        << "at t = " << times[static_cast<std::size_t>(k)];
  }
}

// Two-member chain against the closed form. The daughter's transient -- rising, peaking,
// then decaying -- is where an error in the off-diagonal production term shows up.
TEST(DecayEngine, TwoMemberChainMatchesBateman) {
  const double l0 = 1.0e-3;
  const double l1 = 4.0e-4;
  const NuclearData data = NuclearData::fromArrays(linearChain({l0, l1}));
  const Zai parent{50, 100, 0};
  const Zai daughter{51, 100, 0};

  const std::vector<double> times = {10.0, 500.0, 1500.0, 6000.0};
  const DecayResult result = decay(data, seedOf(parent), times);
  const int ip = resultIndex(result, parent);
  const int id = resultIndex(result, daughter);
  ASSERT_GE(ip, 0);
  ASSERT_GE(id, 0);

  for (int k = 0; k < result.timeCount(); ++k) {
    const double t = times[static_cast<std::size_t>(k)];
    EXPECT_NEAR(result.atomsAt(k)[ip], batemanN0(kSeedAtoms, l0, t),
                batemanN0(kSeedAtoms, l0, t) * 1e-11);
    const double expected = batemanN1(kSeedAtoms, l0, l1, t);
    EXPECT_NEAR(result.atomsAt(k)[id], expected, expected * 1e-10) << "daughter at t = " << t;
  }
}

TEST(DecayEngine, ThreeMemberChainMatchesBateman) {
  const double l0 = 2.0e-3;
  const double l1 = 7.0e-4;
  const double l2 = 3.0e-4;
  const NuclearData data = NuclearData::fromArrays(linearChain({l0, l1, l2}));
  const Zai grandchild{52, 100, 0};

  const std::vector<double> times = {100.0, 1000.0, 4000.0};
  const DecayResult result = decay(data, seedOf({50, 100, 0}), times);
  const int i = resultIndex(result, grandchild);
  ASSERT_GE(i, 0);

  for (int k = 0; k < result.timeCount(); ++k) {
    const double t = times[static_cast<std::size_t>(k)];
    const double expected = batemanN2(kSeedAtoms, l0, l1, l2, t);
    EXPECT_NEAR(result.atomsAt(k)[i], expected, expected * 1e-9) << "at t = " << t;
  }
}

// THE augmented-matrix test. The bottom block of exp([[A,0],[I,0]] t) is claimed to be the
// exact time integral of the inventory; this is the only independent check that it is, and
// there is no substitute for it. Both the parent and the daughter are checked, because the
// parent alone would pass even if the coupling between the blocks were wrong.
TEST(DecayEngine, TimeIntegralMatchesAnalyticBateman) {
  const double l0 = 1.5e-3;
  const double l1 = 5.0e-4;
  const NuclearData data = NuclearData::fromArrays(linearChain({l0, l1}));
  const Zai parent{50, 100, 0};
  const Zai daughter{51, 100, 0};

  const std::vector<double> times = {50.0, 400.0, 2000.0, 9000.0};
  const DecayResult result = decay(data, seedOf(parent), times);
  const int ip = resultIndex(result, parent);
  const int id = resultIndex(result, daughter);

  for (int k = 0; k < result.timeCount(); ++k) {
    const double t = times[static_cast<std::size_t>(k)];
    const double expectedParent = batemanIntegralN0(kSeedAtoms, l0, t);
    EXPECT_NEAR(result.integratedAtomsAt(k)[ip], expectedParent, expectedParent * 1e-10)
        << "parent integral at t = " << t;
    const double expectedDaughter = batemanIntegralN1(kSeedAtoms, l0, l1, t);
    EXPECT_NEAR(result.integratedAtomsAt(k)[id], expectedDaughter, expectedDaughter * 1e-9)
        << "daughter integral at t = " << t;
  }
}

// For a stable nuclide the augmented block is [[0,0],[1,0]] -- a nilpotent Jordan block --
// so G(t) = n0 * t exactly. Passing this asserts that CRAM's rational approximation r(z)
// satisfies both r(0) = 1 and r'(0) = 1. Nothing else in the suite pins the derivative
// condition, and a scheme that got it wrong would still look plausible on decaying species.
TEST(DecayEngine, StableNuclideIntegralIsLinearInTime) {
  StoreArrays arrays;
  arrays.provenance.version = 1;
  arrays.nuclideKey = {Zai{50, 100, 0}.key()};
  arrays.halfLife = {0.0};  // stable
  arrays.modeOffset = {0, 0};
  const NuclearData data = NuclearData::fromArrays(arrays);

  const std::vector<double> times = {1.0, 100.0, 1.0e6};
  const DecayResult result = decay(data, seedOf({50, 100, 0}), times);
  ASSERT_EQ(result.nuclideCount(), 1);

  for (int k = 0; k < result.timeCount(); ++k) {
    const double t = times[static_cast<std::size_t>(k)];
    EXPECT_NEAR(result.atomsAt(k)[0], kSeedAtoms, kSeedAtoms * 1e-13) << "atoms at t = " << t;
    EXPECT_NEAR(result.integratedAtomsAt(k)[0], kSeedAtoms * t, kSeedAtoms * t * 1e-12)
        << "integral at t = " << t;
  }
}

// An interval integral from zero must agree with the cumulative one, and a general interval
// must agree with the difference of the analytic cumulative integrals.
TEST(DecayEngine, IntervalIntegralMatchesAnalytic) {
  const double l0 = 8.0e-4;
  const NuclearData data = NuclearData::fromArrays(linearChain({l0}));
  const Zai parent{50, 100, 0};

  std::vector<std::int64_t> keys;
  const std::vector<double> integral = intervalIntegral(data, seedOf(parent), 500.0, 3000.0, &keys);
  ASSERT_FALSE(keys.empty());
  ASSERT_EQ(keys[0], parent.key());

  const double expected =
      batemanIntegralN0(kSeedAtoms, l0, 3000.0) - batemanIntegralN0(kSeedAtoms, l0, 500.0);
  EXPECT_NEAR(integral[0], expected, expected * 1e-9);
}

// Adjacent intervals must sum to the whole. This catches a restart that carries the wrong
// inventory into the second leg -- an error the cancellation guard could plausibly introduce
// and which no single-interval test would notice.
TEST(DecayEngine, AdjacentIntervalsAreAdditive) {
  const NuclearData data = NuclearData::fromArrays(linearChain({1.2e-3, 6.0e-4}));
  const Inventory seed = seedOf({50, 100, 0});

  const std::vector<double> first = intervalIntegral(data, seed, 0.0, 800.0, nullptr);
  const std::vector<double> second = intervalIntegral(data, seed, 800.0, 2500.0, nullptr);
  const std::vector<double> whole = intervalIntegral(data, seed, 0.0, 2500.0, nullptr);

  ASSERT_EQ(first.size(), whole.size());
  for (std::size_t i = 0; i < whole.size(); ++i) {
    if (whole[i] <= 0.0) {
      continue;
    }
    EXPECT_NEAR(first[i] + second[i], whole[i], whole[i] * 1e-9) << "nuclide index " << i;
  }
}

// The regime the cancellation guard exists for: a one-second window a half-life into the
// decay. G(t1) and G(t2) are both around 6e28 and differ by about 4e19, so the subtraction
// throws away roughly ten of the sixteen available digits.
//
// The reference here is deliberately NOT G(t2) - G(t1). Written that way it is the very
// cancellation under test and cannot serve as its own control -- the first draft of this
// test did exactly that and reported an expected value of zero. The stable form is
//
//     \int_{t1}^{t2} n0 e^{-lambda tau} dtau = (n0/lambda) e^{-lambda t1} (1 - e^{-lambda dt})
//
// with expm1 for the second factor, which is accurate when lambda*dt is tiny.
TEST(DecayEngine, NarrowLateIntervalSurvivesCancellation) {
  const double lambda = 1.0e-9;  // ~22 y half-life
  const NuclearData data = NuclearData::fromArrays(linearChain({lambda}));
  const Inventory seed = seedOf({50, 100, 0});

  const double t1 = 9.46e8;  // ~30 y, a little over one half-life
  const double dt = 1.0;
  const double t2 = t1 + dt;

  const double expected = kSeedAtoms / lambda * std::exp(-lambda * t1) * -std::expm1(-lambda * dt);
  ASSERT_GT(expected, 0.0);

  const std::vector<double> integral = intervalIntegral(data, seed, t1, t2, nullptr);
  EXPECT_NEAR(integral[0], expected, expected * 1e-9);

  // What the unguarded subtraction would have produced, computed from the same public API.
  // It is not garbage at this separation -- it is merely far less accurate -- and asserting
  // the guarded path beats it is what shows the guard is earning its extra solve.
  const DecayResult cumulative = decay(data, seed, std::vector<double>{t1, t2});
  const double naive = cumulative.integratedAtomsAt(1)[0] - cumulative.integratedAtomsAt(0)[0];
  EXPECT_LT(std::abs(integral[0] - expected), std::abs(naive - expected));
}

// Branching must split production in proportion and lose nothing: with two stable daughters
// every atom that leaves the parent has to arrive somewhere.
TEST(DecayEngine, BranchingSplitsInProportionAndConservesAtoms) {
  const double lambda = 1.0e-3;
  const NuclearData data = NuclearData::fromArrays(synth::branchingChain(lambda, 0.7, 0.3));
  const Zai parent{50, 100, 0};
  const Zai up{51, 100, 0};
  const Zai down{49, 100, 0};

  const std::vector<double> times = {5000.0};
  const DecayResult result = decay(data, seedOf(parent), times);
  const int ip = resultIndex(result, parent);
  const int iu = resultIndex(result, up);
  const int id = resultIndex(result, down);
  ASSERT_GE(ip, 0);
  ASSERT_GE(iu, 0);
  ASSERT_GE(id, 0);

  const auto atoms = result.atomsAt(0);
  const double decayed = kSeedAtoms - atoms[ip];
  EXPECT_NEAR(atoms[iu], 0.7 * decayed, decayed * 1e-10);
  EXPECT_NEAR(atoms[id], 0.3 * decayed, decayed * 1e-10);
  EXPECT_NEAR(atoms[ip] + atoms[iu] + atoms[id], kSeedAtoms, kSeedAtoms * 1e-12);
}

// Pruning is claimed to be exact, not an approximation. Solving the full chain and the
// pruned one must agree on every shared nuclide -- otherwise the reachable set is not
// actually closed under production and atoms are being dropped.
TEST(DecayEngine, PrunedAndUnprunedSolvesAgree) {
  // Seed only the middle of a longer chain, so there are unreachable nuclides both before
  // and after the seed's forward closure.
  StoreArrays arrays = linearChain({1.0e-3, 5.0e-4, 2.0e-4});
  const NuclearData data = NuclearData::fromArrays(arrays);
  const Inventory seed = seedOf({51, 100, 0});  // second member, not the first

  const std::vector<double> times = {1000.0, 8000.0};
  DecayOptions pruned;
  pruned.prune = true;
  DecayOptions full;
  full.prune = false;

  const DecayResult a = decay(data, seed, times, pruned);
  const DecayResult b = decay(data, seed, times, full);

  // Pruning must actually have removed something, or the test proves nothing.
  EXPECT_LT(a.nuclideCount(), b.nuclideCount());

  for (int k = 0; k < a.timeCount(); ++k) {
    for (int i = 0; i < a.nuclideCount(); ++i) {
      const int j = resultIndex(b, Zai::fromKey(a.nuclideKeys[static_cast<std::size_t>(i)]));
      ASSERT_GE(j, 0);
      const double expected = b.atomsAt(k)[j];
      EXPECT_NEAR(a.atomsAt(k)[i], expected, std::abs(expected) * 1e-12 + 1e-6);
      const double expectedIntegral = b.integratedAtomsAt(k)[j];
      EXPECT_NEAR(a.integratedAtomsAt(k)[i], expectedIntegral,
                  std::abs(expectedIntegral) * 1e-12 + 1e-6);
    }
  }
}

// The chain in the store is not closed -- its last member's daughter is absent -- so closure
// has to register it, or the atoms decaying into it would vanish and the total would drift.
TEST(DecayEngine, ClosureRegistersUnstagedDaughtersSoAtomsAreConserved) {
  StoreArrays arrays;
  arrays.provenance.version = 1;
  arrays.nuclideKey = {Zai{50, 100, 0}.key()};  // decays, but its daughter is not staged
  arrays.halfLife = {synth::halfLifeFor(1.0e-3)};
  arrays.modeOffset = {0, 1};
  arrays.modeRtyp = {synth::kBetaMinus};
  arrays.modeBranching = {1.0};
  arrays.modeFinalState = {0};
  arrays.modeIsFission = {0};

  const NuclearData data = NuclearData::fromArrays(arrays);
  EXPECT_EQ(data.size(), 2) << "closure should have added the beta-minus daughter";

  const DecayResult result = decay(data, seedOf({50, 100, 0}), std::vector<double>{4000.0});
  double total = 0.0;
  for (const double atoms : result.atomsAt(0)) {
    total += atoms;
  }
  EXPECT_NEAR(total, kSeedAtoms, kSeedAtoms * 1e-12);
}

TEST(DecayEngine, RejectsUnsortedOrNegativeTimes) {
  const NuclearData data = NuclearData::fromArrays(linearChain({1.0e-3}));
  const Inventory seed = seedOf({50, 100, 0});

  EXPECT_THROW(decay(data, seed, std::vector<double>{100.0, 50.0}), InputError);
  EXPECT_THROW(decay(data, seed, std::vector<double>{-1.0}), InputError);
  EXPECT_THROW(decay(data, seed, std::vector<double>{}), InputError);
  EXPECT_THROW(decay(data, seed, std::vector<double>{10.0, 10.0}), InputError);
}

// A nuclide with no evaluated data cannot be decayed, and saying so is much better than
// dropping it: a silently ignored row understates every ranking with no indication why.
TEST(DecayEngine, RejectsInventoryNuclideAbsentFromTheStore) {
  const NuclearData data = NuclearData::fromArrays(linearChain({1.0e-3}));
  Inventory inv;
  inv.add(Zai{96, 244, 0}, 1.0e10);  // not in this chain
  EXPECT_THROW(decay(data, inv, std::vector<double>{100.0}), InputError);
}

TEST(DecayEngine, RejectsBackwardsInterval) {
  const NuclearData data = NuclearData::fromArrays(linearChain({1.0e-3}));
  const Inventory seed = seedOf({50, 100, 0});
  EXPECT_THROW(intervalIntegral(data, seed, 500.0, 500.0, nullptr), InputError);
  EXPECT_THROW(intervalIntegral(data, seed, 900.0, 100.0, nullptr), InputError);
}

// CRAM16 is offered for screening. It should agree with CRAM48 to well within any tolerance
// a screening pass cares about -- if it did not, offering it would be a trap.
TEST(DecayEngine, Cram16AgreesWithCram48ToScreeningAccuracy) {
  const NuclearData data = NuclearData::fromArrays(linearChain({1.0e-3, 4.0e-4}));
  const Inventory seed = seedOf({50, 100, 0});
  const std::vector<double> times = {2000.0};

  DecayOptions fast;
  fast.order = CramOrder::Order16;
  const DecayResult a = decay(data, seed, times, fast);
  const DecayResult b = decay(data, seed, times);  // Order48 default

  for (int i = 0; i < a.nuclideCount(); ++i) {
    const double expected = b.atomsAt(0)[i];
    if (expected > kSeedAtoms * 1e-12) {
      EXPECT_NEAR(a.atomsAt(0)[i], expected, expected * 1e-8);
    }
  }
}

// Threading must not change a single bit. Each time is solved independently and written to its
// own slice, so there is nothing to race on and nothing to accumulate in a different order --
// but that is an argument, and this is the check. A parallel result that differed even in the
// last place would make every golden test depend on the core count of whoever ran it.
TEST(DecayEngine, ThreadingIsBitwiseIdenticalToSerial) {
  const NuclearData data = NuclearData::fromArrays(linearChain({1.0e-3, 5.0e-4, 2.0e-4}));
  const Inventory seed = seedOf({50, 100, 0});

  std::vector<double> times;
  for (int k = 0; k < 40; ++k) {
    times.push_back(10.0 * std::pow(1.3, k));
  }

  DecayOptions serial;
  serial.threads = 1;
  DecayOptions parallel;
  parallel.threads = 8;

  const DecayResult a = decay(data, seed, times, serial);
  const DecayResult b = decay(data, seed, times, parallel);

  ASSERT_EQ(a.atoms.size(), b.atoms.size());
  for (std::size_t i = 0; i < a.atoms.size(); ++i) {
    EXPECT_DOUBLE_EQ(a.atoms[i], b.atoms[i]) << "atoms differ at flat index " << i;
    EXPECT_DOUBLE_EQ(a.integratedAtoms[i], b.integratedAtoms[i])
        << "integrals differ at flat index " << i;
  }
}

// More workers than times must not spawn idle threads or, worse, leave a slice unwritten.
TEST(DecayEngine, MoreThreadsThanTimesIsHarmless) {
  const NuclearData data = NuclearData::fromArrays(linearChain({1.0e-3}));
  const Inventory seed = seedOf({50, 100, 0});
  const std::vector<double> times = {100.0, 200.0};

  DecayOptions many;
  many.threads = 64;
  const DecayResult result = decay(data, seed, times, many);

  ASSERT_EQ(result.timeCount(), 2);
  for (int k = 0; k < result.timeCount(); ++k) {
    EXPECT_GT(result.atomsAt(k)[0], 0.0) << "time slice " << k << " was never written";
  }
}

// An error raised inside a worker has to reach the caller. Letting it escape a thread would
// call std::terminate and lose the message.
TEST(DecayEngine, AnErrorInAWorkerReachesTheCaller) {
  const NuclearData data = NuclearData::fromArrays(linearChain({1.0e-3}));
  Inventory inv;
  inv.add(Zai{96, 244, 0}, 1.0e10);  // not in this chain

  DecayOptions parallel;
  parallel.threads = 4;
  // The seed is validated before any thread starts, so this is the ordinary path -- what
  // matters is that turning threading on does not turn a clean error into a crash.
  EXPECT_THROW(decay(data, inv, std::vector<double>{100.0, 200.0}, parallel), InputError);
}

}  // namespace
}  // namespace nusift
