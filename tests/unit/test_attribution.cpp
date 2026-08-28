#include <gtest/gtest.h>

#include <cmath>
#include <numeric>
#include <string>
#include <vector>

#include "nusift/core/error.hpp"
#include "nusift/engine/adjoint_engine.hpp"
#include "nusift/engine/decay_engine.hpp"
#include "nusift/nucdata/nuclear_data.hpp"
#include "nusift/triage/attribution.hpp"
#include "nusift/triage/ranking.hpp"
#include "nusift/triage/response.hpp"
#include "synthetic_chain.hpp"

namespace nusift {
namespace {

// Z=50 -> 51 -> 52 at mass 100, the last stable. Rates chosen an order apart so the daughter
// both grows and decays over the window the tests use.
constexpr double kL0 = 1.0e-3;
constexpr double kL1 = 4.0e-4;

NuclearData chain() {
  return NuclearData::fromArrays(synth::linearChain({kL0, kL1}));
}

Zai parent() {
  return Zai{50, 100, 0};
}
Zai daughter() {
  return Zai{51, 100, 0};
}
// The chain's terminator. Stable, so an activity metric weights it zero however many atoms of
// it are seeded -- which is the case the inert-seed tests below are about.
Zai terminator() {
  return Zai{52, 100, 0};
}

// The same chain with photon lines on the daughter and an unmodelled continuum beside them, so
// an exposure metric has something to answer with AND something to warn about.
NuclearData emittingChain() {
  StoreArrays s = synth::linearChain({kL0, kL1});
  synth::addLines(s, 1, {6.616e5}, {0.85});
  // A third of the daughter's photon energy left out of the model, comfortably past the 5% the
  // flag trips at. Without it the caveat test would pass on a store that has nothing to say.
  // linearChain leaves the column empty -- it is optional per nuclide -- so it is sized here.
  s.continuumPhotonEv.assign(static_cast<std::size_t>(s.nuclideCount()), 0.0);
  s.continuumPhotonEv[1] = 0.5 * 6.616e5 * 0.85;
  return NuclearData::fromArrays(std::move(s));
}

Inventory seedOf(const std::vector<std::pair<Zai, double>>& rows) {
  Inventory inv;
  for (const auto& [zai, atoms] : rows) {
    inv.add(zai, atoms);
  }
  inv.setProvenance("test seed");
  return inv;
}

double sumOfShares(const SeedAttribution& a) {
  double total = 0.0;
  for (const SeedShare& s : a.shares) {
    total += s.value;
  }
  return total;
}

// THE invariant. Decay is linear, so every seed atom's contribution is exactly its importance
// times its count and the shares partition the response -- they do not approximate it. If this
// ever fails, seed attribution is a sensitivity estimate rather than a decomposition and the
// coveredFraction it reports is meaningless.
TEST(SeedAttribution, SharesPartitionTheTotal) {
  const NuclearData data = chain();
  const Inventory inv = seedOf({{parent(), 1.0e20}, {daughter(), 3.0e19}});

  const SeedAttribution a =
      attributeToSeed(data, inv, 5000.0, ResponseSpec{}, RankRequest{.topN = 0});

  ASSERT_EQ(a.shares.size(), 2u);
  EXPECT_NEAR(sumOfShares(a), a.total, std::abs(a.total) * 1e-12);
  EXPECT_NEAR(a.coveredFraction, 1.0, 1e-12);
  EXPECT_EQ(a.omittedCount, 0);
}

// The two attributions are of ONE number. A forward decay and ranking, and a backward adjoint,
// must agree on the total or one of them is wrong about what the response even is.
TEST(SeedAttribution, TotalAgreesWithTheForwardRanking) {
  const NuclearData data = chain();
  const Inventory inv = seedOf({{parent(), 1.0e20}, {daughter(), 3.0e19}});
  const double t = 5000.0;

  const DecayResult decayed = decay(data, inv, std::vector<double>{t});
  const ResponseTable table = buildResponse(data, decayed, ResponseSpec{});
  const Ranking forward = rank(table, 0, RankRequest{.topN = 0});

  const SeedAttribution a = attributeToSeed(data, inv, t, ResponseSpec{}, RankRequest{.topN = 0});

  EXPECT_NEAR(a.total, forward.total, std::abs(forward.total) * 1e-10);
}

// A seed of the parent alone leaves the daughter's activity attributable only to the parent --
// which is the whole point of the second axis. The ordinary ranking at this time reports both
// nuclides; this reports one, because only one was seeded.
TEST(SeedAttribution, AttributesADaughtersActivityToTheSeededParent) {
  const NuclearData data = chain();
  const Inventory inv = seedOf({{parent(), 1.0e20}});
  const double t = 5000.0;

  const SeedAttribution a = attributeToSeed(data, inv, t, ResponseSpec{}, RankRequest{.topN = 0});
  ASSERT_EQ(a.shares.size(), 1u);
  EXPECT_EQ(a.shares[0].key, parent().key());
  EXPECT_NEAR(a.shares[0].fraction, 1.0, 1e-12);

  // Meanwhile the forward ranking carries the daughter as a contributor in its own right.
  const DecayResult decayed = decay(data, inv, std::vector<double>{t});
  const Ranking forward = rank(buildResponse(data, decayed, ResponseSpec{}), 0, RankRequest{});
  ASSERT_EQ(forward.contributors.size(), 2u);
}

// Checked against Bateman rather than against another solve: the importance of the parent is
// d/dn0 of (l0*n0(t) + l1*n1(t)), and both terms are closed-form for a two-step chain.
TEST(SeedAttribution, ImportanceMatchesTheClosedFormDerivative) {
  const NuclearData data = chain();
  const double t = 5000.0;
  const double n0 = 1.0e20;
  const Inventory inv = seedOf({{parent(), n0}});

  const SeedAttribution a = attributeToSeed(data, inv, t, ResponseSpec{}, RankRequest{.topN = 0});
  ASSERT_EQ(a.shares.size(), 1u);

  // R(n0) is linear in n0, so dR/dn0 is R evaluated at a unit seed.
  const double expected =
      kL0 * synth::batemanN0(1.0, kL0, t) + kL1 * synth::batemanN1(1.0, kL0, kL1, t);
  EXPECT_NEAR(a.shares[0].importance, expected, expected * 1e-9);
  EXPECT_NEAR(a.shares[0].value, expected * n0, expected * n0 * 1e-9);
}

// exp(A^T * 0) is the identity, so at t = 0 the importance of a nuclide is just its own weight
// -- for activity, its decay constant. A shortcut worth pinning: it is the one path that never
// reaches the solver.
TEST(SeedAttribution, ImportanceAtTimeZeroIsTheWeightItself) {
  const NuclearData data = chain();
  const Inventory inv = seedOf({{parent(), 1.0e20}, {daughter(), 3.0e19}});

  const SeedAttribution a = attributeToSeed(data, inv, 0.0, ResponseSpec{}, RankRequest{.topN = 0});
  ASSERT_EQ(a.shares.size(), 2u);
  for (const SeedShare& s : a.shares) {
    const int index = data.indexOf(Zai::fromKey(s.key));
    EXPECT_NEAR(s.importance, data.decayConstant(index), data.decayConstant(index) * 1e-12);
  }
}

// A pin reaches past the cut carrying the rank it really holds, and does not reorder anything
// above it -- the same contract RankRequest::pinned has for an ordinary ranking.
TEST(SeedAttribution, PinsReachPastTheCut) {
  const NuclearData data = chain();
  // The daughter seeded far smaller than the parent, so it cannot make a top-1 cut.
  const Inventory inv = seedOf({{parent(), 1.0e20}, {daughter(), 1.0e12}});

  const SeedAttribution a = attributeToSeed(data, inv, 5000.0, ResponseSpec{},
                                            RankRequest{.topN = 1, .pinned = {daughter().key()}});

  ASSERT_EQ(a.shares.size(), 2u);
  EXPECT_EQ(a.shares[0].key, parent().key());
  EXPECT_FALSE(a.shares[0].pinned);
  EXPECT_EQ(a.shares[1].key, daughter().key());
  EXPECT_TRUE(a.shares[1].pinned);
  EXPECT_EQ(a.shares[1].rank, 2);
  EXPECT_EQ(a.omittedCount, 0);
}

// Truncation must never let the report look complete: coveredFraction is over what was
// returned, and the omitted count says how much was not.
TEST(SeedAttribution, TruncationReportsWhatItLeftOut) {
  const NuclearData data = chain();
  const Inventory inv = seedOf({{parent(), 1.0e20}, {daughter(), 1.0e19}});

  const SeedAttribution a =
      attributeToSeed(data, inv, 5000.0, ResponseSpec{}, RankRequest{.topN = 1});
  ASSERT_EQ(a.shares.size(), 1u);
  EXPECT_EQ(a.omittedCount, 1);
  EXPECT_LT(a.coveredFraction, 1.0);
  EXPECT_GT(a.coveredFraction, 0.0);
}

TEST(SeedAttribution, RefusesADomainOrAggregateItCannotAnswer) {
  const NuclearData data = chain();
  const Inventory inv = seedOf({{parent(), 1.0e20}});

  ResponseSpec interval;
  interval.unit = Unit::Decays;  // a total, not an instantaneous rate
  EXPECT_THROW(attributeToSeed(data, inv, 5000.0, interval, RankRequest{}), InputError);

  ResponseSpec lines;
  lines.aggregate = Aggregate::GammaLine;
  EXPECT_THROW(attributeToSeed(data, inv, 5000.0, lines, RankRequest{}), InputError);

  ResponseSpec chains;
  chains.aggregate = Aggregate::MassChain;
  EXPECT_THROW(attributeToSeed(data, inv, 5000.0, chains, RankRequest{}), InputError);
}

// A pin on a nuclide that was never seeded is refused, rather than answered with a zero row
// that reads as "this seed contributes nothing".
TEST(SeedAttribution, RefusesAPinOnANuclideThatWasNotSeeded) {
  const NuclearData data = chain();
  const Inventory inv = seedOf({{parent(), 1.0e20}});

  EXPECT_EQ(requireSeedPin(data, inv, "Sn-100"), parent().key());
  EXPECT_THROW(requireSeedPin(data, inv, "Sb-100"), InputError);  // in the store, not seeded
  EXPECT_THROW(requireSeedPin(data, inv, "Cs-137"), InputError);  // not in the store at all
}

// The engine's own contract, below the ranking: the weight vector is indexed by the store, and
// a caller that hands over the wrong length is told so rather than reading past the end.
TEST(SeedImportanceEngine, RejectsAMisindexedWeightVector) {
  const NuclearData data = chain();
  const Inventory inv = seedOf({{parent(), 1.0e20}});

  const std::vector<double> tooShort(static_cast<std::size_t>(data.size() - 1), 1.0);
  EXPECT_THROW(seedImportance(data, inv, tooShort, 5000.0), InputError);

  const std::vector<double> ok(static_cast<std::size_t>(data.size()), 1.0);
  EXPECT_THROW(seedImportance(data, inv, ok, -1.0), InputError);
}

// With a weight of 1 on every nuclide the response is a plain atom count, and a pure beta chain
// conserves atoms -- so the total must be the seed itself however far it is decayed. A check on
// the adjoint that does not depend on any decay constant being right.
TEST(SeedImportanceEngine, ConservesAtomsUnderAUnitWeight) {
  const NuclearData data = chain();
  const double n0 = 7.0e19;
  const Inventory inv = seedOf({{parent(), n0}});

  const std::vector<double> unit(static_cast<std::size_t>(data.size()), 1.0);
  for (const double t : {0.0, 100.0, 5000.0, 1.0e6}) {
    const SeedImportance imp = seedImportance(data, inv, unit, t);
    EXPECT_NEAR(imp.response, n0, n0 * 1e-10) << "at t = " << t;
  }
}

// A seed that contributes nothing holds no place in an ordering by value, and is not something
// the reader is missing. Te-100 is the chain's stable terminator: seeding it places real atoms
// that are worth zero becquerel forever. Ranking that row spends a topN slot on it, and
// counting it as omitted reports a gap that could not be closed by showing it.
TEST(SeedAttribution, AnInertSeedIsNeitherRankedNorCountedAsOmitted) {
  const NuclearData data = chain();
  const Inventory inv = seedOf({{parent(), 1.0e20}, {terminator(), 5.0e20}});

  const SeedAttribution untruncated =
      attributeToSeed(data, inv, 5000.0, ResponseSpec{}, RankRequest{});
  ASSERT_EQ(untruncated.shares.size(), 1u);
  EXPECT_EQ(untruncated.shares.front().label, "Sn-100");
  EXPECT_EQ(untruncated.shares.front().rank, 1);
  EXPECT_EQ(untruncated.omittedCount, 0);
  EXPECT_NEAR(untruncated.coveredFraction, 1.0, 1e-12);

  // And the cut does not manufacture one either: the stable seed was never in the ordering the
  // cut applies to.
  const SeedAttribution cut =
      attributeToSeed(data, inv, 5000.0, ResponseSpec{}, RankRequest{.topN = 1});
  ASSERT_EQ(cut.shares.size(), 1u);
  EXPECT_EQ(cut.omittedCount, 0);
}

// Asking about it explicitly is a different question and gets an answer: the row comes back,
// rankless, saying this seed contributes nothing HERE -- which is not the same statement as
// "it was never seeded", the one requireSeedPin refuses to let a pin imply.
TEST(SeedAttribution, AnInertSeedComesBackWhenItIsPinned) {
  const NuclearData data = chain();
  const Inventory inv = seedOf({{parent(), 1.0e20}, {terminator(), 5.0e20}});

  RankRequest request;
  request.pinned.push_back(requireSeedPin(data, inv, "Te-100"));
  const SeedAttribution a = attributeToSeed(data, inv, 5000.0, ResponseSpec{}, request);

  ASSERT_EQ(a.shares.size(), 2u);
  const SeedShare& inert = a.shares.back();
  EXPECT_EQ(inert.label, "Te-100");
  EXPECT_TRUE(inert.pinned);
  EXPECT_EQ(inert.rank, 0);
  EXPECT_EQ(inert.value, 0.0);
  EXPECT_GT(inert.seedAtoms, 0.0);  // seeded, and worth nothing -- not absent
  EXPECT_EQ(a.omittedCount, 0);
}

// The two attributions of one exposure must carry the same reservations. `rank` computes it
// forward and says how far the model was stretched; `attribute` computes it through the adjoint
// and reported the identical figure silently, so the same number looked better characterised
// depending only on which command produced it.
TEST(SeedAttribution, ExposureCarriesTheSameCaveatsTheRankingDoes) {
  const NuclearData data = emittingChain();
  const Inventory inv = seedOf({{parent(), 1.0e20}});

  ResponseSpec spec;
  spec.metric = Metric::Exposure;
  spec.unit = Unit::RoentgenPerHour;
  spec.geometry.distanceM = 100.0;

  const DecayResult forward = decay(data, inv, std::vector<double>{5000.0});
  const Ranking ranking = rank(buildResponse(data, forward, spec), 0, RankRequest{});
  const SeedAttribution a = attributeToSeed(data, inv, 5000.0, spec, RankRequest{});

  EXPECT_NEAR(a.total, ranking.total, ranking.total * 1e-9);
  EXPECT_NEAR(a.unmodeledEnergyFraction, ranking.unmodeledEnergyFraction, 1e-12);
  EXPECT_NEAR(a.meanOpticalDepth, ranking.meanOpticalDepth, 1e-12);
  EXPECT_EQ(a.buildup, ranking.buildup);

  // A caveat of zero would satisfy the equalities above and prove nothing.
  EXPECT_GT(a.unmodeledEnergyFraction, 0.0);
  EXPECT_GT(a.meanOpticalDepth, 0.0);
  EXPECT_EQ(a.unmodeledContinuum, std::vector<std::string>{"Sb-100"});
}

// Activity carries neither caveat, and reporting one would annotate the answer with a
// reservation about a metric it never computed.
TEST(SeedAttribution, ActivityCarriesNoExposureCaveats) {
  const NuclearData data = emittingChain();
  const Inventory inv = seedOf({{parent(), 1.0e20}});

  const SeedAttribution a = attributeToSeed(data, inv, 5000.0, ResponseSpec{}, RankRequest{});
  EXPECT_EQ(a.unmodeledEnergyFraction, 0.0);
  EXPECT_EQ(a.meanOpticalDepth, 0.0);
  EXPECT_TRUE(a.unmodeledContinuum.empty());
}

}  // namespace
}  // namespace nusift
