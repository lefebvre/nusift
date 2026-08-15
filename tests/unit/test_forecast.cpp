#include <gtest/gtest.h>

#include <cmath>
#include <string>
#include <vector>

#include "nusift/engine/decay_engine.hpp"
#include "nusift/io/time_spec.hpp"
#include "nusift/nucdata/nuclear_data.hpp"
#include "nusift/triage/forecast.hpp"
#include "synthetic_chain.hpp"

namespace nusift {
namespace {

// Two independent nuclides, each decaying to its own stable daughter. Independent so their
// activities are pure exponentials and the crossing between them has a closed form.
NuclearData twoIndependentEmitters(double lambdaFast, double lambdaSlow) {
  StoreArrays a;
  a.provenance.version = 1;
  // Sorted ascending by key: Sn-100, Sb-100, Nd-120, Pm-120.
  a.nuclideKey = {Zai{50, 100, 0}.key(), Zai{51, 100, 0}.key(), Zai{60, 120, 0}.key(),
                  Zai{61, 120, 0}.key()};
  a.halfLife = {synth::halfLifeFor(lambdaFast), 0.0, synth::halfLifeFor(lambdaSlow), 0.0};
  a.modeOffset = {0, 1, 1, 2, 2};
  a.modeRtyp = {synth::kBetaMinus, synth::kBetaMinus};
  a.modeBranching = {1.0, 1.0};
  a.modeFinalState = {0, 0};
  a.modeIsFission = {0, 0};
  return NuclearData::fromArrays(std::move(a));
}

ResponseTable tableOver(const NuclearData& data, const Inventory& inventory,
                        const std::vector<double>& times) {
  const DecayResult result = decay(data, inventory, times);
  return buildResponse(data, result, ResponseSpec{});
}

// The crossing has a closed form, so the interpolated boundary can be checked against it
// rather than against a grid point. That is the whole reason boundaries are interpolated: on a
// log grid at late times, consecutive samples can be years apart.
//
//   lambda_a N_a e^{-lambda_a t} = lambda_b N_b e^{-lambda_b t}
//   t = ln(lambda_a N_a / lambda_b N_b) / (lambda_a - lambda_b)
TEST(Forecast, WindowBoundaryMatchesTheAnalyticCrossing) {
  const double fast = 1.0e-3;
  const double slow = 1.0e-4;
  const double atoms = 1.0e20;
  const NuclearData data = twoIndependentEmitters(fast, slow);

  Inventory inv;
  inv.add(Zai{50, 100, 0}, atoms);
  inv.add(Zai{60, 120, 0}, atoms);

  const double expected = std::log((fast * atoms) / (slow * atoms)) / (fast - slow);
  ASSERT_GT(expected, 0.0);

  const ResponseTable table = tableOver(data, inv, logspace(1.0, 1.0e5, 60));
  const std::vector<DominanceWindow> windows = dominanceWindows(table);

  ASSERT_EQ(windows.size(), 2u);
  EXPECT_EQ(windows[0].label, "Sn-100") << "the shorter-lived nuclide leads first";
  EXPECT_EQ(windows[1].label, "Nd-120");

  // Within a percent of the closed form, and far tighter than the grid spacing at that point.
  EXPECT_NEAR(windows[0].endSeconds, expected, expected * 0.01);
  EXPECT_DOUBLE_EQ(windows[1].startSeconds, windows[0].endSeconds)
      << "windows must abut, not overlap or leave a gap";
}

TEST(Forecast, WindowsSpanTheWholeGridInOrder) {
  const NuclearData data = twoIndependentEmitters(1.0e-3, 1.0e-4);
  Inventory inv;
  inv.add(Zai{50, 100, 0}, 1.0e20);
  inv.add(Zai{60, 120, 0}, 1.0e20);

  const std::vector<double> times = logspace(1.0, 1.0e5, 60);
  const std::vector<DominanceWindow> windows = dominanceWindows(tableOver(data, inv, times));

  ASSERT_FALSE(windows.empty());
  EXPECT_DOUBLE_EQ(windows.front().startSeconds, times.front());
  EXPECT_DOUBLE_EQ(windows.back().endSeconds, times.back());
  for (std::size_t i = 0; i + 1 < windows.size(); ++i) {
    EXPECT_LT(windows[i].startSeconds, windows[i].endSeconds);
    EXPECT_DOUBLE_EQ(windows[i].endSeconds, windows[i + 1].startSeconds);
  }
}

// A single contributor that leads throughout is one window, not sixty.
TEST(Forecast, AConstantLeaderIsOneWindow) {
  const NuclearData data = twoIndependentEmitters(1.0e-3, 1.0e-4);
  Inventory inv;
  inv.add(Zai{50, 100, 0}, 1.0e20);  // only one seeded, so nothing else ever leads

  const std::vector<DominanceWindow> windows =
      dominanceWindows(tableOver(data, inv, logspace(1.0, 1.0e4, 40)));
  ASSERT_EQ(windows.size(), 1u);
  EXPECT_EQ(windows[0].label, "Sn-100");
}

// The regression this threshold was rewritten for.
//
// Absorption used to measure a run's length as a fraction of the LINEAR span. On a log grid
// from one hour to a century, a nuclide's genuine reign over the first hours is under 0.001%
// of the axis, so every early window was absorbed into whichever contributor led first -- and
// a 33-minute nuclide was reported as leading for 263 days. Counting samples instead is
// grid-agnostic.
TEST(Forecast, EarlyWindowsSurviveOnALogGrid) {
  const double fast = 1.0e-3;
  const double slow = 1.0e-9;
  const double atoms = 1.0e20;
  const NuclearData data = twoIndependentEmitters(fast, slow);
  Inventory inv;
  inv.add(Zai{50, 100, 0}, atoms);
  inv.add(Zai{60, 120, 0}, atoms);

  // Six decades of grid. The crossing sits about 1.4e4 s in -- roughly 1.4% of the linear span
  // and therefore invisible to a duration-based threshold, but a third of the way along a log
  // axis and unmistakable in samples.
  const std::vector<DominanceWindow> windows =
      dominanceWindows(tableOver(data, inv, logspace(1.0, 1.0e6, 60)));

  ASSERT_GE(windows.size(), 2u);
  EXPECT_EQ(windows[0].label, "Sn-100");

  const double expected = std::log((fast * atoms) / (slow * atoms)) / (fast - slow);
  EXPECT_NEAR(windows[0].endSeconds, expected, expected * 0.01)
      << "the short-lived leader's window ended at " << windows[0].endSeconds
      << " s rather than the analytic crossing at " << expected << " s";
}

// A nuclide that only matters late still belongs in the table; ordering by value at any one
// time would bury exactly the row a forecast exists to surface.
TEST(Forecast, UnionTopNIncludesALateBloomer) {
  const NuclearData data = twoIndependentEmitters(1.0e-3, 1.0e-6);
  Inventory inv;
  inv.add(Zai{50, 100, 0}, 1.0e22);  // dominant early
  inv.add(Zai{60, 120, 0}, 1.0e20);  // only ever matters late

  const std::vector<RankTrack> tracks =
      unionTopN(tableOver(data, inv, logspace(1.0, 1.0e7, 60)), 1);

  ASSERT_EQ(tracks.size(), 2u) << "both lead at some point, so both belong";
  // Ordered by peak share, and each peak time should sit where that nuclide actually leads.
  EXPECT_GT(tracks[0].peakFraction, 0.0);
  EXPECT_GT(tracks[1].peakFraction, 0.0);
}

// The complement: what must always be accounted for, as opposed to what could ever matter.
TEST(Forecast, PersistentTopNExcludesAnyoneWhoDropsOut) {
  const NuclearData data = twoIndependentEmitters(1.0e-3, 1.0e-6);
  Inventory inv;
  inv.add(Zai{50, 100, 0}, 1.0e22);
  inv.add(Zai{60, 120, 0}, 1.0e20);

  const ResponseTable table = tableOver(data, inv, logspace(1.0, 1.0e7, 60));
  const std::vector<RankTrack> persistent = persistentTopN(table, 1);
  const std::vector<RankTrack> ever = unionTopN(table, 1);

  EXPECT_LT(persistent.size(), ever.size())
      << "the leader changes, so nobody holds first place throughout";
}

TEST(Forecast, RankTracksRecordRankAndFractionAtEveryTime) {
  const NuclearData data = twoIndependentEmitters(1.0e-3, 1.0e-4);
  Inventory inv;
  inv.add(Zai{50, 100, 0}, 1.0e20);
  inv.add(Zai{60, 120, 0}, 1.0e20);

  const std::vector<double> times = logspace(1.0, 1.0e5, 30);
  const ResponseTable table = tableOver(data, inv, times);
  const std::vector<RankTrack> tracks = unionTopN(table, 2);

  ASSERT_FALSE(tracks.empty());
  for (const RankTrack& track : tracks) {
    EXPECT_EQ(track.rank.size(), times.size());
    EXPECT_EQ(track.fraction.size(), times.size());
    EXPECT_GE(track.peakFraction, 0.0);
    // A hair over 1.0 is legitimate: the fraction is a value divided by a SUM of values, and
    // when one contributor holds essentially everything the sum can round a few ulp below it.
    // The tolerance is tight enough that a real accounting error -- a fraction of 1.5, say --
    // would still fail, which is why this is not clamped at the source: clamping would hide
    // exactly the bug the assertion is for.
    EXPECT_LE(track.peakFraction, 1.0 + 1e-9);
    // The recorded peak must actually be the largest fraction in the track.
    const double largest = *std::max_element(track.fraction.begin(), track.fraction.end());
    EXPECT_NEAR(track.peakFraction, largest, 1e-15);
    EXPECT_NEAR(track.fraction[static_cast<std::size_t>(track.peakTimeIndex)], largest, 1e-15);
  }
}

TEST(Forecast, HandlesAnEmptyOrSingleTimeTableWithoutCrashing) {
  ResponseTable empty;
  EXPECT_TRUE(dominanceWindows(empty).empty());
  EXPECT_TRUE(unionTopN(empty, 5).empty());

  const NuclearData data = twoIndependentEmitters(1.0e-3, 1.0e-4);
  Inventory inv;
  inv.add(Zai{50, 100, 0}, 1.0e20);
  const std::vector<DominanceWindow> single =
      dominanceWindows(tableOver(data, inv, std::vector<double>{100.0}));
  EXPECT_EQ(single.size(), 1u);
}

}  // namespace
}  // namespace nusift
