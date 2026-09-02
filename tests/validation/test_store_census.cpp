#include <gtest/gtest.h>

#include <algorithm>
#include <string>
#include <vector>

#include "nusift/core/nuclide_name.hpp"
#include "nusift/exposure/air_coefficients.hpp"
#include "nusift/nucdata/nuclear_data.hpp"
#include "nusift/nucdata/photon_lines.hpp"
#include "validation_store.hpp"

namespace nusift::validation {
namespace {

// --- what the shipped evaluation contains ----------------------------------
//
// These are the only checks in the suite that compare the store against itself rather than
// against something published, and they earn their place for a different reason: staging is a
// one-time offline step whose output is committed, so a silent change in what ENDF parsing
// produces would otherwise reach users with nothing to catch it. Every number here is quoted
// in docs/nuclear-data.md and printed by `nusift data info`, so a drift here is a drift in
// documentation that is already published.
//
// They are tripwires on a deliberately-updated file, not golden files for the physics. When
// the store is restaged these counts are expected to move, and validation/README.md documents
// updating them as a step of that process.

struct Census {
  int unstable = 0;
  int withLines = 0;
  int partialContinuum = 0;
  int withClampedLines = 0;
  int clampedLines = 0;
  int noSpectrumAtAll = 0;
  int sfWithoutYields = 0;
  int withWeights = 0;
};

// Deliberately the same traversal `nusift data info` performs (runDataInfo in
// nusift_apps/nusift.cpp). If the two ever disagree, the counts a user is shown and the counts
// CI gates on have diverged, which is worth a failure in its own right.
Census censusOf(const NuclearData& data) {
  Census c;
  for (int i = 0; i < data.size(); ++i) {
    if (data.decayConstant(i) > 0.0) {
      ++c.unstable;
      const LineSpectrum lines = data.lines(i);
      if (!lines.empty()) {
        ++c.withLines;
        if (data.unmodeledPhotonFraction(i) > 0.05) {
          ++c.partialContinuum;
        }
        int clampedHere = 0;
        for (const GammaLine& line : lines) {
          if (exposure::isOutsideTabulatedRange(line.energyEv)) {
            ++clampedHere;
          }
        }
        if (clampedHere > 0) {
          c.clampedLines += clampedHere;
          ++c.withClampedLines;
        }
      } else if (data.emEnergyEv(i) > 0.0) {
        ++c.noSpectrumAtAll;
      }
    }
    if (data.molarMassGPerMol(i) > 0.0) {
      ++c.withWeights;
    }
  }
  c.sfWithoutYields = static_cast<int>(data.spontaneousFissionWithoutYields().size());
  return c;
}

TEST(StoreCensus, ProvenanceIsTheShippedEvaluation) {
  const StoreProvenance& provenance = committedStore().provenance();
  EXPECT_EQ(provenance.version, 1);
  EXPECT_EQ(provenance.library, "ENDF/B-VIII.1");
  EXPECT_EQ(provenance.stagedTapeCount, 3821);
}

// The two counts differ by orders of magnitude in general, and only the first describes what
// the store knows: closure and fission-yield registration inflate the chain with nuclides that
// carry no evaluated data. Pinning both is what stops a future change reporting chain size as
// coverage.
TEST(StoreCensus, StagedAndChainSizesAreUnchanged) {
  EXPECT_EQ(committedStore().stagedCount(), 3828);
  EXPECT_EQ(committedStore().size(), 4012);
}

TEST(StoreCensus, CoverageCountsAreUnchanged) {
  const Census c = censusOf(committedStore());
  EXPECT_EQ(c.unstable, 3562);
  EXPECT_EQ(c.withLines, 1595);
  EXPECT_EQ(c.partialContinuum, 34);
  EXPECT_EQ(c.withWeights, 3576);
}

// The four coverage GAPS, pinned separately because they are four different problems and
// lumping them hid the first behind the second. A nuclide with no evaluated spectrum
// contributes exactly zero to an exposure ranking while genuinely emitting photons; one with a
// continuum tail is present but low; a clamped line is ranked but evaluated with an
// end-of-table coefficient; and a spontaneous-fission branch with no yield set is the one gap
// that costs atoms rather than photons. All four are reported to users, so all four are gated.
TEST(StoreCensus, TheKnownCoverageGapsAreUnchanged) {
  const Census c = censusOf(committedStore());
  EXPECT_EQ(c.noSpectrumAtAll, 1546) << "unstable nuclides emitting photons with no spectrum";
  EXPECT_EQ(c.withClampedLines, 1471) << "nuclides carrying at least one clamped line";
  EXPECT_EQ(c.clampedLines, 5705) << "lines outside the tabulated air-coefficient range";
  EXPECT_EQ(c.sfWithoutYields, 103) << "spontaneous-fission branches with no yield set";
}

// The uncertainty census. How MUCH of an evaluation carries a sigma is what decides whether an
// error budget built on it means anything, so the coverage is pinned like every other coverage
// figure here rather than discovered later by whoever first tries to use it.
//
// ENDF/B-VIII.1 states a half-life uncertainty for 85% of staged nuclides and a branching
// uncertainty for every decay mode it carries. The gap is real and is mostly the short-lived
// exotics, whose half-lives are themselves estimates.
TEST(StoreCensus, HalfLifeUncertaintyCoverageIsUnchanged) {
  const NuclearData& store = committedStore();
  int stated = 0;
  int unstable = 0;
  for (int i = 0; i < store.stagedCount(); ++i) {
    if (!(store.halfLifeSeconds(i) > 0.0)) {
      continue;
    }
    ++unstable;
    if (store.halfLifeUncertainty(i) > 0.0) {
      ++stated;
    }
  }
  EXPECT_EQ(unstable, 3562);
  EXPECT_EQ(stated, 3270) << "nuclides whose evaluated half-life carries a stated uncertainty";
}

// Sanity on the VALUES, not merely on the count. A sigma LARGER than the half-life it qualifies
// would be nonsense and would be the signature of a mispairing; one exactly EQUAL to it is not,
// and three nuclides genuinely carry that.
//
// Kr-100 is evaluated at 7 ms +/- 7 ms -- the two numbers sit side by side in the same MT457
// record -- which is an evaluator saying the value is known to within about a factor of two.
// Mt-266m and Mt-269 say the same thing. Pinned at three so a restage that changed it surfaces
// here rather than in whatever first divides by one of them.
TEST(StoreCensus, NoStatedHalfLifeUncertaintyExceedsItsHalfLife) {
  const NuclearData& store = committedStore();
  int checked = 0;
  int fullyUncertain = 0;
  for (int i = 0; i < store.stagedCount(); ++i) {
    const double halfLife = store.halfLifeSeconds(i);
    const double sigma = store.halfLifeUncertainty(i);
    if (!(halfLife > 0.0) || !(sigma > 0.0)) {
      continue;
    }
    ++checked;
    EXPECT_LE(sigma, halfLife) << formatNuclideName(store.zaiAt(i))
                               << " has an uncertainty larger than its half-life";
    if (sigma >= halfLife) {
      ++fullyUncertain;
    }
  }
  EXPECT_GT(checked, 3000);
  EXPECT_EQ(fullyUncertain, 3) << "nuclides evaluated at 100% uncertainty on their half-life";
}

// Against values published outside NuSIFT, which is what makes this a validation test rather
// than a round-trip. ENDF/B-VIII.1 evaluates Cs-137 at 30.08(9) y and Co-60 at 5.2711(4) y, so
// the relative uncertainties are 0.30% and 0.0076%; a store that had paired the sigma column
// with the wrong nuclide would land nowhere near either.
TEST(StoreCensus, StagedHalfLifeUncertaintiesMatchTheirEvaluations) {
  const NuclearData& store = committedStore();
  struct Expected {
    Zai zai;
    double relative;
  };
  for (const Expected& one :
       {Expected{Zai{55, 137, 0}, 0.0030}, Expected{Zai{27, 60, 0}, 7.6e-5}}) {
    const int index = store.indexOf(one.zai);
    ASSERT_GE(index, 0) << formatNuclideName(one.zai);
    const double relative = store.halfLifeUncertainty(index) / store.halfLifeSeconds(index);
    EXPECT_NEAR(relative, one.relative, one.relative * 0.1) << formatNuclideName(one.zai);
  }
}

// Every one of those is a heavy nuclide, and decay only ever lowers A, so nothing a
// fission-product source produces can reach one. That is the fact behind `data info`
// calling the gap negligible for such a source, and it is checked rather than assumed.
TEST(StoreCensus, TheSpontaneousFissionGapLiesAboveTheFissionProductRange) {
  for (const Zai& zai : committedStore().spontaneousFissionWithoutYields()) {
    EXPECT_GE(zai.a, 180) << formatNuclideName(zai) << " is within reach of a fission product";
  }
}

// Seeding from fission is a headline capability, and it is only available for parents the
// store carries yields for. A membership check rather than a count: the yield sublibrary can
// gain parents without invalidating anything, but losing one of these would remove a
// documented example from the README.
TEST(StoreCensus, TheDocumentedFissileParentsAreAvailable) {
  std::vector<std::string> names;
  for (const Zai& zai : committedStore().fissionYields().parents()) {
    names.push_back(formatNuclideName(zai));
  }
  for (const char* expected : {"U-235", "U-238", "Pu-239", "Pu-241"}) {
    EXPECT_NE(std::find(names.begin(), names.end(), expected), names.end())
        << expected << " should be seedable from fission";
  }
}

}  // namespace
}  // namespace nusift::validation
