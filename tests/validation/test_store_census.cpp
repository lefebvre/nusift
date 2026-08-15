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

// The three coverage GAPS, pinned separately because they are three different problems and
// lumping them hid the first behind the second. A nuclide with no evaluated spectrum
// contributes exactly zero to an exposure ranking while genuinely emitting photons; one with a
// continuum tail is present but low; a clamped line is ranked but evaluated with an
// end-of-table coefficient. All three are reported to users, so all three are gated.
TEST(StoreCensus, TheKnownCoverageGapsAreUnchanged) {
  const Census c = censusOf(committedStore());
  EXPECT_EQ(c.noSpectrumAtAll, 1546) << "unstable nuclides emitting photons with no spectrum";
  EXPECT_EQ(c.withClampedLines, 1471) << "nuclides carrying at least one clamped line";
  EXPECT_EQ(c.clampedLines, 5705) << "lines outside the tabulated air-coefficient range";
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
