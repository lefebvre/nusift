#include <gtest/gtest.h>

#include <algorithm>
#include <unordered_set>
#include <vector>

#include "nusift/core/element_symbols.hpp"
#include "nusift/core/nuclide.hpp"

namespace nusift {
namespace {

// The packed key is the contract between NuSIFT, cram, the HDF5 store, and OpenMC.
// If it ever changes, every staged store silently becomes garbage -- so it is pinned here
// against hand-computed values rather than against the formula that produces it.
TEST(Zai, KeyMatchesTheCanonicalEncoding) {
  EXPECT_EQ((Zai{92, 235, 0}).key(), 922350);    // U-235
  EXPECT_EQ((Zai{55, 137, 0}).key(), 551370);    // Cs-137
  EXPECT_EQ((Zai{56, 137, 1}).key(), 561371);    // Ba-137m
  EXPECT_EQ((Zai{95, 242, 1}).key(), 952421);    // Am-242m
  EXPECT_EQ((Zai{1, 1, 0}).key(), 10010);        // H-1
  EXPECT_EQ((Zai{118, 294, 0}).key(), 1182940);  // Og-294
}

TEST(Zai, FromKeyRoundTrips) {
  const std::vector<Zai> cases = {
      {92, 235, 0}, {55, 137, 0}, {56, 137, 1}, {95, 242, 1}, {1, 1, 0}, {118, 294, 0},
  };
  for (const Zai& zai : cases) {
    EXPECT_EQ(Zai::fromKey(zai.key()), zai) << "round trip failed for key " << zai.key();
  }
}

// A three-digit mass number must not bleed into the Z field. A=294 is the largest mass in
// any evaluation, so this is the real boundary rather than a synthetic one.
TEST(Zai, LargeMassNumberDoesNotOverflowIntoZ) {
  const Zai og{118, 294, 0};
  const Zai back = Zai::fromKey(og.key());
  EXPECT_EQ(back.z, 118);
  EXPECT_EQ(back.a, 294);
  EXPECT_EQ(back.i, 0);
}

TEST(Zai, NeutronCount) {
  EXPECT_EQ((Zai{92, 235, 0}).n(), 143);
  EXPECT_EQ((Zai{1, 1, 0}).n(), 0);
}

// Ordering by key must read as the chart of the nuclides does: by element, then mass, then
// ground state before metastable. The data store relies on this to write a canonical,
// diffable nuclide axis that does not depend on chain construction order.
TEST(Zai, OrdersByElementThenMassThenIsomer) {
  std::vector<Zai> nuclides = {
      {56, 137, 1}, {55, 137, 0}, {56, 137, 0}, {55, 135, 0}, {92, 235, 0},
  };
  std::sort(nuclides.begin(), nuclides.end());

  const std::vector<Zai> expected = {
      {55, 135, 0}, {55, 137, 0}, {56, 137, 0}, {56, 137, 1}, {92, 235, 0},
  };
  EXPECT_EQ(nuclides, expected);
}

TEST(Zai, HashesDistinctlyInAnUnorderedSet) {
  std::unordered_set<Zai, ZaiHash> seen;
  seen.insert({55, 137, 0});
  seen.insert({56, 137, 1});
  seen.insert({56, 137, 0});
  EXPECT_EQ(seen.size(), 3u);
  // The isomer must not collide with its ground state -- they are physically different
  // nuclides with different half-lives and different photon spectra.
  EXPECT_EQ(seen.count(Zai{56, 137, 1}), 1u);
  EXPECT_EQ(seen.count(Zai{56, 137, 2}), 0u);
}

TEST(Zai, MassChainIsTheMassNumber) {
  EXPECT_EQ(massChain(Zai{56, 140, 0}), 140);
  EXPECT_EQ(massChain(Zai{57, 140, 0}), 140);  // Ba-140 and La-140 share a chain
}

TEST(ElementSymbols, KnownSymbols) {
  EXPECT_STREQ(elementSymbol(1), "H");
  EXPECT_STREQ(elementSymbol(55), "Cs");
  EXPECT_STREQ(elementSymbol(92), "U");
  EXPECT_STREQ(elementSymbol(118), "Og");
}

// Out of range returns "?" rather than throwing, so a malformed nuclide shows up as a
// visibly wrong label instead of aborting a long run mid-report.
TEST(ElementSymbols, OutOfRangeIsSentinelNotAnError) {
  EXPECT_STREQ(elementSymbol(0), "?");
  EXPECT_STREQ(elementSymbol(-1), "?");
  EXPECT_STREQ(elementSymbol(119), "?");
}

TEST(ElementSymbols, ReverseLookupIsCaseInsensitive) {
  EXPECT_EQ(atomicNumber("Cs"), 55);
  EXPECT_EQ(atomicNumber("cs"), 55);
  EXPECT_EQ(atomicNumber("CS"), 55);
  EXPECT_EQ(atomicNumber("U"), 92);
  EXPECT_EQ(atomicNumber("u"), 92);
}

TEST(ElementSymbols, ReverseLookupRejectsNonElements) {
  EXPECT_EQ(atomicNumber(""), 0);
  EXPECT_EQ(atomicNumber("Xx"), 0);
  EXPECT_EQ(atomicNumber("Csx"), 0);
}

// Every symbol must map back to its own Z. This catches a duplicated or transposed entry
// in the table, which a spot-check of a handful of elements would miss.
TEST(ElementSymbols, EverySymbolRoundTrips) {
  for (int z = 1; z <= kMaxAtomicNumber; ++z) {
    EXPECT_EQ(atomicNumber(elementSymbol(z)), z) << "symbol table broken at Z=" << z;
  }
}

}  // namespace
}  // namespace nusift
