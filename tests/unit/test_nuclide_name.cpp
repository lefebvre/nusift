#include <gtest/gtest.h>

#include <string>

#include "nusift/core/error.hpp"
#include "nusift/core/nuclide_name.hpp"

namespace nusift {
namespace {

Zai parsed(const std::string& name) {
  const auto zai = parseNuclideName(name);
  EXPECT_TRUE(zai.has_value()) << "failed to parse \"" << name << "\"";
  return zai.value_or(Zai{});
}

// The spellings that actually turn up in user spreadsheets and CLI arguments. All must
// resolve to the same nuclide, because this is the single point where NuSIFT reconciles
// them and a gap here surfaces as "nuclide not in store" much later.
TEST(NuclideName, AcceptsCommonSpellings) {
  const Zai cs137{55, 137, 0};
  EXPECT_EQ(parsed("Cs-137"), cs137);
  EXPECT_EQ(parsed("cs-137"), cs137);
  EXPECT_EQ(parsed("CS-137"), cs137);
  EXPECT_EQ(parsed("Cs137"), cs137);
  EXPECT_EQ(parsed("cs137"), cs137);
  EXPECT_EQ(parsed("  Cs-137  "), cs137);
}

TEST(NuclideName, AcceptsElementPrefixedIaeaForm) {
  EXPECT_EQ(parsed("92-U-235"), (Zai{92, 235, 0}));
  EXPECT_EQ(parsed("092-U-235"), (Zai{92, 235, 0}));
}

// A bare "m" is the first metastable state, matching ENDF's LISO convention and every
// printed chart of the nuclides.
TEST(NuclideName, MetastableSuffix) {
  const Zai am242m{95, 242, 1};
  EXPECT_EQ(parsed("Am-242m"), am242m);
  EXPECT_EQ(parsed("Am-242M"), am242m);
  EXPECT_EQ(parsed("am242m"), am242m);
  EXPECT_EQ(parsed("Am-242m1"), am242m);
  EXPECT_EQ(parsed("Am-242m2"), (Zai{95, 242, 2}));
}

// Ba-137m is the nuclide that actually emits the 662 keV line attributed to Cs-137, so
// confusing it with its ground state would silently move the entire Cs-137 exposure
// contribution. Worth its own case.
TEST(NuclideName, GroundStateAndIsomerAreDistinct) {
  EXPECT_NE(parsed("Ba-137"), parsed("Ba-137m"));
  EXPECT_EQ(parsed("Ba-137").i, 0);
  EXPECT_EQ(parsed("Ba-137m").i, 1);
}

// Store diagnostics and JSON output emit packed keys; pasting one straight back into the
// CLI should work.
TEST(NuclideName, AcceptsRawPackedKey) {
  EXPECT_EQ(parsed("551370"), (Zai{55, 137, 0}));
  EXPECT_EQ(parsed("922350"), (Zai{92, 235, 0}));
  EXPECT_EQ(parsed("561371"), (Zai{56, 137, 1}));
}

TEST(NuclideName, RejectsNonsense) {
  EXPECT_FALSE(parseNuclideName("").has_value());
  EXPECT_FALSE(parseNuclideName("   ").has_value());
  EXPECT_FALSE(parseNuclideName("Xx-137").has_value());  // not an element
  EXPECT_FALSE(parseNuclideName("Cs").has_value());      // no mass number
  EXPECT_FALSE(parseNuclideName("-137").has_value());    // no element
  EXPECT_FALSE(parseNuclideName("hello").has_value());
}

// A mass number below the proton count is not a nuclide. Rejecting it in the parser keeps
// the error at the input boundary, where the offending text is still available to report.
TEST(NuclideName, RejectsPhysicallyImpossibleNuclides) {
  EXPECT_FALSE(parseNuclideName("Cs-0").has_value());
  EXPECT_FALSE(parseNuclideName("U-1").has_value());
  EXPECT_FALSE(parseNuclideName("Cs-54").has_value());  // A < Z
}

// A pathological mass number must be rejected rather than silently wrapping into a
// different, valid-looking nuclide.
TEST(NuclideName, RejectsOverflowingMassNumber) {
  EXPECT_FALSE(parseNuclideName("Cs-99999999999999999999").has_value());
}

TEST(NuclideName, FormatIsCanonical) {
  EXPECT_EQ(formatNuclideName(Zai{55, 137, 0}), "Cs-137");
  EXPECT_EQ(formatNuclideName(Zai{56, 137, 1}), "Ba-137m");
  EXPECT_EQ(formatNuclideName(Zai{95, 242, 2}), "Am-242m2");
  EXPECT_EQ(formatNuclideName(Zai{92, 235, 0}), "U-235");
}

// Formatting is what makes two runs of the same inventory produce byte-identical output
// regardless of how the nuclides were spelled on input.
TEST(NuclideName, ParseFormatRoundTripsToCanonicalForm) {
  const char* inputs[] = {"cs137", "CS-137", "Cs-137", "551370"};
  for (const char* input : inputs) {
    EXPECT_EQ(formatNuclideName(parsed(input)), "Cs-137") << "input was " << input;
  }
  EXPECT_EQ(formatNuclideName(parsed("am242M")), "Am-242m");
}

// The require* form is what input boundaries call, and the message has to contain the text
// the user actually wrote -- otherwise a bad row in a 500-line CSV is unfindable.
TEST(NuclideName, RequireThrowsInputErrorNamingTheToken) {
  try {
    requireNuclideName("Xx-137");
    FAIL() << "expected InputError";
  } catch (const InputError& e) {
    const std::string what = e.what();
    EXPECT_NE(what.find("Xx-137"), std::string::npos) << what;
  }
}

TEST(NuclideName, RequireIncludesContextWhenGiven) {
  try {
    requireNuclideName("bogus", "inventory.csv line 12");
    FAIL() << "expected InputError";
  } catch (const InputError& e) {
    const std::string what = e.what();
    EXPECT_NE(what.find("inventory.csv line 12"), std::string::npos) << what;
    EXPECT_NE(what.find("bogus"), std::string::npos) << what;
  }
}

// A mass number of 1000 or more does not overflow the packed key, it CARRIES: A=1000 lands in
// the Z field, so "U-1000" would be accepted and then persisted and looked up under the key of
// a nuclide that is not it. Nothing physical comes near A=999, which is exactly why the
// validator has to be the thing that catches it.
TEST(NuclideName, RejectsAMassNumberTheKeyCannotHold) {
  EXPECT_FALSE(parseNuclideName("U-1000").has_value());
  EXPECT_FALSE(parseNuclideName("H-99999").has_value());
  // The boundary itself still parses, and still round-trips through the key.
  const auto widest = parseNuclideName("U-999");
  ASSERT_TRUE(widest.has_value());
  EXPECT_EQ(Zai::fromKey(widest->key()), *widest);
}

// The invariant the bound exists to protect: whatever the parser accepts has to survive the
// packed key it will then be persisted and looked up under.
TEST(NuclideName, EverythingAcceptedSurvivesThePackedKey) {
  const char* names[] = {"Cs-137", "am242m", "U-235", "H-1", "U-999", "551370", "092-U-235"};
  for (const char* name : names) {
    const auto zai = parseNuclideName(name);
    ASSERT_TRUE(zai.has_value()) << name;
    EXPECT_EQ(Zai::fromKey(zai->key()), *zai) << name;
  }
}

// InputError must remain catchable as NusiftError and as std::exception: the CLI
// distinguishes them for its exit code, and the Python layer catches only the base.
TEST(NuclideName, InputErrorIsANusiftError) {
  EXPECT_THROW(requireNuclideName("nope"), NusiftError);
  EXPECT_THROW(requireNuclideName("nope"), std::exception);
}

}  // namespace
}  // namespace nusift
