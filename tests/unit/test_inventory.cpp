#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <string>

#include "nusift/core/error.hpp"
#include "nusift/engine/inventory.hpp"
#include "nusift/nucdata/nuclear_data.hpp"
#include "nusift/units.hpp"
#include "synthetic_chain.hpp"

namespace nusift {
namespace {

// A two-nuclide chain with staged atomic weights: one unstable, one stable. That pairing is
// what makes the activity-unit and mass-unit edge cases testable at all.
NuclearData chainWithWeights(double lambda = 1.0e-3) {
  StoreArrays a = synth::linearChain({lambda});
  // Cs-137's actual ENDF atomic weight ratio, giving a molar mass of 136.907 g/mol. A real
  // value rather than a round one, so the test doubles as a reference for what the
  // conversion should produce.
  a.awr = {135.7305, 135.8};
  return NuclearData::fromArrays(std::move(a));
}

const Zai kUnstable{50, 100, 0};
const Zai kStable{51, 100, 0};

TEST(Inventory, MergesRepeatedNuclidesRatherThanReplacingThem) {
  Inventory inv;
  inv.add(kUnstable, 1.0e20);
  inv.add(kUnstable, 5.0e19);
  ASSERT_EQ(inv.size(), 1);
  EXPECT_DOUBLE_EQ(inv.atomsOf(kUnstable), 1.5e20);
}

// Sorted storage is what makes two inventories built in different orders compare and
// serialize identically -- which a report header and a golden test both depend on.
TEST(Inventory, KeepsEntriesSortedByKeyRegardlessOfInsertionOrder) {
  Inventory a;
  a.add(Zai{92, 235, 0}, 1.0);
  a.add(Zai{55, 137, 0}, 2.0);
  a.add(Zai{56, 137, 1}, 3.0);

  Inventory b;
  b.add(Zai{56, 137, 1}, 3.0);
  b.add(Zai{92, 235, 0}, 1.0);
  b.add(Zai{55, 137, 0}, 2.0);

  ASSERT_EQ(a.size(), b.size());
  for (int i = 0; i < a.size(); ++i) {
    EXPECT_EQ(a.entries()[i].zaiKey, b.entries()[i].zaiKey);
    EXPECT_DOUBLE_EQ(a.entries()[i].atoms, b.entries()[i].atoms);
  }
  EXPECT_LT(a.entries()[0].zaiKey, a.entries()[1].zaiKey);
}

// An atom count is non-negative and finite by definition, and the failures the alternatives
// cause are all silent ones: a negative seed decays into negative activities that sort to the
// bottom of every ranking and vanish, and a NaN propagates through CRAM into a column of NaNs
// whose origin is no longer anywhere in the output.
TEST(Inventory, RefusesCountsThatAreNotAtomCounts) {
  const double inf = std::numeric_limits<double>::infinity();
  Inventory inv;
  EXPECT_THROW(inv.add(kUnstable, -1.0), InputError);
  EXPECT_THROW(inv.add(kUnstable, inf), InputError);
  EXPECT_THROW(inv.add(kUnstable, std::nan("")), InputError);
  EXPECT_TRUE(inv.empty()) << "a refused count must not leave a partial entry behind";

  // Zero is a count, if a dull one, and a subsequent bad value must not corrupt a good entry.
  inv.add(kUnstable, 0.0);
  EXPECT_THROW(inv.add(kUnstable, -1.0), InputError);
  EXPECT_EQ(inv.size(), 1);
  EXPECT_DOUBLE_EQ(inv.atomsOf(kUnstable), 0.0);

  try {
    inv.add(kStable, -1.0);
    FAIL() << "expected InputError";
  } catch (const InputError& e) {
    const std::string what = e.what();
    EXPECT_NE(what.find("Sb-100"), std::string::npos)
        << "the message should name the nuclide: " << what;
  }
}

TEST(Inventory, AtomsOfAbsentNuclideIsZero) {
  Inventory inv;
  inv.add(kUnstable, 1.0e20);
  EXPECT_DOUBLE_EQ(inv.atomsOf(Zai{96, 244, 0}), 0.0);
}

TEST(InventoryUnits, AtomsAndMolesConvertWithoutNuclearData) {
  const NuclearData data = chainWithWeights();
  EXPECT_DOUBLE_EQ(toAtoms(1.0e20, Quantity::Atoms, 0, data), 1.0e20);
  EXPECT_DOUBLE_EQ(toAtoms(1.0, Quantity::Moles, 0, data), units::kAvogadro);
}

// grams -> atoms through the staged AWR, not through A. Cs-137's molar mass is 136.907
// g/mol against a mass number of 137 -- a 0.07% systematic error in every mass-specified
// inventory if A were used instead, in the same direction for every row.
TEST(InventoryUnits, MassConvertsThroughTheStagedAtomicWeight) {
  const NuclearData data = chainWithWeights();
  const double molar = data.molarMassGPerMol(0);
  EXPECT_NEAR(molar, 136.907, 0.002);
  EXPECT_GT(std::abs(molar - 137.0), 0.05) << "using A as the molar mass must be detectably wrong";

  const double atoms = toAtoms(1.0, Quantity::Grams, 0, data);
  EXPECT_NEAR(atoms, units::kAvogadro / molar, units::kAvogadro / molar * 1e-12);

  // And the scaled mass units agree with each other.
  EXPECT_NEAR(toAtoms(1.0, Quantity::Kilograms, 0, data), atoms * 1000.0, atoms * 1e-9);
  EXPECT_NEAR(toAtoms(1000.0, Quantity::Milligrams, 0, data), atoms, atoms * 1e-9);
}

TEST(InventoryUnits, ActivityConvertsThroughTheDecayConstant) {
  const double lambda = 1.0e-3;
  const NuclearData data = chainWithWeights(lambda);

  const double atoms = toAtoms(1.0e6, Quantity::Becquerel, 0, data);
  EXPECT_NEAR(atoms, 1.0e6 / lambda, 1.0e6 / lambda * 1e-12);

  // One curie is 3.7e10 Bq exactly, by definition.
  EXPECT_NEAR(toAtoms(1.0, Quantity::Curie, 0, data),
              toAtoms(units::kBqPerCi, Quantity::Becquerel, 0, data), 1e-3);
  EXPECT_NEAR(toAtoms(1.0, Quantity::Millicurie, 0, data),
              toAtoms(1.0, Quantity::Curie, 0, data) * 1e-3, 1e-3);
}

// N = A/lambda is undefined for a stable nuclide. Returning zero or infinity here would put
// a silently wrong number into the inventory; naming the nuclide lets the user fix the row.
TEST(InventoryUnits, ActivityForAStableNuclideThrowsNamingIt) {
  const NuclearData data = chainWithWeights();
  const int stableIndex = data.indexOf(kStable);
  ASSERT_GE(stableIndex, 0);
  ASSERT_DOUBLE_EQ(data.decayConstant(stableIndex), 0.0);

  try {
    toAtoms(1.0e6, Quantity::Becquerel, stableIndex, data);
    FAIL() << "expected InputError";
  } catch (const InputError& e) {
    // The canonical NuSIFT spelling, not a raw (Z, A) triple: the message has to name the
    // nuclide the way the user wrote it in their file, or it is not actionable.
    const std::string what = e.what();
    EXPECT_NE(what.find("Sb-100"), std::string::npos) << what;
    EXPECT_NE(what.find("stable"), std::string::npos) << what;
  }
}

// Without a staged AWR there is no defensible molar mass. Falling back on A would be wrong
// by ~0.1% and unattributable later, so the conversion refuses and says what to do instead.
TEST(InventoryUnits, MassWithoutAtomicWeightThrowsRatherThanApproximating) {
  const NuclearData data = NuclearData::fromArrays(synth::linearChain({1.0e-3}));
  ASSERT_FALSE(data.hasAtomicWeights());

  try {
    toAtoms(1.0, Quantity::Grams, 0, data);
    FAIL() << "expected InputError";
  } catch (const InputError& e) {
    const std::string what = e.what();
    EXPECT_NE(what.find("atomic weight"), std::string::npos) << what;
  }
}

TEST(InventoryUnits, RoundTripsThroughEveryConvertibleUnit) {
  const NuclearData data = chainWithWeights();
  const double atoms = 1.234e20;
  for (const Quantity q :
       {Quantity::Atoms, Quantity::Moles, Quantity::Grams, Quantity::Kilograms,
        Quantity::Milligrams, Quantity::Becquerel, Quantity::Curie, Quantity::Microcurie}) {
    const double value = fromAtoms(atoms, q, 0, data);
    EXPECT_NEAR(toAtoms(value, q, 0, data), atoms, atoms * 1e-10) << "unit " << quantityName(q);
  }
}

TEST(InventoryUnits, ParsesUnitSpellingsCaseInsensitively) {
  Quantity q{};
  EXPECT_TRUE(parseQuantity("g", q));
  EXPECT_EQ(q, Quantity::Grams);
  EXPECT_TRUE(parseQuantity("GRAMS", q));
  EXPECT_EQ(q, Quantity::Grams);
  EXPECT_TRUE(parseQuantity("Ci", q));
  EXPECT_EQ(q, Quantity::Curie);
  EXPECT_TRUE(parseQuantity("mci", q));
  EXPECT_EQ(q, Quantity::Millicurie);
  EXPECT_TRUE(parseQuantity("MBq", q));
  EXPECT_EQ(q, Quantity::Megabecquerel);
  EXPECT_TRUE(parseQuantity("atoms", q));
  EXPECT_EQ(q, Quantity::Atoms);
}

// mCi and MBq differ only in case, and getting them confused is a factor of 27 error in
// the inventory. The parser is case-insensitive by design, so this pins that the two
// spellings remain distinguishable anyway.
TEST(InventoryUnits, MilliCurieAndMegaBecquerelDoNotCollide) {
  Quantity a{};
  Quantity b{};
  ASSERT_TRUE(parseQuantity("mCi", a));
  ASSERT_TRUE(parseQuantity("MBq", b));
  EXPECT_EQ(a, Quantity::Millicurie);
  EXPECT_EQ(b, Quantity::Megabecquerel);
  EXPECT_NE(a, b);
}

TEST(InventoryUnits, RejectsUnknownSpellings) {
  Quantity q{};
  EXPECT_FALSE(parseQuantity("", q));
  EXPECT_FALSE(parseQuantity("furlongs", q));
  EXPECT_FALSE(parseQuantity("becquerels", q));
}

}  // namespace
}  // namespace nusift
