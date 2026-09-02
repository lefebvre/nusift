// Assays taken on different dates, brought to one.
//
// The arithmetic is a decay solve and is checked against the closed form like every other one
// here. What is not arithmetic, and is the reason this file exists, is the direction rule:
// forward is a solve and backward is an inverse problem, and the refusal is tested as carefully
// as the number.
//
#include <gtest/gtest.h>

#include <cmath>
#include <sstream>
#include <vector>

#include "nusift/core/error.hpp"
#include "nusift/engine/reconcile.hpp"
#include "nusift/io/inventory_io.hpp"
#include "nusift/io/number_format.hpp"
#include "nusift/io/time_spec.hpp"
#include "nusift/nucdata/nuclear_data.hpp"
#include "synthetic_chain.hpp"

namespace nusift {
namespace {

constexpr double kLambda = 1.0e-8;  // ~2.2 y half-life, so a year of carry is visible

// Sn-100 decaying to a stable Sb-100. Atoms are conserved along it, which makes the
// atoms-at-assay against atoms-at-epoch column a real check rather than decoration.
NuclearData chain() {
  StoreArrays arrays = synth::linearChain({kLambda});
  arrays.awr = {99.1, 99.2};
  return NuclearData::fromArrays(std::move(arrays));
}

const Zai kParent{50, 100, 0};
const Zai kDaughter{51, 100, 0};

AssayGroup assay(const char* date, const Zai& zai, double atoms, const char* label = "sheet") {
  AssayGroup group;
  group.dateSeconds = parseCalendarDate(date);
  group.inventory.add(zai, atoms);
  group.label = label;
  return group;
}

// --- carrying forward --------------------------------------------------------

TEST(Reconcile, CarriesAnAssayForwardByExactlyTheDecayLaw) {
  const NuclearData data = chain();
  const double atoms = 1.0e20;
  // The second assay exists only to put the epoch a year on, and carries one atom rather than
  // none because an empty assay is refused -- and rightly, since it says nothing about anything.
  std::vector<AssayGroup> groups = {assay("2023-01-01", kParent, atoms),
                                    assay("2024-01-01", kDaughter, 1.0)};
  const double epoch = latestAssayDate(groups);
  const double gap = epoch - groups[0].dateSeconds;
  ASSERT_GT(gap, 0.0);

  const Reconciliation result = reconcile(data, groups, epoch);
  EXPECT_DOUBLE_EQ(result.epochSeconds, epoch);
  EXPECT_DOUBLE_EQ(result.spanSeconds, gap);

  const double expected = atoms * std::exp(-kLambda * gap);
  EXPECT_NEAR(result.inventory.atomsOf(kParent), expected, expected * 1.0e-9);
  // Everything the parent lost went to the stable daughter, plus the one atom assayed there.
  EXPECT_NEAR(result.inventory.atomsOf(kDaughter), atoms - expected + 1.0, atoms * 1.0e-9);
}

// The assay that defines the epoch is not solved for at all, and the observable consequence is
// that its numbers come back BIT EXACT. A solve over a zero interval is the identity in exact
// arithmetic and is not in floating point, so paying for one would put CRAM noise on numbers a
// user supplied by hand.
TEST(Reconcile, LeavesTheEpochsOwnAssayUntouched) {
  const NuclearData data = chain();
  const double atoms = 1.2345678901234e19;
  std::vector<AssayGroup> groups = {assay("2023-01-01", kDaughter, 5.0e18),
                                    assay("2024-01-01", kParent, atoms)};

  const Reconciliation result = reconcile(data, groups, latestAssayDate(groups));
  EXPECT_DOUBLE_EQ(result.inventory.atomsOf(kParent), atoms);
  EXPECT_EQ(result.contributions.back().carriedSeconds, 0.0);
  EXPECT_DOUBLE_EQ(result.contributions.back().atomsAtEpoch,
                   result.contributions.back().atomsAtAssay);
}

// Reconciliation is a sum of independently carried assays, which is exact because the decay
// operator is linear. Splitting one assay in two and reconciling both must give what the whole
// one gives -- the same linearity every other exactness claim in NuSIFT rests on.
TEST(Reconcile, IsLinearInTheAssaysItMerges) {
  const NuclearData data = chain();
  std::vector<AssayGroup> whole = {assay("2022-06-01", kParent, 9.0e19),
                                   assay("2024-01-01", kDaughter, 1.0)};
  std::vector<AssayGroup> split = {assay("2022-06-01", kParent, 4.0e19, "a"),
                                   assay("2022-06-01", kParent, 5.0e19, "b"),
                                   assay("2024-01-01", kDaughter, 1.0)};

  const double epoch = latestAssayDate(whole);
  const Reconciliation one = reconcile(data, whole, epoch);
  const Reconciliation two = reconcile(data, split, epoch);

  EXPECT_NEAR(two.inventory.atomsOf(kParent), one.inventory.atomsOf(kParent),
              one.inventory.atomsOf(kParent) * 1.0e-12);
  // Two assays on the same date are one group in the file reader; passed separately they stay
  // separate here, which is what lets a caller keep two sheets distinguishable in the report.
  EXPECT_EQ(two.contributions.size(), 3u);
}

TEST(Reconcile, AcceptsAnEpochLaterThanEveryAssay) {
  const NuclearData data = chain();
  std::vector<AssayGroup> groups = {assay("2024-01-01", kParent, 1.0e20)};
  const double later = parseCalendarDate("2026-01-01");

  const Reconciliation result = reconcile(data, groups, later);
  const double gap = later - groups[0].dateSeconds;
  const double expected = 1.0e20 * std::exp(-kLambda * gap);
  EXPECT_NEAR(result.inventory.atomsOf(kParent), expected, expected * 1.0e-9);
  EXPECT_DOUBLE_EQ(result.contributions.front().carriedSeconds, gap);
}

// --- the direction rule ------------------------------------------------------

TEST(Reconcile, RefusesToCarryAnAssayBackward) {
  const NuclearData data = chain();
  std::vector<AssayGroup> groups = {assay("2023-01-01", kParent, 1.0e20),
                                    assay("2024-01-01", kDaughter, 1.0e19)};

  // Un-growing a daughter has no unique answer: the daughter measured at assay could have been
  // there from the start or grown in since, and nothing in a one-instant measurement separates
  // them. Refused rather than caveated, because a caveat on a number nobody can check is a
  // disclaimer and not a warning.
  EXPECT_THROW(reconcile(data, groups, parseCalendarDate("2022-01-01")), InputError);
  EXPECT_THROW(reconcile(data, groups, parseCalendarDate("2023-06-01")), InputError);
  // The earliest epoch that carries nothing backward is the latest assay, and it is allowed.
  EXPECT_NO_THROW(reconcile(data, groups, parseCalendarDate("2024-01-01")));
}

TEST(Reconcile, RefusesAnEmptySetAndAnEmptyAssay) {
  const NuclearData data = chain();
  EXPECT_THROW(latestAssayDate({}), InputError);
  EXPECT_THROW(reconcile(data, {}, 0.0), InputError);

  std::vector<AssayGroup> withEmpty = {assay("2024-01-01", kParent, 1.0e20)};
  AssayGroup nothing;
  nothing.dateSeconds = parseCalendarDate("2023-01-01");
  nothing.label = "blank sheet";
  withEmpty.push_back(nothing);
  EXPECT_THROW(reconcile(data, withEmpty, parseCalendarDate("2024-01-01")), InputError);
}

TEST(Reconcile, DefaultsTheEpochToTheLatestAssay) {
  std::vector<AssayGroup> groups = {assay("2024-01-01", kParent, 1.0),
                                    assay("2021-05-05", kParent, 1.0),
                                    assay("2023-11-30", kParent, 1.0)};
  EXPECT_DOUBLE_EQ(latestAssayDate(groups), parseCalendarDate("2024-01-01"));
}

// --- reading dated files -----------------------------------------------------

DatedInventory parseCsv(const std::string& text, const NuclearData& data) {
  std::istringstream in(text);
  return readInventoryDatedCsv(in, data, "test.csv");
}

TEST(DatedInventoryFile, GroupsRowsSharingADateIntoOneAssay) {
  const NuclearData data = chain();
  const DatedInventory dated = parseCsv(
      "nuclide,quantity,unit,assayed\n"
      "Sn-100,1e20,atoms,2024-03-15\n"
      "Sb-100,2e20,atoms,2023-01-10\n"
      "Sn-100,3e19,atoms,2024-03-15\n",
      data);

  EXPECT_TRUE(dated.dated);
  ASSERT_EQ(dated.groups.size(), 2u);
  // Ascending by date, whatever order the rows appeared in.
  EXPECT_DOUBLE_EQ(dated.groups[0].dateSeconds, parseCalendarDate("2023-01-10"));
  EXPECT_DOUBLE_EQ(dated.groups[1].dateSeconds, parseCalendarDate("2024-03-15"));
  // Two rows of one nuclide on one date accumulate, exactly as they do in an undated file.
  EXPECT_DOUBLE_EQ(dated.groups[1].inventory.atomsOf(kParent), 1.3e20);
}

TEST(DatedInventoryFile, AFileWithNoDatesIsOneUndatedAssay) {
  const NuclearData data = chain();
  const DatedInventory dated = parseCsv("Sn-100,1e20,atoms\nSb-100,2e20,atoms\n", data);
  EXPECT_FALSE(dated.dated);
  ASSERT_EQ(dated.groups.size(), 1u);
  EXPECT_EQ(dated.groups[0].inventory.size(), 2);
}

TEST(DatedInventoryFile, RefusesAFileThatDatesSomeRowsAndNotOthers) {
  const NuclearData data = chain();
  // No reading of the undated row is anything but a guess about when it was measured, and the
  // guess would not be visible in any answer, so the file is refused instead.
  EXPECT_THROW(parseCsv("Sn-100,1e20,atoms,2024-03-15\nSb-100,2e20,atoms\n", data), InputError);
  EXPECT_THROW(parseCsv("Sn-100,1e20,atoms\nSb-100,2e20,atoms,2024-03-15\n", data), InputError);
  EXPECT_THROW(parseCsv("Sn-100,1e20,atoms,not-a-date\n", data), InputError);
}

// The whole point of the io half: every command reads inventories through readInventory, so a
// dated file has to arrive already reconciled or the dates would be silently ignored.
TEST(DatedInventoryFile, ReadingADatedFileReconcilesIt) {
  const NuclearData data = chain();
  std::istringstream in(
      "Sn-100,1e20,atoms,2023-01-01\n"
      "Sb-100,1,atoms,2024-01-01\n");
  const Inventory inventory = readInventoryCsv(in, data, "test.csv");

  const double gap = parseCalendarDate("2024-01-01") - parseCalendarDate("2023-01-01");
  const double expected = 1.0e20 * std::exp(-kLambda * gap);
  EXPECT_NEAR(inventory.atomsOf(kParent), expected, expected * 1.0e-9);
  // And the epoch travels in the provenance, because a seed reconciled to one date is a
  // different inventory from the same file reconciled to another.
  EXPECT_NE(inventory.provenance().find("2024-01-01"), std::string::npos);
}

TEST(DatedInventoryFile, ReadingAnUndatedFileSolvesNothing) {
  const NuclearData data = chain();
  const double atoms = 9.87654321e19;
  std::istringstream in("Sn-100," + shortestRoundTrip(atoms) + ",atoms\n");
  const Inventory inventory = readInventoryCsv(in, data, "test.csv");
  // Bit exact: the overwhelmingly common case must stay exactly as cheap and exactly as
  // faithful as it was before dates existed.
  EXPECT_DOUBLE_EQ(inventory.atomsOf(kParent), atoms);
}

TEST(DatedInventoryFile, ReadsDatesFromJsonToo) {
  const NuclearData data = chain();
  std::istringstream in(
      "[{\"nuclide\": \"Sn-100\", \"quantity\": 1e20, \"unit\": \"atoms\", "
      "\"assayed\": \"2023-01-01\"},\n"
      " {\"nuclide\": \"Sb-100\", \"quantity\": 1, \"unit\": \"atoms\", "
      "\"assayed\": \"2024-01-01\"}]\n");
  const DatedInventory dated = readInventoryDatedJson(in, data, "test.json");
  EXPECT_TRUE(dated.dated);
  ASSERT_EQ(dated.groups.size(), 2u);
  EXPECT_DOUBLE_EQ(dated.groups[0].dateSeconds, parseCalendarDate("2023-01-01"));
}

}  // namespace
}  // namespace nusift
