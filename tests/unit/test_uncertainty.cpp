// The error bar the assay puts on the answer.
//
// Two things are worth checking beyond the arithmetic. First that the propagation is EXACT and
// not first-order, which for a single nuclide means a 5% assay gives a 5% answer at every time
// and to every digit -- an approximation would drift. Second the identity that says the per-assay
// adjoints were run at the right times: R formed as the sum over assays of <g_a, n_a> has to
// equal what an ordinary ranking reports for the same spec, and for a DATED file those are
// different solves reaching the same number.
//
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <sstream>
#include <vector>

#include "nusift/core/error.hpp"
#include "nusift/engine/decay_engine.hpp"
#include "nusift/io/inventory_io.hpp"
#include "nusift/io/time_spec.hpp"
#include "nusift/nucdata/nuclear_data.hpp"
#include "nusift/triage/uncertainty.hpp"
#include "synthetic_chain.hpp"

namespace nusift {
namespace {

constexpr double kLambda = 1.0e-6;

NuclearData chain() {
  StoreArrays arrays = synth::linearChain({kLambda, 5.0e-7});
  synth::addLines(arrays, 1, {661657.0}, {0.9});
  return NuclearData::fromArrays(std::move(arrays));
}

const Zai kParent{50, 100, 0};
const Zai kMiddle{51, 100, 0};

AssayGroup assay(const char* date, const Inventory& inventory, const char* label = "sheet") {
  AssayGroup group;
  group.dateSeconds = parseCalendarDate(date);
  group.inventory = inventory;
  group.label = label;
  group.dated = true;
  return group;
}

Inventory one(const Zai& zai, double atoms, double sigma) {
  Inventory inv;
  inv.add(zai, atoms, sigma);
  return inv;
}

// --- exactness ---------------------------------------------------------------

// For a single seeded nuclide R = lambda N0 e^{-lambda T}, so dR/dN0 = lambda e^{-lambda T} and
// sigma_R / R = sigma_N / N0 EXACTLY, at every time. A first-order propagation would agree here
// and a wrong one would not; what makes the test worth running is that it holds to every digit
// rather than to a tolerance chosen to let it pass.
TEST(Uncertainty, IsExactRatherThanFirstOrder) {
  const NuclearData data = chain();
  const double atoms = 1.0e20;
  const double relative = 0.05;

  ResponseSpec spec;  // activity, becquerel
  const std::vector<AssayGroup> assays = {
      assay("2024-01-01", one(kParent, atoms, relative * atoms))};

  // With ONE seeded nuclide the response is proportional to its atom count however long the
  // chain beneath it, so a 5% assay is a 5% answer at every time and to every digit. That is
  // the sharp form of "exact": a first-order propagation would agree at t = 0 and drift, and a
  // tolerance loose enough to hide the drift would not be testing anything.
  for (const double time : {0.0, 1.0e5, 1.0e6, 1.0e7}) {
    const ResponseUncertainty u =
        responseUncertainty(data, assays, parseCalendarDate("2024-01-01"), time, spec);
    EXPECT_NEAR(u.relative, relative, 1.0e-12) << "at t=" << time;
    EXPECT_NEAR(u.sigma, relative * u.response, u.response * 1.0e-12) << "at t=" << time;
  }

  // The absolute figures against the closed form, on a chain short enough to have one: a lone
  // unstable nuclide with a stable terminator holds all the activity itself, so
  // R = lambda N0 e^{-lambda T} exactly.
  StoreArrays lone = synth::linearChain({kLambda});
  const NuclearData single = NuclearData::fromArrays(std::move(lone));
  for (const double time : {0.0, 1.0e5, 1.0e6, 1.0e7}) {
    const ResponseUncertainty u =
        responseUncertainty(single, assays, parseCalendarDate("2024-01-01"), time, spec);
    const double expected = kLambda * atoms * std::exp(-kLambda * time);
    EXPECT_NEAR(u.response, expected, expected * 1.0e-9) << "at t=" << time;
    EXPECT_NEAR(u.sigma, relative * expected, expected * 1.0e-9) << "at t=" << time;
  }
}

// The number a reader is meant to act on. Variance fractions are (g_i sigma_i)^2 normalized, so
// a row with a small share of the answer can dominate its error bar -- which is the whole reason
// the table is ordered by this column and not by share.
TEST(Uncertainty, RanksByVarianceAndNotByShare) {
  const NuclearData data = chain();
  Inventory inv;
  // A large, well-measured row and a small, badly-measured one.
  inv.add(kParent, 1.0e20, 0.01 * 1.0e20);
  inv.add(kMiddle, 1.0e19, 0.50 * 1.0e19);

  ResponseSpec spec;
  const std::vector<AssayGroup> assays = {assay("2024-01-01", inv)};
  const ResponseUncertainty u =
      responseUncertainty(data, assays, parseCalendarDate("2024-01-01"), 0.0, spec);

  ASSERT_EQ(u.seeds.size(), 2u);
  // At t = 0 the parent holds lambda*1e20 = 1e14 and the middle 5e-7*1e19 = 5e12, so the parent
  // is 95% of the answer; but the middle's 50% uncertainty makes it the larger error term.
  EXPECT_GT(u.seeds[0].share, 0.0);
  EXPECT_EQ(u.seeds[0].label, "Sb-100") << "ordered by variance, so the badly measured row leads";
  EXPECT_LT(u.seeds[0].share, u.seeds[1].share) << "and it is the SMALLER share of the answer";

  double fractions = 0.0;
  double variance = 0.0;
  for (const SeedUncertainty& seed : u.seeds) {
    fractions += seed.varianceFraction;
    variance += seed.sigmaContribution * seed.sigmaContribution;
  }
  EXPECT_NEAR(fractions, 1.0, 1.0e-12) << "variance fractions partition, shares of R do not";
  EXPECT_NEAR(std::sqrt(variance), u.sigma, u.sigma * 1.0e-12);
}

// --- the dated identity ------------------------------------------------------

// An assay carried forward by tau and asked about at T must give what the same assay asked about
// at T + tau gives, because the carry and the response interval share one decay matrix and their
// exponentials commute. That identity is the reason a dated file needs no new machinery, so it
// is asserted rather than assumed.
TEST(Uncertainty, CarriesAnAssayByRunningTheAdjointLonger) {
  const NuclearData data = chain();
  const Inventory inv = one(kParent, 1.0e20, 0.07 * 1.0e20);
  const double epoch = parseCalendarDate("2024-01-01");
  const double carry = 90.0 * 86400.0;
  const double time = 30.0 * 86400.0;

  ResponseSpec spec;
  spec.metric = Metric::Exposure;
  spec.unit = Unit::RoentgenPerHour;

  // Dated 90 days before the epoch, asked at 30 days after it.
  AssayGroup old;
  old.dateSeconds = epoch - carry;
  old.inventory = inv;
  old.label = "old sheet";
  old.dated = true;
  const ResponseUncertainty carried =
      responseUncertainty(data, std::vector<AssayGroup>{old}, epoch, time, spec);

  // The same sheet with no carry, asked 120 days out.
  AssayGroup fresh;
  fresh.dateSeconds = epoch;
  fresh.inventory = inv;
  fresh.label = "same sheet";
  fresh.dated = true;
  const ResponseUncertainty direct =
      responseUncertainty(data, std::vector<AssayGroup>{fresh}, epoch, time + carry, spec);

  EXPECT_NEAR(carried.response, direct.response, direct.response * 1.0e-9);
  EXPECT_NEAR(carried.sigma, direct.sigma, direct.sigma * 1.0e-9);
  EXPECT_DOUBLE_EQ(carried.seeds.front().carriedSeconds, carry);
}

// The identity that says the whole construction is sound: R assembled from per-assay adjoints
// equals what the forward path reports for the merged inventory. Two different solves, one
// number, and a disagreement would mean an assay was carried by the wrong amount.
TEST(Uncertainty, ReproducesTheOrdinaryTotalFromSeveralAssays) {
  const NuclearData data = chain();
  const double epoch = parseCalendarDate("2024-06-01");
  const std::vector<AssayGroup> assays = {assay("2023-01-10", one(kParent, 5.0e19, 1.0e19), "a"),
                                          assay("2024-06-01", one(kMiddle, 2.0e19, 1.0e18), "b")};

  ResponseSpec spec;
  const double time = 45.0 * 86400.0;
  const ResponseUncertainty u = responseUncertainty(data, assays, epoch, time, spec);

  // The forward path over the reconciled seed.
  const Reconciliation merged = reconcile(data, assays, epoch);
  const ResponseTable table =
      buildResponse(data, decay(data, merged.inventory, std::vector<double>{time}), spec);
  EXPECT_NEAR(u.response, table.totals.front(), table.totals.front() * 1.0e-9);
}

// --- what it says about itself -----------------------------------------------

TEST(Uncertainty, SaysHowMuchOfTheAnswerCarriesNoUncertaintyAtAll) {
  const NuclearData data = chain();
  Inventory inv;
  inv.add(kParent, 1.0e20, 0.05 * 1.0e20);
  inv.add(kMiddle, 1.0e20, 0.0);  // stated none

  ResponseSpec spec;
  const ResponseUncertainty u =
      responseUncertainty(data, std::vector<AssayGroup>{assay("2024-01-01", inv)},
                          parseCalendarDate("2024-01-01"), 0.0, spec);

  EXPECT_EQ(u.rowsWithSigma, 1);
  EXPECT_EQ(u.rowsWithoutSigma, 1);
  // An error bar propagated from rows holding part of the answer is not an error bar on the
  // answer, and the covered fraction is what lets a reader see that.
  EXPECT_GT(u.coveredFraction, 0.0);
  EXPECT_LT(u.coveredFraction, 1.0);
  // The row that stated nothing contributes nothing to the variance rather than a guess.
  const auto silent = std::find_if(u.seeds.begin(), u.seeds.end(),
                                   [](const SeedUncertainty& s) { return s.label == "Sb-100"; });
  ASSERT_NE(silent, u.seeds.end());
  EXPECT_DOUBLE_EQ(silent->sigmaContribution, 0.0);
  EXPECT_DOUBLE_EQ(silent->varianceFraction, 0.0);
}

// A row measured as zero still states an uncertainty, and it has to count.
//
// This is how a sheet reports a non-detect: 0 +/- MDA, meaning the true quantity is somewhere
// under the detection limit rather than known to be nothing. The response's error bar depends on
// that row through dR/dn0 exactly as any other row's does -- the derivative is a property of the
// decay matrix, not of the seed's own value -- so a propagation that skipped it would quote an
// error bar that silently ignored a measurement the sheet actually made.
//
// Seeded at the MIDDLE of the chain, with the sigma on the parent above it. The parent is not
// forward-reachable from the middle, so nothing but the uncertainty itself puts it in the index
// space: this fails on the pruning as well as on the summation if either drops the row.
TEST(Uncertainty, PropagatesARowMeasuredAsZeroWithAStatedUncertainty) {
  const NuclearData data = chain();
  const double sigma = 1.0e18;
  const double time = 1.0e6;

  Inventory nonDetect;
  nonDetect.add(kParent, 0.0, sigma);  // 0 +/- MDA, and upstream of everything seeded
  nonDetect.add(kMiddle, 1.0e20, 0.0);

  ResponseSpec spec;  // activity, becquerel
  const ResponseUncertainty u =
      responseUncertainty(data, std::vector<AssayGroup>{assay("2024-01-01", nonDetect)},
                          parseCalendarDate("2024-01-01"), time, spec);

  const auto parent = std::find_if(u.seeds.begin(), u.seeds.end(),
                                   [](const SeedUncertainty& s) { return s.label == "Sn-100"; });
  ASSERT_NE(parent, u.seeds.end()) << "a stated sigma has to reach the index space";
  EXPECT_DOUBLE_EQ(parent->seedAtoms, 0.0);
  EXPECT_DOUBLE_EQ(parent->share, 0.0) << "no atoms means no share of R, which stays true";
  EXPECT_GT(parent->sigmaContribution, 0.0) << "but its uncertainty is the whole error bar";
  EXPECT_GT(u.sigma, 0.0);
  EXPECT_EQ(u.rowsWithSigma, 1);

  // And the seed's own value does not decide whether its uncertainty counts. Giving the row one
  // atom instead of none changes R by a part in 1e20 and must leave sigma_R where it was: the
  // importance is the same vector either way, which is the reason the zero row belongs at all.
  Inventory oneAtom;
  oneAtom.add(kParent, 1.0, sigma);
  oneAtom.add(kMiddle, 1.0e20, 0.0);
  const ResponseUncertainty nudged =
      responseUncertainty(data, std::vector<AssayGroup>{assay("2024-01-01", oneAtom)},
                          parseCalendarDate("2024-01-01"), time, spec);
  EXPECT_NEAR(u.sigma, nudged.sigma, 1.0e-12 * nudged.sigma);
}

TEST(Uncertainty, RefusesAQuestionItCannotAnswer) {
  const NuclearData data = chain();
  const std::vector<AssayGroup> assays = {assay("2024-01-01", one(kParent, 1.0e20, 1.0e19))};
  const double epoch = parseCalendarDate("2024-01-01");

  ResponseSpec spec;
  EXPECT_THROW(responseUncertainty(data, {}, epoch, 0.0, spec), InputError);
  EXPECT_THROW(responseUncertainty(data, assays, epoch, -1.0, spec), InputError);
  // Carrying an assay backward is an inverse problem here for the same reason it is in
  // reconcile(), and is refused on the same terms.
  EXPECT_THROW(responseUncertainty(data, assays, parseCalendarDate("2023-01-01"), 0.0, spec),
               InputError);

  // The shares are instantaneous, as attributeToSeed's are: an accrued total needs the
  // integrated adjoint, which is a different solve.
  ResponseSpec interval;
  interval.unit = Unit::Decays;
  EXPECT_THROW(responseUncertainty(data, assays, epoch, 0.0, interval), InputError);

  ResponseSpec lines;
  lines.metric = Metric::Exposure;
  lines.unit = Unit::RoentgenPerHour;
  lines.aggregate = Aggregate::GammaLine;
  EXPECT_THROW(responseUncertainty(data, assays, epoch, 0.0, lines), InputError);
}

// An undated sheet's date is a placeholder zero, not a measurement. Carrying it to a real epoch
// would age it by however far that epoch is from 1970 -- fifty-four years for a 2024 date --
// and nothing in the answer would show why. Refused in both places that carry an assay.
TEST(Uncertainty, WillNotCarryAnUndatedAssayToADate) {
  const NuclearData data = chain();
  AssayGroup undated;
  undated.inventory = one(kParent, 1.0e20, 1.0e19);
  undated.label = "no date given";
  undated.dated = false;
  const std::vector<AssayGroup> assays = {undated};

  ResponseSpec spec;
  // Epoch zero is the only one an undated sheet can mean, and it works.
  EXPECT_NO_THROW(responseUncertainty(data, assays, 0.0, 0.0, spec));
  EXPECT_THROW(responseUncertainty(data, assays, parseCalendarDate("2024-01-01"), 0.0, spec),
               InputError);
  EXPECT_THROW(reconcile(data, assays, parseCalendarDate("2024-01-01")), InputError);
  EXPECT_NO_THROW(reconcile(data, assays, 0.0));
}

// --- reading the column ------------------------------------------------------

DatedInventory parseCsv(const std::string& text, const NuclearData& data) {
  std::istringstream in(text);
  return readInventoryDatedCsv(in, data, "test.csv");
}

TEST(UncertaintyColumn, TakesAnAbsoluteOrARelativeSpelling) {
  const NuclearData data = chain();
  const DatedInventory dated = parseCsv(
      "nuclide,quantity,unit,uncertainty\n"
      "Sn-100,1e20,atoms,5e18\n"
      "Sb-100,2e20,atoms,10%\n",
      data);
  ASSERT_EQ(dated.groups.size(), 1u);
  const Inventory& inv = dated.groups.front().inventory;
  EXPECT_TRUE(inv.hasUncertainties());
  for (const InventoryEntry& entry : inv.entries()) {
    if (entry.zaiKey == kParent.key()) {
      EXPECT_DOUBLE_EQ(entry.sigmaAtoms, 5.0e18);
    } else {
      EXPECT_DOUBLE_EQ(entry.sigmaAtoms, 2.0e19) << "10% of the row's own quantity";
    }
  }
}

// A sigma in becquerel is a sigma on the activity, and toAtoms is value*k with no offset in
// every branch -- so the same call converts both, and the basis needs no separate bookkeeping.
TEST(UncertaintyColumn, ConvertsThroughTheRowsOwnUnit) {
  const NuclearData data = chain();
  const DatedInventory dated =
      parseCsv("nuclide,quantity,unit,uncertainty\nSn-100,1e10,Bq,4%\n", data);
  const Inventory& inv = dated.groups.front().inventory;
  const double atoms = inv.atomsOf(kParent);
  ASSERT_GT(atoms, 0.0);
  EXPECT_NEAR(inv.entries().front().sigmaAtoms, 0.04 * atoms, atoms * 1.0e-12);
}

TEST(UncertaintyColumn, IsFoundByNameWhenAHeaderNamesIt) {
  const NuclearData data = chain();
  // The reason the header mapping exists: `assayed` shipped as the fourth field, so a file
  // wanting an uncertainty and no dates would otherwise need an empty placeholder column.
  const DatedInventory named =
      parseCsv("nuclide,quantity,unit,uncertainty\nSn-100,1e20,atoms,5%\n", data);
  EXPECT_FALSE(named.dated);
  EXPECT_DOUBLE_EQ(named.groups.front().inventory.entries().front().sigmaAtoms, 5.0e18);

  // And the positional form still reads as it always did, placeholder and all.
  const DatedInventory positional = parseCsv("Sn-100,1e20,atoms,,5%\n", data);
  EXPECT_FALSE(positional.dated);
  EXPECT_DOUBLE_EQ(positional.groups.front().inventory.entries().front().sigmaAtoms, 5.0e18);

  // A header that names neither nuclide nor quantity is not describing these columns, and the
  // reader falls back to positional rather than erroring -- a file that worked yesterday has to
  // keep working.
  const DatedInventory unknown = parseCsv("who,how much,in what\nSn-100,1e20,atoms\n", data);
  EXPECT_DOUBLE_EQ(unknown.groups.front().inventory.atomsOf(kParent), 1.0e20);
}

TEST(UncertaintyColumn, MergesTwoRowsOfOneNuclideInQuadrature) {
  const NuclearData data = chain();
  // Two containers, two measurements: the counts add and the variances add.
  const DatedInventory dated = parseCsv("Sn-100,1e20,atoms,,3e18\nSn-100,1e20,atoms,,4e18\n", data);
  const Inventory& inv = dated.groups.front().inventory;
  EXPECT_DOUBLE_EQ(inv.atomsOf(kParent), 2.0e20);
  EXPECT_DOUBLE_EQ(inv.entries().front().sigmaAtoms, 5.0e18) << "3-4-5, in quadrature";
}

TEST(UncertaintyColumn, RefusesAnUncertaintyThatIsNotOne) {
  const NuclearData data = chain();
  EXPECT_THROW(parseCsv("nuclide,quantity,unit,uncertainty\nSn-100,1e20,atoms,later\n", data),
               InputError);
  EXPECT_THROW(parseCsv("nuclide,quantity,unit,uncertainty\nSn-100,1e20,atoms,-5%\n", data),
               InputError);
  // A sigma LARGER than its quantity is accepted: that is what a measurement near a detection
  // limit honestly reports, and refusing it would refuse the rows an error bar is most for.
  EXPECT_NO_THROW(parseCsv("nuclide,quantity,unit,uncertainty\nSn-100,1e20,atoms,300%\n", data));
}

}  // namespace
}  // namespace nusift
