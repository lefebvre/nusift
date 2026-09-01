#include <gtest/gtest.h>

#include <sstream>
#include <string>
#include <vector>

#include "nusift/core/error.hpp"
#include "nusift/core/nuclide_name.hpp"
#include "nusift/engine/inventory.hpp"
#include "nusift/nucdata/coefficient_pack.hpp"
#include "nusift/nucdata/nuclear_data.hpp"
#include "synthetic_chain.hpp"

namespace nusift {
namespace {

// A header with every required field, so a test can vary one thing and leave the rest valid.
std::string headerWith(const std::string& overrides = "") {
  return "# pack: test-pack\n"
         "# version: 1\n"
         "# quantity: a screening index\n"
         "# unit: 1\n"
         "# basis: activity\n"
         "# domain: instant\n"
         "# progeny: excluded\n"
         "# source: invented for this test\n" +
         overrides;
}

CoefficientPack packFrom(const std::string& text) {
  std::istringstream in(text);
  return CoefficientPack::read(in, "test.csv");
}

TEST(CoefficientPack, ReadsTheHeaderAndTheRows) {
  const CoefficientPack pack = packFrom(headerWith() +
                                        "nuclide,coefficient,note\n"
                                        "Sn-100,2.5e-12,\n"
                                        "Sb-100,0,no limit applies\n");
  const PackProvenance& p = pack.provenance();
  EXPECT_EQ(p.name, "test-pack");
  EXPECT_EQ(p.version, "1");
  EXPECT_EQ(p.unit, "1");
  EXPECT_EQ(p.basis, PackBasis::Activity);
  EXPECT_EQ(p.domains, PackDomains::InstantOnly);
  EXPECT_FALSE(p.foldsProgeny);
  EXPECT_EQ(pack.size(), 2);

  const std::int64_t tin = requireNuclideName("Sn-100").key();
  EXPECT_TRUE(pack.covers(tin));
  EXPECT_DOUBLE_EQ(pack.coefficient(tin), 2.5e-12);

  // A zero written in the file and a nuclide the pack does not carry are different statements,
  // and the note is what keeps them apart.
  const std::int64_t antimony = requireNuclideName("Sb-100").key();
  EXPECT_TRUE(pack.covers(antimony));
  EXPECT_DOUBLE_EQ(pack.coefficient(antimony), 0.0);
  EXPECT_EQ(pack.note(antimony), "no limit applies");
  EXPECT_FALSE(pack.covers(requireNuclideName("Cs-137").key()));
  EXPECT_DOUBLE_EQ(pack.coefficient(requireNuclideName("Cs-137").key()), 0.0);
}

// Every header field is part of what an answer from the pack means, so a missing one is refused
// rather than defaulted. A default version would be the worst of them: it would put a number in
// a report under an edition nobody chose.
TEST(CoefficientPack, RefusesAHeaderMissingAnyRequiredField) {
  for (const char* omit :
       {"# pack: test-pack\n", "# version: 1\n", "# quantity: a screening index\n", "# unit: 1\n",
        "# basis: activity\n", "# domain: instant\n", "# progeny: excluded\n",
        "# source: invented for this test\n"}) {
    std::string header = headerWith();
    const std::size_t at = header.find(omit);
    ASSERT_NE(at, std::string::npos);
    header.erase(at, std::string(omit).size());
    EXPECT_THROW(packFrom(header + "nuclide,coefficient\nSn-100,1,\n"), InputError)
        << "omitting " << omit;
  }
}

TEST(CoefficientPack, RefusesRowsItCannotReadAsStated) {
  const std::string columns = "nuclide,coefficient,note\n";
  // A name that is not a nuclide, a coefficient that is not a number, a negative coefficient,
  // and the same nuclide twice -- the last because two coefficients for one nuclide is not a
  // stated quantity and choosing between them would be a guess.
  EXPECT_THROW(packFrom(headerWith() + columns + "Unobtainium-1,1,\n"), InputError);
  EXPECT_THROW(packFrom(headerWith() + columns + "Sn-100,about a picocurie,\n"), InputError);
  EXPECT_THROW(packFrom(headerWith() + columns + "Sn-100,-1,\n"), InputError);
  EXPECT_THROW(packFrom(headerWith() + columns + "Sn-100,1,\nSn-100,2,\n"), InputError);
  EXPECT_THROW(packFrom(headerWith() + "zai,value\n100100,1\n"), InputError);
  EXPECT_THROW(packFrom(headerWith() + columns), InputError) << "a pack with no coefficients";
}

TEST(CoefficientPack, RefusesABasisOrDomainItDoesNotKnow) {
  EXPECT_THROW(packFrom(headerWith("# basis: per curie\n") + "nuclide,coefficient\nSn-100,1\n"),
               InputError);
  EXPECT_THROW(packFrom(headerWith("# domain: whenever\n") + "nuclide,coefficient\nSn-100,1\n"),
               InputError);
}

// --- folded progeny ----------------------------------------------------------

std::string foldedHeader() {
  std::string header = headerWith();
  const std::size_t at = header.find("# progeny: excluded\n");
  header.replace(at, std::string("# progeny: excluded\n").size(), "# progeny: folded\n");
  return header;
}

// The two declarations have to match what the rows actually do, in both directions: a pack that
// says `excluded` while folding, and one that says `folded` while folding nothing, are each
// describing a file that is not there.
TEST(CoefficientPack, TheProgenyDeclarationHasToMatchTheRows) {
  EXPECT_THROW(packFrom(headerWith() + "nuclide,coefficient,folded\nSn-100,1,Sb-100\n"),
               InputError);
  EXPECT_THROW(packFrom(foldedHeader() + "nuclide,coefficient,folded\nSn-100,1,\n"), InputError);
  EXPECT_THROW(packFrom(headerWith("# progeny: rolled up\n") + "nuclide,coefficient\nSn-100,1\n"),
               InputError);
}

// The case the fold list exists for, and the case that makes it inventory-dependent: a daughter
// that ALSO holds a row of its own. Shipped with its parent it is inside the parent's number;
// shipped alone its own applies. SSR-6 does this 36 times.
TEST(CoefficientPack, AFoldedDaughterIsAccountedForByItsParentOnlyWhenTheParentIsSeeded) {
  const CoefficientPack pack = packFrom(foldedHeader() +
                                        "nuclide,coefficient,folded,note\n"
                                        "Sn-100,4e-12,Sb-100,\n"
                                        "Sb-100,8e-12,,\n");
  const std::int64_t parent = requireNuclideName("Sn-100").key();
  const std::int64_t daughter = requireNuclideName("Sb-100").key();
  EXPECT_TRUE(pack.provenance().foldsProgeny);
  ASSERT_EQ(pack.foldedInto(daughter).size(), 1u);
  EXPECT_EQ(pack.foldedInto(daughter).front(), parent);

  Inventory together;
  together.add(Zai::fromKey(parent), 1.0e20);
  Inventory aloneSeed;
  aloneSeed.add(Zai::fromKey(daughter), 1.0e20);

  EXPECT_EQ(pack.coverageOf(daughter, together), PackCoverage::Folded);
  EXPECT_EQ(pack.coverageOf(daughter, aloneSeed), PackCoverage::Own);
  EXPECT_EQ(pack.coverageOf(parent, together), PackCoverage::Own);

  StoreArrays arrays = synth::linearChain({1.0e-3, 5.0e-4});
  const NuclearData data = NuclearData::fromArrays(std::move(arrays));
  const std::vector<double> withParent = pack.weights(data, together);
  const std::vector<double> withoutParent = pack.weights(data, aloneSeed);
  const int index = data.indexOfKey(daughter);
  ASSERT_GE(index, 0);
  EXPECT_DOUBLE_EQ(withParent[static_cast<std::size_t>(index)], 0.0)
      << "weighting it under its parent would count it twice";
  EXPECT_GT(withoutParent[static_cast<std::size_t>(index)], 0.0);
}

// Chains nest, so one daughter can sit under several parents -- SSR-6 folds Tl-208 into Bi-212,
// Pb-212, Ra-224 and Th-228. Any one of them being present accounts for it.
TEST(CoefficientPack, ADaughterMayBeFoldedIntoSeveralParents) {
  const CoefficientPack pack = packFrom(foldedHeader() +
                                        "nuclide,coefficient,folded\n"
                                        "Sn-100,1e-12,Sb-100\n"
                                        "Sn-102,2e-12,Sb-100\n");
  const std::int64_t daughter = requireNuclideName("Sb-100").key();
  EXPECT_EQ(pack.foldedInto(daughter).size(), 2u);

  Inventory viaSecond;
  viaSecond.add(requireNuclideName("Sn-102"), 1.0e20);
  EXPECT_EQ(pack.coverageOf(daughter, viaSecond), PackCoverage::Folded);

  Inventory neither;
  neither.add(requireNuclideName("Cs-137"), 1.0e20);
  EXPECT_EQ(pack.coverageOf(daughter, neither), PackCoverage::None)
      << "no parent present and no row of its own";
}

// A fold naming a parent the pack does not carry cannot account for anything, and saying it did
// would overstate the coverage figure by exactly the rows that are missing.
TEST(CoefficientPack, AFoldIntoAParentThePackLacksDoesNotCover) {
  const CoefficientPack pack = packFrom(foldedHeader() +
                                        "nuclide,coefficient,folded\n"
                                        "Sn-100,1e-12,Sb-100\n");
  Inventory seed;
  seed.add(requireNuclideName("Sn-100"), 1.0e20);
  EXPECT_EQ(pack.coverageOf(requireNuclideName("Sb-100").key(), seed), PackCoverage::Folded);

  const CoefficientPack orphan = packFrom(foldedHeader() +
                                          "nuclide,coefficient,folded\n"
                                          "Sn-102,1e-12,Sb-100\n");
  Inventory other;
  other.add(requireNuclideName("Sn-100"), 1.0e20);
  EXPECT_EQ(orphan.coverageOf(requireNuclideName("Sb-100").key(), other), PackCoverage::None);
}

// --- concentration ------------------------------------------------------------

std::string concentrationHeader(const char* per = "m3") {
  std::string header = headerWith();
  const std::size_t at = header.find("# basis: activity\n");
  header.replace(at, std::string("# basis: activity\n").size(),
                 std::string("# basis: concentration\n# per: ") + per + "\n");
  return header;
}

// An inventory is atoms; a concentration is atoms over an extent, and nothing in the material
// says what that extent is. So a concentration pack requires it, and every other basis refuses
// it -- a number per becquerel does not become a different number when told how much space the
// becquerels occupy.
TEST(CoefficientPack, AConcentrationPackNeedsTheExtentAndTheRestRefuseIt) {
  const CoefficientPack cloud =
      packFrom(concentrationHeader() + "nuclide,coefficient\nSn-100,2e-14\n");
  EXPECT_EQ(cloud.provenance().basis, PackBasis::Concentration);
  EXPECT_EQ(cloud.provenance().per, "m3");

  StoreArrays arrays = synth::linearChain({1.0e-3, 5.0e-4});
  const NuclearData data = NuclearData::fromArrays(std::move(arrays));
  const Inventory empty;
  const int index = data.indexOfKey(requireNuclideName("Sn-100").key());
  ASSERT_GE(index, 0);

  EXPECT_THROW(cloud.weights(data, empty), InputError) << "no extent given";
  EXPECT_THROW(cloud.weights(data, empty, PackExtent{1.0e6, "m2"}), InputError)
      << "a volume pack given an area";

  // The extent divides: the same inventory in twice the volume is half the concentration.
  const double small =
      cloud.weights(data, empty, PackExtent{1.0e6, "m3"})[static_cast<std::size_t>(index)];
  const double large =
      cloud.weights(data, empty, PackExtent{2.0e6, "m3"})[static_cast<std::size_t>(index)];
  EXPECT_DOUBLE_EQ(small, 2.0e-14 * data.decayConstant(index) / 1.0e6);
  EXPECT_DOUBLE_EQ(large, small / 2.0);

  const CoefficientPack perBq = packFrom(headerWith() + "nuclide,coefficient\nSn-100,1\n");
  EXPECT_THROW(perBq.weights(data, empty, PackExtent{1.0e6, "m3"}), InputError);
}

TEST(CoefficientPack, TheConcentrationDenominatorIsDeclaredAndChecked) {
  // Deposition per square metre and a cloud per cubic metre are different questions whose
  // coefficients look alike, so `per` is required and is not free text.
  EXPECT_THROW(
      packFrom(headerWith("# basis: concentration\n") + "nuclide,coefficient\nSn-100,1e-14\n"),
      InputError)
      << "concentration with no per";
  EXPECT_THROW(packFrom(concentrationHeader("furlong") + "nuclide,coefficient\nSn-100,1e-14\n"),
               InputError);
  EXPECT_THROW(packFrom(headerWith("# per: m3\n") + "nuclide,coefficient\nSn-100,1\n"), InputError)
      << "per given for an activity-basis pack";

  for (const char* per : {"m2", "m3", "kg"}) {
    const CoefficientPack pack =
        packFrom(concentrationHeader(per) + "nuclide,coefficient\nSn-100,1e-14\n");
    EXPECT_EQ(pack.provenance().per, per);
  }
}

// --- the basis ---------------------------------------------------------------

// What the coefficient multiplies is declared rather than inferred, and getting it wrong scales
// every answer by a decay constant -- so the two bases have to produce visibly different
// weights from the same number.
TEST(CoefficientPack, TheBasisDecidesWhatTheCoefficientMultiplies) {
  StoreArrays arrays = synth::linearChain({1.0e-3, 5.0e-4});
  const NuclearData data = NuclearData::fromArrays(std::move(arrays));
  const Inventory empty;
  const int index = data.indexOfKey(requireNuclideName("Sn-100").key());
  ASSERT_GE(index, 0);

  const CoefficientPack perBq = packFrom(headerWith() + "nuclide,coefficient\nSn-100,1\n");
  const CoefficientPack perAtom =
      packFrom(headerWith("# basis: atoms\n") + "nuclide,coefficient\nSn-100,1\n");
  EXPECT_DOUBLE_EQ(perAtom.weights(data, empty)[static_cast<std::size_t>(index)], 1.0);
  EXPECT_DOUBLE_EQ(perBq.weights(data, empty)[static_cast<std::size_t>(index)],
                   data.decayConstant(index));
}

}  // namespace
}  // namespace nusift
