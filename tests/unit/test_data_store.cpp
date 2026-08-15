#include <gtest/gtest.h>
#include <hdf5.h>

#include <filesystem>
#include <string>

#include "nusift/core/error.hpp"
#include "nusift/nucdata/data_store.hpp"
#include "nusift/nucdata/nuclear_data.hpp"
#include "synthetic_chain.hpp"

namespace nusift {
namespace {

namespace fs = std::filesystem;

constexpr const char* kVersionAttr = "nusift_store_version";

// Overwrite the version attribute of an existing store. There is no production path that
// writes anything but the current version, so forging one here is the only way to exercise
// the reader's forward-compatibility guard at all.
void stampVersionAttribute(const std::string& path, int version) {
  const hid_t file = H5Fopen(path.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  ASSERT_GE(file, 0);
  H5Adelete(file, kVersionAttr);
  const hid_t space = H5Screate(H5S_SCALAR);
  const hid_t attr =
      H5Acreate2(file, kVersionAttr, H5T_NATIVE_INT, space, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attr, H5T_NATIVE_INT, &version);
  H5Aclose(attr);
  H5Sclose(space);
  H5Fclose(file);
}

void removeVersionAttribute(const std::string& path) {
  const hid_t file = H5Fopen(path.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  ASSERT_GE(file, 0);
  H5Adelete(file, kVersionAttr);
  H5Fclose(file);
}

// A unique path under the system temp directory, removed on destruction so a failing test
// cannot leave a stale store behind for the next run to read.
class TempStore {
public:
  explicit TempStore(const std::string& name)
      : path_(fs::temp_directory_path() / ("nusift_test_" + name + ".h5")) {
    fs::remove(path_);
  }
  ~TempStore() {
    std::error_code ignored;
    fs::remove(path_, ignored);
  }
  TempStore(const TempStore&) = delete;
  TempStore& operator=(const TempStore&) = delete;

  std::string str() const { return path_.string(); }

private:
  fs::path path_;
};

StoreArrays sampleStore() {
  StoreArrays a = synth::linearChain({1.0e-3, 4.0e-4});
  a.provenance.library = "ENDF/B-VIII.1";
  a.provenance.createdUtc = "2026-08-13T00:00:00Z";
  a.provenance.nusiftVersion = "0.1.0";
  a.provenance.stagedTapeCount = 3;
  a.provenance.decaySource = DataSource::Endf;
  a.provenance.linesSource = DataSource::Endf;
  a.provenance.yieldsSource = DataSource::None;

  a.awr = {99.1, 99.2, 99.3};
  a.emEnergyEv = {6.6e5, 1.2e6, 0.0};
  a.continuumPhotonEv = {0.0, 3.0e5, 0.0};

  // Two lines on the first nuclide, one on the second, none on the stable terminator --
  // enough for the CSR offsets to be wrong in a detectable way if they are wrong at all.
  a.lineOffset = {0, 2, 3, 3};
  a.lineEnergyEv = {661657.0, 32000.0, 1173228.0};
  a.lineIntensity = {0.851, 0.056, 0.9985};
  a.lineStyp = {0, 9, 0};

  return a;
}

TEST(DataStore, RoundTripsEveryField) {
  const TempStore path("roundtrip");
  const StoreArrays original = sampleStore();
  writeStore(path.str(), original);

  const StoreArrays read = readStore(path.str());

  EXPECT_EQ(read.provenance.version, kStoreVersion);
  EXPECT_EQ(read.provenance.library, "ENDF/B-VIII.1");
  EXPECT_EQ(read.provenance.createdUtc, "2026-08-13T00:00:00Z");
  EXPECT_EQ(read.provenance.nusiftVersion, "0.1.0");
  EXPECT_EQ(read.provenance.stagedTapeCount, 3);
  EXPECT_EQ(read.provenance.decaySource, DataSource::Endf);
  EXPECT_EQ(read.provenance.linesSource, DataSource::Endf);
  EXPECT_EQ(read.provenance.yieldsSource, DataSource::None);

  EXPECT_EQ(read.nuclideKey, original.nuclideKey);
  EXPECT_EQ(read.halfLife, original.halfLife);
  EXPECT_EQ(read.awr, original.awr);
  EXPECT_EQ(read.emEnergyEv, original.emEnergyEv);
  EXPECT_EQ(read.continuumPhotonEv, original.continuumPhotonEv);
  EXPECT_EQ(read.modeOffset, original.modeOffset);
  EXPECT_EQ(read.modeRtyp, original.modeRtyp);
  EXPECT_EQ(read.modeBranching, original.modeBranching);
}

// The line CSR is the schema's headline addition, and an off-by-one in the offsets would
// attribute one nuclide's photons to another -- a wrong answer that still looks plausible.
TEST(DataStore, RoundTripsThePhotonLineCsr) {
  const TempStore path("lines");
  writeStore(path.str(), sampleStore());
  const StoreArrays read = readStore(path.str());

  EXPECT_EQ(read.lineOffset, (std::vector<int>{0, 2, 3, 3}));
  ASSERT_EQ(read.lineEnergyEv.size(), 3u);
  EXPECT_DOUBLE_EQ(read.lineEnergyEv[0], 661657.0);
  EXPECT_DOUBLE_EQ(read.lineIntensity[0], 0.851);
  EXPECT_EQ(read.lineStyp[1], 9);

  // And the spans NuclearData hands out have to land on the right nuclide.
  const NuclearData data = NuclearData::fromArrays(read);
  const LineSpectrum first = data.lines(0);
  ASSERT_EQ(first.size(), 2u);
  EXPECT_DOUBLE_EQ(first[0].energyEv, 661657.0);
  EXPECT_EQ(first[1].type, SpectrumType::XrayOrAnnih);
  EXPECT_EQ(data.lines(1).size(), 1u);
  EXPECT_EQ(data.lines(2).size(), 0u);
}

// An absent optional dataset means "never staged", which is a different statement from
// "staged as zero" and has to survive the round trip as such.
TEST(DataStore, OmitsAbsentOptionalFieldsRatherThanZeroFillingThem) {
  const TempStore path("minimal");
  StoreArrays minimal = synth::linearChain({1.0e-3});
  writeStore(path.str(), minimal);

  const StoreArrays read = readStore(path.str());
  EXPECT_TRUE(read.awr.empty());
  EXPECT_TRUE(read.lineEnergyEv.empty());
  EXPECT_TRUE(read.emEnergyEv.empty());

  const NuclearData data = NuclearData::fromArrays(read);
  EXPECT_FALSE(data.hasAtomicWeights());
  EXPECT_FALSE(data.hasPhotonLines());
  EXPECT_EQ(data.lines(0).size(), 0u);
}

TEST(DataStore, RejectsAFileThatIsNotAStore) {
  EXPECT_THROW(readStore("definitely_does_not_exist_12345.h5"), NusiftError);
}

// Reading a newer store partially would be worse than refusing it: the fields this build
// does not know are exactly the ones that would change the answer. Written by stamping a
// future version onto a real store, since writeStore will only ever emit the current one.
TEST(DataStore, RejectsANewerSchemaVersion) {
  const TempStore path("future");
  writeStore(path.str(), sampleStore());
  ASSERT_NO_THROW(readStore(path.str()));

  stampVersionAttribute(path.str(), kStoreVersion + 1);

  try {
    readStore(path.str());
    FAIL() << "expected a newer store to be refused";
  } catch (const NusiftError& e) {
    const std::string what = e.what();
    EXPECT_NE(what.find(std::to_string(kStoreVersion + 1)), std::string::npos) << what;
    EXPECT_NE(what.find("restage"), std::string::npos) << what;
  }
}

// A file that is valid HDF5 but carries no version attribute is somebody else's file, and
// saying so is much more useful than reporting a missing dataset three reads later.
TEST(DataStore, RejectsAnHdf5FileThatIsNotANusiftStore) {
  const TempStore path("foreign");
  writeStore(path.str(), sampleStore());
  removeVersionAttribute(path.str());

  try {
    readStore(path.str());
    FAIL() << "expected a versionless file to be refused";
  } catch (const NusiftError& e) {
    EXPECT_NE(std::string(e.what()).find("not a NuSIFT data store"), std::string::npos);
  }
}

// validateStoreArrays is the single contract both a file and a hand-built test chain are
// held to, so these malformed cases must fail identically whichever way they arrive.
TEST(StoreArraysValidation, RejectsUnsortedNuclideAxis) {
  StoreArrays a;
  a.nuclideKey = {551370, 501000};  // descending
  a.halfLife = {1.0, 2.0};
  a.modeOffset = {0, 0, 0};
  EXPECT_THROW(validateStoreArrays(a), NusiftError);
}

TEST(StoreArraysValidation, RejectsDuplicateNuclides) {
  StoreArrays a;
  a.nuclideKey = {551370, 551370};
  a.halfLife = {1.0, 2.0};
  a.modeOffset = {0, 0, 0};
  EXPECT_THROW(validateStoreArrays(a), NusiftError);
}

TEST(StoreArraysValidation, RejectsMismatchedArrayLengths) {
  StoreArrays a;
  a.nuclideKey = {501000, 511000};
  a.halfLife = {1.0};  // one short
  a.modeOffset = {0, 0, 0};
  EXPECT_THROW(validateStoreArrays(a), NusiftError);
}

TEST(StoreArraysValidation, RejectsBrokenCsrOffsets) {
  StoreArrays a;
  a.nuclideKey = {501000, 511000};
  a.halfLife = {1.0, 2.0};
  a.modeOffset = {0, 1, 0};  // decreasing
  a.modeRtyp = {1.0};
  a.modeBranching = {1.0};
  a.modeFinalState = {0};
  a.modeIsFission = {0};
  EXPECT_THROW(validateStoreArrays(a), NusiftError);
}

TEST(StoreArraysValidation, RejectsCsrOffsetsThatDoNotSpanTheValues) {
  StoreArrays a;
  a.nuclideKey = {501000, 511000};
  a.halfLife = {1.0, 2.0};
  a.modeOffset = {0, 1, 1};  // claims one mode
  a.modeRtyp = {1.0, 2.0};   // but two are present
  a.modeBranching = {1.0, 1.0};
  a.modeFinalState = {0, 0};
  a.modeIsFission = {0, 0};
  EXPECT_THROW(validateStoreArrays(a), NusiftError);
}

// The continuum field is what keeps "this nuclide emits no photons" distinguishable from
// "this nuclide's photons are in a spectrum NuSIFT does not model". Both would otherwise
// present as a zero exposure contribution, and only one of them is correct.
TEST(NuclearDataPhotons, ReportsTheUnmodeledContinuumFraction) {
  const NuclearData data = NuclearData::fromArrays(sampleStore());

  // First nuclide: lines only, no continuum.
  EXPECT_DOUBLE_EQ(data.continuumPhotonEv(0), 0.0);
  EXPECT_DOUBLE_EQ(data.unmodeledPhotonFraction(0), 0.0);

  // Second: one 1.173 MeV line at ~1.0 intensity plus 3.0e5 eV of continuum, so roughly a
  // fifth of its photon energy is invisible to the model.
  EXPECT_GT(data.continuumPhotonEv(1), 0.0);
  const double fraction = data.unmodeledPhotonFraction(1);
  EXPECT_GT(fraction, 0.15);
  EXPECT_LT(fraction, 0.30);

  // Third has neither, and must read as 0 rather than as a division by zero.
  EXPECT_DOUBLE_EQ(data.unmodeledPhotonFraction(2), 0.0);
}

}  // namespace
}  // namespace nusift
