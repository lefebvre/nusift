#include <gtest/gtest.h>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

#include "nusift/nucdata/store_locator.hpp"

namespace nusift {
namespace {

namespace fs = std::filesystem;

// A scratch install prefix -- <root>/bin/nusift beside <root>/share/nusift -- holding whatever
// .h5 names a test asks for. The locator looks only at names and existence, so an empty file
// stands in for a store.
class ScratchPrefix {
public:
  ScratchPrefix() : root_(fs::temp_directory_path() / "nusift_store_locator_test") {
    fs::remove_all(root_);
    fs::create_directories(root_ / "share" / "nusift");
  }
  ~ScratchPrefix() {
    std::error_code ec;
    fs::remove_all(root_, ec);
  }
  ScratchPrefix(const ScratchPrefix&) = delete;
  ScratchPrefix& operator=(const ScratchPrefix&) = delete;

  std::string store(const std::string& name) const {
    const fs::path path = root_ / "share" / "nusift" / name;
    std::ofstream(path).put('\n');
    return path.string();
  }

  StoreSearch search() const {
    StoreSearch s;
    s.executablePath = (root_ / "bin" / "nusift").string();
    return s;
  }

private:
  fs::path root_;
};

// $NUSIFT_DATA_STORE outranks the prefix, so a developer with one set would be testing their
// own environment here rather than the locator.
void requireNoEnvironmentStore() {
  if (const char* env = std::getenv("NUSIFT_DATA_STORE"); env != nullptr && *env != '\0') {
    GTEST_SKIP() << "NUSIFT_DATA_STORE is set";
  }
}

// A directory holding two stores is ambiguous. The rule is fixed -- first by name -- but a fixed
// rule silently applied is exactly how a source tree with a staged fixture beside the committed
// store ends up answering from the fixture. The choice has to be said, and so does what it
// passed over.
TEST(StoreLocator, SaysWhichStoreItChoseWhenADirectoryHoldsSeveral) {
  requireNoEnvironmentStore();
  const ScratchPrefix scratch;
  const std::string first = scratch.store("a_fixture.h5");
  const std::string second = scratch.store("b_evaluation.h5");

  std::ostringstream warnings;
  StoreSearch search = scratch.search();
  search.warnings = &warnings;

  EXPECT_EQ(locateStore(search), first);
  const std::string said = warnings.str();
  EXPECT_NE(said.find("2 stores"), std::string::npos) << said;
  EXPECT_NE(said.find("first by name"), std::string::npos) << said;
  EXPECT_NE(said.find(second), std::string::npos) << "the store passed over is named: " << said;
  EXPECT_NE(said.find("--store"), std::string::npos) << "and the way to choose: " << said;
}

TEST(StoreLocator, IsSilentWhenTheChoiceIsUnambiguous) {
  requireNoEnvironmentStore();
  const ScratchPrefix scratch;
  const std::string only = scratch.store("evaluation.h5");

  std::ostringstream warnings;
  StoreSearch search = scratch.search();
  search.warnings = &warnings;

  EXPECT_EQ(locateStore(search), only);
  EXPECT_TRUE(warnings.str().empty()) << warnings.str();
}

// The library default: nobody listening, and the same store chosen.
TEST(StoreLocator, ChoosesTheSameStoreWhenNobodyIsListening) {
  requireNoEnvironmentStore();
  const ScratchPrefix scratch;
  const std::string first = scratch.store("a.h5");
  scratch.store("b.h5");
  EXPECT_EQ(locateStore(scratch.search()), first);
}

// An explicit path is a statement of intent, and there is nothing to warn about: the user
// chose, so the store beside the one they named is not a candidate at all.
TEST(StoreLocator, AnExplicitPathIsNeverSecondGuessed) {
  const ScratchPrefix scratch;
  scratch.store("a.h5");
  const std::string chosen = scratch.store("b.h5");

  std::ostringstream warnings;
  StoreSearch search = scratch.search();
  search.explicitPath = chosen;
  search.warnings = &warnings;

  EXPECT_EQ(locateStore(search), chosen);
  EXPECT_TRUE(warnings.str().empty()) << warnings.str();
}

}  // namespace
}  // namespace nusift
