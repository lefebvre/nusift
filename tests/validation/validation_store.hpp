#pragma once
//
// Opens the committed nuclear-data store for the validation suite.
//
// The unit suite is built the other way round: it needs no data files at all and constructs
// every chain from StoreArrays, so it passes on a fresh clone with nothing staged. Validation
// asks a different question -- does the evaluation NuSIFT actually ships reproduce published
// physics -- and that question is only meaningful against a specific file.
//
// The path is compiled in rather than passed through the environment so the binary can be run
// directly, which is what anyone debugging a failure does. NUSIFT_VALIDATION_STORE overrides
// it, deliberately NOT the general NUSIFT_DATA_STORE: a developer with that set to some other
// evaluation would otherwise silently swap the file underneath census counts whose entire
// meaning is "the store this repository commits".
//
// A missing store FAILS rather than skips. It is a committed file, so its absence means a
// broken checkout, and a suite that quietly skipped would report success for having checked
// nothing -- the exact failure a validation suite exists to prevent. Running without it is
// what the CTest label is for: `ctest -LE validation`.
//
#include <cstdlib>
#include <filesystem>
#include <string>

#include "nusift/core/error.hpp"
#include "nusift/nucdata/nuclear_data.hpp"

namespace nusift::validation {

inline std::string validationStorePath() {
  if (const char* override = std::getenv("NUSIFT_VALIDATION_STORE")) {
    return override;
  }
  return NUSIFT_COMMITTED_STORE;
}

// One store for the whole binary. Opening it costs an HDF5 read and a chain closure, and
// nothing any test does can modify it.
//
// A missing store throws rather than aborting: GoogleTest turns the exception into an ordinary
// failure, so every test reports the same explanation and the run finishes with a nonzero exit
// instead of a crash on whichever test happened to be scheduled first.
inline const NuclearData& committedStore() {
  static const NuclearData data = [] {
    const std::string path = validationStorePath();
    if (!std::filesystem::exists(path)) {
      throw InputError(
          "the committed nuclear-data store is missing: " + path +
          "\nThis file is tracked in git; its absence means an incomplete checkout."
          "\nSet NUSIFT_VALIDATION_STORE to point elsewhere, or exclude these tests with"
          " `ctest -LE validation`.");
    }
    return NuclearData::open(path);
  }();
  return data;
}

}  // namespace nusift::validation
