#pragma once
/**
 * @file
 * @brief Finding the nuclear-data store.
 * @ingroup nucdata
 */
//
// One search order, in one place, used by every entry point. Duplicating it -- once in the
// CLI and again in a language binding -- is how two front ends end up quietly reading
// different evaluations and reporting different answers for the same question.
//
// The order, first match wins:
//   1. an explicit path from --store or an equivalent argument
//   2. $NUSIFT_DATA_STORE
//   3. caller-supplied extra paths (a language binding passes its packaged resource here)
//   4. <install prefix>/share/nusift/*.h5, derived from the running executable
//   5. ./data/*.h5, for working in a source tree
//
// When nothing is found the error names every place that was searched and gives the command
// that would produce a store, because "no data store found" on its own leaves a new user
// with nowhere to go.
//
#include <string>
#include <vector>

namespace nusift {

struct StoreSearch {
  std::string explicitPath;             // from --store; empty if not given
  std::vector<std::string> extraPaths;  // caller-injected candidates
  std::string executablePath;           // argv[0], for deriving the install prefix
};

// Returns the first store found. Throws InputError listing everywhere it looked when there
// is none.
std::string locateStore(const StoreSearch& search);

// The candidate paths, in order, without touching the filesystem. Exposed so the error
// message and the search can never disagree about what was tried.
std::vector<std::string> storeSearchPaths(const StoreSearch& search);

}  // namespace nusift
