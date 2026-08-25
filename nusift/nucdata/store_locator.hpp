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
// A directory in step 4 or 5 that holds more than one store is ambiguous by nature. The first
// by name is used, and the choice is echoed in every report header -- but a source tree with
// locally staged fixtures beside the committed store makes that the wrong one, silently, so
// the search also says what it passed over when given somewhere to say it.
//
#include <iosfwd>
#include <string>
#include <vector>

namespace nusift {

struct StoreSearch {
  std::string explicitPath;             // from --store; empty if not given
  std::vector<std::string> extraPaths;  // caller-injected candidates
  std::string executablePath;           // argv[0], for deriving the install prefix
  // Where to report a directory that held more than one store, naming the one used and the
  // ones passed over. Null means nobody is listening, which is the library default; the CLI
  // hands in stderr and the Python binding routes it through the warnings module.
  std::ostream* warnings = nullptr;
};

// Returns the first store found. Throws InputError listing everywhere it looked when there
// is none.
std::string locateStore(const StoreSearch& search);

// The candidate paths, in order, without touching the filesystem. Exposed so the error
// message and the search can never disagree about what was tried.
std::vector<std::string> storeSearchPaths(const StoreSearch& search);

}  // namespace nusift
