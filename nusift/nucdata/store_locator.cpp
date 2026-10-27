#include "nusift/nucdata/store_locator.hpp"

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <ostream>
#include <string>
#include <vector>

#include "nusift/core/error.hpp"

namespace nusift {
namespace {

namespace fs = std::filesystem;

constexpr const char* kModule = "nucdata store";
constexpr const char* kEnvVar = "NUSIFT_DATA_STORE";

// One place the search would look. `directory` is set when the path came from scanning a
// directory rather than from being named, which is the only case in which the search has a
// choice to make and therefore the only case worth reporting.
struct Candidate {
  std::string path;
  std::string directory;
};

// Any .h5 in `directory`, in a fixed order so the choice is deterministic: the stores NAMED as
// NuSIFT stores first -- nusift_*.h5, which is what staging and the install write -- and every
// other .h5 after them, each group sorted by name. The first is taken.
//
// Name order alone put a fixture first: a source tree holds the committed nusift_b8.1.h5 beside
// locally staged ones like am242m.h5 and test_store.h5, and "a" sorts before "n", so every report
// run from the tree -- and the Python package's own README example -- answered from a
// three-nuclide fixture with a warning nobody reads. Between two stores that are both named as
// NuSIFT's, the rule still cannot know which is meant, and still says what it passed over.
bool namedAsStore(const std::string& path) {
  return fs::path(path).filename().string().rfind("nusift_", 0) == 0;
}

void storesIn(const fs::path& directory, std::vector<Candidate>& out) {
  std::vector<std::string> found;
  std::error_code ec;
  if (!fs::is_directory(directory, ec)) {
    return;
  }
  for (const fs::directory_entry& entry : fs::directory_iterator(directory, ec)) {
    if (entry.is_regular_file(ec) && entry.path().extension() == ".h5") {
      found.push_back(entry.path().string());
    }
  }
  std::sort(found.begin(), found.end(), [](const std::string& a, const std::string& b) {
    const bool storeA = namedAsStore(a);
    const bool storeB = namedAsStore(b);
    if (storeA != storeB) {
      return storeA;
    }
    return a < b;
  });
  for (std::string& path : found) {
    out.push_back(Candidate{std::move(path), directory.string()});
  }
}

std::vector<Candidate> candidatesFor(const StoreSearch& search) {
  std::vector<Candidate> candidates;

  if (!search.explicitPath.empty()) {
    candidates.push_back(Candidate{search.explicitPath, {}});
    // An explicit path is a statement of intent: if it is wrong the user wants to hear that,
    // not to have a different store silently substituted.
    return candidates;
  }

  if (const char* fromEnv = std::getenv(kEnvVar); fromEnv != nullptr && *fromEnv != '\0') {
    candidates.push_back(Candidate{fromEnv, {}});
  }

  for (const std::string& extra : search.extraPaths) {
    candidates.push_back(Candidate{extra, {}});
  }

  // <prefix>/bin/nusift -> <prefix>/share/nusift. Derived from the executable rather than a
  // path baked in at build time, so a relocated install still finds its own data.
  if (!search.executablePath.empty()) {
    std::error_code ec;
    const fs::path exe = fs::absolute(search.executablePath, ec);
    if (!ec) {
      const fs::path prefix = exe.parent_path().parent_path();
      storesIn(prefix / "share" / "nusift", candidates);
    }
  }

  storesIn(fs::path("data"), candidates);

  return candidates;
}

}  // namespace

std::vector<std::string> storeSearchPaths(const StoreSearch& search) {
  std::vector<std::string> paths;
  for (const Candidate& candidate : candidatesFor(search)) {
    paths.push_back(candidate.path);
  }
  return paths;
}

std::string locateStore(const StoreSearch& search) {
  const std::vector<Candidate> candidates = candidatesFor(search);
  std::error_code ec;
  for (const Candidate& candidate : candidates) {
    if (!fs::is_regular_file(candidate.path, ec)) {
      continue;
    }
    if (!candidate.directory.empty() && search.warnings != nullptr) {
      std::vector<std::string> passedOver;
      for (const Candidate& other : candidates) {
        if (other.directory == candidate.directory && other.path != candidate.path) {
          passedOver.push_back(other.path);
        }
      }
      if (!passedOver.empty()) {
        *search.warnings << "nusift: " << candidate.directory << " holds "
                         << (passedOver.size() + 1) << " stores; using " << candidate.path
                         << ", the first by name (nusift_*.h5 ahead of the rest). Name a store "
                         << "explicitly (--store, or a path to NuclearData.open) or set " << kEnvVar
                         << " to choose. Passed over:";
        for (const std::string& path : passedOver) {
          *search.warnings << ' ' << path;
        }
        *search.warnings << '\n';
      }
    }
    return candidate.path;
  }

  std::string message = "no nuclear-data store found. Looked in:\n";
  if (candidates.empty()) {
    message += "  (nowhere -- no --store, no " + std::string(kEnvVar) +
               ", and no store beside the executable or under ./data)\n";
  } else {
    for (const Candidate& candidate : candidates) {
      message += "  " + candidate.path + "\n";
    }
  }
  message += "Point --store at one, set " + std::string(kEnvVar) +
             ", or build one from ENDF tapes:\n"
             "  nusift_stage_data --decay-dir <endf-decay-dir> -o data/nusift.h5";
  throw InputError(tagged(kModule, message));
}

}  // namespace nusift
