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

// Any .h5 in `directory`, sorted by name so the choice is deterministic. Which one is taken is
// the first: not because the first is likelier to be right, but because no ordering rule can
// know whether `nusift_b9.0.h5` or `am242m.h5` is the one meant, and a fixed rule that is
// echoed and warned about beats a clever one that is neither.
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
  std::sort(found.begin(), found.end());
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
                         << ", the first by name. Pass --store or set " << kEnvVar
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
