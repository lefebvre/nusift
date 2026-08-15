#include "nusift/nucdata/store_locator.hpp"

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <string>
#include <vector>

#include "nusift/core/error.hpp"

namespace nusift {
namespace {

namespace fs = std::filesystem;

constexpr const char* kModule = "nucdata store";
constexpr const char* kEnvVar = "NUSIFT_DATA_STORE";

// Any .h5 in `directory`, newest name last so the choice is at least deterministic. A
// directory holding two evaluations is ambiguous by nature; the resolved path is echoed in
// every report header so which one was used is never in doubt.
std::vector<std::string> storesIn(const fs::path& directory) {
  std::vector<std::string> found;
  std::error_code ec;
  if (!fs::is_directory(directory, ec)) {
    return found;
  }
  for (const fs::directory_entry& entry : fs::directory_iterator(directory, ec)) {
    if (entry.is_regular_file(ec) && entry.path().extension() == ".h5") {
      found.push_back(entry.path().string());
    }
  }
  std::sort(found.begin(), found.end());
  return found;
}

}  // namespace

std::vector<std::string> storeSearchPaths(const StoreSearch& search) {
  std::vector<std::string> candidates;

  if (!search.explicitPath.empty()) {
    candidates.push_back(search.explicitPath);
    // An explicit path is a statement of intent: if it is wrong the user wants to hear that,
    // not to have a different store silently substituted.
    return candidates;
  }

  if (const char* fromEnv = std::getenv(kEnvVar); fromEnv != nullptr && *fromEnv != '\0') {
    candidates.emplace_back(fromEnv);
  }

  for (const std::string& extra : search.extraPaths) {
    candidates.push_back(extra);
  }

  // <prefix>/bin/nusift -> <prefix>/share/nusift. Derived from the executable rather than a
  // path baked in at build time, so a relocated install still finds its own data.
  if (!search.executablePath.empty()) {
    std::error_code ec;
    const fs::path exe = fs::absolute(search.executablePath, ec);
    if (!ec) {
      const fs::path prefix = exe.parent_path().parent_path();
      for (const std::string& path : storesIn(prefix / "share" / "nusift")) {
        candidates.push_back(path);
      }
    }
  }

  for (const std::string& path : storesIn(fs::path("data"))) {
    candidates.push_back(path);
  }

  return candidates;
}

std::string locateStore(const StoreSearch& search) {
  const std::vector<std::string> candidates = storeSearchPaths(search);
  std::error_code ec;
  for (const std::string& candidate : candidates) {
    if (fs::is_regular_file(candidate, ec)) {
      return candidate;
    }
  }

  std::string message = "no nuclear-data store found. Looked in:\n";
  if (candidates.empty()) {
    message += "  (nowhere -- no --store, no " + std::string(kEnvVar) +
               ", and no store beside the executable or under ./data)\n";
  } else {
    for (const std::string& candidate : candidates) {
      message += "  " + candidate + "\n";
    }
  }
  message += "Point --store at one, set " + std::string(kEnvVar) +
             ", or build one from ENDF tapes:\n"
             "  nusift_stage_data --decay-dir <endf-decay-dir> -o data/nusift.h5";
  throw InputError(tagged(kModule, message));
}

}  // namespace nusift
