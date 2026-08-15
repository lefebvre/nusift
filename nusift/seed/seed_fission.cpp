#include "nusift/seed/seed_fission.hpp"

#include <algorithm>
#include <cctype>
#include <charconv>
#include <cmath>
#include <cstdio>
#include <string>

#include "nusift/core/error.hpp"
#include "nusift/core/nuclide_name.hpp"
#include "nusift/nucdata/nuclear_data.hpp"

namespace nusift::seed {
namespace {

constexpr const char* kModule = "fission seed";

}  // namespace

double fissionsFromKt(double kt, double meVPerFission) {
  return kt * kJoulesPerKt / (meVPerFission * kJoulesPerMeV);
}

double fissionsFromEnergyJ(double joules, double meVPerFission) {
  return joules / (meVPerFission * kJoulesPerMeV);
}

double ktFromFissions(double fissions, double meVPerFission) {
  return fissions * meVPerFission * kJoulesPerMeV / kJoulesPerKt;
}

bool parseMeVPerFission(std::string_view text, double& meVPerFission) {
  std::string key(text);
  std::transform(key.begin(), key.end(), key.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  if (key == "explosive") {
    meVPerFission = kMeVPerFissionExplosiveYield;
    return true;
  }
  if (key == "recoverable") {
    meVPerFission = kMeVPerFissionRecoverable;
    return true;
  }
  double value = 0.0;
  const auto result = std::from_chars(key.data(), key.data() + key.size(), value);
  if (result.ec == std::errc{} && result.ptr == key.data() + key.size() && value > 0.0) {
    meVPerFission = value;
    return true;
  }
  return false;
}

Inventory seedFromFission(const NuclearData& data, const FissionSeed& seed) {
  if (!(seed.fissions > 0.0)) {
    throw InputError(tagged(kModule, "the number of fissions must be positive"));
  }

  const FissionYieldTable& table = data.fissionYields();
  const FissionYieldSet* set = table.nearest(seed.fissile, seed.incidentEnergyEv);
  if (set == nullptr) {
    // The set of fissionable nuclides in an evaluation is small, so listing it beats leaving
    // the user to guess which spelling or which nuclide the store actually knows.
    std::string message = "this data store has no fission yields for " +
                          formatNuclideName(seed.fissile) + ". It carries yields for: ";
    const std::vector<Zai> parents = table.parents();
    if (parents.empty()) {
      message += "(none at all -- stage with --nfy-dir and --sfy-dir)";
    } else {
      for (std::size_t i = 0; i < parents.size(); ++i) {
        message += (i == 0 ? "" : ", ") + formatNuclideName(parents[i]);
      }
    }
    throw InputError(tagged(kModule, message));
  }

  Inventory inventory;
  int missing = 0;
  for (const FissionProduct& product : set->products) {
    if (product.yield <= 0.0) {
      continue;
    }
    // A product the chain does not know would have nowhere to decay to. In practice closure
    // registers every yield product, so this is a belt-and-braces count rather than an
    // expected condition -- but a silent drop here would quietly lose inventory.
    if (data.indexOf(product.nuclide) < 0) {
      ++missing;
      continue;
    }
    inventory.add(product.nuclide, seed.fissions * product.yield);
  }

  if (inventory.empty()) {
    throw InputError(tagged(kModule, "the yield set for " + formatNuclideName(seed.fissile) +
                                         " produced no usable products"));
  }

  char buffer[320];
  const double totalYield = set->totalYield();
  std::snprintf(buffer, sizeof(buffer),
                "fission: %s at %s (%.4g eV), %.4g fissions, %zu products, sum Y_indep = %.4g",
                formatNuclideName(seed.fissile).c_str(), describeIncidentEnergy(set->energyEv),
                set->energyEv, seed.fissions, inventory.entries().size(), totalYield);
  std::string provenance = buffer;

  // Independent yields sum to about 2.0 because fission makes two fragments. A materially
  // different sum means the wrong ENDF section was staged -- cumulative yields sum to far more
  // -- and every number downstream is wrong by that factor. Cheap to check, and invisible
  // otherwise until someone compares against a reference.
  if (totalYield < 1.8 || totalYield > 2.2) {
    provenance += "  [WARNING: independent yields should sum to about 2.0]";
  }
  if (missing > 0) {
    std::snprintf(buffer, sizeof(buffer), "  [%d products absent from the chain, skipped]",
                  missing);
    provenance += buffer;
  }
  inventory.setProvenance(provenance);
  return inventory;
}

}  // namespace nusift::seed
