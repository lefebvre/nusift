#include "nusift/nucdata/fission_yield.hpp"

#include <algorithm>
#include <cctype>
#include <charconv>
#include <cmath>
#include <string>

namespace nusift {
namespace {

std::string lower(std::string_view text) {
  std::string out(text);
  std::transform(out.begin(), out.end(), out.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return out;
}

}  // namespace

double FissionYieldSet::totalYield() const {
  double total = 0.0;
  for (const FissionProduct& product : products) {
    total += product.yield;
  }
  return total;
}

void FissionYieldTable::add(FissionYieldSet set) {
  sets_.push_back(std::move(set));
}

std::vector<Zai> FissionYieldTable::parents() const {
  std::vector<Zai> found;
  found.reserve(sets_.size());
  for (const FissionYieldSet& set : sets_) {
    found.push_back(set.parent);
  }
  std::sort(found.begin(), found.end());
  found.erase(std::unique(found.begin(), found.end()), found.end());
  return found;
}

const FissionYieldSet* FissionYieldTable::nearest(const Zai& parent, double energyEv) const {
  const FissionYieldSet* best = nullptr;
  double bestDistance = 0.0;
  for (const FissionYieldSet& set : sets_) {
    if (!(set.parent == parent)) {
      continue;
    }
    const double distance = std::abs(set.energyEv - energyEv);
    if (best == nullptr || distance < bestDistance) {
      best = &set;
      bestDistance = distance;
    }
  }
  return best;
}

bool parseIncidentEnergy(std::string_view text, double& energyEv) {
  const std::string key = lower(text);
  if (key == "spontaneous" || key == "sf") {
    energyEv = kSpontaneousEv;
    return true;
  }
  if (key == "thermal") {
    energyEv = kThermalEv;
    return true;
  }
  if (key == "fast") {
    energyEv = kFastEv;
    return true;
  }
  if (key == "14mev" || key == "fusion") {
    energyEv = kFusionEv;
    return true;
  }
  double value = 0.0;
  const auto result = std::from_chars(key.data(), key.data() + key.size(), value);
  if (result.ec == std::errc{} && result.ptr == key.data() + key.size() && value >= 0.0) {
    energyEv = value;
    return true;
  }
  return false;
}

const char* describeIncidentEnergy(double energyEv) {
  // Tolerances are wide because the point is to name what the evaluation tabulated, and
  // evaluations differ slightly in where they put "fast".
  if (energyEv <= 0.0) {
    return "spontaneous";
  }
  if (energyEv < 1.0) {
    return "thermal";
  }
  if (energyEv < 5.0e6) {
    return "fast";
  }
  return "14 MeV";
}

}  // namespace nusift
