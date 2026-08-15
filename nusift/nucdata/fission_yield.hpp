#pragma once
/**
 * @file
 * @brief Independent fission yields, and lookup by fissioning nuclide and incident energy.
 * @ingroup nucdata
 */
//
// INDEPENDENT yields (ENDF MF8/MT454), never cumulative (MT459). The distinction is not a
// detail: independent yields are what fission produces directly, before any decay. NuSIFT
// decays the chain explicitly, so seeding with cumulative yields -- which already include
// everything a precursor decays into -- would count every precursor decay twice. The staging
// tool reads MT454 for this reason and the reader has no option to do otherwise.
//
// A set is tabulated at a specific incident neutron energy. ENDF evaluations typically give
// thermal (0.0253 eV), fast (around 0.5 MeV), and 14 MeV; spontaneous fission is tabulated at
// energy 0. Lookup snaps to the nearest tabulated energy rather than interpolating, because
// yields between tabulated points are not a smooth function anyone has evaluated -- and
// silently interpolating would invent data the evaluation does not support.
//
#include <cstdint>
#include <span>
#include <string_view>
#include <vector>

#include "nusift/core/nuclide.hpp"

namespace nusift {

struct FissionProduct {
  Zai nuclide;
  double yield = 0.0;  // atoms produced per fission
};

struct FissionYieldSet {
  Zai parent;
  double energyEv = 0.0;  // incident neutron energy; 0 means spontaneous fission
  std::vector<FissionProduct> products;

  // Sum of the independent yields. For a binary fission this is close to 2.0, because each
  // fission makes two fragments. A set summing to about 1.0 or about 4.0 is the signature of
  // the wrong MT having been staged, so it is worth reporting rather than assuming.
  double totalYield() const;
};

class FissionYieldTable {
public:
  void add(FissionYieldSet set);

  bool empty() const { return sets_.empty(); }
  std::span<const FissionYieldSet> sets() const { return sets_; }

  // Every nuclide with tabulated yields, ascending by ZAI key. What `nusift data info` lists
  // so a user can find out what they are allowed to seed from.
  std::vector<Zai> parents() const;

  // The set for `parent` whose tabulated energy is nearest `energyEv`, or null if that nuclide
  // has no yields at all. Returns the nearest rather than requiring an exact match, since a
  // user asking for "fast" means whatever the evaluation tabulated near there.
  const FissionYieldSet* nearest(const Zai& parent, double energyEv) const;

private:
  std::vector<FissionYieldSet> sets_;
};

// Named incident energies, matching what evaluations actually tabulate.
inline constexpr double kSpontaneousEv = 0.0;
inline constexpr double kThermalEv = 0.0253;
inline constexpr double kFastEv = 5.0e5;
inline constexpr double kFusionEv = 1.4e7;

// Parse "thermal", "fast", "14mev", "spontaneous", or a bare number in eV. Returns false if
// unrecognised.
bool parseIncidentEnergy(std::string_view text, double& energyEv);

// The nearest named energy, for reporting back what a lookup actually landed on.
const char* describeIncidentEnergy(double energyEv);

}  // namespace nusift
