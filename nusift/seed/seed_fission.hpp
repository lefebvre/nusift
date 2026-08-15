#pragma once
/**
 * @file
 * @brief Building an inventory from fission.
 * @ingroup seed
 */
//
// "Born from fission" is half of what NuSIFT exists to triage, and this is where that
// inventory comes from: N fissions of a nuclide at an incident energy, distributed over the
// fission products by their independent yields.
//
// The result is the inventory at the instant fission stops. Everything after that -- decay,
// ranking, exposure -- is the ordinary pipeline, which is why this returns a plain Inventory
// rather than anything fission-specific.
//
#include <string>

#include "nusift/core/nuclide.hpp"
#include "nusift/engine/inventory.hpp"
#include "nusift/nucdata/fission_yield.hpp"
#include "nusift/seed/fission_energy.hpp"

namespace nusift {

class NuclearData;

namespace seed {

struct FissionSeed {
  Zai fissile;
  double incidentEnergyEv = kThermalEv;
  double fissions = 0.0;

  // Recorded only so the provenance line can say how the fission count was arrived at when the
  // user gave kilotons or joules rather than a count.
  double meVPerFission = kMeVPerFissionExplosiveYield;
};

// Distribute `seed.fissions` fissions over the products of the yield set nearest the requested
// incident energy.
//
// Throws InputError when the store carries no yields for that nuclide, naming what it does
// carry -- the set of fissionable nuclides in an evaluation is small and worth listing rather
// than leaving the user to guess.
//
// The provenance records the resolved fission count, the energy the lookup actually landed on,
// and the sum of the independent yields. That last one is a cheap integrity check: independent
// yields sum to about 2.0 because fission makes two fragments, so a sum near 1.0 or 4.0 means
// the wrong ENDF section was staged and every number downstream is wrong by that factor.
Inventory seedFromFission(const NuclearData& data, const FissionSeed& seed);

}  // namespace seed
}  // namespace nusift
