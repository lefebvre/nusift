#pragma once
/**
 * @file
 * @brief Physical constants and unit conversions shared across NuSIFT.
 * @ingroup core
 */
//
// One place for every constant that appears in more than one module, so a value can never
// be spelled two ways in two files. Values are CODATA 2018 / SI where a definition exists.
//
// Naming convention throughout NuSIFT: a quantity's units are part of its name
// (energyEv, halfLifeSeconds, distanceM, sigmaBarn). This file is where the conversions
// between those spellings live.
//
#include <numbers>

namespace nusift::units {

// --- fundamental ------------------------------------------------------------
inline constexpr double kAvogadro = 6.02214076e23;        // /mol (exact, SI 2019)
inline constexpr double kEvToJ = 1.602176634e-19;         // J/eV (exact, SI 2019)
inline constexpr double kNeutronMassAmu = 1.00866491595;  // u -- converts ENDF AWR to g/mol
inline constexpr double kLn2 = std::numbers::ln2;

// --- time -------------------------------------------------------------------
// A year is the Julian year, 365.25 d. Stated here and in `--help` because the choice is
// arbitrary but its consequences are not: over a 100 y decay the difference against a
// 365 d year is about three months of Cs-137 decay.
inline constexpr double kSecondsPerMinute = 60.0;
inline constexpr double kSecondsPerHour = 3600.0;
inline constexpr double kSecondsPerDay = 86400.0;
inline constexpr double kSecondsPerYear = 365.25 * kSecondsPerDay;

// --- activity ---------------------------------------------------------------
inline constexpr double kBqPerCi = 3.7e10;  // exact by definition of the curie

// --- exposure and dose ------------------------------------------------------
// Air kerma per unit exposure. The roentgen is defined as 2.58e-4 C/kg, which with a mean
// ionization energy of 33.97 J/C gives 8.76e-3 Gy per R. NuSIFT reports one physical
// quantity -- photon exposure -- and converts to Sv with a radiation weighting factor of 1;
// this is NOT an ICRP-74 fluence-to-H*(10) conversion, and the headers that use it say so.
inline constexpr double kGyPerR = 0.00876;
inline constexpr double kRadiationWeightingPhoton = 1.0;

// Dry air at ~20 C, 101.325 kPa. Scales the mass attenuation coefficient to a linear one;
// override it for a site at elevation.
inline constexpr double kAirDensityKgM3 = 1.205;

// --- cross sections ---------------------------------------------------------
inline constexpr double kBarnToCm2 = 1.0e-24;

// --- derived helpers --------------------------------------------------------
// Decay constant from half-life. Returns 0 for a non-positive half-life, which is how
// NuSIFT (and the data store) encode "stable" -- a stable nuclide is a terminator with no
// removal rate, not an error.
constexpr double decayConstant(double halfLifeSeconds) {
  return halfLifeSeconds > 0.0 ? kLn2 / halfLifeSeconds : 0.0;
}

// Molar mass from the ENDF atomic weight ratio (mass relative to the neutron). Staged per
// nuclide so gram/atom conversions do not fall back on using A as the molar mass, which is
// wrong by ~0.1% for mid-A nuclides and considerably worse for light ones.
constexpr double molarMassFromAwr(double awr) {
  return awr * kNeutronMassAmu;
}

}  // namespace nusift::units
