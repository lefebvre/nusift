#pragma once
/**
 * @file
 * @brief Converting a fission source between fission count, energy release, and kilotons.
 * @ingroup seed
 */
//
// A fission source can be specified three ways, and they all reduce to a count of fissions:
//
//   fissions      the physical quantity everything downstream uses
//   kilotons TNT  how an explosive yield is quoted
//   joules        how a reactor energy release is quoted
//
// THE CONSTANT MATTERS, AND THE RIGHT ONE DEPENDS ON THE QUESTION. Fission releases about
// 200 MeV in total, but only about 180 MeV of that appears promptly. The remainder arrives as
// delayed beta and gamma energy over the following hours and days -- which is precisely the
// decay NuSIFT is being asked to model, so counting it as part of the driving yield would
// double-count it.
//
//   180 MeV  explosive yield. Reproduces the canonical 1.45e23 fissions per kiloton.
//   200 MeV  total recoverable energy. The reactor convention, giving 1.31e23 per kiloton.
//
// The two differ by 11% in every downstream number. Neither is wrong; picking silently would
// be, which is why the choice is a named parameter and why the resolved fission count is
// reported in the provenance line of every run.
//
#include <string_view>

namespace nusift::seed {

inline constexpr double kJoulesPerKt = 4.184e12;  // exact: 1 kt TNT is defined as 1e12 calories
inline constexpr double kJoulesPerMeV = 1.602176634e-13;

inline constexpr double kMeVPerFissionExplosiveYield = 180.0;
inline constexpr double kMeVPerFissionRecoverable = 200.0;

// Fissions releasing `kt` kilotons at the given energy per fission.
double fissionsFromKt(double kt, double meVPerFission = kMeVPerFissionExplosiveYield);

// Fissions releasing `joules`.
double fissionsFromEnergyJ(double joules, double meVPerFission = kMeVPerFissionExplosiveYield);

// The inverse, so a run specified in fissions can still report its kiloton equivalent.
double ktFromFissions(double fissions, double meVPerFission = kMeVPerFissionExplosiveYield);

// Parse "explosive", "recoverable", or a bare number in MeV. Returns false if unrecognised.
bool parseMeVPerFission(std::string_view text, double& meVPerFission);

}  // namespace nusift::seed
