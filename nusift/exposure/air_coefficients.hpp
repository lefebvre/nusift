#pragma once
/**
 * @file
 * @brief Photon interaction coefficients for dry air.
 * @ingroup exposure
 */
//
// Two energy-dependent coefficients, both from the NIST X-Ray Mass Attenuation Tables
// (Hubbell & Seltzer, NISTIR 5632), tabulated for dry air near sea level:
//
//   mu/rho      mass attenuation -- how much a beam is removed traversing air
//   mu_en/rho   mass energy-absorption -- how much energy air actually absorbs
//
// Exposure needs both, and they are not interchangeable. Attenuation governs what reaches the
// point of interest; energy absorption governs what that fluence deposits once it arrives.
//
// The energy dependence is the reason NuSIFT persists photon lines rather than a single
// per-nuclide constant. Both coefficients vary by more than two orders of magnitude across the
// tabulated range, so a spectrum cannot be collapsed to one effective energy without fixing
// the geometry at the same time.
//
namespace nusift::exposure {

// Lowest and highest tabulated energies. Outside this range the coefficients are CLAMPED to
// the end values rather than extrapolated: below 10 keV photoelectric absorption rises so
// steeply that a power-law extrapolation is badly wrong, and such photons are absorbed by any
// realistic source encapsulation long before reaching air. Callers that care can compare
// against these bounds.
inline constexpr double kMinTabulatedEv = 1.0e4;  // 10 keV
inline constexpr double kMaxTabulatedEv = 1.0e7;  // 10 MeV

// Mass energy-absorption coefficient (mu_en/rho) for dry air [m^2/kg], log-log interpolated.
double airMassEnergyAbsorption(double energyEv);

// Mass attenuation coefficient (mu/rho) for dry air [m^2/kg], log-log interpolated.
double airMassAttenuation(double energyEv);

// True when `energyEv` falls outside the tabulated range, so the value returned above is a
// clamped end point rather than an interpolation. Reported by `nusift data info` as a store
// coverage gap and marked against the individual line by `nusift data nuclide`, rather than
// treated as an error -- a soft X-ray line is real data, it is just outside where this model
// can say anything useful.
bool isOutsideTabulatedRange(double energyEv);

}  // namespace nusift::exposure
