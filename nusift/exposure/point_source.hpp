#pragma once
/**
 * @file
 * @brief Photon exposure from an unshielded point source at a distance.
 * @ingroup exposure
 */
//
// The model. An unshielded point source in air. Photons spread over 4*pi*d^2, are attenuated
// along the path by air, and deposit energy in air at the point of interest according to the
// mass energy-absorption coefficient:
//
//   mu_air(E) = (mu/rho)_air(E) * rho_air                                            [1/m]
//   k(E)      = exp(-mu_air(E) * d) / (4 pi d^2)      <- geometry and path attenuation
//               * E * 1.602176634e-19                 <- eV to J, giving energy fluence
//               * (mu_en/rho)_air(E)                  <- air kerma, assuming CPE
//               * 3600 / 0.00876                      <- kerma to exposure, per hour
//               * buildup
//                                        [R/h per photon/s emitted at energy E]
//
//   Xdot = A[Bq] * SUM_j y_j k(E_j)                                                  [R/h]
//
// WHY THE LINES CANNOT BE COLLAPSED. mu_air is energy-dependent, so exp(-mu_air(E) d) sits
// inside the sum over lines and cannot be factored out of it. There is therefore no single
// per-nuclide constant that is correct at more than one distance -- which is precisely why the
// data store persists every discrete line instead of one number per nuclide. gammaConstant()
// below exists only for the vacuum case, where the exponential is 1 and the sum does factor.
//
// WHAT IS NOT MODELLED, and would each raise the answer:
//   * scattered photons (buildup) -- the default factor is 1.0, i.e. uncollided fluence only.
//     For a bare source at a metre in air this is a small correction; through any shielding it
//     is not, which is why the field is exposed rather than hidden.
//   * source self-absorption -- a point source has no volume to absorb its own photons.
//   * bremsstrahlung and any continuous photon spectrum. NuSIFT models discrete lines only;
//     what is missing is recorded per nuclide and reported alongside the ranking.
//   * beta, alpha, and neutron dose entirely.
//
// UNITS. The physical quantity computed is photon exposure in air, in roentgen. Sieverts are
// obtained by the air-kerma conversion with a radiation weighting factor of 1. That is NOT an
// ICRP-74 fluence-to-H*(10) operational quantity, and it is not an effective dose to a person.
//
#include <span>

#include "nusift/nucdata/photon_lines.hpp"

namespace nusift::exposure {

struct PointSourceGeometry {
  double distanceM = 1.0;

  // Scales the mass attenuation coefficient to a linear one. The default is dry air at about
  // 20 C and one atmosphere; a site at elevation is meaningfully thinner.
  double airDensityKgM3 = 1.205;

  // When false the path attenuation term is dropped and the result is pure inverse-square.
  // Useful for comparing against published gamma constants, which are vacuum quantities.
  bool airAttenuation = true;

  // Multiplicative scatter buildup. 1.0 means uncollided photons only, which understates a
  // real measurement. Exposed rather than assumed so the assumption is visible in the caller.
  double buildup = 1.0;
};

// Exposure-rate coefficient for a single photon energy, in R/h per photon/s emitted.
double pointExposureCoeff(double energyEv, const PointSourceGeometry& geometry);

// The specific gamma-ray constant: exposure rate per unit activity at unit distance, in
// vacuum, [R*m^2/(h*Bq)]. Distance-independent by construction, since with no attenuation the
// only distance dependence is the 1/d^2 that has been divided out.
//
// This is the quantity published tables give (usually as R*cm^2/(h*mCi)), so it is what to
// compare against a reference. It is NOT what NuSIFT uses to compute an exposure rate --
// exposureRate() sums over lines with attenuation inside the sum.
double gammaConstant(LineSpectrum lines);

// Exposure rate from `activityBq` of a nuclide with this line spectrum, in R/h.
double exposureRate(LineSpectrum lines, double activityBq, const PointSourceGeometry& geometry);

// Exposure rate per becquerel, in R/h per Bq. The per-nuclide weight the response layer
// multiplies activity by; separated out because it depends only on the spectrum and the
// geometry, so it is computed once per nuclide rather than once per nuclide per time.
double exposureRatePerBecquerel(LineSpectrum lines, const PointSourceGeometry& geometry);

// Roentgen to absorbed dose in air, and to equivalent dose with a photon radiation weighting
// factor of 1. Both are the same number; they are named separately because the quantities are
// different and conflating them in a report is how a dose gets misread.
constexpr double kGrayPerRoentgen = 0.00876;
constexpr double roentgenToGray(double roentgen) {
  return roentgen * kGrayPerRoentgen;
}
constexpr double roentgenToSievert(double roentgen) {
  return roentgen * kGrayPerRoentgen;
}

}  // namespace nusift::exposure
