#pragma once
/**
 * @file
 * @brief Photon exposure and fluence from an unshielded point source at a distance.
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
// WHAT IS NOT MODELED, and would each raise the answer:
//   * scattered photons (buildup) -- the default factor is 1.0, i.e. uncollided fluence only.
//     For a bare source at a meter in air this is a small correction; through any shielding it
//     is not, which is why the field is exposed rather than hidden.
//   * source self-absorption -- a point source has no volume to absorb its own photons.
//   * bremsstrahlung and any continuous photon spectrum. NuSIFT models discrete lines only;
//     what is missing is recorded per nuclide and reported alongside the ranking.
//   * beta, alpha, and neutron dose entirely.
//
// UNITS. Two quantities come out of this file, and which one a caller gets is decided by the
// unit it asks for rather than displayed by it. Roentgen and gray are photon exposure and air
// kerma: what the field would deposit in air. A sievert is ICRP 116 effective dose to a person,
// computed down its own path -- fluence at the point times the fluence-to-dose coefficient for
// the photon's energy and the stated irradiation geometry -- and it is not a multiple of an air
// kerma. It is still not H*(10), the operational quantity a survey meter reads, which needs
// ICRP 74 and is not here.
//
#include <span>

#include "nusift/exposure/dose_coefficients.hpp"
#include "nusift/nucdata/photon_lines.hpp"
#include "nusift/units.hpp"

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

  // How the body is oriented in the field, for the effective-dose units and them alone -- an
  // exposure in roentgen and an air kerma in gray have no phantom in them to orient.
  //
  // It sits here, beside the distance, rather than in the response spec, because it is the same
  // KIND of fact: not an input that can be forgotten once the number exists, but half of what
  // the number means. A sievert with no irradiation geometry named is as incomplete as an
  // exposure with no distance, and keeping the two together is what makes it hard to report one
  // without the other -- every path that carries a geometry into a table, an adjoint, a
  // counterfactual or a report carries this with it and cannot silently disagree about it.
  //
  // AP by default: facing the source is the most exposing orientation for the organs carrying
  // the largest tissue weights, so it is the conservative choice and the one regulatory
  // screening reaches for.
  Irradiation irradiation = Irradiation::AP;
};

// Exposure-rate coefficient for a single photon energy, in R/h per photon/s emitted.
double pointExposureCoeff(double energyEv, const PointSourceGeometry& geometry);

// Fluence-rate coefficient for a single photon energy: photons arriving at the point, per
// (m^2 * s), per photon/s emitted at that energy. This is pointExposureCoeff with the energy
// deposition factor removed -- inverse-square spreading and path attenuation only, no kerma
// factor -- which makes it the native quantity for a detector or a transport-code source term:
// it says how many photons get there, and leaves what they do once they arrive to the consumer.
//
// Like the exposure coefficient it carries the uncollided assumption: scattered photons are
// counted only through the explicit buildup factor, which is why the same geometry (and the
// same honest report of it) rides along.
double pointFluenceCoeff(double energyEv, const PointSourceGeometry& geometry);

// Effective-dose-rate coefficient for a single photon energy, in Sv/h per photon/s emitted.
//
// The fluence coefficient with ICRP 116's fluence-to-effective-dose conversion applied, which
// is the honest route to a sievert: how many photons arrive, times what a photon of that energy
// does to a person standing that way. Note which factor is NOT here -- air's mass energy
// absorption. Effective dose does not care what air would have absorbed; the phantom
// calculation already carries what a body absorbs, and multiplying by both would be counting
// the interaction twice.
double pointEffectiveDoseCoeff(double energyEv, const PointSourceGeometry& geometry);

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

// Fluence rate at the point per becquerel, in photons/(m^2 * s * Bq). The photon-metric
// counterpart of the above, for the units that carry a point geometry.
double fluenceRatePerBecquerel(LineSpectrum lines, const PointSourceGeometry& geometry);

// Effective dose rate per becquerel, in Sv/h per Bq, summed over the lines with the conversion
// inside the sum. Inside, because the conversion varies by three orders of magnitude across a
// spectrum and applying an average of it to a total would be a different quantity -- the same
// reason attenuation lives inside the exposure sum.
double effectiveDoseRatePerBecquerel(LineSpectrum lines, const PointSourceGeometry& geometry);

// The optical depth of the air path in mean free paths, mu_air(E) * d, averaged over the lines
// with each weighted by the exposure it actually delivers at the point of interest.
//
// Buildup is a function of exactly this quantity, so it is the number that says whether
// leaving the factor at 1.0 is a small omission or a large one. At a meter it is a few
// hundredths and the scattered photons are a percent-level correction; past about half a mean
// free path they are tens of percent of the uncollided value, and beyond one they exceed it.
// Zero with attenuation off, where there is no path to be thick, and zero for a spectrum that
// delivers nothing.
double meanOpticalDepth(LineSpectrum lines, const PointSourceGeometry& geometry);

// Roentgen to absorbed dose in air. The constant is units::kGyPerR -- one spelling of the
// roentgen, in the one file that owns every conversion.
//
// There is deliberately no roentgenToSievert beside it any more. It existed to produce the
// sievert column by multiplying air kerma by a photon radiation weighting factor of 1, which
// names a quantity nobody wants: not air kerma, which Gy already says, and not dose to a
// person, which is what a reader takes a sievert to mean. A sievert now comes from
// effectiveDoseRatePerBecquerel() and from nowhere else.
constexpr double roentgenToGray(double roentgen) {
  return roentgen * units::kGyPerR;
}

}  // namespace nusift::exposure
