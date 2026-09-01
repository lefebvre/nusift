#pragma once
/**
 * @file
 * @brief Fluence-to-effective-dose conversion coefficients for photons.
 * @ingroup exposure
 */
//
// ICRP Publication 116, Table A.1: effective dose per unit fluence, for monoenergetic photons
// incident on the reference computational phantoms in six irradiation geometries.
//
// This is what makes a sievert in NuSIFT a sievert. Air kerma answers "how much energy would
// air absorb here"; effective dose answers "what does this field do to a person", and the two
// are not proportional. They agree within a couple of percent near 1 MeV, which is a
// coincidence of energy rather than a relationship, and diverge by a factor of five by 60 keV,
// where air keeps absorbing strongly and a body's organs are shielded by everything in front
// of them. A tool that reported the first while writing Sv would overstate the hazard from
// every soft emitter it ranked, in the direction that looks conservative and is simply wrong.
//
// THE GEOMETRY IS PART OF THE ANSWER. Effective dose is defined for a person in a field, so it
// depends on how the person stands in it: facing the source (AP) is the most exposing for the
// organs that carry the largest tissue weights, and by 30 keV the same fluence gives four times
// the effective dose in AP that it does in PA. Which geometry a number was computed in is
// therefore not a setting, it is half of what the number means, and every report that prints a
// sievert prints the geometry beside it.
//
// WHAT THIS IS NOT. Effective dose is a protection quantity: it is defined on reference
// phantoms with sex-averaged, tissue-weighted organ doses, for setting and checking limits.
// It is not a measurable quantity, it is not the dose to any particular person, and it is not
// what a survey meter reads -- that is ambient dose equivalent H*(10), which comes from ICRP 74
// and is not in this publication or in NuSIFT.
//
#include <string_view>

namespace nusift::exposure {

// How the body is oriented in the field. ICRP 116's own names, because a report that renamed
// them would be harder to check against the publication it took its numbers from.
enum class Irradiation {
  AP,    // antero-posterior: front to back, the conservative default
  PA,    // postero-anterior: back to front
  LLAT,  // left lateral
  RLAT,  // right lateral
  ROT,   // rotational: the body turning in a horizontal beam
  ISO,   // isotropic: the field arriving equally from every direction
};

const char* irradiationName(Irradiation irradiation);

// Parse a spelling case-insensitively against the names irradiationName() prints.
bool parseIrradiation(std::string_view text, Irradiation& out);

// Lowest and highest tabulated energies. The table runs from 10 keV to 10 GeV, so the high end
// is far above anything a decay photon reaches and only the low end is ever met in practice.
// Below 10 keV the coefficient is CLAMPED rather than extrapolated, for the reason the air
// coefficients are: the curve is turning over steeply there and any extrapolation of it would
// be invention. Am-241 is the nuclide this bites -- a real part of its emission sits below the
// table -- so the clamp is reported rather than left silent.
inline constexpr double kMinTabulatedDoseEv = 1.0e4;   // 10 keV
inline constexpr double kMaxTabulatedDoseEv = 1.0e13;  // 10 GeV

// Effective dose per unit fluence at `energyEv`, in Sv*m^2, log-log interpolated.
//
// ICRP prints these in pSv*cm^2 and this returns SI, which is the printed value times 1e-16.
// The conversion happens here rather than in the table so that every number in the table is
// literally the number on the page.
double effectiveDosePerFluence(double energyEv, Irradiation irradiation);

// True when `energyEv` falls outside the tabulated range, so the value above is a clamped end
// point rather than an interpolation.
bool isOutsideTabulatedDoseRange(double energyEv);

}  // namespace nusift::exposure
