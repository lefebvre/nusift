#pragma once
/**
 * @file
 * @brief A discrete photon line, and a non-owning view of one nuclide's line spectrum.
 * @ingroup nucdata
 */
//
// NuSIFT persists the full per-nuclide discrete photon spectrum from ENDF MF8/MT457 rather
// than collapsing it to a single exposure-rate constant at staging time. The reason is
// physical, not just informational: the air attenuation factor exp(-mu_air(E) d) is
// energy-dependent, so it cannot be pulled out of the sum over lines. There is no single
// constant that is correct at more than one distance, which means the lines have to survive
// to evaluation time. Persisting them additionally buys ranking of individual gamma LINES,
// runtime geometry changes with no restage, and a path to shielding.
//
// Lines are stored CSR-packed across all nuclides in the data store, so a LineSpectrum is a
// span into that flat array -- no per-nuclide allocation, and no copy when a nuclide's
// spectrum is handed to the exposure evaluator.
//
#include <cstddef>
#include <span>

namespace nusift {

// ENDF MT457 spectrum types NuSIFT models. Both are photons and both deposit energy in air,
// so both are kept, but they are distinguished in the store because the 511 keV annihilation
// line is reported under STYP 9 while ordinary gammas are STYP 0 -- and a nuclide can emit
// at 511 keV under both, which is a double-counting hazard. The staging tool checks for it and
// warns; it does not drop either entry, since which one is the duplicate is a judgement about
// the evaluation rather than something the tool can know.
enum class SpectrumType : int {
  Gamma = 0,        // ENDF STYP 0
  XrayOrAnnih = 9,  // ENDF STYP 9: X-rays and annihilation radiation
};

// One discrete photon line: energy, and absolute intensity in photons per decay.
//
// `intensity` is ABSOLUTE (FD * RI in ENDF terms -- the discrete normalisation factor times
// the relative intensity), not relative to the strongest line. Storing it absolute is what
// lets the exposure sum be written as a plain dot product against the line array.
struct GammaLine {
  double energyEv = 0.0;
  double intensity = 0.0;  // photons per decay
  SpectrumType type = SpectrumType::Gamma;
};

// A non-owning view of one nuclide's discrete photon lines. Valid only while the
// NuclearData it came from is alive.
using LineSpectrum = std::span<const GammaLine>;

// Total photon energy carried by the discrete lines, in eV per decay. Compared against the
// staged average electromagnetic decay energy to detect how much of a nuclide's photon
// output lives in a continuum NuSIFT does not model -- see NuclearData::continuumPhotonEv.
inline double discretePhotonEnergyEv(LineSpectrum lines) {
  double total = 0.0;
  for (const GammaLine& line : lines) {
    total += line.energyEv * line.intensity;
  }
  return total;
}

}  // namespace nusift
