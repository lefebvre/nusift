#include "nusift/exposure/air_coefficients.hpp"

#include <cmath>
#include <cstddef>

namespace nusift::exposure {
namespace {

// Dry air, from the NIST X-Ray Mass Attenuation Tables (Hubbell & Seltzer, NISTIR 5632).
// NIST tabulates in cm^2/g; these are SI, m^2/kg, which is cm^2/g divided by 10.
//
// The grid is NIST's own, not a resampling. Keeping their energy points means the
// interpolation error is only what log-log interpolation introduces between tabulated values,
// with nothing added by a prior regridding step.
struct AirCoefficient {
  double energyEv;
  double massAttenuation;       // mu/rho    [m^2/kg]
  double massEnergyAbsorption;  // mu_en/rho [m^2/kg]
};

constexpr AirCoefficient kAir[] = {
    {1.0e4, 0.5120, 0.4742},     {1.5e4, 0.1614, 0.1334},     {2.0e4, 0.07779, 0.05389},
    {3.0e4, 0.03538, 0.01537},   {4.0e4, 0.02485, 0.006833},  {5.0e4, 0.02080, 0.004098},
    {6.0e4, 0.01875, 0.003041},  {8.0e4, 0.01662, 0.002407},  {1.0e5, 0.01541, 0.002325},
    {1.5e5, 0.01356, 0.002496},  {2.0e5, 0.01233, 0.002672},  {3.0e5, 0.01067, 0.002872},
    {4.0e5, 0.009549, 0.002949}, {5.0e5, 0.008712, 0.002966}, {6.0e5, 0.008055, 0.002953},
    {8.0e5, 0.007074, 0.002882}, {1.0e6, 0.006358, 0.002789}, {1.25e6, 0.005687, 0.002666},
    {1.5e6, 0.005175, 0.002547}, {2.0e6, 0.004447, 0.002345}, {3.0e6, 0.003581, 0.002057},
    {4.0e6, 0.003079, 0.001870}, {5.0e6, 0.002751, 0.001740}, {6.0e6, 0.002522, 0.001647},
    {8.0e6, 0.002225, 0.001525}, {1.0e7, 0.002045, 0.001450},
};
constexpr std::size_t kCount = sizeof(kAir) / sizeof(kAir[0]);

// Log-log interpolation, clamped at both ends.
//
// Log-log rather than linear because both coefficients are close to power laws in energy over
// each interval; interpolating them linearly on a grid this coarse would misplace the values
// by percent-level amounts in exactly the few-hundred-keV region where most decay photons sit.
double interpolate(double energyEv, bool absorption) {
  const auto value = [absorption](std::size_t i) {
    return absorption ? kAir[i].massEnergyAbsorption : kAir[i].massAttenuation;
  };
  if (energyEv <= kAir[0].energyEv) {
    return value(0);
  }
  if (energyEv >= kAir[kCount - 1].energyEv) {
    return value(kCount - 1);
  }
  std::size_t hi = 1;
  while (hi < kCount && kAir[hi].energyEv < energyEv) {
    ++hi;
  }
  const std::size_t lo = hi - 1;
  const double t = (std::log(energyEv) - std::log(kAir[lo].energyEv)) /
                   (std::log(kAir[hi].energyEv) - std::log(kAir[lo].energyEv));
  return std::exp(std::log(value(lo)) + t * (std::log(value(hi)) - std::log(value(lo))));
}

}  // namespace

double airMassEnergyAbsorption(double energyEv) {
  return interpolate(energyEv, /*absorption=*/true);
}

double airMassAttenuation(double energyEv) {
  return interpolate(energyEv, /*absorption=*/false);
}

bool isOutsideTabulatedRange(double energyEv) {
  return energyEv < kMinTabulatedEv || energyEv > kMaxTabulatedEv;
}

}  // namespace nusift::exposure
