#include "nusift/exposure/point_source.hpp"

#include <cmath>
#include <numbers>

#include "nusift/core/error.hpp"
#include "nusift/exposure/air_coefficients.hpp"
#include "nusift/units.hpp"

namespace nusift::exposure {
namespace {

constexpr const char* kModule = "exposure";

// Everything in the coefficient except the geometry: energy fluence to air kerma, then kerma
// to exposure per hour. Factored out because gammaConstant() needs exactly this without the
// distance terms.
double kermaToExposurePerDecayEnergy(double energyEv) {
  return energyEv * units::kEvToJ * airMassEnergyAbsorption(energyEv) * units::kSecondsPerHour /
         kGrayPerRoentgen;
}

void requireUsableGeometry(const PointSourceGeometry& geometry) {
  if (!(geometry.distanceM > 0.0)) {
    throw InputError(tagged(kModule,
                            "distance must be positive; a point source has no "
                            "exposure rate defined at zero distance"));
  }
  if (!(geometry.airDensityKgM3 >= 0.0)) {
    throw InputError(tagged(kModule, "air density cannot be negative"));
  }
  if (!(geometry.buildup > 0.0)) {
    throw InputError(tagged(kModule, "buildup factor must be positive"));
  }
}

}  // namespace

double pointExposureCoeff(double energyEv, const PointSourceGeometry& geometry) {
  requireUsableGeometry(geometry);
  if (!(energyEv > 0.0)) {
    return 0.0;
  }

  const double d = geometry.distanceM;
  const double geometric = 1.0 / (4.0 * std::numbers::pi * d * d);

  double attenuation = 1.0;
  if (geometry.airAttenuation) {
    const double linear = airMassAttenuation(energyEv) * geometry.airDensityKgM3;  // 1/m
    attenuation = std::exp(-linear * d);
  }

  return geometric * attenuation * kermaToExposurePerDecayEnergy(energyEv) * geometry.buildup;
}

double gammaConstant(LineSpectrum lines) {
  // Vacuum, so the attenuation exponential is 1 and the 1/(4 pi d^2) factors out to leave a
  // distance-independent constant. This is the only configuration in which a per-nuclide
  // scalar is a complete description of the spectrum.
  double total = 0.0;
  for (const GammaLine& line : lines) {
    if (line.energyEv > 0.0 && line.intensity > 0.0) {
      total += line.intensity * kermaToExposurePerDecayEnergy(line.energyEv);
    }
  }
  return total / (4.0 * std::numbers::pi);
}

double exposureRatePerBecquerel(LineSpectrum lines, const PointSourceGeometry& geometry) {
  requireUsableGeometry(geometry);
  double total = 0.0;
  for (const GammaLine& line : lines) {
    total += line.intensity * pointExposureCoeff(line.energyEv, geometry);
  }
  return total;
}

double exposureRate(LineSpectrum lines, double activityBq, const PointSourceGeometry& geometry) {
  return activityBq * exposureRatePerBecquerel(lines, geometry);
}

}  // namespace nusift::exposure
