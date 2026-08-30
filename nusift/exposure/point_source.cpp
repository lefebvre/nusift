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
         units::kGyPerR;
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

// The part of the point kernel that is geometry and path only: spreading over 4*pi*d^2,
// attenuation along the path, and the explicit buildup factor. Shared by the exposure and the
// fluence coefficients so the two can never disagree about where a photon gets to.
double spreadingAndAttenuation(double energyEv, const PointSourceGeometry& geometry) {
  const double d = geometry.distanceM;
  const double geometric = 1.0 / (4.0 * std::numbers::pi * d * d);

  double attenuation = 1.0;
  if (geometry.airAttenuation) {
    const double linear = airMassAttenuation(energyEv) * geometry.airDensityKgM3;  // 1/m
    attenuation = std::exp(-linear * d);
  }

  return geometric * attenuation * geometry.buildup;
}

}  // namespace

double pointExposureCoeff(double energyEv, const PointSourceGeometry& geometry) {
  requireUsableGeometry(geometry);
  if (!(energyEv > 0.0)) {
    return 0.0;
  }
  return spreadingAndAttenuation(energyEv, geometry) * kermaToExposurePerDecayEnergy(energyEv);
}

double pointFluenceCoeff(double energyEv, const PointSourceGeometry& geometry) {
  requireUsableGeometry(geometry);
  if (!(energyEv > 0.0)) {
    return 0.0;
  }
  return spreadingAndAttenuation(energyEv, geometry);
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

double fluenceRatePerBecquerel(LineSpectrum lines, const PointSourceGeometry& geometry) {
  requireUsableGeometry(geometry);
  double total = 0.0;
  for (const GammaLine& line : lines) {
    total += line.intensity * pointFluenceCoeff(line.energyEv, geometry);
  }
  return total;
}

double exposureRate(LineSpectrum lines, double activityBq, const PointSourceGeometry& geometry) {
  return activityBq * exposureRatePerBecquerel(lines, geometry);
}

double meanOpticalDepth(LineSpectrum lines, const PointSourceGeometry& geometry) {
  requireUsableGeometry(geometry);
  if (!geometry.airAttenuation) {
    return 0.0;
  }
  // Weighted by what arrives rather than by what is emitted: a soft line that the path has
  // already absorbed should not pull the average toward its own thickness, because it is not
  // carrying the exposure the buildup question is about.
  double weighted = 0.0;
  double total = 0.0;
  for (const GammaLine& line : lines) {
    const double share = line.intensity * pointExposureCoeff(line.energyEv, geometry);
    if (share <= 0.0) {
      continue;
    }
    const double depth =
        airMassAttenuation(line.energyEv) * geometry.airDensityKgM3 * geometry.distanceM;
    weighted += share * depth;
    total += share;
  }
  return total > 0.0 ? weighted / total : 0.0;
}

}  // namespace nusift::exposure
