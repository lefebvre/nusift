#include "nusift/triage/spectrum.hpp"

#include <algorithm>
#include <cmath>
#include <set>
#include <string>

#include "nusift/core/error.hpp"
#include "nusift/core/nuclide.hpp"
#include "nusift/core/nuclide_name.hpp"
#include "nusift/nucdata/nuclear_data.hpp"
#include "nusift/nucdata/photon_lines.hpp"

namespace nusift {
namespace {

constexpr const char* kModule = "spectrum";

// The same threshold response.cpp flags a ranking row at, for the same reason: below a percent
// the continuum is smaller than the spread between evaluations and does not earn a name.
constexpr double kUnmodeledContinuumFlag = 0.01;

[[noreturn]] void fail(const std::string& what) {
  throw InputError(tagged(kModule, what));
}

// Locate a value on an ascending edge grid. Returns the bin index, or -1 below the grid and
// nBins at or above the top edge.
//
// The top edge is CLOSED -- a line exactly at maxEv lands in the last bin rather than falling
// off the end -- because the default grid is built from the hardest line present, so the
// half-open convention would put that line out of range in every default invocation. Every
// other boundary stays half-open, which is what a histogram means.
std::ptrdiff_t binOf(std::span<const double> edges, double energyEv) {
  const std::size_t nBins = edges.size() - 1;
  if (energyEv < edges.front()) {
    return -1;
  }
  if (energyEv > edges.back()) {
    return static_cast<std::ptrdiff_t>(nBins);
  }
  if (energyEv == edges.back()) {
    return static_cast<std::ptrdiff_t>(nBins) - 1;
  }
  const auto it = std::upper_bound(edges.begin(), edges.end(), energyEv);
  return std::distance(edges.begin(), it) - 1;
}

// The softest and hardest evaluated line over the nuclides present with a positive amount.
// Bounded by what is actually THERE, not by what the store carries: an inventory of Cs-137
// should not get a grid stretched to Tl-208's 2.6 MeV because some other nuclide in the
// closure could have emitted it.
struct LineRange {
  double lowestEv = 0.0;
  double highestEv = 0.0;
  bool any = false;
};

LineRange lineRange(const NuclearData& data, std::span<const std::int64_t> keys,
                    std::span<const double> amounts) {
  LineRange range;
  for (std::size_t i = 0; i < keys.size(); ++i) {
    if (!(amounts[i] > 0.0)) {
      continue;
    }
    const int index = data.indexOfKey(keys[i]);
    if (index < 0 || !(data.decayConstant(index) > 0.0)) {
      continue;
    }
    for (const GammaLine& line : data.lines(index)) {
      if (!(line.energyEv > 0.0 && line.intensity > 0.0)) {
        continue;
      }
      if (!range.any) {
        range.lowestEv = line.energyEv;
        range.highestEv = line.energyEv;
        range.any = true;
        continue;
      }
      range.lowestEv = std::min(range.lowestEv, line.energyEv);
      range.highestEv = std::max(range.highestEv, line.energyEv);
    }
  }
  return range;
}

// The body both entry points share. `amounts` is atoms for an instant and atom-seconds for an
// interval, and multiplying by lambda turns either into the right thing without a second code
// path: lambda*n is photons per second, lambda times an integral is a count of them.
BinnedSpectrum build(const NuclearData& data, std::span<const std::int64_t> keys,
                     std::span<const double> amounts, const BinningSpec& binning) {
  if (keys.size() != amounts.size()) {
    fail("keys and amounts are different lengths");
  }

  const LineRange range = lineRange(data, keys, amounts);
  // An inventory with no evaluated photon lines at all still gets a grid, because the caller
  // asked for one and an empty histogram is a truthful answer: this material emits no photons
  // NuSIFT models. Explicit edges are honoured; a generated grid needs a range to generate on,
  // and 1 keV to 10 MeV is the span decay photons occupy.
  const double lowest = range.any ? range.lowestEv : 1.0e3;
  const double highest = range.any ? range.highestEv : 1.0e7;

  BinnedSpectrum spectrum;
  spectrum.edgesEv = resolveBinEdges(binning, lowest, highest);
  spectrum.values.assign(spectrum.edgesEv.size() - 1, 0.0);

  std::set<std::string> continuum;
  std::set<std::int64_t> emitters;
  double modelledEnergy = 0.0;
  double unmodeledEnergy = 0.0;

  for (std::size_t i = 0; i < keys.size(); ++i) {
    const double amount = amounts[i];
    if (!(amount > 0.0)) {
      continue;
    }
    const int index = data.indexOfKey(keys[i]);
    if (index < 0) {
      continue;
    }
    const double lambda = data.decayConstant(index);
    if (!(lambda > 0.0)) {
      continue;
    }
    const double decays = lambda * amount;

    // The continuum accounting is over the whole inventory and not only over emitters that
    // reached a bin, for the reason contextFor() gives on the ranking side: a nuclide whose
    // photon output is ENTIRELY continuum places no line anywhere, so a scan over what landed
    // in the histogram would drop precisely the worst case.
    //
    // Modeled energy is the DISCRETE line sum rather than the staged average electromagnetic
    // energy, which matches unmodeledEnergyFractions() in response.cpp exactly. The two are
    // close but not equal, and the whole value of this figure is that a source deck and a
    // ranking of the same inventory cannot report different amounts of missing photon power.
    modelledEnergy += decays * discretePhotonEnergyEv(data.lines(index));
    unmodeledEnergy += decays * data.continuumPhotonEv(index);
    if (data.unmodeledPhotonFraction(index) > kUnmodeledContinuumFlag) {
      continuum.insert(formatNuclideName(Zai::fromKey(keys[i])));
    }

    for (const GammaLine& line : data.lines(index)) {
      if (!(line.energyEv > 0.0 && line.intensity > 0.0)) {
        continue;
      }
      const double emission = decays * line.intensity;
      spectrum.total += emission;
      const std::ptrdiff_t bin = binOf(spectrum.edgesEv, line.energyEv);
      if (bin < 0) {
        spectrum.belowRange += emission;
      } else if (bin >= static_cast<std::ptrdiff_t>(spectrum.values.size())) {
        spectrum.aboveRange += emission;
      } else {
        spectrum.values[static_cast<std::size_t>(bin)] += emission;
        ++spectrum.lineCount;
        emitters.insert(keys[i]);
      }
    }
  }

  spectrum.emitterCount = static_cast<int>(emitters.size());
  const double photonEnergy = modelledEnergy + unmodeledEnergy;
  spectrum.unmodeledEnergyFraction = photonEnergy > 0.0 ? unmodeledEnergy / photonEnergy : 0.0;
  spectrum.unmodeledContinuum.assign(continuum.begin(), continuum.end());
  return spectrum;
}

}  // namespace

bool parseBinScale(std::string_view text, BinScale& out) {
  if (text == "linear" || text == "lin") {
    out = BinScale::Linear;
    return true;
  }
  if (text == "logarithmic" || text == "log") {
    out = BinScale::Logarithmic;
    return true;
  }
  return false;
}

std::vector<double> resolveBinEdges(const BinningSpec& binning, double lowestEv, double highestEv) {
  if (!binning.edgesEv.empty()) {
    // Refused rather than silently overridden: a caller who gave both a grid and a shape for
    // one meant one of them, and guessing which is how a deck ends up on edges nobody chose.
    if (binning.minEv > 0.0 || binning.maxEv > 0.0) {
      fail(
          "explicit edges already fix the range, so a minimum or maximum energy cannot also "
          "apply. Give one or the other");
    }
    if (binning.edgesEv.size() < 2) {
      fail("explicit edges need at least two values, which is one bin");
    }
    for (std::size_t i = 1; i < binning.edgesEv.size(); ++i) {
      if (!(binning.edgesEv[i] > binning.edgesEv[i - 1])) {
        fail("explicit edges must be strictly increasing");
      }
    }
    if (!(binning.edgesEv.front() >= 0.0)) {
      fail("explicit edges must be non-negative energies in eV");
    }
    return binning.edgesEv;
  }

  if (binning.count < 1) {
    fail("bin count must be at least 1, got " + std::to_string(binning.count));
  }

  const double low = binning.minEv > 0.0 ? binning.minEv : lowestEv;
  const double high = binning.maxEv > 0.0 ? binning.maxEv : highestEv;
  if (!(high > low)) {
    fail("the energy range must end above where it starts, got " + std::to_string(low) + " to " +
         std::to_string(high) + " eV");
  }
  if (binning.scale == BinScale::Logarithmic && !(low > 0.0)) {
    fail("logarithmic bins cannot start at zero; give a positive minimum energy");
  }

  std::vector<double> edges(static_cast<std::size_t>(binning.count) + 1);
  if (binning.scale == BinScale::Logarithmic) {
    const double logLow = std::log(low);
    const double step = (std::log(high) - logLow) / binning.count;
    for (int k = 0; k <= binning.count; ++k) {
      edges[static_cast<std::size_t>(k)] = std::exp(logLow + step * k);
    }
  } else {
    const double step = (high - low) / binning.count;
    for (int k = 0; k <= binning.count; ++k) {
      edges[static_cast<std::size_t>(k)] = low + step * k;
    }
  }
  // Rounding in the loop can leave the last edge a few ulps off the range the caller named,
  // which would put the hardest line just outside a grid built to contain it.
  edges.front() = low;
  edges.back() = high;
  return edges;
}

BinnedSpectrum binnedSpectrum(const NuclearData& data, const DecayResult& result, int timeIndex,
                              const BinningSpec& binning) {
  if (timeIndex < 0 || timeIndex >= result.timeCount()) {
    fail("time index " + std::to_string(timeIndex) + " is outside the solve");
  }
  BinnedSpectrum spectrum = build(data, result.nuclideKeys, result.atomsAt(timeIndex), binning);
  spectrum.domain = Domain::Instant;
  spectrum.unit = Unit::PhotonsPerSecond;
  spectrum.timeSeconds = result.times[static_cast<std::size_t>(timeIndex)];
  spectrum.timeEndSeconds = spectrum.timeSeconds;
  return spectrum;
}

BinnedSpectrum binnedIntervalSpectrum(const NuclearData& data, std::span<const std::int64_t> keys,
                                      std::span<const double> integral, double t1, double t2,
                                      const BinningSpec& binning) {
  if (!(t2 > t1)) {
    fail("an interval must end after it starts");
  }
  BinnedSpectrum spectrum = build(data, keys, integral, binning);
  spectrum.domain = Domain::Interval;
  spectrum.unit = Unit::Photons;
  spectrum.timeSeconds = t1;
  spectrum.timeEndSeconds = t2;
  return spectrum;
}

}  // namespace nusift
