#pragma once
/**
 * @file
 * @brief The inventory's photon emission binned into an energy histogram.
 * @ingroup triage
 */
//
// The other half of the photon metric. `Metric::Photon` answers how MANY photons come out --
// one number, over the whole spectrum -- and ranking by gamma line answers which single lines
// carry it. Neither is a source definition: a transport code wants the emission resolved in
// ENERGY and summed over every emitter, which is the axis both of those collapse.
//
// The quantity is emission, not field. There is no geometry here and there deliberately is
// none: what leaves the material is a property of the material, and where it goes is the
// transport code's question, not NuSIFT's. That is why a bin value is photons per second (or
// photons over an interval) and never a fluence -- handing MCNP a number that already had an
// inverse-square in it would be handing it the answer twice.
//
// The histogram is built from the SAME per-line data the exposure sum walks, with no floor.
// The ranking view in response.cpp drops lines below a fraction of their emitter because a
// table with ten thousand rows is not a table anyone reads; a source has no such excuse, since
// binning collapses the rows anyway. So every evaluated line is placed, and
//
//     sum(values) + belowRange + aboveRange == total
//
// holds to floating point. That identity is the whole audit: a deck whose bins do not span the
// spectrum is not wrong here, it is short by a number this file reports.
//
#include <cstdint>
#include <span>
#include <string>
#include <vector>

#include "nusift/engine/decay_result.hpp"
#include "nusift/triage/response.hpp"

namespace nusift {

class NuclearData;

// How generated bin edges are spaced. A decay spectrum runs from X-rays at a few keV to a
// couple of MeV, which is two and a half decades: linear bins resolve the MeV lines that carry
// the dose, logarithmic ones resolve the X-ray end that carries the counts. Neither is right
// for both, which is why the choice is the caller's.
enum class BinScale {
  Linear,
  Logarithmic,
};

bool parseBinScale(std::string_view text, BinScale& out);

// How to lay the bins out. `edgesEv` wins when it is non-empty: explicit edges are how a user
// matches a group structure their transport model already has, and a generated grid that
// silently disagreed with it would be worse than no grid at all.
struct BinningSpec {
  BinScale scale = BinScale::Linear;
  int count = 100;

  // The range to cover, in eV. Zero means "take it from the lines actually present", which is
  // the only default that cannot silently drop photons: a fixed floor of 10 keV would discard
  // every X-ray in the inventory without saying so.
  double minEv = 0.0;
  double maxEv = 0.0;

  // Explicit edges, ascending, at least two. Overrides everything above.
  std::vector<double> edgesEv;
};

// The edges `binning` asks for, given the softest and hardest line present. Throws InputError
// for a spec that cannot produce a grid -- fewer than two explicit edges, edges out of order,
// a non-positive count, an inverted or empty range, or a logarithmic scale reaching zero.
//
// Exposed rather than kept inside the builder because a caller that wants to bin several times
// onto ONE grid -- a source that changes with cooling time is only comparable on fixed edges --
// has to be able to resolve the grid once and hand the same edges back in.
std::vector<double> resolveBinEdges(const BinningSpec& binning, double lowestEv, double highestEv);

// A photon emission spectrum, binned.
struct BinnedSpectrum {
  Domain domain = Domain::Instant;
  // PhotonsPerSecond for an instant, Photons for an interval. Not a free choice: the value in
  // a bin is a rate or a count for the same reason every other quantity here is, and the two
  // fluence units cannot appear because there is no geometry to make a fluence with.
  Unit unit = Unit::PhotonsPerSecond;

  double timeSeconds = 0.0;
  double timeEndSeconds = 0.0;  // interval domain only

  std::vector<double> edgesEv;  // [nBins + 1], ascending
  std::vector<double> values;   // [nBins], in `unit`

  // Every photon the discrete lines emit, whether or not a bin caught it, and the two ways one
  // can miss. Reported rather than folded in, because a source deck built from a grid that does
  // not reach the 2.6 MeV Tl-208 line is missing exactly the photons a shielding calculation
  // was run to find, and the deck itself cannot say so.
  double total = 0.0;
  double belowRange = 0.0;
  double aboveRange = 0.0;

  int lineCount = 0;     // evaluated lines placed in a bin, counted per emitter
  int emitterCount = 0;  // nuclides contributing at least one of them

  // The share of emitted photon ENERGY sitting in continua NuSIFT does not model, over the
  // whole inventory. It belongs in a source report more than anywhere else in the tool: every
  // other metric that carries this flag is understated by it, but a SOURCE that is understated
  // propagates the shortfall through someone else's transport run, where nothing downstream
  // knows to doubt it. An energy fraction, not a photon fraction -- see ResponseTable.
  double unmodeledEnergyFraction = 0.0;
  // The emitters carrying it, sorted, so the note names names.
  std::vector<std::string> unmodeledContinuum;

  int binCount() const { return static_cast<int>(values.size()); }
};

// The spectrum emitted at one time of a solve, in photons per second.
BinnedSpectrum binnedSpectrum(const NuclearData& data, const DecayResult& result, int timeIndex,
                              const BinningSpec& binning);

// The spectrum emitted over one interval, in photons -- a COUNT, integrated exactly, which is
// the right source term for a job that runs while the inventory decays under it. `integral` is
// per-nuclide atom-seconds in the index space of `keys`, exactly as intervalIntegral produces.
BinnedSpectrum binnedIntervalSpectrum(const NuclearData& data, std::span<const std::int64_t> keys,
                                      std::span<const double> integral, double t1, double t2,
                                      const BinningSpec& binning);

}  // namespace nusift
