// The binned photon source: does it conserve, does it agree with the metric it is a view of,
// and do the two decks say what the codes reading them will believe?
//
// The load-bearing test in this file is SourceAgreesWithPhotonMetric. A source term that did
// not sum to the same number `rank --metric photon` prints would mean NuSIFT has two answers
// for how many photons an inventory emits, and the one that leaves the building in an SDEF is
// the one nobody would check.
//
#include <gtest/gtest.h>

#include <cmath>
#include <sstream>
#include <string>
#include <vector>

#include "nusift/core/error.hpp"
#include "nusift/engine/decay_engine.hpp"
#include "nusift/io/source_report.hpp"
#include "nusift/nucdata/nuclear_data.hpp"
#include "nusift/triage/response.hpp"
#include "nusift/triage/spectrum.hpp"
#include "synthetic_chain.hpp"

namespace nusift {
namespace {

// Two emitters with lines far apart in energy, so a bin can hold one and not the other. The
// seeded parent emits nothing, which keeps "has activity" and "has photons" distinguishable.
NuclearData twoEmitterChain() {
  StoreArrays arrays = synth::linearChain({1.0e-3, 5.0e-4});
  synth::addLines(arrays, 1, {50000.0, 661657.0}, {0.4, 0.9});
  return NuclearData::fromArrays(std::move(arrays));
}

Inventory seeded() {
  Inventory inv;
  inv.add(Zai{50, 100, 0}, 1.0e20);
  return inv;
}

// --- bin edges ---------------------------------------------------------------

TEST(SpectrumEdges, LinearSpacingHitsBothEnds) {
  BinningSpec binning;
  binning.count = 4;
  binning.minEv = 0.0;
  binning.maxEv = 0.0;
  const std::vector<double> edges = resolveBinEdges(binning, 100.0, 500.0);
  ASSERT_EQ(edges.size(), 5u);
  EXPECT_DOUBLE_EQ(edges.front(), 100.0);
  EXPECT_DOUBLE_EQ(edges.back(), 500.0);
  EXPECT_DOUBLE_EQ(edges[2], 300.0);
}

TEST(SpectrumEdges, LogarithmicSpacingIsGeometric) {
  BinningSpec binning;
  binning.scale = BinScale::Logarithmic;
  binning.count = 3;
  const std::vector<double> edges = resolveBinEdges(binning, 1.0e3, 1.0e6);
  ASSERT_EQ(edges.size(), 4u);
  EXPECT_DOUBLE_EQ(edges.front(), 1.0e3);
  EXPECT_DOUBLE_EQ(edges.back(), 1.0e6);
  EXPECT_NEAR(edges[1], 1.0e4, 1.0e-6);
  EXPECT_NEAR(edges[2], 1.0e5, 1.0e-5);
}

TEST(SpectrumEdges, ExplicitEdgesPassThrough) {
  BinningSpec binning;
  binning.edgesEv = {1.0, 2.0, 8.0};
  // The generated range is ignored entirely, which is the point of giving edges: a group
  // structure a model already has must not be reshaped to fit the lines that happen to exist.
  const std::vector<double> edges = resolveBinEdges(binning, 1.0e3, 1.0e6);
  EXPECT_EQ(edges, binning.edgesEv);
}

TEST(SpectrumEdges, RefusesGridsThatCannotMeanAnything) {
  BinningSpec zeroBins;
  zeroBins.count = 0;
  EXPECT_THROW(resolveBinEdges(zeroBins, 1.0, 2.0), InputError);

  BinningSpec inverted;
  inverted.minEv = 900.0;
  inverted.maxEv = 100.0;
  EXPECT_THROW(resolveBinEdges(inverted, 1.0, 2.0), InputError);

  BinningSpec logFromZero;
  logFromZero.scale = BinScale::Logarithmic;
  EXPECT_THROW(resolveBinEdges(logFromZero, 0.0, 100.0), InputError);

  BinningSpec oneEdge;
  oneEdge.edgesEv = {5.0};
  EXPECT_THROW(resolveBinEdges(oneEdge, 1.0, 2.0), InputError);

  BinningSpec unsorted;
  unsorted.edgesEv = {5.0, 4.0, 9.0};
  EXPECT_THROW(resolveBinEdges(unsorted, 1.0, 2.0), InputError);

  // Edges and a range are alternatives, and silently preferring one would put a deck on a grid
  // its author did not choose.
  BinningSpec both;
  both.edgesEv = {1.0, 2.0};
  both.minEv = 10.0;
  EXPECT_THROW(resolveBinEdges(both, 1.0, 2.0), InputError);
}

// --- what lands where --------------------------------------------------------

TEST(Spectrum, PlacesEachLineInItsOwnBin) {
  const NuclearData data = twoEmitterChain();
  const DecayResult result = decay(data, seeded(), std::vector<double>{1000.0});

  BinningSpec binning;
  binning.edgesEv = {0.0, 1.0e5, 1.0e6};
  const BinnedSpectrum spectrum = binnedSpectrum(data, result, 0, binning);

  ASSERT_EQ(spectrum.binCount(), 2);
  // Intensities 0.4 and 0.9 on the same emitter, so the two bins stand in that ratio exactly:
  // one decay rate, two lines, no geometry to disturb it.
  EXPECT_GT(spectrum.values[0], 0.0);
  EXPECT_NEAR(spectrum.values[1] / spectrum.values[0], 0.9 / 0.4, 1.0e-12);
  EXPECT_EQ(spectrum.lineCount, 2);
  EXPECT_EQ(spectrum.emitterCount, 1);
  EXPECT_EQ(spectrum.unit, Unit::PhotonsPerSecond);
  EXPECT_EQ(spectrum.domain, Domain::Instant);
}

TEST(Spectrum, ConservesEveryPhotonAcrossTheGrid) {
  const NuclearData data = twoEmitterChain();
  const DecayResult result = decay(data, seeded(), std::vector<double>{1000.0});

  // A grid that deliberately excludes both lines on one side and one on the other.
  BinningSpec binning;
  binning.edgesEv = {1.0e5, 3.0e5};
  const BinnedSpectrum spectrum = binnedSpectrum(data, result, 0, binning);

  double binned = 0.0;
  for (const double value : spectrum.values) {
    binned += value;
  }
  EXPECT_NEAR(binned + spectrum.belowRange + spectrum.aboveRange, spectrum.total,
              spectrum.total * 1.0e-12);
  // The 50 keV line is below the grid and the 661 keV line above it, so nothing is in a bin
  // and the histogram is honestly empty rather than quietly wrong.
  EXPECT_DOUBLE_EQ(binned, 0.0);
  EXPECT_GT(spectrum.belowRange, 0.0);
  EXPECT_GT(spectrum.aboveRange, 0.0);
  EXPECT_EQ(spectrum.lineCount, 0);
}

TEST(Spectrum, ClosesTheTopEdgeSoTheHardestLineIsInside) {
  const NuclearData data = twoEmitterChain();
  const DecayResult result = decay(data, seeded(), std::vector<double>{1000.0});

  // The default grid is built FROM the lines present, so its top edge is the hardest line. A
  // half-open top would put that line out of range in every default invocation.
  const BinnedSpectrum spectrum = binnedSpectrum(data, result, 0, BinningSpec{});
  EXPECT_DOUBLE_EQ(spectrum.aboveRange, 0.0);
  EXPECT_DOUBLE_EQ(spectrum.belowRange, 0.0);
  EXPECT_DOUBLE_EQ(spectrum.edgesEv.back(), 661657.0);
  EXPECT_DOUBLE_EQ(spectrum.edgesEv.front(), 50000.0);
}

// --- agreement with the metric it is a view of -------------------------------

TEST(Spectrum, SourceAgreesWithPhotonMetric) {
  const NuclearData data = twoEmitterChain();
  const DecayResult result = decay(data, seeded(), std::vector<double>{1000.0});

  ResponseSpec spec;
  spec.metric = Metric::Photon;
  spec.unit = Unit::PhotonsPerSecond;
  const ResponseTable table = buildResponse(data, result, spec);

  const BinnedSpectrum spectrum = binnedSpectrum(data, result, 0, BinningSpec{});
  // Two ways of asking how many photons leave the material, and they are the same question.
  EXPECT_NEAR(spectrum.total, table.totals[0], table.totals[0] * 1.0e-12);
}

TEST(Spectrum, IntervalSourceIsACountAndAgreesToo) {
  const NuclearData data = twoEmitterChain();
  const Inventory inv = seeded();
  std::vector<std::int64_t> keys;
  const std::vector<double> integral = intervalIntegral(data, inv, 0.0, 1000.0, &keys);

  ResponseSpec spec;
  spec.metric = Metric::Photon;
  spec.unit = Unit::Photons;
  const ResponseTable table = buildIntervalResponse(data, keys, integral, 0.0, 1000.0, spec);

  const BinnedSpectrum spectrum =
      binnedIntervalSpectrum(data, keys, integral, 0.0, 1000.0, BinningSpec{});
  EXPECT_EQ(spectrum.unit, Unit::Photons);
  EXPECT_EQ(spectrum.domain, Domain::Interval);
  EXPECT_DOUBLE_EQ(spectrum.timeSeconds, 0.0);
  EXPECT_DOUBLE_EQ(spectrum.timeEndSeconds, 1000.0);
  EXPECT_NEAR(spectrum.total, table.totals[0], table.totals[0] * 1.0e-12);
}

// A source built with the ranking view's per-emitter floor would be short by the lines it
// drops. This is the case that would catch that: one line four orders of magnitude below its
// emitter's strongest, which assembleLines() carries but a coarser floor would not.
TEST(Spectrum, KeepsTraceLinesTheRankingWouldNotShow) {
  StoreArrays arrays = synth::linearChain({1.0e-3, 5.0e-4});
  synth::addLines(arrays, 1, {1.0e5, 6.0e5}, {1.0, 1.0e-9});
  const NuclearData data = NuclearData::fromArrays(std::move(arrays));
  const DecayResult result = decay(data, seeded(), std::vector<double>{1000.0});

  BinningSpec binning;
  binning.edgesEv = {5.0e5, 7.0e5};
  const BinnedSpectrum spectrum = binnedSpectrum(data, result, 0, binning);
  EXPECT_GT(spectrum.values[0], 0.0);
  EXPECT_EQ(spectrum.lineCount, 1);
}

// --- the decks ---------------------------------------------------------------

std::string deck(const BinnedSpectrum& spectrum, SourceFormat format) {
  std::ostringstream out;
  writeSourceSpectrum(out, spectrum, ReportContext{}, format);
  return out.str();
}

// Count whitespace-separated numbers following `card` up to the next line that starts a new
// card. Both decks continue a list over several lines, and a writer that miscounted would
// produce a file the reading code rejects rather than one a test would notice by eye.
int countValuesAfter(const std::string& text, const std::string& card) {
  const std::size_t start = text.find(card);
  if (start == std::string::npos) {
    return -1;
  }
  std::istringstream lines(text.substr(start + card.size()));
  std::string line;
  int count = 0;
  while (std::getline(lines, line)) {
    if (line.empty() || (line[0] != ' ' && line[0] != '\t')) {
      break;
    }
    std::istringstream fields(line);
    std::string field;
    while (fields >> field) {
      ++count;
    }
  }
  return count;
}

TEST(SourceDeck, McnpHistogramHasOneMoreBoundaryThanBinAndALeadingZero) {
  const NuclearData data = twoEmitterChain();
  const DecayResult result = decay(data, seeded(), std::vector<double>{1000.0});
  BinningSpec binning;
  binning.count = 7;
  const BinnedSpectrum spectrum = binnedSpectrum(data, result, 0, binning);

  const std::string text = deck(spectrum, SourceFormat::McnpSdef);
  EXPECT_NE(text.find("SDEF PAR=P ERG=D1"), std::string::npos);
  // SI1 carries the 8 boundaries of 7 bins; SP1 carries 8 probabilities, the first standing
  // for the bin below the lowest boundary, which is empty by construction.
  EXPECT_EQ(countValuesAfter(text, "SI1 H\n"), 8);
  EXPECT_EQ(countValuesAfter(text, "SP1 D\n"), 8);
  const std::size_t sp = text.find("SP1 D\n");
  ASSERT_NE(sp, std::string::npos);
  EXPECT_NE(text.find("0.00000E+00", sp), std::string::npos);

  // MeV on the card, eV in the library. A deck in the wrong units is arithmetically fine and
  // physically nonsense, and nothing downstream would flag it.
  EXPECT_NE(text.find("6.61657E-01"), std::string::npos);
}

TEST(SourceDeck, McnpSaysTheCardCarriesNoStrength) {
  const NuclearData data = twoEmitterChain();
  const DecayResult result = decay(data, seeded(), std::vector<double>{1000.0});
  const std::string text =
      deck(binnedSpectrum(data, result, 0, BinningSpec{}), SourceFormat::McnpSdef);
  // The one thing a reader must not assume. MCNP normalizes SP, so a deck that did not say so
  // would silently be off by fourteen orders of magnitude.
  EXPECT_NE(text.find("SHAPE, NOT A STRENGTH"), std::string::npos);
  EXPECT_NE(text.find("photons/s"), std::string::npos);
}

TEST(SourceDeck, OpenmcEmissionIsOneShorterThanEdges) {
  const NuclearData data = twoEmitterChain();
  const DecayResult result = decay(data, seeded(), std::vector<double>{1000.0});
  BinningSpec binning;
  binning.count = 5;
  const BinnedSpectrum spectrum = binnedSpectrum(data, result, 0, binning);

  const std::string text = deck(spectrum, SourceFormat::OpenmcPython);
  EXPECT_EQ(countValuesAfter(text, "edges = np.array([\n"), 6);
  EXPECT_EQ(countValuesAfter(text, "emission = np.array([\n"), 5);
  EXPECT_NE(text.find("interpolation=\"histogram\""), std::string::npos);
  EXPECT_NE(text.find("particle=\"photon\""), std::string::npos);
  EXPECT_NE(text.find("strength="), std::string::npos);
}

// A truncated grid is the failure this whole design exists to make visible, so both decks have
// to carry the warning -- not just the report a person reads.
TEST(SourceDeck, WarnsInsideTheDeckWhenEmissionFallsOutsideTheBins) {
  const NuclearData data = twoEmitterChain();
  const DecayResult result = decay(data, seeded(), std::vector<double>{1000.0});
  BinningSpec binning;
  binning.edgesEv = {1.0e5, 3.0e5};
  const BinnedSpectrum spectrum = binnedSpectrum(data, result, 0, binning);

  EXPECT_NE(deck(spectrum, SourceFormat::McnpSdef).find("falls outside the bins"),
            std::string::npos);
  EXPECT_NE(deck(spectrum, SourceFormat::OpenmcPython).find("falls outside the bins"),
            std::string::npos);
  EXPECT_NE(deck(spectrum, SourceFormat::Text).find("falls outside the bins"), std::string::npos);
}

TEST(SourceDeck, CsvCarriesTheOutOfRangeTalliesWithTheBins) {
  const NuclearData data = twoEmitterChain();
  const DecayResult result = decay(data, seeded(), std::vector<double>{1000.0});
  BinningSpec binning;
  binning.edgesEv = {1.0e5, 3.0e5};
  const std::string text = deck(binnedSpectrum(data, result, 0, binning), SourceFormat::Csv);
  // A header a spreadsheet drops would take the audit with it, so the three tallies are rows.
  EXPECT_NE(text.find("\nbelow,"), std::string::npos);
  EXPECT_NE(text.find("\nabove,"), std::string::npos);
  EXPECT_NE(text.find("\ntotal,"), std::string::npos);
}

TEST(SourceFormats, ParseTheSpellingsTheFrontEndsOffer) {
  SourceFormat format = SourceFormat::Text;
  EXPECT_TRUE(parseSourceFormat("mcnp", format));
  EXPECT_EQ(format, SourceFormat::McnpSdef);
  EXPECT_TRUE(parseSourceFormat("sdef", format));
  EXPECT_EQ(format, SourceFormat::McnpSdef);
  EXPECT_TRUE(parseSourceFormat("openmc", format));
  EXPECT_EQ(format, SourceFormat::OpenmcPython);
  EXPECT_FALSE(parseSourceFormat("serpent", format));

  BinScale scale = BinScale::Linear;
  EXPECT_TRUE(parseBinScale("log", scale));
  EXPECT_EQ(scale, BinScale::Logarithmic);
  EXPECT_FALSE(parseBinScale("quadratic", scale));
}

}  // namespace
}  // namespace nusift
