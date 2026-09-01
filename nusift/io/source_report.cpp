#include "nusift/io/source_report.hpp"

#include <cmath>
#include <cstdio>
#include <ostream>
#include <string>
#include <vector>

#include "nusift/io/number_format.hpp"
#include "nusift/io/time_spec.hpp"
#include "nusift/version.hpp"

namespace nusift {
namespace {

std::string sci(double value) {
  char buffer[32];
  std::snprintf(buffer, sizeof(buffer), "%.4e", value);
  return buffer;
}

std::string jsonNumber(double value) {
  return std::isfinite(value) ? shortestRoundTrip(value) : std::string("null");
}

// keV, because that is the unit a decay photon is spoken in and the one the line ranking
// already prints. The decks convert on the way out, each to what its code reads.
std::string keV(double energyEv) {
  char buffer[32];
  std::snprintf(buffer, sizeof(buffer), "%.4g", energyEv / 1000.0);
  return buffer;
}

const char* unitText(const BinnedSpectrum& spectrum) {
  return spectrum.domain == Domain::Interval ? "photons" : "photons/s";
}

// "at 30 d" or "over 1 h to 30 d". The distinction is the difference between a rate and a
// count, and it is the first thing anyone reading a source term needs to know.
std::string whenText(const BinnedSpectrum& spectrum) {
  if (spectrum.domain == Domain::Interval) {
    return "over " + formatDuration(spectrum.timeSeconds) + " to " +
           formatDuration(spectrum.timeEndSeconds);
  }
  return "at " + formatDuration(spectrum.timeSeconds);
}

// The sentence both decks need and neither code can infer: what the tally has to be multiplied
// by, and why the deck cannot do it. MCNP normalises SP; OpenMC scores per source particle.
// Different mechanisms, same arithmetic left to the user, so it is said once and worded twice.
std::string strengthNote(const BinnedSpectrum& spectrum) {
  return sci(spectrum.total) + " " + unitText(spectrum);
}

// What the histogram does not carry. Returns an empty vector when there is nothing to say,
// which is the common case for a clean inventory on a default grid.
std::vector<std::string> shortfalls(const BinnedSpectrum& spectrum) {
  std::vector<std::string> notes;
  const double outside = spectrum.belowRange + spectrum.aboveRange;
  if (outside > 0.0 && spectrum.total > 0.0) {
    char share[32];
    std::snprintf(share, sizeof(share), "%.3g%%", 100.0 * outside / spectrum.total);
    std::string note = std::string(share) + " of the emission falls outside the bins (";
    note += sci(spectrum.belowRange) + " below " + keV(spectrum.edgesEv.front()) + " keV, ";
    note += sci(spectrum.aboveRange) + " above " + keV(spectrum.edgesEv.back()) +
            " keV) and is NOT in this source. Widen the range to include it.";
    notes.push_back(note);
  }
  if (spectrum.unmodeledEnergyFraction > 0.0) {
    char share[32];
    std::snprintf(share, sizeof(share), "%.3g%%", 100.0 * spectrum.unmodeledEnergyFraction);
    std::string note = std::string(share) +
                       " of the emitted photon ENERGY is in continua NuSIFT does not model "
                       "(bremsstrahlung, internal conversion continua), so this source is "
                       "understated by about that much and softer than it looks.";
    if (!spectrum.unmodeledContinuum.empty()) {
      note += " Carried by: ";
      for (std::size_t i = 0; i < spectrum.unmodeledContinuum.size(); ++i) {
        note += (i == 0 ? "" : ", ") + spectrum.unmodeledContinuum[i];
      }
      note += ".";
    }
    notes.push_back(note);
  }
  return notes;
}

// Wrap `text` to `width` columns. The two prefixes differ so a marked note -- "c  ! " on a
// comment card, "  ! " in a report -- keeps its marker on the first line and stays aligned
// under it on the rest. Comment cards in both decks are read by people, and a 400-character
// line is read by nobody.
void writeWrapped(std::ostream& out, const std::string& firstPrefix, const std::string& contPrefix,
                  const std::string& text, std::size_t width) {
  std::size_t start = 0;
  bool first = true;
  while (start < text.size()) {
    std::size_t take = std::min(width, text.size() - start);
    if (start + take < text.size()) {
      const std::size_t space = text.rfind(' ', start + take);
      if (space != std::string::npos && space > start) {
        take = space - start;
      }
    }
    out << (first ? firstPrefix : contPrefix) << text.substr(start, take) << "\n";
    first = false;
    start += take;
    while (start < text.size() && text[start] == ' ') {
      ++start;
    }
  }
}

// The common case: one prefix on every line.
void writeWrapped(std::ostream& out, const std::string& prefix, const std::string& text,
                  std::size_t width) {
  writeWrapped(out, prefix, prefix, text, width);
}

void writeProvenance(std::ostream& out, const BinnedSpectrum& spectrum,
                     const ReportContext& context, const std::string& prefix) {
  out << prefix << "Photon emission source, written by NuSIFT " << NUSIFT_VERSION_STRING << "\n";
  if (!context.seedProvenance.empty()) {
    writeWrapped(out, prefix, "inventory: " + context.seedProvenance, 72);
  }
  if (!context.storePath.empty()) {
    std::string store = "store: " + context.storePath;
    if (!context.storeLibrary.empty()) {
      store += " (" + context.storeLibrary;
      if (!context.storeCreatedUtc.empty()) {
        store += ", staged " + context.storeCreatedUtc;
      }
      store += ")";
    }
    writeWrapped(out, prefix, store, 72);
  }
  out << prefix << "emission " << whenText(spectrum) << ": " << strengthNote(spectrum) << " over "
      << spectrum.binCount() << " bins, " << keV(spectrum.edgesEv.front()) << " keV to "
      << keV(spectrum.edgesEv.back()) << " keV\n";
  out << prefix << spectrum.lineCount << " evaluated lines from " << spectrum.emitterCount
      << " emitters\n";
}

// --- text -------------------------------------------------------------------

void writeText(std::ostream& out, const BinnedSpectrum& spectrum, const ReportContext& context) {
  out << "NuSIFT photon emission source\n";
  out << "  " << whenText(spectrum) << "    total = " << sci(spectrum.total) << " "
      << unitText(spectrum) << "\n";
  out << "  grid:  " << spectrum.binCount() << " " << keV(spectrum.edgesEv.front()) << "-"
      << keV(spectrum.edgesEv.back()) << " keV bins over " << spectrum.lineCount << " lines from "
      << spectrum.emitterCount << " emitters\n";
  // Said in the header rather than in a footnote, because it is the one thing about this
  // quantity a reader is most likely to assume wrongly: every other number this tool prints
  // with an energy on it is a field somewhere.
  out << "  model: emission from the whole inventory -- no geometry, no distance, no "
         "self-absorption\n";
  if (!context.seedProvenance.empty()) {
    out << "  seed:  " << context.seedProvenance << "\n";
  }
  if (!context.storeLibrary.empty()) {
    out << "  store: " << context.storeLibrary << " (" << context.storeNuclideCount
        << " nuclides, staged " << context.storeCreatedUtc << ")\n";
  }
  out << "\n";

  out << "  bin      E_low [keV]    E_high [keV]     " << unitText(spectrum) << "      share\n";
  out << "  ------------------------------------------------------------------------\n";
  const double binned = spectrum.total - spectrum.belowRange - spectrum.aboveRange;
  for (int b = 0; b < spectrum.binCount(); ++b) {
    const double value = spectrum.values[static_cast<std::size_t>(b)];
    char row[160];
    std::snprintf(row, sizeof(row), "  %4d  %14.4f  %14.4f  %14s  %7.3f%%\n", b + 1,
                  spectrum.edgesEv[static_cast<std::size_t>(b)] / 1000.0,
                  spectrum.edgesEv[static_cast<std::size_t>(b) + 1] / 1000.0, sci(value).c_str(),
                  binned > 0.0 ? 100.0 * value / binned : 0.0);
    out << row;
  }
  out << "  ------------------------------------------------------------------------\n";
  out << "  in bins  " << sci(binned) << " " << unitText(spectrum) << "\n";
  out << "  total    " << sci(spectrum.total) << " " << unitText(spectrum)
      << "   (every evaluated line, binned or not)\n";

  const std::vector<std::string> notes = shortfalls(spectrum);
  if (!notes.empty()) {
    out << "\n";
    for (const std::string& note : notes) {
      writeWrapped(out, "  ! ", "    ", note, 86);
    }
  }
}

// --- csv --------------------------------------------------------------------

void writeCsv(std::ostream& out, const BinnedSpectrum& spectrum) {
  out << "bin,energy_low_ev,energy_high_ev,emission,unit\n";
  for (int b = 0; b < spectrum.binCount(); ++b) {
    out << (b + 1) << "," << shortestRoundTrip(spectrum.edgesEv[static_cast<std::size_t>(b)]) << ","
        << shortestRoundTrip(spectrum.edgesEv[static_cast<std::size_t>(b) + 1]) << ","
        << shortestRoundTrip(spectrum.values[static_cast<std::size_t>(b)]) << ","
        << unitText(spectrum) << "\n";
  }
  // The three tallies that make the table auditable travel with it rather than in a header a
  // spreadsheet would drop: their sum is the first column's sum plus what missed the grid.
  out << "below," << shortestRoundTrip(0.0) << "," << shortestRoundTrip(spectrum.edgesEv.front())
      << "," << shortestRoundTrip(spectrum.belowRange) << "," << unitText(spectrum) << "\n";
  out << "above," << shortestRoundTrip(spectrum.edgesEv.back()) << ",,"
      << shortestRoundTrip(spectrum.aboveRange) << "," << unitText(spectrum) << "\n";
  out << "total,,," << shortestRoundTrip(spectrum.total) << "," << unitText(spectrum) << "\n";
}

// --- json -------------------------------------------------------------------

void writeJson(std::ostream& out, const BinnedSpectrum& spectrum, const ReportContext& context) {
  out << "{\n";
  out << "  \"quantity\": \"photon emission\",\n";
  out << "  \"unit\": \"" << unitText(spectrum) << "\",\n";
  out << "  \"domain\": \"" << (spectrum.domain == Domain::Interval ? "interval" : "instant")
      << "\",\n";
  out << "  \"time_seconds\": " << jsonNumber(spectrum.timeSeconds) << ",\n";
  if (spectrum.domain == Domain::Interval) {
    out << "  \"time_end_seconds\": " << jsonNumber(spectrum.timeEndSeconds) << ",\n";
  }
  out << "  \"store\": \"" << context.storePath << "\",\n";
  out << "  \"library\": \"" << context.storeLibrary << "\",\n";

  out << "  \"edges_ev\": [";
  for (std::size_t i = 0; i < spectrum.edgesEv.size(); ++i) {
    out << (i == 0 ? "" : ", ") << jsonNumber(spectrum.edgesEv[i]);
  }
  out << "],\n";
  out << "  \"emission\": [";
  for (std::size_t i = 0; i < spectrum.values.size(); ++i) {
    out << (i == 0 ? "" : ", ") << jsonNumber(spectrum.values[i]);
  }
  out << "],\n";

  out << "  \"total\": " << jsonNumber(spectrum.total) << ",\n";
  out << "  \"below_range\": " << jsonNumber(spectrum.belowRange) << ",\n";
  out << "  \"above_range\": " << jsonNumber(spectrum.aboveRange) << ",\n";
  out << "  \"line_count\": " << spectrum.lineCount << ",\n";
  out << "  \"emitter_count\": " << spectrum.emitterCount << ",\n";
  out << "  \"unmodeled_energy_fraction\": " << jsonNumber(spectrum.unmodeledEnergyFraction)
      << ",\n";
  out << "  \"unmodeled_continuum\": [";
  for (std::size_t i = 0; i < spectrum.unmodeledContinuum.size(); ++i) {
    out << (i == 0 ? "" : ", ") << "\"" << spectrum.unmodeledContinuum[i] << "\"";
  }
  out << "]\n";
  out << "}\n";
}

// --- MCNP -------------------------------------------------------------------

// Values per continuation line. MCNP's classic card is 80 columns and a continuation must be
// indented; five values at %.5E fit inside that with room to spare.
constexpr int kMcnpPerLine = 5;

void writeMcnpValues(std::ostream& out, const std::vector<double>& values, double scale) {
  for (std::size_t i = 0; i < values.size(); ++i) {
    if (i % kMcnpPerLine == 0) {
      out << (i == 0 ? "     " : "\n     ");
    }
    char buffer[24];
    std::snprintf(buffer, sizeof(buffer), " %13.5E", values[i] * scale);
    out << buffer;
  }
  out << "\n";
}

void writeMcnp(std::ostream& out, const BinnedSpectrum& spectrum, const ReportContext& context) {
  out << "c ==========================================================================\n";
  writeProvenance(out, spectrum, context, "c  ");
  out << "c --------------------------------------------------------------------------\n";
  out << "c  THE SP CARD BELOW IS A SHAPE, NOT A STRENGTH. MCNP normalises SP to unity\n";
  out << "c  and scores every tally per source particle, so multiply your tally by\n";
  out << "c\n";
  out << "c      " << strengthNote(spectrum) << "\n";
  out << "c\n";
  writeWrapped(out, "c  ",
               "to get an absolute result. An FM card is the usual place for it. This deck "
               "cannot apply it for you: nothing in an SDEF carries an absolute rate.",
               70);
  out << "c\n";
  writeWrapped(out, "c  ",
               "The source is EMISSION from the material and carries no geometry -- no "
               "distance, no self-absorption, no container. POS, CEL, RAD and AXS are yours "
               "to set; with none given MCNP puts a point source at the origin.",
               70);
  const std::vector<std::string> notes = shortfalls(spectrum);
  if (!notes.empty()) {
    out << "c\n";
    for (const std::string& note : notes) {
      writeWrapped(out, "c  ! ", "c    ", note, 66);
    }
  }
  out << "c\n";
  writeWrapped(out, "c  ",
               "Binning is lossy in one direction only: bin totals are exact, but the discrete "
               "lines inside a bin cannot be recovered from it. For a detector-resolution "
               "question use the line list, not this.",
               70);
  out << "c  PAR=P is MCNP6; write PAR=2 for MCNP5.\n";
  out << "c ==========================================================================\n";
  out << "SDEF PAR=P ERG=D1\n";
  out << "c  SI1: bin boundaries in MeV. SP1: the first entry is the bin below the\n";
  out << "c  lowest boundary, which is empty by construction.\n";
  out << "SI1 H\n";
  writeMcnpValues(out, spectrum.edgesEv, 1.0e-6);
  out << "SP1 D\n";
  std::vector<double> probabilities;
  probabilities.reserve(spectrum.values.size() + 1);
  probabilities.push_back(0.0);
  for (const double value : spectrum.values) {
    probabilities.push_back(value);
  }
  writeMcnpValues(out, probabilities, 1.0);
}

// --- OpenMC -----------------------------------------------------------------

void writePythonValues(std::ostream& out, const std::vector<double>& values, double scale) {
  for (std::size_t i = 0; i < values.size(); ++i) {
    if (i % 4 == 0) {
      out << (i == 0 ? "    " : "\n    ");
    }
    out << shortestRoundTrip(values[i] * scale) << ",";
    if ((i + 1) % 4 != 0 && i + 1 != values.size()) {
      out << " ";
    }
  }
  out << "\n";
}

void writeOpenmc(std::ostream& out, const BinnedSpectrum& spectrum, const ReportContext& context) {
  out << "# ==========================================================================\n";
  writeProvenance(out, spectrum, context, "#  ");
  out << "# --------------------------------------------------------------------------\n";
  writeWrapped(out, "#  ",
               "OpenMC scores every tally per source particle, so `strength` below records "
               "the absolute rate but does not apply it: multiply your tally by "
               "source.strength to get an absolute result.",
               70);
  out << "#\n";
  writeWrapped(out, "#  ",
               "The source is EMISSION from the material and carries no geometry -- no "
               "distance, no self-absorption, no container. The `space` distribution below is "
               "a placeholder point at the origin; replace it with your geometry.",
               70);
  const std::vector<std::string> notes = shortfalls(spectrum);
  if (!notes.empty()) {
    out << "#\n";
    for (const std::string& note : notes) {
      writeWrapped(out, "#  ! ", "#    ", note, 66);
    }
  }
  out << "#\n";
  writeWrapped(out, "#  ",
               "Binning is lossy in one direction only: bin totals are exact, but the discrete "
               "lines inside a bin cannot be recovered from it. For a detector-resolution "
               "question use the line list, not this.",
               70);
  out << "# ==========================================================================\n";
  out << "import numpy as np\n";
  out << "import openmc\n\n";

  out << "# Bin boundaries in eV, ascending.\n";
  out << "edges = np.array([\n";
  writePythonValues(out, spectrum.edgesEv, 1.0);
  out << "], dtype=float)\n\n";

  out << "# Emission in each bin, in " << unitText(spectrum)
      << ". Histogram interpolation takes one\n";
  out << "# value per bin, so this is exactly one shorter than `edges`.\n";
  out << "emission = np.array([\n";
  writePythonValues(out, spectrum.values, 1.0);
  out << "], dtype=float)\n\n";

  out << "energy = openmc.stats.Tabular(edges, emission, interpolation=\"histogram\")\n\n";
  out << "source = openmc.IndependentSource(\n";
  out << "    space=openmc.stats.Point((0.0, 0.0, 0.0)),  # placeholder -- your geometry here\n";
  out << "    angle=openmc.stats.Isotropic(),\n";
  out << "    energy=energy,\n";
  out << "    particle=\"photon\",\n";
  out << "    strength=" << shortestRoundTrip(spectrum.total) << ",  # " << unitText(spectrum)
      << "\n";
  out << ")\n";
}

}  // namespace

bool parseSourceFormat(std::string_view text, SourceFormat& out) {
  if (text == "text") {
    out = SourceFormat::Text;
    return true;
  }
  if (text == "csv") {
    out = SourceFormat::Csv;
    return true;
  }
  if (text == "json") {
    out = SourceFormat::Json;
    return true;
  }
  if (text == "mcnp" || text == "sdef") {
    out = SourceFormat::McnpSdef;
    return true;
  }
  if (text == "openmc") {
    out = SourceFormat::OpenmcPython;
    return true;
  }
  return false;
}

void writeSourceSpectrum(std::ostream& out, const BinnedSpectrum& spectrum,
                         const ReportContext& context, SourceFormat format) {
  switch (format) {
    case SourceFormat::Text:
      writeText(out, spectrum, context);
      return;
    case SourceFormat::Csv:
      writeCsv(out, spectrum);
      return;
    case SourceFormat::Json:
      writeJson(out, spectrum, context);
      return;
    case SourceFormat::McnpSdef:
      writeMcnp(out, spectrum, context);
      return;
    case SourceFormat::OpenmcPython:
      writeOpenmc(out, spectrum, context);
      return;
  }
}

}  // namespace nusift
