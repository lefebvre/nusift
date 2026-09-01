#include "nusift/io/report.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <ostream>
#include <string>

#include "nusift/core/error.hpp"
#include "nusift/io/number_format.hpp"
#include "nusift/io/time_spec.hpp"
#include "nusift/triage/response.hpp"

namespace nusift {
namespace {

constexpr const char* kModule = "report";

// The text format's %.4e is a display choice for a terminal and stays one; CSV and JSON are
// written with shortestRoundTrip(), because they exist to be parsed again.
std::string sci(double value) {
  char buffer[32];
  std::snprintf(buffer, sizeof(buffer), "%.4e", value);
  return buffer;
}

// JSON has no infinity and no NaN, so a non-finite value written as a bare `inf` produces a
// document that fails to parse -- in whatever consumes the report rather than here. `null` is
// what every parser accepts for "this is not a number".
std::string jsonNumber(double value) {
  return std::isfinite(value) ? shortestRoundTrip(value) : std::string("null");
}

// Half a mean free path is where scattered photons stop being a percent-level correction:
// there they are tens of percent of the uncollided value, and at one mean free path they
// exceed it. Below this the omission is smaller than the spread between published constants
// and does not earn a paragraph.
constexpr double kThickAirPathMfp = 0.5;

// The buildup caveat, for an exposure -- or a photon fluence -- whose air path is thick and
// whose buildup was left at 1.0. A caller who set a factor has made their own assumption about
// scatter and is not told again; a caller who set none has, by default, left the scattered
// photons out, and past this optical depth that is the largest thing the number is missing.
// Both metrics share the note because they share the uncollided assumption: it says nothing
// about what the photons DO at the point, only how many the model counted on the way. (The
// photon-strength metric never reaches the depth check: its tables carry no optical depth.)
void writeThickAirPathNote(std::ostream& out, Metric metric, double opticalDepth, double buildup) {
  if ((metric != Metric::Exposure && metric != Metric::Photon) || buildup != 1.0 ||
      !(opticalDepth > kThickAirPathMfp)) {
    return;
  }
  char depth[32];
  std::snprintf(depth, sizeof(depth), "%.2g", opticalDepth);
  out << "  ! the air path is " << depth
      << " mean free paths thick at the energies carrying this answer, and buildup\n"
      << "    is 1.0, so scattered photons are left out. Past about half a mean free path they\n"
      << "    add tens of percent to the uncollided value, and beyond one they exceed it. Set\n"
      << "    --buildup to include them.\n";
}

// RFC 4180 quoting. Labels today are nuclide names and energies and carry no commas, but that
// is a property of the data rather than of the format: one label with a comma in it shifts
// every column right of it by one, silently, in a file nobody re-reads by eye.
std::string csvField(const std::string& text) {
  if (text.find_first_of(",\"\r\n") == std::string::npos) {
    return text;
  }
  std::string out = "\"";
  for (const char c : text) {
    if (c == '"') {
      out += '"';
    }
    out += c;
  }
  out += '"';
  return out;
}

std::string percent(double fraction) {
  char buffer[32];
  const double value = fraction * 100.0;
  // A contributor at 0.011% renders as "0.0%" under one decimal place, which reads as
  // "nothing" when it is really "small but present". Two significant figures below 0.1%
  // keeps that distinction without widening the column for the common case.
  if (value > 0.0 && value < 0.1) {
    std::snprintf(buffer, sizeof(buffer), "%.2g%%", value);
  } else {
    std::snprintf(buffer, sizeof(buffer), "%.1f%%", value);
  }
  return buffer;
}

// Quotes and backslashes are the obvious cases; control characters are the ones that actually
// occur. A provenance string is a file path or a command line the user supplied, and a newline
// or a tab in one produces a report no JSON parser will read -- a failure that surfaces in
// whatever consumes the report rather than here, where it was caused.
std::string escapeJson(const std::string& text) {
  std::string out;
  out.reserve(text.size());
  for (const char c : text) {
    switch (c) {
      case '"':
        out += "\\\"";
        break;
      case '\\':
        out += "\\\\";
        break;
      case '\b':
        out += "\\b";
        break;
      case '\f':
        out += "\\f";
        break;
      case '\n':
        out += "\\n";
        break;
      case '\r':
        out += "\\r";
        break;
      case '\t':
        out += "\\t";
        break;
      default: {
        const auto byte = static_cast<unsigned char>(c);
        if (byte < 0x20) {
          char buffer[8];
          std::snprintf(buffer, sizeof(buffer), "\\u%04x", static_cast<unsigned>(byte));
          out += buffer;
        } else {
          // Everything at 0x20 and above passes through unchanged, including the bytes of a
          // multi-byte UTF-8 sequence: JSON strings are UTF-8, so escaping them would corrupt
          // a nuclide name or a path that is already correct.
          out += c;
        }
        break;
      }
    }
  }
  return out;
}

// The time a ranking describes: a point, or a window.
std::string whenOf(const Ranking& ranking) {
  if (ranking.domain == Domain::Interval) {
    return formatDuration(ranking.time) + " to " + formatDuration(ranking.timeEnd);
  }
  return formatDuration(ranking.time);
}

// The unit as it should be printed: a pack's own spelling when it has one, and the enum's name
// otherwise. Written once because a report that spelled the unit differently in the header and
// in the column would be describing two quantities.
std::string unitOf(const Ranking& ranking) {
  return ranking.unitLabel.empty() ? std::string(unitName(ranking.unit)) : ranking.unitLabel;
}

// The same question asked of a spec rather than of a ranking, for the writers that carry
// criteria instead of tables. Unit::PackDefined has no spelling of its own, so printing the
// enum's name here would put "pack-defined" in a column beside a real number.
std::string unitOf(const InterventionStudy& study) {
  return study.unitLabel.empty() ? std::string(unitName(study.unit)) : study.unitLabel;
}

std::string unitOf(const ResponseSpec& spec) {
  if (spec.metric == Metric::Pack && spec.pack != nullptr && spec.pack->pack != nullptr) {
    return spec.pack->pack->provenance().unit;
  }
  return std::string(unitName(spec.unit));
}

void writeTextHeader(std::ostream& out, const Ranking& ranking, const ReportContext& context) {
  out << "NuSIFT " << metricName(ranking.metric) << " ranking by "
      << aggregateName(ranking.aggregate) << '\n';
  out << "  t = " << whenOf(ranking) << "    total = " << sci(ranking.total) << ' '
      << unitOf(ranking) << '\n';
  if (!context.pack.empty()) {
    out << "  pack:  " << context.pack << '\n';
  }
  if (!context.geometry.empty()) {
    out << "  model: " << context.geometry << '\n';
  }
  if (!context.seedProvenance.empty()) {
    out << "  seed:  " << context.seedProvenance << '\n';
  }
  if (!context.storeLibrary.empty() || context.storeNuclideCount > 0) {
    out << "  store: ";
    if (!context.storeLibrary.empty()) {
      out << context.storeLibrary;
    }
    if (context.storeNuclideCount > 0) {
      out << " (" << context.storeNuclideCount << " nuclides";
      if (!context.storeCreatedUtc.empty()) {
        out << ", staged " << context.storeCreatedUtc;
      }
      out << ')';
    }
    out << '\n';
  }
  out << '\n';
}

void writeTextRows(std::ostream& out, const Ranking& ranking) {
  if (ranking.contributors.empty()) {
    out << "  (nothing contributes; the inventory has no " << metricName(ranking.metric) << ")\n";
    return;
  }

  // Wide enough for the heading too, or the header row shifts right of the values it labels.
  static constexpr std::string_view kLabelHeading = "contributor";
  std::size_t labelWidth = kLabelHeading.size();
  // And wide enough for the largest rank shown. Four holds every ranking that is only a top-N,
  // but a pinned row reaches past the cut and reports where it truly stands -- 12358th, in a
  // gamma-line table -- and a rank that overruns its column shifts the whole row right of it.
  int rankWidth = 4;
  for (const Contributor& c : ranking.contributors) {
    labelWidth = std::max(labelWidth, c.label.size());
    rankWidth = std::max(rankWidth, static_cast<int>(std::to_string(c.rank).size()));
  }

  out << std::right << std::setw(rankWidth) << "#" << "  " << std::left
      << std::setw(static_cast<int>(labelWidth)) << kLabelHeading << std::right << std::setw(13)
      << unitOf(ranking) << std::setw(9) << "frac" << std::setw(9) << "cum" << '\n';

  bool separated = false;
  for (const Contributor& c : ranking.contributors) {
    // The pinned rows are a tail, not a continuation of the ranking: their ranks jump, and run
    // together with the prefix above they would read as one list with numbers missing from it.
    // The heading is what says the rows below were asked for rather than reached.
    if (c.pinned && !separated) {
      out << "  pinned:\n";
      separated = true;
    }

    // A contributor with no rank contributes nothing at this time and holds no place in the
    // ordering; printing a 0 there would look like one. Its cumulative is meaningless for the
    // same reason -- there is no "everything down to it" -- while its own share, zero, is not.
    if (c.rank > 0) {
      out << std::right << std::setw(rankWidth) << c.rank;
    } else {
      out << std::right << std::setw(rankWidth) << "-";
    }
    out << "  " << std::left << std::setw(static_cast<int>(labelWidth)) << c.label << std::right
        << std::setw(13) << sci(c.value) << std::setw(9) << percent(c.fraction) << std::setw(9)
        << (c.rank > 0 ? percent(c.cumulativeFraction) : std::string("-"));
    if ((c.flags & kFlagUnmodeledContinuum) != 0) {
      out << "  !";
    }
    out << '\n';
  }
}

// The photon-coverage caveat: how much of the emitted photon energy sits outside the model, and
// which emitters carry it. Shared by the ranking footer and the attribution footer, because the
// two report the SAME figure and understate it by the same amount -- worded differently they
// would read as two separate reservations about one number.
void writeUnmodeledEnergyNote(std::ostream& out, Metric metric, double unmodeledEnergyFraction,
                              const std::vector<std::string>& named) {
  // The quantity the note understates, named the way the table above names it. The exposure
  // understatement follows mu_en/rho per missing energy; the photon-COUNT understatement
  // follows the photon number per missing energy. Both climb where the missing spectrum is
  // soft, so the "of that order, and larger when softer" statement holds for both.
  const char* singular = metric == Metric::Photon ? "photon count" : "exposure";
  const char* plural = metric == Metric::Photon ? "photon counts" : "exposures";

  // The magnitude first, because it is what decides whether the count matters at all. A
  // hundred flagged nuclides contributing 0.01% of the photon output is a footnote; three
  // contributing 30% is a reason not to trust the number above.
  //
  // Stated as the order of the understatement rather than as its size. The fraction is one of
  // emitted ENERGY, and the metric per unit of missing energy climbs steeply for soft spectra
  // -- mu_en/rho for exposure, photon number for a count -- so a continuum softer than the
  // lines, as bremsstrahlung usually is, understates more than its share of the energy says.
  if (unmodeledEnergyFraction > 0.0) {
    out << "  ! " << percent(unmodeledEnergyFraction)
        << " of the emitted photon energy is in spectra NuSIFT does not model. The " << singular
        << "\n"
        << "    understatement is of that order, and larger where the missing spectrum is\n"
        << "    softer than the lines, as bremsstrahlung usually is";
    // Terminated here unless the named list below continues the sentence. Left open, the line
    // runs into whatever is written next -- and with several --at times that is the blank line
    // writeRankings lays between rankings, which then disappears.
    if (named.empty()) {
      out << '\n';
    }
  }

  if (!named.empty()) {
    const std::size_t total = named.size();
    const bool one = total == 1;
    // Indented under the magnitude line when there is one, since it is the detail behind it.
    if (unmodeledEnergyFraction > 0.0) {
      out << " (" << total << " nuclide" << (one ? "" : "s") << ")";
    } else {
      out << "  ! " << total << " contributor" << (one ? "" : "s") << (one ? " carries" : " carry")
          << " photon energy NuSIFT does not model, so " << (one ? "its " : "their ")
          << (one ? singular : plural) << (one ? " is" : " are") << " understated";
    }

    // Naming every one of them is what a real evaluation turns this into: a full store flags
    // several hundred, overwhelmingly short-lived species that contribute nothing, and an
    // unbounded list buries the answer it was meant to annotate. The count is the signal; a
    // handful of names makes it concrete.
    constexpr std::size_t kMaxNamed = 8;
    const std::size_t show = std::min(total, kMaxNamed);
    out << ":\n    ";
    for (std::size_t i = 0; i < show; ++i) {
      out << (i == 0 ? "" : ", ") << named[i];
    }
    if (total > show) {
      out << ", and " << (total - show) << " more";
    }
    out << "\n    (see `nusift data info` for the store's photon coverage)\n";
  }
}

// What share of the inventory the pack could speak for. A sum of fractions over 97% of the
// activity is a screening index; the same number over 23% of it is arithmetic on whichever
// nuclides happened to be listed, and the two are indistinguishable without this line.
//
// Stated whenever it is not complete, and stated as a shortfall rather than as a coverage
// figure alone: "covers 23.4%" invites being read as a quality score, while naming what is
// missing says what the number is short by.
void writePackCoverageNote(std::ostream& out, const Ranking& ranking) {
  if (ranking.metric != Metric::Pack || ranking.packCoverage >= 0.9995) {
    return;
  }
  out << "  ! this pack carries a coefficient for " << percent(ranking.packCoverage)
      << " of the inventory, measured in\n"
      << "    the quantity its coefficients multiply. The rest is in nuclides it does not\n"
      << "    list, so the total above is a lower bound and the ranking is over what the\n"
      << "    pack covers rather than over what is present\n";
}

void writeTextFooter(std::ostream& out, const Ranking& ranking, const ReportContext& context) {
  // The honesty line. Without it a top-10 worth 40% and one worth 99% look identical.
  if (ranking.omittedCount > 0) {
    // Never let rounding claim the whole total while something is still omitted: "cover
    // 100.0% ... 1 further contributor omitted" reads as a contradiction even though it is
    // only a display artefact of a contributor at 0.0004%.
    std::string covered = percent(ranking.coveredFraction);
    if (covered == "100.0%") {
      covered = ">99.9%";
    }
    out << '\n'
        << "  shown rows cover " << covered << " of the total; " << ranking.omittedCount
        << " further contributor" << (ranking.omittedCount == 1 ? "" : "s") << " omitted\n";
  } else if (!ranking.contributors.empty()) {
    out << '\n' << "  shown rows cover the entire total\n";
  }

  // Only a pinned row can be rankless, and a dash in a column of numbers deserves one line of
  // explanation. It is also a real answer worth stating plainly: a pure beta emitter pinned in
  // an exposure ranking is not missing from the table, it contributes nothing to the metric.
  const bool anyRankless = std::any_of(ranking.contributors.begin(), ranking.contributors.end(),
                                       [](const Contributor& c) { return c.rank == 0; });
  if (anyRankless) {
    out << "  a pinned row with no rank contributes nothing to this " << metricName(ranking.metric)
        << " at this time\n";
  }

  writeUnmodeledEnergyNote(out, ranking.metric, ranking.unmodeledEnergyFraction,
                           context.unmodeledContinuum);

  writeThickAirPathNote(out, ranking.metric, ranking.meanOpticalDepth, ranking.buildup);

  writePackCoverageNote(out, ranking);
}

// The best place a contributor holds anywhere on the grid, or 0 if it never holds one at all.
int bestRankOf(const RankTrack& track) {
  int best = 0;
  for (const int rank : track.rank) {
    if (rank > 0 && (best == 0 || rank < best)) {
      best = rank;
    }
  }
  return best;
}

void writeForecastText(std::ostream& out, const std::vector<DominanceWindow>& windows,
                       const std::vector<RankTrack>& tracks, const ResponseTable& table,
                       const ReportContext& context) {
  out << "NuSIFT " << metricName(table.metric) << " forecast by " << aggregateName(table.aggregate)
      << '\n';
  if (table.timeCount() > 0) {
    out << "  " << formatDuration(table.times.front()) << " to "
        << formatDuration(table.times.back()) << ", " << table.timeCount() << " points\n";
  }
  if (!context.geometry.empty()) {
    out << "  model: " << context.geometry << '\n';
  }
  if (!context.seedProvenance.empty()) {
    out << "  seed:  " << context.seedProvenance << '\n';
  }
  out << '\n';

  if (windows.empty()) {
    out << "  (nothing contributes over this window)\n";
    return;
  }

  std::size_t labelWidth = 8;
  for (const DominanceWindow& window : windows) {
    labelWidth = std::max(labelWidth, window.label.size());
  }

  out << "  leads:\n";
  for (const DominanceWindow& window : windows) {
    out << "    " << std::left << std::setw(static_cast<int>(labelWidth)) << window.label
        << std::right << "  " << std::setw(10) << formatDuration(window.startSeconds) << " to "
        << std::setw(10) << formatDuration(window.endSeconds) << "   peak "
        << percent(window.peakFraction) << '\n';
  }

  if (tracks.empty()) {
    return;
  }

  // Every contributor that reaches the top at any point, with where it peaks. A nuclide that
  // matters only at thirty years belongs here; ordering by its value at any one time would
  // bury exactly the row a forecast exists to surface.
  std::size_t trackWidth = 8;
  for (const RankTrack& track : tracks) {
    trackWidth = std::max(trackWidth, track.label.size());
  }

  // The pinned tracks are separated for the same reason a pinned ranking row is: they are here
  // because someone asked after them, and listing them among contributors the forecast found
  // would say they came close when the whole point may be that they never did.
  for (const bool pinned : {false, true}) {
    const bool any = std::any_of(tracks.begin(), tracks.end(), [pinned](const RankTrack& track) {
      return track.pinned == pinned;
    });
    if (!any) {
      continue;
    }
    out << (pinned ? "\n  pinned:\n" : "\n  ever near the top:\n");
    for (const RankTrack& track : tracks) {
      if (track.pinned != pinned) {
        continue;
      }
      out << "    " << std::left << std::setw(static_cast<int>(trackWidth)) << track.label
          << std::right;

      // Best place held ANYWHERE on the grid, which need not be where the contributor peaks:
      // a share is measured against the total, and a shrinking total can lift a rank while the
      // share falls. Said as "anywhere" for that reason -- read as a property of the peak it
      // would be two different times reported as one.
      //
      // Not worth stating for a contributor the forecast surfaced, since it reached the top by
      // definition, but for a pinned one it is the number the reader came for: 3rd at best is a
      // different situation from 40th at best, and both peak somewhere.
      const int best = bestRankOf(track);
      if (best == 0) {
        out << "  contributes nothing over this grid\n";
        continue;
      }
      out << "  peaks at " << std::setw(10)
          << formatDuration(table.times[static_cast<std::size_t>(track.peakTimeIndex)]) << "  ("
          << percent(track.peakFraction) << " of the total)";
      if (pinned) {
        out << ", best rank anywhere " << best;
      }
      out << '\n';
    }
  }
}

void writeForecastJson(std::ostream& out, const std::vector<DominanceWindow>& windows,
                       const std::vector<RankTrack>& tracks, const ResponseTable& table) {
  out << "{\n";
  out << "  \"metric\": \"" << metricName(table.metric) << "\",\n";
  out << "  \"aggregate\": \"" << aggregateName(table.aggregate) << "\",\n";
  out << "  \"unit\": \"" << unitName(table.unit) << "\",\n";
  out << "  \"windows\": [\n";
  for (std::size_t i = 0; i < windows.size(); ++i) {
    const DominanceWindow& window = windows[i];
    out << "    {\"label\": \"" << escapeJson(window.label) << "\", \"key\": " << window.id.key
        << ", \"start_s\": " << jsonNumber(window.startSeconds)
        << ", \"end_s\": " << jsonNumber(window.endSeconds)
        << ", \"peak_fraction\": " << jsonNumber(window.peakFraction) << "}"
        << (i + 1 < windows.size() ? ",\n" : "\n");
  }
  out << "  ],\n  \"tracks\": [\n";
  for (std::size_t i = 0; i < tracks.size(); ++i) {
    const RankTrack& track = tracks[i];
    out << "    {\"label\": \"" << escapeJson(track.label) << "\", \"key\": " << track.id.key
        << ", \"peak_fraction\": " << jsonNumber(track.peakFraction) << ", \"peak_time_s\": "
        << jsonNumber(table.times[static_cast<std::size_t>(track.peakTimeIndex)])
        << ", \"best_rank\": " << bestRankOf(track)
        << ", \"pinned\": " << (track.pinned ? "true" : "false") << "}"
        << (i + 1 < tracks.size() ? ",\n" : "\n");
  }
  out << "  ]\n}\n";
}

void writeForecastCsv(std::ostream& out, const std::vector<DominanceWindow>& windows) {
  out << "label,key,start_s,end_s,peak_fraction\n";
  for (const DominanceWindow& window : windows) {
    out << csvField(window.label) << ',' << window.id.key << ','
        << shortestRoundTrip(window.startSeconds) << ',' << shortestRoundTrip(window.endSeconds)
        << ',' << shortestRoundTrip(window.peakFraction) << '\n';
  }
}

void writeCsvRows(std::ostream& out, const Ranking& ranking, bool withHeader) {
  if (withHeader) {
    out << "time_s,time_end_s,rank,contributor,key,value,unit,fraction,cumulative_fraction,"
           "flags,pinned\n";
  }
  for (const Contributor& c : ranking.contributors) {
    out << shortestRoundTrip(ranking.time) << ',';
    if (ranking.domain == Domain::Interval) {
      out << shortestRoundTrip(ranking.timeEnd);
    }
    // The unit is not quoted: it comes from a closed enum of spellings that contain no comma,
    // so unlike a label it cannot acquire one.
    //
    // `pinned` is last so that adding it did not renumber the columns anyone already reads by
    // position, and it is here at all because without it a loaded table cannot tell a row that
    // placed from one that was fetched from below the cut -- which is the difference between a
    // top-N and a top-N plus an aside.
    out << ',' << c.rank << ',' << csvField(c.label) << ',' << c.id.key << ','
        << shortestRoundTrip(c.value) << ',' << csvField(unitOf(ranking)) << ','
        << shortestRoundTrip(c.fraction) << ',' << shortestRoundTrip(c.cumulativeFraction) << ','
        << c.flags << ',' << (c.pinned ? 1 : 0) << '\n';
  }
}

void writeJsonRanking(std::ostream& out, const Ranking& ranking, const ReportContext& context,
                      int indent) {
  const std::string pad(static_cast<std::size_t>(indent), ' ');
  out << pad << "{\n";
  out << pad << "  \"metric\": \"" << metricName(ranking.metric) << "\",\n";
  out << pad << "  \"aggregate\": \"" << aggregateName(ranking.aggregate) << "\",\n";
  out << pad << "  \"unit\": \"" << escapeJson(unitOf(ranking)) << "\",\n";
  out << pad << "  \"time_s\": " << jsonNumber(ranking.time) << ",\n";
  if (ranking.domain == Domain::Interval) {
    out << pad << "  \"time_end_s\": " << jsonNumber(ranking.timeEnd) << ",\n";
  }
  out << pad << "  \"total\": " << jsonNumber(ranking.total) << ",\n";
  out << pad << "  \"covered_fraction\": " << jsonNumber(ranking.coveredFraction) << ",\n";
  out << pad << "  \"omitted_count\": " << ranking.omittedCount << ",\n";
  // The two caveats the text footer states, as numbers a script can act on: how much of the
  // photon energy the model does not carry, and how thick the air path was left uncorrected.
  // Mean optical depth reads as zero for a photon-strength answer, which is the truth of it:
  // that number used no path at all.
  if (ranking.metric == Metric::Exposure || ranking.metric == Metric::Photon) {
    out << pad << "  \"unmodeled_energy_fraction\": " << jsonNumber(ranking.unmodeledEnergyFraction)
        << ",\n";
    out << pad << "  \"mean_optical_depth\": " << jsonNumber(ranking.meanOpticalDepth) << ",\n";
    out << pad << "  \"buildup\": " << jsonNumber(ranking.buildup) << ",\n";
  }
  // The pack's caveat in the same form: a consumer that reads the total without this cannot
  // tell an index over the whole inventory from one over a quarter of it.
  if (ranking.metric == Metric::Pack) {
    out << pad << "  \"pack_coverage\": " << jsonNumber(ranking.packCoverage) << ",\n";
    if (!context.pack.empty()) {
      out << pad << "  \"pack\": \"" << escapeJson(context.pack) << "\",\n";
    }
  }
  if (!context.seedProvenance.empty()) {
    out << pad << "  \"seed\": \"" << escapeJson(context.seedProvenance) << "\",\n";
  }
  if (!context.storeLibrary.empty()) {
    out << pad << "  \"library\": \"" << escapeJson(context.storeLibrary) << "\",\n";
  }
  if (!context.geometry.empty()) {
    out << pad << "  \"model\": \"" << escapeJson(context.geometry) << "\",\n";
  }
  out << pad << "  \"contributors\": [\n";
  for (std::size_t i = 0; i < ranking.contributors.size(); ++i) {
    const Contributor& c = ranking.contributors[i];
    out << pad << "    {\"rank\": " << c.rank << ", \"label\": \"" << escapeJson(c.label)
        << "\", \"key\": " << c.id.key << ", \"value\": " << jsonNumber(c.value)
        << ", \"fraction\": " << jsonNumber(c.fraction)
        << ", \"cumulative_fraction\": " << jsonNumber(c.cumulativeFraction)
        << ", \"flags\": " << c.flags << ", \"pinned\": " << (c.pinned ? "true" : "false") << "}";
    out << (i + 1 < ranking.contributors.size() ? ",\n" : "\n");
  }
  out << pad << "  ]\n";
  out << pad << "}";
}

// --- located events ----------------------------------------------------------

// The provenance of an instant: the grid interval that observed it, how tightly it was placed,
// and whether real evaluations did the placing. Printed on every event because without it a
// crossing time is indistinguishable from a grid artefact.
std::string bracketNote(const TrajectoryEvent& event) {
  std::string note = "bracket " + formatDuration(event.bracketStartSeconds) + " to " +
                     formatDuration(event.bracketEndSeconds) + ", within " +
                     formatDuration(event.locatedToSeconds);
  note += event.refined ? "" : " (interpolated)";
  if (!event.converged) {
    note += " (refinement did not reach its tolerance)";
  }
  return note;
}

bool isTurn(const TrajectoryEvent& event) {
  return event.kind == EventKind::Maximum || event.kind == EventKind::Minimum;
}

// Said once, wherever a search came back with less than the reader might expect. The grid is
// what decides which events exist to be found, and a reader who does not know that will read
// an empty list as "it never happens" rather than "the sampling did not resolve it".
void writeGridNote(std::ostream& out) {
  out << "\n  Events are found only where consecutive samples straddle them. An excursion that\n"
         "  rises and falls back between two samples leaves no sign change and is not found; a\n"
         "  denser grid is what resolves one.\n";
}

void writeEventsText(std::ostream& out, const EventReport& report, const ReportContext& context) {
  out << "NuSIFT " << report.metric << " events on " << report.curve << '\n';
  if (report.gridPoints > 0) {
    out << "  " << formatDuration(report.gridStartSeconds) << " to "
        << formatDuration(report.gridEndSeconds) << ", " << report.gridPoints << " points\n";
  }
  if (!context.geometry.empty()) {
    out << "  model: " << context.geometry << '\n';
  }
  if (!context.seedProvenance.empty()) {
    out << "  seed:  " << context.seedProvenance << '\n';
  }
  if (report.hasLevel) {
    out << "  level: " << sci(report.level);
    if (!report.unit.empty()) {
      out << ' ' << report.unit;
    }
    out << '\n';
  }
  out << '\n';

  const bool anyCrossing = std::any_of(report.events.begin(), report.events.end(),
                                       [](const TrajectoryEvent& event) { return !isTurn(event); });
  const bool anyTurn = std::any_of(report.events.begin(), report.events.end(), isTurn);

  if (report.hasLevel) {
    if (anyCrossing) {
      out << "  crossings:\n";
      for (const TrajectoryEvent& event : report.events) {
        if (isTurn(event)) {
          continue;
        }
        out << "    " << std::left << std::setw(8) << eventKindName(event.kind) << std::right
            << "  at " << std::setw(10) << formatDuration(event.timeSeconds) << "   "
            << bracketNote(event) << '\n';
      }
    } else {
      out << "  crossings: none -- the grid never observed this curve cross the level\n";
    }

    if (!report.windows.empty()) {
      const char* side = report.windowsBelowLevel ? "below" : "above";
      out << "\n  " << side << " the level:\n";
      for (const LevelWindow& window : report.windows) {
        out << "    " << std::setw(10) << formatDuration(window.startSeconds) << " to "
            << std::setw(10) << formatDuration(window.endSeconds);
        // An edge the grid never observed is a bound, not a crossing, and the difference is
        // the whole reason the flags exist.
        if (!window.entryObserved && !window.exitObserved) {
          out << "   (open at both ends: never observed to rise above or fall below)";
        } else if (!window.entryObserved) {
          out << "   (already " << side << " when the grid started)";
        } else if (!window.exitObserved) {
          out << "   (still " << side << " when the grid ended)";
        }
        out << '\n';
      }
    }
  }

  if (anyTurn) {
    out << (report.hasLevel ? "\n  turns:\n" : "  turns:\n");
    for (const TrajectoryEvent& event : report.events) {
      if (!isTurn(event)) {
        continue;
      }
      out << "    " << std::left << std::setw(8) << eventKindName(event.kind) << std::right
          << "  at " << std::setw(10) << formatDuration(event.timeSeconds) << "   value "
          << sci(event.value);
      if (!report.unit.empty()) {
        out << ' ' << report.unit;
      }
      out << "   " << bracketNote(event) << '\n';
    }
  } else if (!report.hasLevel) {
    out << "  turns: none -- the grid resolved no maximum or minimum on this curve\n";
  }

  writeGridNote(out);
}

void writeEventsCsv(std::ostream& out, const EventReport& report) {
  // Events and windows in one table under a `kind` column: two tables would not be a CSV, and
  // a reader that got only the crossings would lose which grid edges were never observed.
  out << "kind,time_s,end_s,value,bracket_start_s,bracket_end_s,located_to_s,refined,converged,"
         "entry_observed,exit_observed,side\n";
  for (const TrajectoryEvent& event : report.events) {
    out << eventKindName(event.kind) << ',' << shortestRoundTrip(event.timeSeconds) << ",,"
        << shortestRoundTrip(event.value) << ',' << shortestRoundTrip(event.bracketStartSeconds)
        << ',' << shortestRoundTrip(event.bracketEndSeconds) << ','
        << shortestRoundTrip(event.locatedToSeconds) << ',' << (event.refined ? "true" : "false")
        << ',' << (event.converged ? "true" : "false") << ",,,\n";
  }
  for (const LevelWindow& window : report.windows) {
    out << "window," << shortestRoundTrip(window.startSeconds) << ','
        << shortestRoundTrip(window.endSeconds) << ",,,,,,,"
        << (window.entryObserved ? "true" : "false") << ','
        << (window.exitObserved ? "true" : "false") << ','
        << (report.windowsBelowLevel ? "below" : "above") << '\n';
  }
}

void writeEventsJson(std::ostream& out, const EventReport& report) {
  out << "{\n";
  out << "  \"metric\": \"" << escapeJson(report.metric) << "\",\n";
  out << "  \"curve\": \"" << escapeJson(report.curve) << "\",\n";
  out << "  \"unit\": " << (report.unit.empty() ? "null" : "\"" + escapeJson(report.unit) + "\"")
      << ",\n";
  out << "  \"level\": " << (report.hasLevel ? jsonNumber(report.level) : "null") << ",\n";
  out << "  \"grid\": {\"start_s\": " << jsonNumber(report.gridStartSeconds)
      << ", \"end_s\": " << jsonNumber(report.gridEndSeconds)
      << ", \"points\": " << report.gridPoints << "},\n";

  out << "  \"events\": [\n";
  for (std::size_t i = 0; i < report.events.size(); ++i) {
    const TrajectoryEvent& event = report.events[i];
    out << "    {\"kind\": \"" << eventKindName(event.kind)
        << "\", \"time_s\": " << jsonNumber(event.timeSeconds)
        << ", \"value\": " << jsonNumber(event.value)
        << ", \"bracket_start_s\": " << jsonNumber(event.bracketStartSeconds)
        << ", \"bracket_end_s\": " << jsonNumber(event.bracketEndSeconds)
        << ", \"located_to_s\": " << jsonNumber(event.locatedToSeconds)
        << ", \"refined\": " << (event.refined ? "true" : "false")
        << ", \"converged\": " << (event.converged ? "true" : "false") << "}"
        << (i + 1 < report.events.size() ? ",\n" : "\n");
  }
  out << "  ],\n  \"windows_side\": \"" << (report.windowsBelowLevel ? "below" : "above")
      << "\",\n  \"windows\": [\n";
  for (std::size_t i = 0; i < report.windows.size(); ++i) {
    const LevelWindow& window = report.windows[i];
    out << "    {\"start_s\": " << jsonNumber(window.startSeconds)
        << ", \"end_s\": " << jsonNumber(window.endSeconds)
        << ", \"entry_observed\": " << (window.entryObserved ? "true" : "false")
        << ", \"exit_observed\": " << (window.exitObserved ? "true" : "false") << "}"
        << (i + 1 < report.windows.size() ? ",\n" : "\n");
  }
  out << "  ]\n}\n";
}

// --- maximum allowable scale --------------------------------------------------

std::string limitingSummary(const AllowableScale& at) {
  std::string text;
  for (std::size_t i = 0; i < at.limiting.size(); ++i) {
    if (i > 0) {
      text += ", ";
    }
    text += at.limiting[i].label + " " + percent(at.limiting[i].fraction);
  }
  return text;
}

void writeAllowableText(std::ostream& out, const std::vector<AllowableScale>& scaled,
                        std::span<const Criterion> criteria, const ReportContext& context) {
  out << "NuSIFT maximum allowable scale\n";
  if (!scaled.empty()) {
    out << "  " << formatDuration(scaled.front().timeSeconds) << " to "
        << formatDuration(scaled.back().timeSeconds) << ", " << scaled.size() << " points\n";
  }
  if (!context.geometry.empty()) {
    out << "  model: " << context.geometry << '\n';
  }
  if (!context.seedProvenance.empty()) {
    out << "  seed:  " << context.seedProvenance << '\n';
  }

  std::size_t nameWidth = 8;
  for (const Criterion& criterion : criteria) {
    nameWidth = std::max(nameWidth, criterion.name.size());
  }

  out << "\n  criteria:\n";
  for (const Criterion& criterion : criteria) {
    out << "    " << std::left << std::setw(static_cast<int>(nameWidth)) << criterion.name
        << std::right << "  " << sci(criterion.limit) << ' ' << unitOf(criterion.spec) << '\n';
  }

  out << "\n  " << std::setw(12) << "time" << "  " << std::setw(11) << "scale" << "  " << std::left
      << std::setw(static_cast<int>(nameWidth)) << "binding" << std::right << "  driven by\n";
  for (const AllowableScale& at : scaled) {
    out << "  " << std::setw(12) << formatDuration(at.timeSeconds) << "  ";
    if (!at.bounded) {
      // Not a very large number. Nothing constrains the inventory here at all, and the two
      // statements are different.
      out << std::setw(11) << "unbounded" << "  " << std::left
          << std::setw(static_cast<int>(nameWidth)) << "--" << std::right << '\n';
      continue;
    }
    out << std::setw(11) << sci(at.scale) << "  " << std::left
        << std::setw(static_cast<int>(nameWidth))
        << criteria[static_cast<std::size_t>(at.bindingIndex)].name << std::right << "  "
        << limitingSummary(at) << '\n';
  }

  const bool anyUnbounded = std::any_of(scaled.begin(), scaled.end(),
                                        [](const AllowableScale& at) { return !at.bounded; });
  if (anyUnbounded) {
    out << "\n  \"unbounded\" means no criterion has anything to constrain at that time -- not\n"
           "  that a very large multiple is permitted.\n";
  }
}

void writeAllowableCsv(std::ostream& out, const std::vector<AllowableScale>& scaled,
                       std::span<const Criterion> criteria) {
  // Long format, one row per time per criterion: every criterion's headroom is in the table
  // rather than only the binding one's, which is what makes "how close was the runner-up"
  // answerable from the file.
  out << "time_s,criterion,unit,response,limit,fraction,criterion_scale,binding,unbounded,"
         "allowed_scale\n";
  for (const AllowableScale& at : scaled) {
    for (std::size_t q = 0; q < at.criteria.size(); ++q) {
      const CriterionHeadroom& headroom = at.criteria[q];
      out << shortestRoundTrip(at.timeSeconds) << ',' << csvField(headroom.name) << ','
          << csvField(unitOf(criteria[q].spec)) << ',' << shortestRoundTrip(headroom.response)
          << ',' << shortestRoundTrip(headroom.limit) << ',';
      if (headroom.unbounded) {
        out << ",,";
      } else {
        out << shortestRoundTrip(headroom.fraction) << ',' << shortestRoundTrip(headroom.scale)
            << ',';
      }
      out << (headroom.binding ? "true" : "false") << ',' << (headroom.unbounded ? "true" : "false")
          << ',';
      if (at.bounded) {
        out << shortestRoundTrip(at.scale);
      }
      out << '\n';
    }
  }
}

void writeAllowableJson(std::ostream& out, const std::vector<AllowableScale>& scaled,
                        std::span<const Criterion> criteria) {
  out << "{\n  \"criteria\": [\n";
  for (std::size_t q = 0; q < criteria.size(); ++q) {
    out << "    {\"name\": \"" << escapeJson(criteria[q].name) << "\", \"unit\": \""
        << escapeJson(unitOf(criteria[q].spec))
        << "\", \"limit\": " << jsonNumber(criteria[q].limit) << "}"
        << (q + 1 < criteria.size() ? ",\n" : "\n");
  }
  out << "  ],\n  \"times\": [\n";
  for (std::size_t k = 0; k < scaled.size(); ++k) {
    const AllowableScale& at = scaled[k];
    out << "    {\"time_s\": " << jsonNumber(at.timeSeconds) << ",\n";
    // null rather than a number: the only encoding a parser cannot mistake for a bound of zero.
    out << "     \"scale\": " << (at.bounded ? jsonNumber(at.scale) : "null") << ",\n";
    out << "     \"binding\": "
        << (at.bounded
                ? "\"" + escapeJson(criteria[static_cast<std::size_t>(at.bindingIndex)].name) + "\""
                : "null")
        << ",\n";
    out << "     \"criteria\": [";
    for (std::size_t q = 0; q < at.criteria.size(); ++q) {
      const CriterionHeadroom& headroom = at.criteria[q];
      out << (q == 0 ? "\n" : ",\n") << "       {\"name\": \"" << escapeJson(headroom.name)
          << "\", \"response\": " << jsonNumber(headroom.response)
          << ", \"fraction\": " << (headroom.unbounded ? "null" : jsonNumber(headroom.fraction))
          << ", \"scale\": " << (headroom.unbounded ? "null" : jsonNumber(headroom.scale))
          << ", \"unbounded\": " << (headroom.unbounded ? "true" : "false") << "}";
    }
    out << (at.criteria.empty() ? "" : "\n     ") << "],\n";
    out << "     \"limiting\": [";
    for (std::size_t i = 0; i < at.limiting.size(); ++i) {
      out << (i == 0 ? "\n" : ",\n") << "       {\"label\": \"" << escapeJson(at.limiting[i].label)
          << "\", \"key\": " << at.limiting[i].id.key
          << ", \"fraction\": " << jsonNumber(at.limiting[i].fraction) << "}";
    }
    out << (at.limiting.empty() ? "" : "\n     ") << "]}" << (k + 1 < scaled.size() ? ",\n" : "\n");
  }
  out << "  ]\n}\n";
}

// --- counterfactual interventions ----------------------------------------------

std::string drivenBy(const InterventionEffect& effect, std::size_t most) {
  std::string text;
  for (std::size_t i = 0; i < effect.contributors.size() && i < most; ++i) {
    if (i > 0) {
      text += ", ";
    }
    text += effect.contributors[i].label + " " + percent(effect.contributors[i].fraction);
  }
  return text;
}

void writeInterventionsText(std::ostream& out, const InterventionStudy& study,
                            const ReportContext& context) {
  out << "NuSIFT " << metricName(study.metric) << " interventions\n";
  out << "  remove at " << formatDuration(study.interventionTimeSeconds) << ", response at "
      << formatDuration(study.responseTimeSeconds) << '\n';
  if (!context.geometry.empty()) {
    out << "  model: " << context.geometry << '\n';
  }
  if (!study.seedProvenance.empty()) {
    out << "  seed:  " << study.seedProvenance << '\n';
  }
  out << "\n  baseline: " << sci(study.baseline) << ' ' << unitOf(study)
      << " with nothing removed\n\n";

  std::size_t nameWidth = 12;
  for (const InterventionEffect& effect : study.effects) {
    nameWidth = std::max(nameWidth, effect.name.size());
  }

  out << "  " << std::left << std::setw(static_cast<int>(nameWidth)) << "intervention" << std::right
      << "  " << std::setw(12) << "after" << "  " << std::setw(12) << "removed"
      << "  " << std::setw(8) << "of base" << "  driven by\n";
  for (const InterventionEffect& effect : study.effects) {
    out << "  " << std::left << std::setw(static_cast<int>(nameWidth)) << effect.name << std::right
        << "  " << std::setw(12) << sci(effect.response) << "  " << std::setw(12)
        << sci(effect.removed) << "  " << std::setw(8) << percent(effect.removedFraction) << "  "
        << drivenBy(effect, 3) << '\n';
  }

  // Alternatives, not a sequence. Summing two rows would describe a schedule nobody computed,
  // and the numbers look perfectly addable, so the report has to say it.
  out << "\n  Each row is an alternative measured against the same baseline, not a step in a\n"
         "  sequence -- two rows do not add. A schedule of removals on different dates is a\n"
         "  separate run per date.\n";
}

void writeInterventionsCsv(std::ostream& out, const InterventionStudy& study) {
  // Long format, one row per intervention per nuclide, with the intervention's own totals
  // repeated. A summary-only table would drop what the benefit was actually made of.
  out << "intervention,response,removed,removed_fraction,unit,nuclide,key,atoms_removed,value,"
         "fraction\n";
  for (const InterventionEffect& effect : study.effects) {
    if (effect.contributors.empty()) {
      out << csvField(effect.name) << ',' << shortestRoundTrip(effect.response) << ','
          << shortestRoundTrip(effect.removed) << ',' << shortestRoundTrip(effect.removedFraction)
          << ',' << unitOf(study) << ",,,,,\n";
      continue;
    }
    for (const RemovedContributor& contributor : effect.contributors) {
      out << csvField(effect.name) << ',' << shortestRoundTrip(effect.response) << ','
          << shortestRoundTrip(effect.removed) << ',' << shortestRoundTrip(effect.removedFraction)
          << ',' << unitOf(study) << ',' << csvField(contributor.label) << ',' << contributor.key
          << ',' << shortestRoundTrip(contributor.atomsRemoved) << ','
          << shortestRoundTrip(contributor.value) << ',' << shortestRoundTrip(contributor.fraction)
          << '\n';
    }
  }
}

void writeInterventionsJson(std::ostream& out, const InterventionStudy& study) {
  out << "{\n";
  out << "  \"metric\": \"" << metricName(study.metric) << "\",\n";
  out << "  \"unit\": \"" << unitOf(study) << "\",\n";
  out << "  \"intervention_time_s\": " << jsonNumber(study.interventionTimeSeconds) << ",\n";
  out << "  \"response_time_s\": " << jsonNumber(study.responseTimeSeconds) << ",\n";
  out << "  \"baseline\": " << jsonNumber(study.baseline) << ",\n";
  out << "  \"interventions\": [\n";
  for (std::size_t q = 0; q < study.effects.size(); ++q) {
    const InterventionEffect& effect = study.effects[q];
    out << "    {\"name\": \"" << escapeJson(effect.name)
        << "\", \"response\": " << jsonNumber(effect.response)
        << ", \"removed\": " << jsonNumber(effect.removed)
        << ", \"removed_fraction\": " << jsonNumber(effect.removedFraction)
        << ", \"contributors\": [";
    for (std::size_t i = 0; i < effect.contributors.size(); ++i) {
      const RemovedContributor& contributor = effect.contributors[i];
      out << (i == 0 ? "\n" : ",\n") << "      {\"label\": \"" << escapeJson(contributor.label)
          << "\", \"key\": " << contributor.key
          << ", \"atoms_removed\": " << jsonNumber(contributor.atomsRemoved)
          << ", \"value\": " << jsonNumber(contributor.value)
          << ", \"fraction\": " << jsonNumber(contributor.fraction) << "}";
    }
    out << (effect.contributors.empty() ? "" : "\n    ") << "]}"
        << (q + 1 < study.effects.size() ? ",\n" : "\n");
  }
  out << "  ]\n}\n";
}

}  // namespace

bool parseReportFormat(std::string_view text, ReportFormat& out) {
  if (text == "text") {
    out = ReportFormat::Text;
    return true;
  }
  if (text == "csv") {
    out = ReportFormat::Csv;
    return true;
  }
  if (text == "json") {
    out = ReportFormat::Json;
    return true;
  }
  return false;
}

void writeForecast(std::ostream& out, const std::vector<DominanceWindow>& windows,
                   const std::vector<RankTrack>& tracks, const ResponseTable& table,
                   const ReportContext& context, ReportFormat format) {
  switch (format) {
    case ReportFormat::Text: {
      writeForecastText(out, windows, tracks, table, context);
      // Judged at the thickest point on the grid: the spectrum hardens and softens as the
      // leaders turn over, and a forecast that is fine at one end and not at the other should
      // say so rather than average the caveat away.
      double thickest = 0.0;
      for (const double depth : table.meanOpticalDepth) {
        thickest = std::max(thickest, depth);
      }
      writeThickAirPathNote(out, table.metric, thickest, table.geometry.buildup);
      break;
    }
    case ReportFormat::Csv:
      writeForecastCsv(out, windows);
      break;
    case ReportFormat::Json:
      writeForecastJson(out, windows, tracks, table);
      break;
  }
}

void writeAttribution(std::ostream& out, const SeedAttribution& a, const ReportContext& context,
                      ReportFormat format) {
  if (format == ReportFormat::Csv) {
    // shortestRoundTrip, not sci(): this table exists to be loaded again, and %.4e is a
    // terminal's five digits. Xe-140 and Cs-140 importances differ in the eleventh, so at
    // four they serialize identically -- and `value` no longer equals seed_atoms x importance
    // in the reloaded frame. See the note on sci() at the top of this file.
    //
    // The time and the unit ride on every row, as they do in writeCsvRows: a share is a
    // number of becquerel or of R/h at one instant, and a table stating neither cannot be
    // interpreted at all once it has left the terminal it was printed in.
    out << "time_s,rank,seed,key,seed_atoms,importance,value,unit,fraction,cumulative,pinned\n";
    for (const SeedShare& s : a.shares) {
      out << shortestRoundTrip(a.time) << ',' << s.rank << ',' << csvField(s.label) << ',' << s.key
          << ',' << shortestRoundTrip(s.seedAtoms) << ',' << shortestRoundTrip(s.importance) << ','
          << shortestRoundTrip(s.value) << ',' << unitName(a.unit) << ','
          << shortestRoundTrip(s.fraction) << ',' << shortestRoundTrip(s.cumulativeFraction) << ','
          << (s.pinned ? 1 : 0) << '\n';
    }
    return;
  }
  if (format == ReportFormat::Json) {
    out << "{\"metric\":\"" << metricName(a.metric) << "\",\"attribution\":\"seed\"";
    out << ",\"unit\":\"" << unitName(a.unit) << "\"";
    out << ",\"time_s\":" << jsonNumber(a.time);
    out << ",\"total\":" << jsonNumber(a.total);
    out << ",\"covered_fraction\":" << jsonNumber(a.coveredFraction);
    out << ",\"omitted\":" << a.omittedCount;
    // The same keys writeJsonRanking emits, spelled the same way. A consumer reading an
    // exposure -- or a photon figure -- out of one document and out of the other should not
    // have to know which command produced it to find out how far the model was stretched to
    // get it.
    if (a.metric == Metric::Exposure || a.metric == Metric::Photon) {
      out << ",\"unmodeled_energy_fraction\":" << jsonNumber(a.unmodeledEnergyFraction);
      out << ",\"mean_optical_depth\":" << jsonNumber(a.meanOpticalDepth);
      out << ",\"buildup\":" << jsonNumber(a.buildup);
    }
    if (!context.seedProvenance.empty()) {
      out << ",\"seed_provenance\":\"" << escapeJson(context.seedProvenance) << "\"";
    }
    if (!context.storeLibrary.empty()) {
      out << ",\"library\":\"" << escapeJson(context.storeLibrary) << "\"";
    }
    if (!context.geometry.empty()) {
      out << ",\"model\":\"" << escapeJson(context.geometry) << "\"";
    }
    out << ",\"shares\":[";
    for (std::size_t i = 0; i < a.shares.size(); ++i) {
      const SeedShare& s = a.shares[i];
      if (i > 0) {
        out << ',';
      }
      out << "{\"rank\":" << s.rank << ",\"seed\":\"" << escapeJson(s.label) << "\""
          << ",\"key\":" << s.key << ",\"seed_atoms\":" << jsonNumber(s.seedAtoms)
          << ",\"importance\":" << jsonNumber(s.importance) << ",\"value\":" << jsonNumber(s.value)
          << ",\"fraction\":" << jsonNumber(s.fraction)
          << ",\"cumulative_fraction\":" << jsonNumber(s.cumulativeFraction)
          << ",\"pinned\":" << (s.pinned ? "true" : "false") << '}';
    }
    out << "]}\n";
    return;
  }

  out << "NuSIFT " << metricName(a.metric) << " attributed to the seed\n";
  out << "  t = " << formatDuration(a.time) << "    total = " << sci(a.total) << ' '
      << unitName(a.unit) << '\n';
  if (!context.geometry.empty()) {
    out << "  model: " << context.geometry << '\n';
  }
  if (!context.seedProvenance.empty()) {
    out << "  seed:  " << context.seedProvenance << '\n';
  }
  out << '\n';

  if (a.shares.empty()) {
    out << "  (nothing was seeded that contributes)\n";
    return;
  }

  static constexpr std::string_view kHeading = "seed";
  std::size_t labelWidth = kHeading.size();
  int rankWidth = 4;
  for (const SeedShare& s : a.shares) {
    labelWidth = std::max(labelWidth, s.label.size());
    rankWidth = std::max(rankWidth, static_cast<int>(std::to_string(s.rank).size()));
  }

  out << std::right << std::setw(rankWidth) << "#" << "  " << std::left
      << std::setw(static_cast<int>(labelWidth)) << kHeading << std::right << std::setw(13)
      << "seed atoms" << std::setw(13) << unitName(a.unit) << std::setw(9) << "frac" << std::setw(9)
      << "cum" << '\n';

  bool separated = false;
  for (const SeedShare& s : a.shares) {
    if (s.pinned && !separated) {
      out << "  pinned:\n";
      separated = true;
    }
    if (s.rank > 0) {
      out << std::right << std::setw(rankWidth) << s.rank;
    } else {
      out << std::right << std::setw(rankWidth) << "-";
    }
    out << "  " << std::left << std::setw(static_cast<int>(labelWidth)) << s.label << std::right
        << std::setw(13) << sci(s.seedAtoms) << std::setw(13) << sci(s.value) << std::setw(9)
        << percent(s.fraction) << std::setw(9)
        << (s.rank > 0 ? percent(s.cumulativeFraction) : std::string("-")) << '\n';
  }

  if (a.omittedCount > 0) {
    std::string covered = percent(a.coveredFraction);
    if (covered == "100.0%") {
      covered = ">99.9%";
    }
    out << "\n  shown rows cover " << covered << " of the total; " << a.omittedCount
        << " further seed" << (a.omittedCount == 1 ? "" : "s") << " omitted\n";
  } else if (!a.shares.empty()) {
    out << "\n  shown rows cover the entire total\n";
  }

  // Only an inert pinned row is rankless here, and a dash in a column of numbers earns one
  // line of explanation. It is also a real answer: a stable seed pinned into an attribution is
  // not missing from the table, it contributes nothing to this metric.
  const bool anyRankless =
      std::any_of(a.shares.begin(), a.shares.end(), [](const SeedShare& s) { return s.rank == 0; });
  if (anyRankless) {
    out << "  a pinned seed with no rank contributes nothing to this " << metricName(a.metric)
        << " at this time\n";
  }

  // The caveats the FORWARD ranking of this same number prints. An exposure -- or a photon
  // answer -- does not become better characterised by being decomposed, and a reader who ran
  // `attribute` instead of `rank` has asked a different question about an identically
  // uncertain figure.
  writeUnmodeledEnergyNote(out, a.metric, a.unmodeledEnergyFraction, a.unmodeledContinuum);
  writeThickAirPathNote(out, a.metric, a.meanOpticalDepth, a.buildup);
}

void writeRanking(std::ostream& out, const Ranking& ranking, const ReportContext& context,
                  ReportFormat format) {
  switch (format) {
    case ReportFormat::Text:
      writeTextHeader(out, ranking, context);
      writeTextRows(out, ranking);
      writeTextFooter(out, ranking, context);
      break;
    case ReportFormat::Csv:
      writeCsvRows(out, ranking, /*withHeader=*/true);
      break;
    case ReportFormat::Json:
      writeJsonRanking(out, ranking, context, 0);
      out << '\n';
      break;
  }
}

void writeRankings(std::ostream& out, const std::vector<Ranking>& rankings,
                   const ReportContext& context, ReportFormat format) {
  writeRankings(out, rankings, std::vector<ReportContext>(rankings.size(), context), format);
}

void writeRankings(std::ostream& out, const std::vector<Ranking>& rankings,
                   const std::vector<ReportContext>& contexts, ReportFormat format) {
  if (rankings.empty()) {
    return;
  }
  if (contexts.size() != rankings.size()) {
    throw NusiftError(tagged(kModule, "a report needs one context per ranking"));
  }
  switch (format) {
    case ReportFormat::Text:
      for (std::size_t k = 0; k < rankings.size(); ++k) {
        if (k > 0) {
          out << "\n";
        }
        writeTextHeader(out, rankings[k], contexts[k]);
        writeTextRows(out, rankings[k]);
        writeTextFooter(out, rankings[k], contexts[k]);
      }
      break;
    case ReportFormat::Csv:
      // One flat table with a time column rather than a section per time, because that is
      // what loads into a spreadsheet or a dataframe without further work.
      for (std::size_t k = 0; k < rankings.size(); ++k) {
        writeCsvRows(out, rankings[k], /*withHeader=*/k == 0);
      }
      break;
    case ReportFormat::Json:
      out << "[\n";
      for (std::size_t k = 0; k < rankings.size(); ++k) {
        writeJsonRanking(out, rankings[k], contexts[k], 2);
        out << (k + 1 < rankings.size() ? ",\n" : "\n");
      }
      out << "]\n";
      break;
  }
}

void writeEvents(std::ostream& out, const EventReport& report, const ReportContext& context,
                 ReportFormat format) {
  switch (format) {
    case ReportFormat::Text:
      writeEventsText(out, report, context);
      return;
    case ReportFormat::Csv:
      writeEventsCsv(out, report);
      return;
    case ReportFormat::Json:
      writeEventsJson(out, report);
      return;
  }
}

namespace {

// The one sentence a stay-time report exists to make unambiguous, worded the same way in every
// format that has room for words.
std::string neverSpentNote(const StayTime& stay, const std::string& unit) {
  char share[32];
  std::snprintf(share, sizeof(share), "%.1f%%",
                stay.budget > 0.0 ? 100.0 * stay.accruedAtMax / stay.budget : 0.0);
  return "not spent within " + formatDuration(stay.maxDurationSeconds) +
         " -- a stay that long "
         "accrues " +
         sci(stay.accruedAtMax) + (unit.empty() ? "" : " " + unit) + ", " + share +
         " of the budget";
}

void writeStayText(std::ostream& out, const StayReport& report, const ReportContext& context) {
  out << "NuSIFT " << report.metric << " stay time\n";
  out << "  budget: " << sci(report.budget);
  if (!report.unit.empty()) {
    out << ' ' << report.unit;
  }
  out << ", over stays of up to " << formatDuration(report.maxDurationSeconds) << '\n';
  if (!context.geometry.empty()) {
    out << "  model: " << context.geometry << '\n';
  }
  if (!context.seedProvenance.empty()) {
    out << "  seed:  " << context.seedProvenance << '\n';
  }
  out << '\n';

  out << "   going in at        may stay        located to      integrals\n";
  for (const StayTime& stay : report.stays) {
    out << "   " << std::left << std::setw(16) << formatDuration(stay.startSeconds) << std::right;
    if (stay.bounded) {
      out << std::setw(12) << formatDuration(stay.durationSeconds) << "  " << std::setw(14)
          << (stay.locatedToSeconds > 0.0 ? formatDuration(stay.locatedToSeconds) : "exactly")
          << "  " << std::setw(9) << stay.samples;
      if (!stay.converged) {
        out << "   (search stopped before the tolerance was reached)";
      }
      out << '\n';
    } else {
      // Deliberately not a number in the duration column. A very large figure there would be
      // read as a duration, and the honest statement is that no duration was located.
      out << std::setw(12) << "--" << "  " << std::setw(14) << "--" << "  " << std::setw(9)
          << stay.samples << '\n';
      out << "     " << neverSpentNote(stay, report.unit) << '\n';
    }
  }

  out << "\n  A stay time is the exact interval integral inverted: the accrued total is closed\n"
         "  form within the decay model, so what is approximated is only where the root sits,\n"
         "  and the located-to column is that width rather than a tolerance that was asked for.\n";
}

void writeStayCsv(std::ostream& out, const StayReport& report) {
  out << "start_s,budget,unit,max_duration_s,bounded,duration_s,accrued_at_max,located_to_s,"
         "converged,samples\n";
  for (const StayTime& stay : report.stays) {
    out << shortestRoundTrip(stay.startSeconds) << ',' << shortestRoundTrip(stay.budget) << ','
        << csvField(report.unit) << ',' << shortestRoundTrip(stay.maxDurationSeconds) << ','
        << (stay.bounded ? "true" : "false") << ',';
    // Empty rather than zero: a spreadsheet column of durations with a 0 in it reads as "leave
    // immediately", which is the opposite of what an unbounded stay means.
    if (stay.bounded) {
      out << shortestRoundTrip(stay.durationSeconds);
    }
    out << ',' << shortestRoundTrip(stay.accruedAtMax) << ','
        << shortestRoundTrip(stay.locatedToSeconds) << ',' << (stay.converged ? "true" : "false")
        << ',' << stay.samples << '\n';
  }
}

void writeStayJson(std::ostream& out, const StayReport& report) {
  out << "{\n";
  out << "  \"metric\": \"" << report.metric << "\",\n";
  out << "  \"unit\": \"" << report.unit << "\",\n";
  out << "  \"budget\": " << jsonNumber(report.budget) << ",\n";
  out << "  \"max_duration_s\": " << jsonNumber(report.maxDurationSeconds) << ",\n";
  out << "  \"stays\": [\n";
  for (std::size_t i = 0; i < report.stays.size(); ++i) {
    const StayTime& stay = report.stays[i];
    out << "    {\"start_s\": " << jsonNumber(stay.startSeconds)
        << ", \"bounded\": " << (stay.bounded ? "true" : "false") << ", \"duration_s\": "
        << (stay.bounded ? jsonNumber(stay.durationSeconds) : std::string("null"))
        << ", \"accrued_at_max\": " << jsonNumber(stay.accruedAtMax)
        << ", \"located_to_s\": " << jsonNumber(stay.locatedToSeconds)
        << ", \"converged\": " << (stay.converged ? "true" : "false")
        << ", \"samples\": " << stay.samples << "}";
    out << (i + 1 == report.stays.size() ? "\n" : ",\n");
  }
  out << "  ]\n}\n";
}

}  // namespace

void writeStayTimes(std::ostream& out, const StayReport& report, const ReportContext& context,
                    ReportFormat format) {
  switch (format) {
    case ReportFormat::Text:
      writeStayText(out, report, context);
      return;
    case ReportFormat::Csv:
      writeStayCsv(out, report);
      return;
    case ReportFormat::Json:
      writeStayJson(out, report);
      return;
  }
}

void writeAllowable(std::ostream& out, const std::vector<AllowableScale>& scaled,
                    std::span<const Criterion> criteria, const ReportContext& context,
                    ReportFormat format) {
  switch (format) {
    case ReportFormat::Text:
      writeAllowableText(out, scaled, criteria, context);
      return;
    case ReportFormat::Csv:
      writeAllowableCsv(out, scaled, criteria);
      return;
    case ReportFormat::Json:
      writeAllowableJson(out, scaled, criteria);
      return;
  }
}

void writeInterventions(std::ostream& out, const InterventionStudy& study,
                        const ReportContext& context, ReportFormat format) {
  switch (format) {
    case ReportFormat::Text:
      writeInterventionsText(out, study, context);
      return;
    case ReportFormat::Csv:
      writeInterventionsCsv(out, study);
      return;
    case ReportFormat::Json:
      writeInterventionsJson(out, study);
      return;
  }
}

}  // namespace nusift
