#include "nusift/io/report.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <ostream>
#include <string>

#include "nusift/core/error.hpp"
#include "nusift/core/nuclide_name.hpp"
#include "nusift/io/number_format.hpp"
#include "nusift/io/time_spec.hpp"
#include "nusift/nucdata/nuclear_data.hpp"
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

// Wrap a long note under a "  ! " marker, aligned under it on continuation lines. Most notes
// here are short enough to hand-break; the ones that are not should not be re-flowed by hand
// every time a word changes.
void writeWrappedNote(std::ostream& out, const std::string& text, std::size_t width = 84) {
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
    out << (first ? "" : "    ") << text.substr(start, take) << "\n";
    first = false;
    start += take;
    while (start < text.size() && text[start] == ' ') {
      ++start;
    }
  }
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
      << "  " << std::setw(12) << "after" << "  " << std::setw(12) << "removed" << "  "
      << std::setw(8) << "of base" << "  driven by\n";
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

void writeDecaySensitivitiesText(std::ostream& out, const DecaySensitivities& s,
                                 const ReportContext& context) {
  out << "NuSIFT decay-constant sensitivity\n";
  out << "  t = " << formatDuration(s.time) << "    R = " << sci(s.response) << '\n';
  if (!context.geometry.empty()) {
    out << "  model: " << context.geometry << '\n';
  }
  if (!context.seedProvenance.empty()) {
    out << "  seed:  " << context.seedProvenance << '\n';
  }
  out << "  quadrature: " << s.endRefinements << " end refinements for a "
      << formatDuration(s.shortestRemovalSeconds) << " removal time, " << s.solves << " solves\n";
  out << '\n';

  out << "   nuclide        elasticity     T half sigma      of sigma_R   implicit+explicit\n";
  const std::size_t shown = std::min<std::size_t>(s.nuclides.size(), 12);
  for (std::size_t i = 0; i < shown; ++i) {
    const DecaySensitivity& one = s.nuclides[i];
    char terms[64];
    std::snprintf(terms, sizeof(terms), "%9.2e %+9.2e", one.implicit, one.explicitWeight);
    out << "   " << std::left << std::setw(12) << one.label << std::right << std::setw(13)
        << sci(one.elasticity) << std::setw(14)
        << (one.relativeUncertainty > 0.0 ? percent(one.relativeUncertainty) : std::string("none"))
        << std::setw(15)
        << (one.relativeUncertainty > 0.0 ? percent(one.sigmaContribution) : std::string("--"))
        << "   " << terms << '\n';
  }
  if (s.nuclides.size() > shown) {
    out << "   ... and " << (s.nuclides.size() - shown) << " more\n";
  }

  out << "\n  half-lives, root-sum-square: " << percent(s.relativeNorm) << " of R\n";

  // The branching half, and the one correlation that is derivable. Both figures are shown
  // because the DIFFERENCE between them is the result: it says what treating a nuclide's
  // branchings as independent would have cost, and the constraint that removes it is a fact
  // about the data model rather than an assumption about the evaluation.
  if (!s.branchingBlocks.empty()) {
    out << "  branchings:                  " << percent(s.branchingNorm) << " of R   ("
        << percent(s.branchingNormDiagonal) << " if taken independent)\n";
    out << "    over " << s.branchingBlocks.size()
        << " nuclide(s) with more than one evaluated branch; a single-mode nuclide has b = 1 by\n"
           "    construction and contributes nothing, whatever sigma is stated for it\n";
  } else {
    out << "  branchings:                  nothing to contribute -- every reachable nuclide is\n"
           "    single-mode, so its branching is 1 by construction\n";
  }

  out << "  ";
  writeWrappedNote(out,
                   "The half-life figure takes Sigma diagonal and is a sensitivity NORM rather "
                   "than an error budget: ENDF carries no covariance for decay data, and "
                   "evaluated half-lives are not independent of the yields fitted alongside "
                   "them. The branching figure does NOT assume independence -- the modes of one "
                   "nuclide sum to one, and that constraint is derived rather than imported. "
                   "Fission yields are neither, and are not in either figure.");
  if (s.withoutUncertainty > 0) {
    out << "  ! ";
    writeWrappedNote(out, percent(s.coveredFraction) +
                              " of the total sensitivity sits on nuclides whose half-life "
                              "uncertainty IS evaluated; " +
                              std::to_string(s.withoutUncertainty) +
                              " carry none and contribute nothing to the figure above.");
  }
  if (s.refinementCapped) {
    out << "  ! ";
    writeWrappedNote(out,
                     "The refinement cap bound before the smallest quadrature piece reached the "
                     "shortest removal time, so this is UNDER-REFINED by cram's own criterion "
                     "and may be wrong by tens of percent. Nothing in the numbers would show it. "
                     "Raise the cap, or ask about a shorter time.");
  }
  out << "  ";
  writeWrappedNote(out,
                   "The two terms are printed because their SUM is what cancels: at secular "
                   "equilibrium they annihilate across orders of magnitude, and a small total "
                   "beside two large terms is a pinned nuclide rather than an unimportant one.");
}

void writeDecaySensitivitiesCsv(std::ostream& out, const DecaySensitivities& s) {
  out << "nuclide,half_life_s,decay_constant,implicit,explicit,basis,total,elasticity,"
         "relative_uncertainty,sigma_contribution\n";
  for (const DecaySensitivity& one : s.nuclides) {
    out << csvField(one.label) << ',' << shortestRoundTrip(one.halfLifeSeconds) << ','
        << shortestRoundTrip(one.decayConstant) << ',' << shortestRoundTrip(one.implicit) << ','
        << shortestRoundTrip(one.explicitWeight) << ',' << shortestRoundTrip(one.basis) << ','
        << shortestRoundTrip(one.total) << ',' << shortestRoundTrip(one.elasticity) << ',';
    // Empty rather than zero: an evaluation that stated no uncertainty has said nothing, which
    // is not the same claim as an uncertainty of zero.
    if (one.relativeUncertainty > 0.0) {
      out << shortestRoundTrip(one.relativeUncertainty) << ','
          << shortestRoundTrip(one.sigmaContribution);
    } else {
      out << ',';
    }
    out << '\n';
  }
}

void writeDecaySensitivitiesJson(std::ostream& out, const DecaySensitivities& s) {
  out << "{\n";
  out << "  \"time_s\": " << jsonNumber(s.time) << ",\n";
  out << "  \"response\": " << jsonNumber(s.response) << ",\n";
  out << "  \"relative_norm\": " << jsonNumber(s.relativeNorm) << ",\n";
  out << "  \"covered_fraction\": " << jsonNumber(s.coveredFraction) << ",\n";
  out << "  \"with_uncertainty\": " << s.withUncertainty << ",\n";
  out << "  \"without_uncertainty\": " << s.withoutUncertainty << ",\n";
  out << "  \"end_refinements\": " << s.endRefinements << ",\n";
  out << "  \"refinement_capped\": " << (s.refinementCapped ? "true" : "false") << ",\n";
  out << "  \"solves\": " << s.solves << ",\n";
  out << "  \"nuclides\": [\n";
  for (std::size_t i = 0; i < s.nuclides.size(); ++i) {
    const DecaySensitivity& one = s.nuclides[i];
    out << "    {\"nuclide\": \"" << one.label
        << "\", \"half_life_s\": " << jsonNumber(one.halfLifeSeconds)
        << ", \"implicit\": " << jsonNumber(one.implicit)
        << ", \"explicit\": " << jsonNumber(one.explicitWeight)
        << ", \"basis\": " << jsonNumber(one.basis) << ", \"total\": " << jsonNumber(one.total)
        << ", \"elasticity\": " << jsonNumber(one.elasticity) << ", \"relative_uncertainty\": "
        << (one.relativeUncertainty > 0.0 ? jsonNumber(one.relativeUncertainty)
                                          : std::string("null"))
        << ", \"sigma_contribution\": "
        << (one.relativeUncertainty > 0.0 ? jsonNumber(one.sigmaContribution) : std::string("null"))
        << "}" << (i + 1 == s.nuclides.size() ? "\n" : ",\n");
  }
  out << "  ]\n}\n";
}

}  // namespace

void writeDecaySensitivities(std::ostream& out, const DecaySensitivities& sensitivities,
                             const ReportContext& context, ReportFormat format) {
  switch (format) {
    case ReportFormat::Text:
      writeDecaySensitivitiesText(out, sensitivities, context);
      return;
    case ReportFormat::Csv:
      writeDecaySensitivitiesCsv(out, sensitivities);
      return;
    case ReportFormat::Json:
      writeDecaySensitivitiesJson(out, sensitivities);
      return;
  }
}

namespace {

void writeUncertaintyText(std::ostream& out, const ResponseUncertainty& u,
                          const ReportContext& context) {
  out << "NuSIFT " << metricName(u.metric) << " uncertainty from the assay\n";
  out << "  t = " << formatDuration(u.time) << "    " << sci(u.response) << " +/- " << sci(u.sigma)
      << ' ' << unitName(u.unit);
  if (u.response > 0.0) {
    out << "   (" << percent(u.relative) << ")";
  }
  out << '\n';
  if (!context.geometry.empty()) {
    out << "  model: " << context.geometry << '\n';
  }
  if (!context.seedProvenance.empty()) {
    out << "  seed:  " << context.seedProvenance << '\n';
  }
  out << '\n';

  out << "   seed          of variance    1-sigma on R      row sigma     share of R\n";
  const std::size_t shown = std::min<std::size_t>(u.seeds.size(), 12);
  for (std::size_t i = 0; i < shown; ++i) {
    const SeedUncertainty& seed = u.seeds[i];
    out << "   " << std::left << std::setw(12) << seed.label << std::right << std::setw(11)
        << (seed.sigmaAtoms > 0.0 ? percent(seed.varianceFraction) : std::string("--"))
        << std::setw(16)
        << (seed.sigmaAtoms > 0.0 ? sci(seed.sigmaContribution) : std::string("--"))
        << std::setw(15)
        << (seed.sigmaAtoms > 0.0 ? sci(seed.sigmaAtoms) : std::string("none stated"))
        << std::setw(15) << percent(u.response > 0.0 ? seed.share / u.response : 0.0) << '\n';
  }
  if (u.seeds.size() > shown) {
    out << "   ... and " << (u.seeds.size() - shown) << " more seeds\n";
  }

  out << "\n  ";
  writeWrappedNote(out,
                   "This is exact rather than first order: R is linear in the seed, so "
                   "sigma_R^2 = g^T Sigma g needs no expansion and no estimated derivative. "
                   "What is assumed is the SHAPE of Sigma, below.");

  // The figure that decides whether the error bar means anything at all. An uncertainty
  // propagated from rows carrying half the answer is not an uncertainty on the answer.
  if (u.rowsWithoutSigma > 0) {
    out << "  ! ";
    char covered[32];
    std::snprintf(covered, sizeof(covered), "%.1f%%", 100.0 * u.coveredFraction);
    writeWrappedNote(out, std::string(covered) +
                              " of the response comes from rows that stated an uncertainty; the "
                              "other rows are treated as contributing none, which understates "
                              "sigma_R by an amount nothing here can bound. " +
                              std::to_string(u.rowsWithoutSigma) +
                              " seeded rows stated no uncertainty.");
  }
  out << "  ! ";
  writeWrappedNote(out,
                   "A per-row sigma is the DIAGONAL of Sigma, which asserts the assay errors are "
                   "independent. Aliquots counted on one detector against one standard are not, "
                   "and rows fitted to a total are not; off-diagonal terms move sigma_R in "
                   "either direction. Nothing here can detect that.");
  out << "  ! ";
  writeWrappedNote(out,
                   "The nuclear data is taken as exact. Half-life, branching and yield "
                   "uncertainties are a separate parameter class this does not touch: the first "
                   "two are what `sensitivity` reports, and the third is what `uncertainty "
                   "--seed-fission --yield-covariance` reports.");
}

void writeUncertaintyCsv(std::ostream& out, const ResponseUncertainty& u) {
  out << "seed,assay,carried_s,seed_atoms,sigma_atoms,importance,share,sigma_contribution,"
         "variance_fraction\n";
  for (const SeedUncertainty& seed : u.seeds) {
    out << csvField(seed.label) << ',' << csvField(seed.assay) << ','
        << shortestRoundTrip(seed.carriedSeconds) << ',' << shortestRoundTrip(seed.seedAtoms)
        << ',';
    // Empty rather than zero: "stated no uncertainty" and "stated an uncertainty of zero" are
    // different claims, and only one of them is ever true of a measurement.
    if (seed.sigmaAtoms > 0.0) {
      out << shortestRoundTrip(seed.sigmaAtoms);
    }
    out << ',' << shortestRoundTrip(seed.importance) << ',' << shortestRoundTrip(seed.share) << ','
        << shortestRoundTrip(seed.sigmaContribution) << ','
        << shortestRoundTrip(seed.varianceFraction) << '\n';
  }
}

void writeUncertaintyJson(std::ostream& out, const ResponseUncertainty& u) {
  out << "{\n";
  out << "  \"metric\": \"" << metricName(u.metric) << "\",\n";
  out << "  \"unit\": \"" << unitName(u.unit) << "\",\n";
  out << "  \"time_s\": " << jsonNumber(u.time) << ",\n";
  out << "  \"response\": " << jsonNumber(u.response) << ",\n";
  out << "  \"sigma\": " << jsonNumber(u.sigma) << ",\n";
  out << "  \"relative\": " << jsonNumber(u.relative) << ",\n";
  out << "  \"covered_fraction\": " << jsonNumber(u.coveredFraction) << ",\n";
  out << "  \"rows_with_sigma\": " << u.rowsWithSigma << ",\n";
  out << "  \"rows_without_sigma\": " << u.rowsWithoutSigma << ",\n";
  out << "  \"seeds\": [\n";
  for (std::size_t i = 0; i < u.seeds.size(); ++i) {
    const SeedUncertainty& seed = u.seeds[i];
    out << "    {\"seed\": \"" << seed.label << "\", \"assay\": \"" << seed.assay
        << "\", \"carried_s\": " << jsonNumber(seed.carriedSeconds)
        << ", \"seed_atoms\": " << jsonNumber(seed.seedAtoms) << ", \"sigma_atoms\": "
        << (seed.sigmaAtoms > 0.0 ? jsonNumber(seed.sigmaAtoms) : std::string("null"))
        << ", \"importance\": " << jsonNumber(seed.importance)
        << ", \"share\": " << jsonNumber(seed.share)
        << ", \"sigma_contribution\": " << jsonNumber(seed.sigmaContribution)
        << ", \"variance_fraction\": " << jsonNumber(seed.varianceFraction) << "}"
        << (i + 1 == u.seeds.size() ? "\n" : ",\n");
  }
  out << "  ]\n}\n";
}

}  // namespace

void writeUncertainty(std::ostream& out, const ResponseUncertainty& uncertainty,
                      const ReportContext& context, ReportFormat format) {
  switch (format) {
    case ReportFormat::Text:
      writeUncertaintyText(out, uncertainty, context);
      return;
    case ReportFormat::Csv:
      writeUncertaintyCsv(out, uncertainty);
      return;
    case ReportFormat::Json:
      writeUncertaintyJson(out, uncertainty);
      return;
  }
}

namespace {

void writeYieldUncertaintyText(std::ostream& out, const YieldUncertainty& u,
                               const ReportContext& context) {
  out << "NuSIFT " << metricName(u.metric) << " uncertainty from the evaluated fission yields\n";
  out << "  t = " << formatDuration(u.time) << "    " << sci(u.response) << ' ' << unitName(u.unit)
      << "    " << sci(u.fissions) << " fissions\n";

  // The two figures together, because the comparison IS the result. Both run over the products
  // the matrix carries, so their whole difference is the off-diagonal.
  out << "  diagonal    +/- " << sci(u.sigmaDiagonalMatched) << "   ("
      << percent(u.relativeDiagonalMatched) << ")\n";
  if (u.varianceNegative) {
    out << "  correlated  INDEFINITE: g^T Sigma g = " << sci(u.varianceCorrelated)
        << ", which has no square root\n";
  } else {
    out << "  correlated  +/- " << sci(u.sigmaCorrelated) << "   (" << percent(u.relativeCorrelated)
        << ")   " << std::fixed << std::setprecision(2) << u.varianceRatio
        << "x the diagonal variance\n"
        << std::defaultfloat << std::setprecision(6);
  }
  if (!context.geometry.empty()) {
    out << "  model: " << context.geometry << '\n';
  }
  if (!context.seedProvenance.empty()) {
    out << "  seed:  " << context.seedProvenance << '\n';
  }
  out << '\n';

  out << "   product       of variance    1-sigma on R     yield sigma     share of R\n";
  const std::size_t shown = std::min<std::size_t>(u.products.size(), 12);
  for (std::size_t i = 0; i < shown; ++i) {
    const YieldContribution& product = u.products[i];
    const bool stated = product.sigmaYield > 0.0;
    out << "   " << (product.correlated ? ' ' : '*') << std::left << std::setw(12) << product.label
        << std::right << std::setw(11)
        << (stated ? percent(product.varianceFraction) : std::string("--")) << std::setw(16)
        << (stated ? sci(product.sigmaContribution) : std::string("--")) << std::setw(16)
        << (stated ? sci(product.sigmaYield) : std::string("none stated")) << std::setw(15)
        << percent(u.response > 0.0 ? product.share / u.response : 0.0) << '\n';
  }
  if (u.products.size() > shown) {
    out << "   ... and " << (u.products.size() - shown) << " more products\n";
  }
  if (u.productsUnmatched > 0) {
    out << "   * not carried by the correlation matrix; diagonal term only\n";
  }

  out << "\n  ";
  writeWrappedNote(out,
                   "A fission seed is n0 = N_f Y with N_f exact, so this is the same exact "
                   "propagation the assay error bar uses -- sigma_R^2 = g^T (N_f^2 Sigma_Y) g -- "
                   "with the off-diagonal of Sigma_Y imported rather than assumed away.");

  out << "  ! ";
  writeWrappedNote(out, "Correlation: " + u.provenance.system + " from " + u.provenance.library +
                            ", paired with the sigma_Y this store stages from " + u.storeLibrary +
                            ". A correlation is dimensionless and structural, which is the "
                            "argument that it carries across an edition better than a variance "
                            "would; it is an argument and not a proof.");
  out << "    " << u.provenance.citation << '\n';

  // Indefiniteness is a property of the published product, reported wherever it shows up rather
  // than only when it bites. A reader who sees a correlated figure without this line would have
  // no way to know the matrix is not a covariance in the mathematical sense.
  if (u.smallestEigenvalue < 0.0) {
    char eigen[64];
    std::snprintf(eigen, sizeof(eigen), "%.4g", u.smallestEigenvalue);
    out << "  ! ";
    writeWrappedNote(out,
                     std::string("The imported matrix is not positive semi-definite on these "
                                 "products: its smallest eigenvalue is ") +
                         eigen +
                         ". That is expected of a stochastic estimate and is NOT repaired here -- "
                         "projecting onto the nearest PSD matrix would replace a published "
                         "correlation with one no evaluation stands behind.");
  }
  if (u.varianceNegative) {
    out << "  ! ";
    writeWrappedNote(out,
                     "The contraction came out negative, which the indefiniteness above permits. "
                     "No standard deviation is reported for it. The diagonal figure still holds.");
  }
  if (u.productsUnmatched > 0 || u.coveredFraction < 1.0) {
    char covered[32];
    std::snprintf(covered, sizeof(covered), "%.1f%%", 100.0 * u.coveredFraction);
    out << "  ! ";
    writeWrappedNote(out, std::string(covered) +
                              " of the response sits on products the matrix carries; " +
                              std::to_string(u.productsUnmatched) +
                              " seeded products are outside it and keep their diagonal term "
                              "only. The correlated figure is over the covered set.");
  }
  if (u.productsWithoutSigma > 0) {
    char covered[32];
    std::snprintf(covered, sizeof(covered), "%.1f%%", 100.0 * u.sigmaCoveredFraction);
    out << "  ! ";
    writeWrappedNote(out, std::string(covered) +
                              " of the response comes from products whose yield states an "
                              "uncertainty; " +
                              std::to_string(u.productsWithoutSigma) +
                              " products state none and are treated as contributing no variance, "
                              "which understates sigma_R by an amount nothing here can bound.");
  }
  out << "  ! ";
  writeWrappedNote(out,
                   "This is ONE parameter class. Combining it with the half-life and branching "
                   "figures needs the cross-covariance between the three, which no evaluation "
                   "publishes and this resampling did not estimate, so there is no total here.");
  out << "  ! ";
  writeWrappedNote(out,
                   "The number of fissions is taken as exact. An uncertainty on the source term "
                   "itself scales R and this sigma together and is the user's to carry.");
}

void writeYieldUncertaintyCsv(std::ostream& out, const YieldUncertainty& u) {
  out << "product,yield,sigma_yield,seed_atoms,importance,share,sigma_contribution,"
         "variance_fraction,correlated\n";
  for (const YieldContribution& product : u.products) {
    out << csvField(product.label) << ',' << shortestRoundTrip(product.yield) << ',';
    // Empty rather than zero, as the assay report does it: "states no uncertainty" and "states
    // an uncertainty of zero" are different claims about an evaluation.
    if (product.sigmaYield > 0.0) {
      out << shortestRoundTrip(product.sigmaYield);
    }
    out << ',' << shortestRoundTrip(product.seedAtoms) << ','
        << shortestRoundTrip(product.importance) << ',' << shortestRoundTrip(product.share) << ','
        << shortestRoundTrip(product.sigmaContribution) << ','
        << shortestRoundTrip(product.varianceFraction) << ',' << (product.correlated ? 1 : 0)
        << '\n';
  }
}

void writeYieldUncertaintyJson(std::ostream& out, const YieldUncertainty& u) {
  out << "{\n";
  out << "  \"metric\": \"" << escapeJson(metricName(u.metric)) << "\",\n";
  out << "  \"unit\": \"" << escapeJson(unitName(u.unit)) << "\",\n";
  out << "  \"time_s\": " << jsonNumber(u.time) << ",\n";
  out << "  \"fissions\": " << jsonNumber(u.fissions) << ",\n";
  out << "  \"incident_energy_ev\": " << jsonNumber(u.incidentEnergyEv) << ",\n";
  out << "  \"response\": " << jsonNumber(u.response) << ",\n";
  out << "  \"sigma_diagonal\": " << jsonNumber(u.sigmaDiagonal) << ",\n";
  out << "  \"relative_diagonal\": " << jsonNumber(u.relativeDiagonal) << ",\n";
  out << "  \"sigma_diagonal_matched\": " << jsonNumber(u.sigmaDiagonalMatched) << ",\n";
  out << "  \"relative_diagonal_matched\": " << jsonNumber(u.relativeDiagonalMatched) << ",\n";
  // Null rather than zero when the form is negative: there is no standard deviation, and a zero
  // would read as a vanishing one.
  out << "  \"sigma_correlated\": "
      << (u.varianceNegative ? std::string("null") : jsonNumber(u.sigmaCorrelated)) << ",\n";
  out << "  \"relative_correlated\": "
      << (u.varianceNegative ? std::string("null") : jsonNumber(u.relativeCorrelated)) << ",\n";
  out << "  \"variance_correlated\": " << jsonNumber(u.varianceCorrelated) << ",\n";
  out << "  \"variance_ratio\": "
      << (u.varianceNegative ? std::string("null") : jsonNumber(u.varianceRatio)) << ",\n";
  out << "  \"smallest_eigenvalue\": " << jsonNumber(u.smallestEigenvalue) << ",\n";
  out << "  \"covered_fraction\": " << jsonNumber(u.coveredFraction) << ",\n";
  out << "  \"sigma_covered_fraction\": " << jsonNumber(u.sigmaCoveredFraction) << ",\n";
  out << "  \"products_matched\": " << u.productsMatched << ",\n";
  out << "  \"products_unmatched\": " << u.productsUnmatched << ",\n";
  out << "  \"products_with_sigma\": " << u.productsWithSigma << ",\n";
  out << "  \"products_without_sigma\": " << u.productsWithoutSigma << ",\n";
  out << "  \"products_unused_in_matrix\": " << u.productsUnusedInMatrix << ",\n";
  // Every one of these is free text from a file path, a CLI flag or a store attribute, so each
  // goes through escapeJson rather than straight into the quotes.
  out << "  \"correlation\": {\"system\": \"" << escapeJson(u.provenance.system)
      << "\", \"library\": \"" << escapeJson(u.provenance.library) << "\", \"path\": \""
      << escapeJson(u.provenance.path) << "\", \"citation\": \""
      << escapeJson(u.provenance.citation) << "\"},\n";
  out << "  \"store_library\": \"" << escapeJson(u.storeLibrary) << "\",\n";
  out << "  \"products\": [\n";
  for (std::size_t i = 0; i < u.products.size(); ++i) {
    const YieldContribution& product = u.products[i];
    out << "    {\"product\": \"" << escapeJson(product.label)
        << "\", \"yield\": " << jsonNumber(product.yield) << ", \"sigma_yield\": "
        << (product.sigmaYield > 0.0 ? jsonNumber(product.sigmaYield) : std::string("null"))
        << ", \"seed_atoms\": " << jsonNumber(product.seedAtoms)
        << ", \"importance\": " << jsonNumber(product.importance)
        << ", \"share\": " << jsonNumber(product.share)
        << ", \"sigma_contribution\": " << jsonNumber(product.sigmaContribution)
        << ", \"variance_fraction\": " << jsonNumber(product.varianceFraction)
        << ", \"correlated\": " << (product.correlated ? "true" : "false") << "}"
        << (i + 1 == u.products.size() ? "\n" : ",\n");
  }
  out << "  ]\n}\n";
}

}  // namespace

void writeYieldUncertainty(std::ostream& out, const YieldUncertainty& uncertainty,
                           const ReportContext& context, ReportFormat format) {
  switch (format) {
    case ReportFormat::Text:
      writeYieldUncertaintyText(out, uncertainty, context);
      return;
    case ReportFormat::Csv:
      writeYieldUncertaintyCsv(out, uncertainty);
      return;
    case ReportFormat::Json:
      writeYieldUncertaintyJson(out, uncertainty);
      return;
  }
}

namespace {

void writeTaskPlanText(std::ostream& out, const TaskPlan& plan, const ReportContext& context) {
  out << "NuSIFT task plan\n";
  out << "  starting " << formatDuration(plan.startSeconds) << ", " << plan.legs.size()
      << " legs over " << formatDuration(plan.elapsedSeconds) << " ("
      << formatDuration(plan.exposedSeconds) << " in the field)\n";
  out << "  total: " << sci(plan.total) << ' ' << unitName(plan.unit) << '\n';
  if (!context.geometry.empty()) {
    out << "  model: " << context.geometry << '\n';
  }
  if (!context.seedProvenance.empty()) {
    out << "  seed:  " << context.seedProvenance << '\n';
  }
  out << '\n';

  // "for 3 m" beside "at 4 m" is three minutes beside four meters, and the two columns cannot
  // be told apart by eye. The distance is a bare number under a unit-bearing header instead.
  out << "   leg                 starts   duration   dist/m     occ           "
      << unitName(plan.unit) << "     frac       per h\n";
  for (const LegResult& leg : plan.legs) {
    out << "   " << std::left << std::setw(18) << leg.name << std::right << std::setw(9)
        << formatDuration(leg.startSeconds) << std::setw(11)
        << formatDuration(leg.endSeconds - leg.startSeconds);

    if (leg.isBreak) {
      // A break has no distance and no rate. Printing zeros in those columns would read as
      // "standing at the source, receiving nothing", which is a different and untrue claim.
      out << std::setw(9) << "--" << std::setw(8) << "break" << std::setw(14) << "--"
          << std::setw(9) << "--" << std::setw(12) << "--" << '\n';
      continue;
    }
    char distance[16];
    std::snprintf(distance, sizeof(distance), "%.4g", leg.distanceM);
    out << std::setw(9) << distance << std::setw(8) << percent(leg.occupancy) << std::setw(14)
        << sci(leg.accrued) << std::setw(9) << percent(leg.fraction) << std::setw(12)
        << sci(leg.meanRate * 3600.0) << '\n';
  }

  if (plan.budget > 0.0) {
    out << "\n  budget: " << sci(plan.budget) << ' ' << unitName(plan.unit);
    if (!plan.budgetSpent) {
      out << "  -- the plan fits, with " << sci(plan.budget - plan.total) << " to spare\n";
    } else if (plan.spentInLeg >= 0) {
      const LegResult& leg = plan.legs[static_cast<std::size_t>(plan.spentInLeg)];
      out << "  -- SPENT during \"" << leg.name << "\", "
          << formatDuration(plan.spentAtSeconds - leg.startSeconds) << " into it\n";
      out << "    the plan as written accrues " << sci(plan.total) << ", which is "
          << percent(plan.total / plan.budget) << " of it\n";
    }
  }

  out << "\n  ";
  writeWrappedNote(out,
                   "Each leg is an exact interval integral with its own geometry, summed. No leg "
                   "is shielded -- every one is an unshielded point source in air -- so a leg "
                   "behind a wall is overstated by an amount this cannot know.");
  if (std::any_of(plan.legs.begin(), plan.legs.end(),
                  [](const LegResult& leg) { return !leg.isBreak && leg.occupancy < 1.0; })) {
    out << "  ";
    writeWrappedNote(out,
                     "An occupancy below 1 is applied as a factor, which assumes the presence is "
                     "spread evenly over the leg. Exact for a leg short against the decay; for a "
                     "long one, split it into the stretches actually spent there.");
  }
}

void writeTaskPlanCsv(std::ostream& out, const TaskPlan& plan) {
  out << "leg,start_s,end_s,distance_m,occupancy,accrued,unit,fraction,cumulative,mean_rate_per_s,"
         "is_break\n";
  for (const LegResult& leg : plan.legs) {
    out << csvField(leg.name) << ',' << shortestRoundTrip(leg.startSeconds) << ','
        << shortestRoundTrip(leg.endSeconds) << ',';
    // Empty rather than zero for a break: a distance of 0 m is a real and very different claim.
    if (!leg.isBreak) {
      out << shortestRoundTrip(leg.distanceM);
    }
    out << ',' << shortestRoundTrip(leg.occupancy) << ',' << shortestRoundTrip(leg.accrued) << ','
        << unitName(plan.unit) << ',' << shortestRoundTrip(leg.fraction) << ','
        << shortestRoundTrip(leg.cumulative) << ',' << shortestRoundTrip(leg.meanRate) << ','
        << (leg.isBreak ? "true" : "false") << '\n';
  }
}

void writeTaskPlanJson(std::ostream& out, const TaskPlan& plan) {
  out << "{\n";
  out << "  \"unit\": \"" << unitName(plan.unit) << "\",\n";
  out << "  \"start_s\": " << jsonNumber(plan.startSeconds) << ",\n";
  out << "  \"end_s\": " << jsonNumber(plan.endSeconds) << ",\n";
  out << "  \"total\": " << jsonNumber(plan.total) << ",\n";
  out << "  \"elapsed_s\": " << jsonNumber(plan.elapsedSeconds) << ",\n";
  out << "  \"exposed_s\": " << jsonNumber(plan.exposedSeconds) << ",\n";
  if (plan.budget > 0.0) {
    out << "  \"budget\": " << jsonNumber(plan.budget) << ",\n";
    out << "  \"budget_spent\": " << (plan.budgetSpent ? "true" : "false") << ",\n";
    // Null rather than -1 and 0: a parser reading an index of -1 as a leg would be reading the
    // last one, and a time of zero as an instant would be reading the start of the run.
    out << "  \"spent_in_leg\": "
        << (plan.budgetSpent ? std::to_string(plan.spentInLeg) : std::string("null")) << ",\n";
    out << "  \"spent_at_s\": "
        << (plan.budgetSpent ? jsonNumber(plan.spentAtSeconds) : std::string("null")) << ",\n";
  }
  out << "  \"legs\": [\n";
  for (std::size_t i = 0; i < plan.legs.size(); ++i) {
    const LegResult& leg = plan.legs[i];
    out << "    {\"name\": \"" << leg.name << "\", \"start_s\": " << jsonNumber(leg.startSeconds)
        << ", \"end_s\": " << jsonNumber(leg.endSeconds)
        << ", \"distance_m\": " << (leg.isBreak ? std::string("null") : jsonNumber(leg.distanceM))
        << ", \"occupancy\": " << jsonNumber(leg.occupancy)
        << ", \"accrued\": " << jsonNumber(leg.accrued)
        << ", \"fraction\": " << jsonNumber(leg.fraction)
        << ", \"cumulative\": " << jsonNumber(leg.cumulative)
        << ", \"mean_rate_per_s\": " << jsonNumber(leg.meanRate)
        << ", \"is_break\": " << (leg.isBreak ? "true" : "false") << "}"
        << (i + 1 == plan.legs.size() ? "\n" : ",\n");
  }
  out << "  ]\n}\n";
}

}  // namespace

void writeTaskPlan(std::ostream& out, const TaskPlan& plan, const ReportContext& context,
                   ReportFormat format) {
  switch (format) {
    case ReportFormat::Text:
      writeTaskPlanText(out, plan, context);
      return;
    case ReportFormat::Csv:
      writeTaskPlanCsv(out, plan);
      return;
    case ReportFormat::Json:
      writeTaskPlanJson(out, plan);
      return;
  }
}

namespace {

void writeTriageSetText(std::ostream& out, const TriageSet& set,
                        std::span<const CoverageRequirement> requirements,
                        const ReportContext& context) {
  out << "NuSIFT robust triage set\n";
  out << "  " << set.members.size() << " of " << set.candidateCount << " contributors, holding "
      << percent(set.binding.required) << " of every requirement at every time\n";
  out << "  requirements:\n";
  for (const CoverageRequirement& requirement : requirements) {
    out << "    " << std::left << std::setw(28) << requirement.label << std::right
        << percent(requirement.fraction) << " at each of " << requirement.table->timeCount()
        << " times\n";
  }
  if (!context.geometry.empty()) {
    out << "  model: " << context.geometry << '\n';
  }
  if (!context.seedProvenance.empty()) {
    out << "  seed:  " << context.seedProvenance << '\n';
  }
  out << '\n';

  // First, because it is the number that says whether the answer is comfortable. A reader who
  // saw only the list would have no way to tell a set with room to spare from one that clears
  // its floor by a thousandth at a single instant.
  out << "  closest to failing: " << set.binding.requirement << " at "
      << formatDuration(set.binding.timeSeconds) << " -- holds " << percent(set.binding.achieved)
      << " against " << percent(set.binding.required) << '\n';
  out << '\n';

  out << "   #  contributor        peak share   closed\n";
  for (const SetMember& member : set.members) {
    out << "  " << std::setw(2) << member.order << "  " << std::left << std::setw(18)
        << member.label << std::right << std::setw(10) << percent(member.peakFraction) << "   "
        << sci(member.closedShortfall) << '\n';
  }

  out << "\n  The order is not a ranking. The second member is whichever most improved the\n"
         "  constraints still unmet GIVEN the first, which is usually not the second largest\n"
         "  contributor to anything -- a member with a small peak share is there to hold up one\n"
         "  particular instant, and dropping it is what the floor forbids.\n";

  if (!set.shortfalls.empty()) {
    out << "\n  ! " << set.shortfalls.size()
        << " constraints are not met by ANY set these tables can form:\n";
    const std::size_t shown = std::min<std::size_t>(set.shortfalls.size(), 5);
    for (std::size_t i = 0; i < shown; ++i) {
      out << "      " << set.shortfalls[i].requirement << " at "
          << formatDuration(set.shortfalls[i].timeSeconds) << ": "
          << percent(set.shortfalls[i].achieved) << " of a required "
          << percent(set.shortfalls[i].required) << '\n';
    }
    if (set.shortfalls.size() > shown) {
      out << "      ... and " << (set.shortfalls.size() - shown) << " more\n";
    }
    out << "    ";
    writeWrappedNote(out,
                     "Every column together falls short, so the floor is above what the table can "
                     "express -- a gamma-line total counts lines below the column threshold, "
                     "which no column carries.");
  }

  out << "\n  ";
  writeWrappedNote(out,
                   "Set cover is NP-hard and this is a deterministic greedy: a small set that "
                   "meets the floor, not a proof that none smaller exists.");
}

void writeTriageSetCsv(std::ostream& out, const TriageSet& set) {
  out << "kind,order,contributor,peak_fraction,closed_shortfall,requirement,time_s,achieved,"
         "required\n";
  for (const SetMember& member : set.members) {
    out << "member," << member.order << ',' << csvField(member.label) << ','
        << shortestRoundTrip(member.peakFraction) << ','
        << shortestRoundTrip(member.closedShortfall) << ",,,,\n";
  }
  out << "binding,,,,," << csvField(set.binding.requirement) << ','
      << shortestRoundTrip(set.binding.timeSeconds) << ','
      << shortestRoundTrip(set.binding.achieved) << ',' << shortestRoundTrip(set.binding.required)
      << '\n';
  for (const CoveragePoint& point : set.shortfalls) {
    out << "shortfall,,,,," << csvField(point.requirement) << ','
        << shortestRoundTrip(point.timeSeconds) << ',' << shortestRoundTrip(point.achieved) << ','
        << shortestRoundTrip(point.required) << '\n';
  }
}

void writeTriageSetJson(std::ostream& out, const TriageSet& set) {
  out << "{\n";
  out << "  \"candidates\": " << set.candidateCount << ",\n";
  out << "  \"constraints\": " << set.constraintCount << ",\n";
  out << "  \"members\": [\n";
  for (std::size_t i = 0; i < set.members.size(); ++i) {
    const SetMember& member = set.members[i];
    out << "    {\"order\": " << member.order << ", \"contributor\": \"" << member.label
        << "\", \"peak_fraction\": " << jsonNumber(member.peakFraction)
        << ", \"closed_shortfall\": " << jsonNumber(member.closedShortfall) << "}"
        << (i + 1 == set.members.size() ? "\n" : ",\n");
  }
  out << "  ],\n";
  out << "  \"binding\": {\"requirement\": \"" << set.binding.requirement
      << "\", \"time_s\": " << jsonNumber(set.binding.timeSeconds)
      << ", \"achieved\": " << jsonNumber(set.binding.achieved)
      << ", \"required\": " << jsonNumber(set.binding.required) << "},\n";
  out << "  \"shortfalls\": [\n";
  for (std::size_t i = 0; i < set.shortfalls.size(); ++i) {
    const CoveragePoint& point = set.shortfalls[i];
    out << "    {\"requirement\": \"" << point.requirement
        << "\", \"time_s\": " << jsonNumber(point.timeSeconds)
        << ", \"achieved\": " << jsonNumber(point.achieved)
        << ", \"required\": " << jsonNumber(point.required) << "}"
        << (i + 1 == set.shortfalls.size() ? "\n" : ",\n");
  }
  out << "  ]\n}\n";
}

}  // namespace

void writeTriageSet(std::ostream& out, const TriageSet& set,
                    std::span<const CoverageRequirement> requirements, const ReportContext& context,
                    ReportFormat format) {
  switch (format) {
    case ReportFormat::Text:
      writeTriageSetText(out, set, requirements, context);
      return;
    case ReportFormat::Csv:
      writeTriageSetCsv(out, set);
      return;
    case ReportFormat::Json:
      writeTriageSetJson(out, set);
      return;
  }
}

namespace {

// Said in every format that has room for words, because it is the one thing about a merged
// inventory that no column in it can show.
constexpr const char* kCarryNote =
    "Reconciliation propagates what a sheet MEASURED. A daughter that grew in during the carry "
    "is modeled and appears here; a daughter that was present at an assay and not written down "
    "is not recovered by anything. So the merged rows mix measured and modeled amounts, and "
    "which a row is depends on how far its assay was carried.";

void writeReconciliationText(std::ostream& out, const Reconciliation& reconciled,
                             const NuclearData& data, const ReportContext& context) {
  out << "NuSIFT inventory reconciliation\n";
  out << "  epoch: " << formatCalendarDate(reconciled.epochSeconds)
      << "   -- every assay carried forward to it, none carried back\n";
  if (reconciled.spanSeconds > 0.0) {
    out << "  span:  " << formatDuration(reconciled.spanSeconds)
        << " between the earliest and latest assay\n";
  }
  if (!context.seedProvenance.empty()) {
    out << "  seed:  " << context.seedProvenance << '\n';
  }
  if (!context.storeLibrary.empty()) {
    out << "  store: " << context.storeLibrary << " (" << context.storeNuclideCount
        << " nuclides, staged " << context.storeCreatedUtc << ")\n";
  }
  out << '\n';

  out << "   assayed              carried   nuclides    atoms at assay    atoms at epoch\n";
  for (const AssayContribution& one : reconciled.contributions) {
    out << "   " << std::left << std::setw(20) << formatCalendarDate(one.dateSeconds) << std::right
        << std::setw(9)
        << (one.carriedSeconds > 0.0 ? formatDuration(one.carriedSeconds) : std::string("--"))
        << "   " << std::setw(8) << one.nuclides << "    " << std::setw(14) << sci(one.atomsAtAssay)
        << "    " << std::setw(14) << sci(one.atomsAtEpoch) << '\n';
  }
  out << '\n';
  out << "  merged: " << reconciled.inventory.size() << " nuclides, "
      << sci(reconciled.inventory.totalAtoms()) << " atoms\n";

  // The largest rows of the result, so the reader sees what they are about to seed a run with
  // rather than only how it was assembled.
  std::vector<InventoryEntry> rows(reconciled.inventory.entries().begin(),
                                   reconciled.inventory.entries().end());
  std::sort(rows.begin(), rows.end(),
            [](const InventoryEntry& a, const InventoryEntry& b) { return a.atoms > b.atoms; });
  const std::size_t shown = std::min<std::size_t>(rows.size(), 10);
  out << '\n';
  out << "   nuclide             atoms            Bq\n";
  for (std::size_t i = 0; i < shown; ++i) {
    const Zai zai = Zai::fromKey(rows[i].zaiKey);
    const int index = data.indexOf(zai);
    const double activity = index >= 0 ? data.decayConstant(index) * rows[i].atoms : 0.0;
    out << "   " << std::left << std::setw(12) << formatNuclideName(zai) << std::right
        << std::setw(14) << sci(rows[i].atoms) << "    " << std::setw(10)
        << (activity > 0.0 ? sci(activity) : std::string("stable")) << '\n';
  }
  if (rows.size() > shown) {
    out << "   ... and " << (rows.size() - shown) << " more\n";
  }

  out << "\n  ! ";
  writeWrappedNote(out, kCarryNote);
}

void writeReconciliationCsv(std::ostream& out, const Reconciliation& reconciled,
                            const NuclearData& data) {
  // Two tables would not be a CSV, so the assays and the merged rows share one under a `kind`
  // column -- the same shape writeEventsCsv() uses for events and windows.
  out << "kind,label,assayed,carried_s,nuclide,atoms,atoms_at_assay,atoms_at_epoch\n";
  for (const AssayContribution& one : reconciled.contributions) {
    out << "assay," << csvField(one.label) << ',' << formatCalendarDate(one.dateSeconds) << ','
        << shortestRoundTrip(one.carriedSeconds) << ",," << ','
        << shortestRoundTrip(one.atomsAtAssay) << ',' << shortestRoundTrip(one.atomsAtEpoch)
        << '\n';
  }
  for (const InventoryEntry& entry : reconciled.inventory.entries()) {
    out << "merged,," << formatCalendarDate(reconciled.epochSeconds) << ",,"
        << formatNuclideName(Zai::fromKey(entry.zaiKey)) << ',' << shortestRoundTrip(entry.atoms)
        << ",,\n";
  }
  (void)data;
}

void writeReconciliationJson(std::ostream& out, const Reconciliation& reconciled) {
  out << "{\n";
  out << "  \"epoch\": \"" << formatCalendarDate(reconciled.epochSeconds) << "\",\n";
  out << "  \"epoch_s\": " << jsonNumber(reconciled.epochSeconds) << ",\n";
  out << "  \"span_s\": " << jsonNumber(reconciled.spanSeconds) << ",\n";
  out << "  \"note\": \"" << kCarryNote << "\",\n";
  out << "  \"assays\": [\n";
  for (std::size_t i = 0; i < reconciled.contributions.size(); ++i) {
    const AssayContribution& one = reconciled.contributions[i];
    out << "    {\"label\": \"" << one.label << "\", \"assayed\": \""
        << formatCalendarDate(one.dateSeconds)
        << "\", \"carried_s\": " << jsonNumber(one.carriedSeconds)
        << ", \"nuclides\": " << one.nuclides
        << ", \"atoms_at_assay\": " << jsonNumber(one.atomsAtAssay)
        << ", \"atoms_at_epoch\": " << jsonNumber(one.atomsAtEpoch) << "}"
        << (i + 1 == reconciled.contributions.size() ? "\n" : ",\n");
  }
  out << "  ],\n";
  out << "  \"merged\": [\n";
  const std::span<const InventoryEntry> entries = reconciled.inventory.entries();
  for (std::size_t i = 0; i < entries.size(); ++i) {
    out << "    {\"nuclide\": \"" << formatNuclideName(Zai::fromKey(entries[i].zaiKey))
        << "\", \"atoms\": " << jsonNumber(entries[i].atoms) << "}"
        << (i + 1 == entries.size() ? "\n" : ",\n");
  }
  out << "  ]\n}\n";
}

}  // namespace

void writeReconciliation(std::ostream& out, const Reconciliation& reconciled,
                         const NuclearData& data, const ReportContext& context,
                         ReportFormat format) {
  switch (format) {
    case ReportFormat::Text:
      writeReconciliationText(out, reconciled, data, context);
      return;
    case ReportFormat::Csv:
      writeReconciliationCsv(out, reconciled, data);
      return;
    case ReportFormat::Json:
      writeReconciliationJson(out, reconciled);
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
