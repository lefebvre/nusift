#include "nusift/io/report.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <ostream>
#include <string>

#include "nusift/core/error.hpp"
#include "nusift/io/time_spec.hpp"
#include "nusift/triage/response.hpp"

namespace nusift {
namespace {

constexpr const char* kModule = "report";

std::string sci(double value) {
  char buffer[32];
  std::snprintf(buffer, sizeof(buffer), "%.4e", value);
  return buffer;
}

// The shortest decimal that reads back as the same double.
//
// CSV and JSON exist to be parsed again, so surviving the round trip is the requirement, not a
// nicety: `--at 1.23456789y` printed at six significant digits comes back as a different time
// than the one the report describes. The text format's %.4e is a display choice for a terminal
// and stays one.
//
// The ladder is what keeps the output readable. %.17g always round-trips but renders 0.99999
// as 0.99999000000000005; trying the shorter forms first prints the digits the value actually
// has and falls back only when they are not enough.
std::string exact(double value) {
  char buffer[40];
  for (const int digits : {15, 16, 17}) {
    std::snprintf(buffer, sizeof(buffer), "%.*g", digits, value);
    if (std::strtod(buffer, nullptr) == value) {
      break;
    }
  }
  return buffer;
}

// JSON has no infinity and no NaN, so a non-finite value written as a bare `inf` produces a
// document that fails to parse -- in whatever consumes the report rather than here. `null` is
// what every parser accepts for "this is not a number".
std::string jsonNumber(double value) {
  return std::isfinite(value) ? exact(value) : std::string("null");
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

void writeTextHeader(std::ostream& out, const Ranking& ranking, const ReportContext& context) {
  out << "NuSIFT " << metricName(ranking.metric) << " ranking by "
      << aggregateName(ranking.aggregate) << '\n';
  out << "  t = " << whenOf(ranking) << "    total = " << sci(ranking.total) << ' '
      << unitName(ranking.unit) << '\n';
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
      << unitName(ranking.unit) << std::setw(9) << "frac" << std::setw(9) << "cum" << '\n';

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

  // The magnitude first, because it is what decides whether the count matters at all. A
  // hundred flagged nuclides contributing 0.01% of the photon output is a footnote; three
  // contributing 30% is a reason not to trust the number above.
  if (ranking.unmodeledEnergyFraction > 0.0) {
    out << "  ! " << percent(ranking.unmodeledEnergyFraction)
        << " of the emitted photon energy is in spectra NuSIFT does not model,\n"
        << "    so these exposures are understated by roughly that much";
    // Terminated here unless the named list below continues the sentence. Left open, the line
    // runs into whatever is written next -- and with several --at times that is the blank line
    // writeRankings lays between rankings, which then disappears.
    if (context.unmodeledContinuum.empty()) {
      out << '\n';
    }
  }

  if (!context.unmodeledContinuum.empty()) {
    const std::size_t total = context.unmodeledContinuum.size();
    const bool one = total == 1;
    // Indented under the magnitude line when there is one, since it is the detail behind it.
    if (ranking.unmodeledEnergyFraction > 0.0) {
      out << " (" << total << " nuclide" << (one ? "" : "s") << ")";
    } else {
      out << "  ! " << total << " contributor" << (one ? "" : "s") << (one ? " carries" : " carry")
          << " photon energy NuSIFT does not model, so " << (one ? "its" : "their")
          << (one ? " exposure is" : " exposures are") << " understated";
    }

    // Naming every one of them is what a real evaluation turns this into: a full store flags
    // several hundred, overwhelmingly short-lived species that contribute nothing, and an
    // unbounded list buries the answer it was meant to annotate. The count is the signal; a
    // handful of names makes it concrete.
    constexpr std::size_t kMaxNamed = 8;
    const std::size_t named = std::min(total, kMaxNamed);
    out << ":\n    ";
    for (std::size_t i = 0; i < named; ++i) {
      out << (i == 0 ? "" : ", ") << context.unmodeledContinuum[i];
    }
    if (total > named) {
      out << ", and " << (total - named) << " more";
    }
    out << "\n    (see `nusift data info` for the store's photon coverage)\n";
  }
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
    out << csvField(window.label) << ',' << window.id.key << ',' << exact(window.startSeconds)
        << ',' << exact(window.endSeconds) << ',' << exact(window.peakFraction) << '\n';
  }
}

void writeCsvRows(std::ostream& out, const Ranking& ranking, bool withHeader) {
  if (withHeader) {
    out << "time_s,time_end_s,rank,contributor,key,value,unit,fraction,cumulative_fraction,"
           "flags,pinned\n";
  }
  for (const Contributor& c : ranking.contributors) {
    out << exact(ranking.time) << ',';
    if (ranking.domain == Domain::Interval) {
      out << exact(ranking.timeEnd);
    }
    // The unit is not quoted: it comes from a closed enum of spellings that contain no comma,
    // so unlike a label it cannot acquire one.
    //
    // `pinned` is last so that adding it did not renumber the columns anyone already reads by
    // position, and it is here at all because without it a loaded table cannot tell a row that
    // placed from one that was fetched from below the cut -- which is the difference between a
    // top-N and a top-N plus an aside.
    out << ',' << c.rank << ',' << csvField(c.label) << ',' << c.id.key << ',' << exact(c.value)
        << ',' << unitName(ranking.unit) << ',' << exact(c.fraction) << ','
        << exact(c.cumulativeFraction) << ',' << c.flags << ',' << (c.pinned ? 1 : 0) << '\n';
  }
}

void writeJsonRanking(std::ostream& out, const Ranking& ranking, const ReportContext& context,
                      int indent) {
  const std::string pad(static_cast<std::size_t>(indent), ' ');
  out << pad << "{\n";
  out << pad << "  \"metric\": \"" << metricName(ranking.metric) << "\",\n";
  out << pad << "  \"aggregate\": \"" << aggregateName(ranking.aggregate) << "\",\n";
  out << pad << "  \"unit\": \"" << unitName(ranking.unit) << "\",\n";
  out << pad << "  \"time_s\": " << jsonNumber(ranking.time) << ",\n";
  if (ranking.domain == Domain::Interval) {
    out << pad << "  \"time_end_s\": " << jsonNumber(ranking.timeEnd) << ",\n";
  }
  out << pad << "  \"total\": " << jsonNumber(ranking.total) << ",\n";
  out << pad << "  \"covered_fraction\": " << jsonNumber(ranking.coveredFraction) << ",\n";
  out << pad << "  \"omitted_count\": " << ranking.omittedCount << ",\n";
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
    case ReportFormat::Text:
      writeForecastText(out, windows, tracks, table, context);
      break;
    case ReportFormat::Csv:
      writeForecastCsv(out, windows);
      break;
    case ReportFormat::Json:
      writeForecastJson(out, windows, tracks, table);
      break;
  }
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

}  // namespace nusift
