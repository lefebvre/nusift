#pragma once
/**
 * @file
 * @brief Rendering rankings as a human table, CSV, or JSON.
 * @ingroup io
 */
//
// Every report carries its provenance in the header -- which store, which evaluation, which
// inventory, what coverage the shown rows account for. For a triage tool that is not
// decoration: the answer to "which isotopes dominate" is only as good as the evaluation it
// came from, and a table of bare numbers with no attribution cannot be checked by anyone
// later.
//
#include <iosfwd>
#include <string>
#include <vector>

#include "nusift/triage/forecast.hpp"
#include "nusift/triage/ranking.hpp"

namespace nusift {

class NuclearData;

enum class ReportFormat {
  Text,  // aligned table for a terminal
  Csv,   // for a spreadsheet
  Json,  // for a script
};

bool parseReportFormat(std::string_view text, ReportFormat& out);

struct ReportContext {
  std::string storePath;
  std::string storeLibrary;  // e.g. "ENDF/B-VIII.1"
  std::string storeCreatedUtc;
  int storeNuclideCount = 0;
  std::string seedProvenance;  // where the inventory came from
  // How the exposure was computed, when the metric is exposure. An exposure figure with no
  // stated distance is not interpretable, so this rides in the header rather than being left
  // to the reader to remember from the command line.
  std::string geometry;
  // Contributors whose photon output is partly in a continuum NuSIFT does not model. Named
  // in a footnote so an understated row is visible rather than merely flagged in a column
  // nobody reads.
  std::vector<std::string> unmodeledContinuum;
};

// Render one ranking.
void writeRanking(std::ostream& out, const Ranking& ranking, const ReportContext& context,
                  ReportFormat format);

// Render a dominance forecast: who leads over which windows, and the contributors that reach
// the top at any point. `tracks` may be empty, in which case only the windows are shown.
void writeForecast(std::ostream& out, const std::vector<DominanceWindow>& windows,
                   const std::vector<RankTrack>& tracks, const ResponseTable& table,
                   const ReportContext& context, ReportFormat format);

// Render a ranking per time, as produced by rankAll. Text output separates them with
// headings; CSV and JSON emit one flat table with a time column, which is what a consumer
// wants to load.
void writeRankings(std::ostream& out, const std::vector<Ranking>& rankings,
                   const ReportContext& context, ReportFormat format);

// The same, when the rankings do NOT share a context. Integrating several intervals solves
// each one separately, so each has its own set of flagged emitters; one context for all of
// them footnotes every ranking with the last interval's list, naming nuclides that need not
// appear in the ranking above at all. `contexts` must be the same length as `rankings`.
void writeRankings(std::ostream& out, const std::vector<Ranking>& rankings,
                   const std::vector<ReportContext>& contexts, ReportFormat format);

}  // namespace nusift
