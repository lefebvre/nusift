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
#include <span>
#include <string>
#include <vector>

#include "nusift/triage/allowable.hpp"
#include "nusift/triage/attribution.hpp"
#include "nusift/triage/events.hpp"
#include "nusift/triage/forecast.hpp"
#include "nusift/triage/intervention.hpp"
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
  // How the exposure -- or a photon fluence -- was computed, when the metric uses the point
  // geometry. A figure with no stated distance is not interpretable, so this rides in the
  // header rather than being left to the reader to remember from the command line. Empty for
  // activity and for photon strength, which uses no geometry at all.
  std::string geometry;
  // Contributors whose photon output is partly in a continuum NuSIFT does not model. Named
  // in a footnote so an understated row is visible rather than merely flagged in a column
  // nobody reads.
  std::vector<std::string> unmodeledContinuum;
};

// Render one ranking.
void writeRanking(std::ostream& out, const Ranking& ranking, const ReportContext& context,
                  ReportFormat format);

// Render a seed attribution: which seeded nuclides the response is riding on. The importance
// column is what distinguishes it from an ordinary ranking -- a large share can come from a
// large seed or from a potent one, and only the two columns together say which.
void writeAttribution(std::ostream& out, const SeedAttribution& attribution,
                      const ReportContext& context, ReportFormat format);

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

// What a located-event report is about: which curve was searched, in what unit, and against
// what level. Bundled rather than passed as six parameters because every one of them is needed
// to read the numbers -- a crossing time with no curve and no level named is not an answer.
struct EventReport {
  std::string curve;  // "total", a contributor's label, or "A / B" for a ratio
  std::string unit;   // the curve's unit; empty for a ratio, which is dimensionless
  std::string metric;
  double gridStartSeconds = 0.0;
  double gridEndSeconds = 0.0;
  int gridPoints = 0;

  // Whether a level was asked about at all. Without one only the turns are reported, and
  // printing a level of zero would look like one that was.
  bool hasLevel = false;
  double level = 0.0;

  std::vector<TrajectoryEvent> events;
  std::vector<LevelWindow> windows;
};

// Render located events: crossings and turns, each with the bracket that found it, and the
// windows the crossings pair into.
//
// The bracket travels into every format because it is the honest error bar on the instant --
// a consumer that loads a crossing time without knowing how well it is placed has lost the
// only thing distinguishing a refined answer from a grid artefact. CSV carries events and
// windows in one table under a `kind` column, which is what a spreadsheet can actually load.
void writeEvents(std::ostream& out, const EventReport& report, const ReportContext& context,
                 ReportFormat format);

// Render the maximum allowable scale over time: what binds, what it permits, and which
// contributors drive it.
//
// A time where nothing constrains the inventory is printed as such rather than as a very large
// number, in every format -- JSON gives it a null scale and a null binding criterion, which is
// the only encoding a parser cannot mistake for a bound of zero.
void writeAllowable(std::ostream& out, const std::vector<AllowableScale>& scaled,
                    std::span<const Criterion> criteria, const ReportContext& context,
                    ReportFormat format);

// Render a set of counterfactual interventions: what each one is worth against the same
// baseline, and which nuclides the benefit came from.
//
// The interventions are ALTERNATIVES compared against one baseline, not a sequence, and the
// text form says so -- a reader who added two rows together would be describing a schedule
// nobody computed.
void writeInterventions(std::ostream& out, const InterventionStudy& study,
                        const ReportContext& context, ReportFormat format);

}  // namespace nusift
