#pragma once
/**
 * @file
 * @brief Rendering a binned photon source: as a report, or as a transport code's source card.
 * @ingroup io
 */
//
// Two of the five formats here are not reports at all -- they are input decks for someone
// else's code -- and that is what makes this file its own rather than three more cases inside
// report.cpp. A ranking is read by a person who can weigh a footnote; a source deck is read by
// MCNP, which cannot, and whatever the deck does not say is simply not true of the run.
//
// So the decks carry their caveats as comment cards, at the top, in the imperative: what the
// number is, what has to be multiplied by what, and what the source is missing. A deck that
// merely got the arithmetic right and left those out would be the more dangerous artefact of
// the two, because it would look finished.
//
#include <iosfwd>
#include <string_view>

#include "nusift/io/report.hpp"
#include "nusift/triage/spectrum.hpp"

namespace nusift {

enum class SourceFormat {
  Text,  // aligned histogram for a terminal
  Csv,   // one row per bin
  Json,
  McnpSdef,      // an SDEF card with a histogram ERG distribution
  OpenmcPython,  // a Python snippet building an openmc.IndependentSource
};

bool parseSourceFormat(std::string_view text, SourceFormat& out);

// Render `spectrum` in `format`. `context` supplies the provenance every format carries: which
// store, which evaluation, which inventory. The deck formats carry it as comments, because a
// source card that cannot be traced back to the inventory and evaluation it came from is not
// reproducible by anyone, including the person who wrote it.
void writeSourceSpectrum(std::ostream& out, const BinnedSpectrum& spectrum,
                         const ReportContext& context, SourceFormat format);

}  // namespace nusift
