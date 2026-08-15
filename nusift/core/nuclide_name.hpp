#pragma once
/**
 * @file
 * @brief Parsing and formatting of human-written nuclide names ("Cs-137", "cs137", "Am-242m").
 * @ingroup core
 */
//
// Every nuclide a user types -- in an inventory CSV, a --seed-fission argument, a run
// config -- arrives as a string, and the spellings in the wild are inconsistent. This is
// the single place that reconciles them, so the CLI, the config reader, the CSV reader,
// and the Python bindings cannot disagree about what "Am242M" means.
//
// Accepted forms, case-insensitive, with optional surrounding whitespace:
//
//   Cs-137   cs137   CS-137        ground state
//   Am-242m  am242M  Am242m1       first metastable state
//   Am-242m2                       second metastable state
//   U-235    92-U-235  922350      element-prefixed and raw-key forms
//
// The raw-key form exists because store diagnostics and JSON output emit keys, and being
// able to paste one straight back into the CLI is worth the eight lines it costs.
//
// cram grows its own parseNuclideName() on the burnup branch. The semantics here are
// deliberately a superset, so adopting it means forwarding and handling whatever cram's
// parser rejects, rather than reworking call sites.
//
#include <optional>
#include <string>
#include <string_view>

#include "nusift/core/nuclide.hpp"

namespace nusift {

// Parse a nuclide name. Returns nullopt for anything unrecognised -- including a
// well-formed name with an impossible Z or A, such as "Xx-137" or "Cs-0".
std::optional<Zai> parseNuclideName(std::string_view name);

// As above, but raises InputError naming the offending token instead of returning nullopt.
// Preferred at input boundaries, where a bad name is a user error that must be reported
// with the text the user actually wrote. `context` is prepended when non-empty, e.g.
// "inventory.csv line 12".
Zai requireNuclideName(std::string_view name, std::string_view context = {});

// Canonical name: "Cs-137", "Am-242m", "Am-242m2". This is the spelling NuSIFT emits
// everywhere -- reports, JSON, labels -- regardless of how the nuclide was written on input,
// so two runs of the same inventory produce byte-identical output.
std::string formatNuclideName(const Zai& zai);

}  // namespace nusift
