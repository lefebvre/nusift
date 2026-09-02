#pragma once
/**
 * @file
 * @brief Reading and writing isotopic inventories.
 * @ingroup io
 */
//
// CSV is the canonical inventory format, because the people who have inventories have them
// in spreadsheets. The format is three columns with a tolerant reader:
//
//   # comments, blank lines, and a UTF-8 BOM are all accepted
//   nuclide, quantity, unit
//   Cs-137,  1.2e14,   Bq
//   Sr-90,   3.5,      g
//   Pu-239,  0.8,      Ci
//
// A header row naming the columns is optional and detected rather than required. JSON is the
// secondary format: it round-trips provenance, which CSV cannot carry, and it is what the
// eventual Python layer will hand back and forth.
//
// A FOURTH COLUMN CARRIES THE DATE THE ROW WAS MEASURED, which is what a real assay sheet has
// and what makes several sheets combinable:
//
//   nuclide, quantity, unit, assayed
//   Cs-137,  1.2e14,   Bq,   2024-03-15
//   Sr-90,   3.5,      g,    2023-11-02
//
// Rows sharing a date are one assay. Reading such a file carries every assay forward to the
// latest date and merges them -- see reconcile.hpp for why forward and only forward -- so what
// comes back is an ordinary Inventory describing the material at one instant, which is the only
// thing the rest of NuSIFT can be handed.
//
// A file either dates every row or dates none. A MIXED file is refused rather than defaulted,
// because there is no reading of an undated row among dated ones that is not a guess: treating
// it as measured at the epoch silently ages it by however long the others were carried, and
// treating it as measured at the earliest date silently does the opposite. Neither is visible in
// the answer, so neither is offered.
//
#include <iosfwd>
#include <string>
#include <vector>

#include "nusift/engine/inventory.hpp"
#include "nusift/engine/reconcile.hpp"

namespace nusift {

class NuclearData;

struct InventoryReadOptions {
  // Report and skip rows naming a nuclide the data store does not carry, rather than
  // failing. Off by default: a silently dropped row understates every ranking that follows,
  // with nothing in the output to say so.
  bool ignoreUnknown = false;
  // Where skipped rows are reported. Null suppresses the warnings entirely.
  std::ostream* warnings = nullptr;
};

// The assays a file describes, before anything is carried anywhere. Exactly one group, with
// `dated` false, for a file that names no dates.
//
// This is the form to read when the reconciliation itself is the answer -- which sheet
// contributed what, and how far each was carried. Everything else wants readInventory().
struct DatedInventory {
  std::vector<AssayGroup> groups;  // ascending by date
  bool dated = false;
};

// Read an inventory, choosing the parser by file extension (.json for JSON, anything else
// CSV). Throws InputError naming the file and line for a malformed row.
//
// A dated file is reconciled to its latest assay date on the way out, with DEFAULT solver
// options -- deliberately not the caller's. A file has to mean the same inventory whichever
// command reads it, and a seed that shifted with the CRAM order the caller happened to pass
// would be a different inventory from the same sheet.
Inventory readInventory(const std::string& path, const NuclearData& data,
                        const InventoryReadOptions& options = {});

// The same file, read as assays rather than as one inventory.
DatedInventory readInventoryDated(const std::string& path, const NuclearData& data,
                                  const InventoryReadOptions& options = {});
DatedInventory readInventoryDatedCsv(std::istream& in, const NuclearData& data,
                                     const std::string& sourceName,
                                     const InventoryReadOptions& options = {});
DatedInventory readInventoryDatedJson(std::istream& in, const NuclearData& data,
                                      const std::string& sourceName,
                                      const InventoryReadOptions& options = {});

// Parse CSV from an already-open stream. `sourceName` appears in error messages, so it
// should be the file path when there is one.
Inventory readInventoryCsv(std::istream& in, const NuclearData& data, const std::string& sourceName,
                           const InventoryReadOptions& options = {});

Inventory readInventoryJson(std::istream& in, const NuclearData& data,
                            const std::string& sourceName,
                            const InventoryReadOptions& options = {});

// Write in `unit`. Nuclides that cannot be expressed in it -- an activity unit for a stable
// nuclide, a mass unit with no staged atomic weight -- fall back to atoms for that row and
// are marked as such, rather than aborting the write or emitting a wrong number.
void writeInventoryCsv(std::ostream& out, const Inventory& inventory, const NuclearData& data,
                       Quantity unit = Quantity::Atoms);

void writeInventoryJson(std::ostream& out, const Inventory& inventory, const NuclearData& data,
                        Quantity unit = Quantity::Atoms);

}  // namespace nusift
