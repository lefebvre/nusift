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
#include <iosfwd>
#include <string>

#include "nusift/engine/inventory.hpp"

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

// Read an inventory, choosing the parser by file extension (.json for JSON, anything else
// CSV). Throws InputError naming the file and line for a malformed row.
Inventory readInventory(const std::string& path, const NuclearData& data,
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
