#pragma once
/**
 * @file
 * @brief Assays taken on different dates, carried to one epoch and merged.
 * @ingroup engine
 */
//
// An Inventory is atoms at ONE instant. That is not a limitation of the type, it is what makes
// every downstream answer well posed: a ranking, an integral and an adjoint share all assume the
// vector they were handed describes the material at a single moment. So an inventory assembled
// from several assays taken on different dates cannot simply be added up. The rows have to be
// brought to a common instant first, and bringing them there is a decay solve rather than a
// bookkeeping step.
//
// That is the whole of this file:
//
//     n(epoch) = sum over assays a of  decay(n_a, epoch - date_a)
//
// Each assay is carried forward by its own gap and the results are summed, which is exact
// because the decay operator is linear -- the same linearity every other exactness claim in
// NuSIFT rests on. Nothing is approximated here that is not already approximated in a solve.
//
// FORWARD ONLY, AND THIS IS THE LOAD-BEARING RULE. Carrying an assay forward is well posed: the
// atoms present later are determined by the atoms present now. Carrying one BACKWARD is not.
// The daughters measured at assay have two indistinguishable histories -- present from the
// start, or grown in since -- and no measurement of the mixture at one instant separates them.
// Inverting the decay operator produces a vector that reproduces the measurement and is not the
// composition that existed, and it will contain negative atom counts as soon as the data has any
// noise in it at all. An epoch earlier than any assay is therefore REFUSED rather than caveated:
// a caveat on a number nobody can check is not a warning, it is a disclaimer.
//
#include <span>
#include <string>
#include <vector>

#include "nusift/engine/decay_engine.hpp"
#include "nusift/engine/inventory.hpp"

namespace nusift {

class NuclearData;

// One assay: what was measured, and when it was measured.
//
// The date is seconds on a shared timeline -- parseCalendarDate() produces one from an ISO date,
// and only differences between assays reach any answer, so which origin the timeline uses does
// not matter as long as one file uses one origin.
struct AssayGroup {
  double dateSeconds = 0.0;
  Inventory inventory;
  // Where this assay came from, for the report. A reconciliation that cannot say which sheet
  // contributed what has merged away the only thing making it auditable.
  std::string label;

  // Whether the file actually stated a date. An undated sheet gets dateSeconds = 0, which is a
  // PLACEHOLDER and not a measurement -- so carrying it to any epoch but zero would be carrying
  // it from a date nobody supplied. Without this flag that mistake is silent and large: an
  // undated sheet reconciled to 2024 would be aged fifty-four years, and every number after it
  // would be wrong with nothing to show why.
  bool dated = false;
};

// What one assay contributed to the merged inventory.
struct AssayContribution {
  std::string label;
  double dateSeconds = 0.0;
  // How far this assay was carried to reach the epoch. Zero for the assay that defines it.
  double carriedSeconds = 0.0;

  int nuclides = 0;  // rows in the assay as measured
  double atomsAtAssay = 0.0;
  // Atoms this assay accounts for at the epoch. NOT equal to atomsAtAssay: decay moves atoms
  // between nuclides but a chain ending in a stable nuclide conserves their number, while one
  // whose closure was pruned does not. Both are reported so the difference is visible rather
  // than absorbed.
  double atomsAtEpoch = 0.0;
};

struct Reconciliation {
  // Everything, carried to `epochSeconds` and summed. This is the inventory every other part of
  // NuSIFT consumes, and it is an ordinary one -- the dates are spent in producing it and do not
  // travel any further, because at the epoch there is only one date left.
  Inventory inventory;
  double epochSeconds = 0.0;

  // Per assay, in date order.
  std::vector<AssayContribution> contributions;

  // The span the assays were taken over. A reconciliation across a wide span is not wrong, but
  // it is carrying an old measurement a long way on nothing but the decay model, and how far is
  // something a reader should be told rather than have to compute from the rows.
  double spanSeconds = 0.0;
};

// The epoch reconcile() uses when none is named: the LATEST assay date, which is the only choice
// that carries every assay forward and none backward. Throws InputError on an empty span.
double latestAssayDate(std::span<const AssayGroup> groups);

// Carry every assay to `epochSeconds` and merge.
//
// Throws InputError for an empty group list, an assay with no rows, or an epoch earlier than any
// assay date -- see the refusal above. An epoch equal to an assay's date costs that assay no
// solve at all, which is the ordinary case for the latest one.
Reconciliation reconcile(const NuclearData& data, std::span<const AssayGroup> groups,
                         double epochSeconds, const DecayOptions& options = {});

}  // namespace nusift
