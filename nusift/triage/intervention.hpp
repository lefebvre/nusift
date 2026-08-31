#pragma once
/**
 * @file
 * @brief What removing something on a given date is worth to a later response.
 * @ingroup triage
 */
//
// Attribution says which seeded nuclide a response is riding on. This asks the question a
// process engineer arrives with instead: if I take something OUT on a given date, what is that
// worth later. What does a Cs/Sr separation before storage actually buy against the dose rate
// at thirty years -- in the tool's own currency, rather than as a rule of thumb.
//
// THE WHOLE THING IS ONE ADJOINT SOLVE.
//
// Decay is linear, so R(T) = <g, n(t0)> where g_i = dR(T)/dn_i(t0) is the importance of an atom
// of nuclide i present at the intervention date. The adjoint delivers every g_i from a single
// solve over the window [t0, T] -- see engine/adjoint_engine.hpp for the duality identity that
// makes that true. Once g is in hand, removing a fraction f of nuclide i costs the response
// exactly f * n_i(t0) * g_i, and every intervention on the list is a dot product against the
// same vector. Comparing a dozen processing schedules is a dozen dot products, not a dozen
// solves.
//
// It is EXACT in the sense the methodology docs use, for the same reason the attribution shares
// are: no search, no perturbation, no finite difference, and no tolerance. What is not exact is
// the premise -- that the removal is instantaneous and complete to the stated fraction.
//
// WHAT REMOVING A PARENT DOES, AND DOES NOT DO.
//
// Removing a nuclide at t0 removes the atoms of THAT nuclide and, with them, everything they
// would have gone on to produce after t0. That second part is not an approximation bolted on:
// the importance g_i already carries the whole forward evolution from t0 to T, so a parent's
// value includes the daughters it would have fed.
//
// What it does not do is remove the daughters already present at t0. Take out caesium and the
// barium standing in the drum at that instant stays, because it is barium. That is exactly what
// a chemical separation does, and stating it matters: the benefit of removing Cs-137 an hour
// before the response is nearly nothing, since the Ba-137m doing the emitting is already there
// and is not caesium.
//
// NOT MODELLED HERE: separation efficiency as a function of chemistry, where the removed
// material goes, the dose incurred DURING the processing, and any interval response -- the
// shares are of an instantaneous R(T), for the same reason attribution refuses intervals.
//
#include <cstdint>
#include <span>
#include <string>
#include <vector>

#include "nusift/engine/decay_engine.hpp"
#include "nusift/engine/inventory.hpp"
#include "nusift/triage/response.hpp"

namespace nusift {

class NuclearData;

// One thing taken out, and how much of it.
struct Removal {
  // What to remove, named the way a user writes it:
  //
  //   "Cs-137"  one nuclide
  //   "Cs"      every isotope of an element -- what a chemical separation actually removes
  //   "Z=55"    the same, by atomic number
  //   "A=137"   every nuclide on a mass chain
  //
  // An element form is the physically meaningful one for processing: a separation cannot pick
  // one isotope out of another, so "remove Cs-137 but leave Cs-134" is not a process that
  // exists. Both are allowed because the nuclide form is still the right way to ask what a
  // single nuclide is worth.
  std::string selector;

  // The fraction taken out, 0 to 1. One is a complete separation, which is an upper bound on
  // any real one rather than a description of one.
  double fraction = 1.0;
};

// One counterfactual: a set of removals applied at the same instant.
struct Intervention {
  std::string name;  // "Cs/Sr separation" -- named because naming the winner is the answer
  std::vector<Removal> removals;
};

// What one removal cost the response, per nuclide, so the benefit can be read rather than
// merely totalled.
struct RemovedContributor {
  std::int64_t key = 0;
  std::string label;
  double atomsRemoved = 0.0;
  // The response this nuclide's removed atoms were carrying: atomsRemoved * dR(T)/dn_i(t0).
  double value = 0.0;
  double fraction = 0.0;  // of everything this intervention removed
};

struct InterventionEffect {
  std::string name;

  // R(T) once the intervention has been applied.
  double response = 0.0;
  // What it bought: baseline - response. Never negative -- removing atoms from a non-negative
  // weight cannot raise the response.
  double removed = 0.0;
  double removedFraction = 0.0;  // of the baseline

  // The nuclides the benefit came from, most first.
  std::vector<RemovedContributor> contributors;
};

struct InterventionStudy {
  Metric metric = Metric::Activity;
  Unit unit = Unit::Becquerel;

  double interventionTimeSeconds = 0.0;
  double responseTimeSeconds = 0.0;

  // R(T) with nothing removed. Equal, to rounding, to the total an ordinary ranking reports for
  // the same spec and time -- the duality identity guarantees it and a test asserts it.
  double baseline = 0.0;

  std::string seedProvenance;
  std::vector<InterventionEffect> effects;
};

// What each intervention is worth to R at `responseTime`, all from one adjoint solve over
// [`interventionTime`, `responseTime`].
//
// Every intervention is applied to the SAME inventory at t0 -- they are alternatives being
// compared, not a sequence being accumulated. A schedule of removals at different dates is
// enumeration over several calls, which is what the linearity makes cheap.
//
// Throws InputError for a response time before the intervention time, a negative or non-finite
// time, an unnamed intervention, a fraction outside [0, 1], a selector that names nothing the
// inventory's chain reaches, and for a selector named twice within one intervention -- removing
// half of something twice is not a stated quantity, and guessing which of 75% or 100% was meant
// would be worse than refusing.
//
// A selector that resolves but finds no atoms at t0 is NOT an error: "there is no caesium left
// to remove by then" is a real answer to a question someone actually asked.
InterventionStudy compareInterventions(const NuclearData& data, const Inventory& inventory,
                                       double interventionTime, double responseTime,
                                       const ResponseSpec& spec,
                                       std::span<const Intervention> interventions,
                                       const DecayOptions& options = {});

}  // namespace nusift
