#pragma once
/**
 * @file
 * @brief Importance of each seeded nuclide to a weighted response, from one adjoint solve.
 * @ingroup engine
 */
//
// THE ADJOINT, AND WHY IT IS A PARTITION RATHER THAN A SENSITIVITY.
//
// Every NuSIFT metric is R = <w, n(T)> for a fixed per-nuclide weight w -- activity takes
// w = lambda, exposure takes lambda times a photon sum (see triage/response.hpp). Decay is
// linear, so
//
//     R = <w, exp(AT) n0> = sum_i n0_i * <w, exp(AT) e_i>
//
// and the inner product in that sum is exactly dR/dn0_i, the importance of seed nuclide i.
// The adjoint delivers all of them at once, because <w, exp(AT) n0> = <exp(A^T T) w, n0>:
// ONE solve on the transposed matrix gives a vector whose dot product with the seed is R and
// whose entries are the per-nuclide derivatives.
//
// The consequence is the point. `importance[i] * seedAtoms[i]` is not an estimate of nuclide
// i's influence, it is nuclide i's exact share of R, and those shares SUM TO R. Seed
// attribution is therefore a ranking with a total and a covered fraction -- the same shape as
// every other ranking in the tool -- and not an uncertainty analysis that would need evaluated
// sigmas to mean anything.
//
// What it adds over the ordinary ranking is a second attribution axis. The ordinary ranking
// attributes a response to the nuclide that EMITS it, which is why a Cs-137 source's exposure
// correctly lands on Ba-137m. This attributes the same response to the seeded nuclide it came
// FROM. Neither is derivable from the other, and for a fission seed they name largely disjoint
// sets: the emitters are what is present now, the seeds are the short-lived precursors that
// produced them.
//
// COST. One adjoint solve, against one forward solve for the ordinary ranking -- it is the
// cheaper of the two, since the response falls out of the same vector rather than needing a
// second pass. There is no quadrature and no convergence parameter anywhere in this path.
//
#include <cstdint>
#include <span>
#include <vector>

#include "nusift/engine/decay_engine.hpp"
#include "nusift/engine/inventory.hpp"

namespace nusift {

class NuclearData;

struct SeedImportance {
  // The pruned index space, identical to what decay() would produce for this seed -- the two
  // engines share one prepare() precisely so these line up.
  std::vector<std::int64_t> nuclideKeys;

  // dR/dn0_i: how much the response rises per additional atom of nuclide i in the SEED, in
  // whatever unit `weight` was expressed in, per atom.
  std::vector<double> importance;

  // n0_i over the same index space, so a caller can form the share without re-deriving it.
  std::vector<double> seedAtoms;

  // R = <importance, seedAtoms>. Equal, to rounding, to the total the forward path reports for
  // the same weight and time; the duality identity is what guarantees it, and a test asserts
  // it rather than trusting it.
  double response = 0.0;
};

// Importance of every seeded nuclide to R = <weight, n(time)>.
//
// `weight` is indexed by the DATA STORE's index space and must have data.size() entries --
// the same vector triage/response.cpp forms for a metric, which is why responseWeights()
// exists to build it. Entries for nuclides outside the pruned set are ignored, since they are
// identically zero for all time.
//
// Throws InputError for a negative or non-finite time, a weight vector of the wrong length, or
// an inventory naming a nuclide the store does not carry.
SeedImportance seedImportance(const NuclearData& data, const Inventory& inventory,
                              std::span<const double> weight, double time,
                              const DecayOptions& options = {});

}  // namespace nusift
