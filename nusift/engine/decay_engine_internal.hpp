#pragma once
/**
 * @file
 * @brief The pruned, restricted decay problem shared by the forward and adjoint engines.
 * @ingroup engine
 */
//
// Like nucdata/nuclear_data_internal.hpp, this exposes cram and Eigen types and is therefore
// NOT installed. It exists because the forward solve and the adjoint solve must agree on the
// index space down to the ordering: an importance vector indexed differently from the
// inventory it multiplies is not an error any assertion would catch, it is a silently
// permuted answer. One prepare() shared between them is what makes that impossible rather
// than merely unlikely.
//
#include <Eigen/SparseCore>
#include <cstdint>
#include <vector>

#include "nusift/engine/decay_engine.hpp"
#include "nusift/engine/inventory.hpp"

namespace nusift {

class NuclearData;

namespace engine_internal {

// The seed, the nuclides forward-reachable from it, and the decay matrix restricted to them.
struct Prepared {
  std::vector<int> keep;                // chain indices retained, ascending
  std::vector<std::int64_t> keys;       // their ZAI keys, same order
  Eigen::SparseMatrix<double> reduced;  // A over `keep`
  Eigen::VectorXd seed;                 // n0 over `keep`
};

// Seed, prune, and restrict. Throws InputError for a nuclide the store does not carry and for
// an inventory that reaches nothing.
Prepared prepare(const NuclearData& data, const Inventory& inventory, const DecayOptions& options);

// [[A, 0], [I, 0]] -- the augmented generator whose exponential carries the inventory in its
// top block and the exact cumulative integral in its bottom.
Eigen::SparseMatrix<double> augment(const Eigen::SparseMatrix<double>& a);

}  // namespace engine_internal
}  // namespace nusift
