#include "nusift/engine/adjoint_engine.hpp"

#include <Eigen/SparseCore>
#include <cmath>
#include <string>
#include <vector>

#include "cram/adjoint.hpp"
#include "nusift/core/error.hpp"
#include "nusift/engine/decay_engine_internal.hpp"
#include "nusift/nucdata/nuclear_data.hpp"

// The third translation unit that sees cram and Eigen, after nucdata/nuclear_data.cpp and
// engine/decay_engine.cpp, and the only one that includes cram/adjoint.hpp. Keeping the
// adjoint header out of decay_engine.cpp is deliberate: adjoint.hpp pulls in cram's burnup
// API (deplete.hpp, integrator.hpp, reaction.hpp), and the decay engine has no business
// depending on a depletion system it never builds.

namespace nusift {
namespace {

constexpr const char* kModule = "adjoint engine";

cram::CramOrder toCramOrder(CramOrder order) {
  return order == CramOrder::Order16 ? cram::CramOrder::CRAM16 : cram::CramOrder::CRAM48;
}

}  // namespace

SeedImportance seedImportance(const NuclearData& data, const Inventory& inventory,
                              std::span<const double> weight, double time,
                              const DecayOptions& options) {
  if (!std::isfinite(time) || time < 0.0) {
    throw InputError(tagged(
        kModule, "the response time must be non-negative and finite; got " + std::to_string(time)));
  }
  if (static_cast<int>(weight.size()) != data.size()) {
    throw InputError(tagged(kModule, "the weight vector has " + std::to_string(weight.size()) +
                                         " entries but the store carries " +
                                         std::to_string(data.size()) +
                                         " nuclides; it must be indexed by store index"));
  }

  const engine_internal::Prepared prepared = engine_internal::prepare(data, inventory, options);
  const int m = static_cast<int>(prepared.keys.size());

  Eigen::VectorXd w = Eigen::VectorXd::Zero(m);
  for (int k = 0; k < m; ++k) {
    const int global = prepared.keep[static_cast<std::size_t>(k)];
    w(k) = weight[static_cast<std::size_t>(global)];
  }

  SeedImportance out;
  out.nuclideKeys = prepared.keys;
  out.importance.assign(static_cast<std::size_t>(m), 0.0);
  out.seedAtoms.assign(static_cast<std::size_t>(m), 0.0);

  // exp(A^T * 0) = I, so at t = 0 the importance IS the weight. Solving would be harmless but
  // wasteful, and it is the same shortcut decay() takes on its forward path.
  const Eigen::VectorXd nStar =
      time > 0.0 ? cram::cramSolveAdjoint(prepared.reduced, w, time, toCramOrder(options.order))
                 : w;

  for (int k = 0; k < m; ++k) {
    out.importance[static_cast<std::size_t>(k)] = nStar(k);
    out.seedAtoms[static_cast<std::size_t>(k)] = prepared.seed(k);
  }
  out.response = nStar.dot(prepared.seed);
  return out;
}

}  // namespace nusift
