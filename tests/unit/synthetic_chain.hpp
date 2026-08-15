#pragma once
//
// Builders for synthetic nuclear data, and the closed-form Bateman solutions to check the
// engine against.
//
// The engine is verified against analytic solutions rather than against a reference
// implementation or a stored baseline. A golden file would only prove the code still does
// what it did yesterday; the Bateman equations prove it does what the physics says. That
// matters most for the augmented-matrix time integral, which has no independent
// implementation anywhere to compare against.
//
// Every chain here is built through StoreArrays, the same path a real HDF5 store takes, so
// these tests exercise the production construction code rather than a test-only shortcut.
//
#include <cmath>
#include <cstdint>
#include <vector>

#include "nusift/nucdata/store_arrays.hpp"
#include "nusift/units.hpp"

// Namespace is `synth`, not `testing`: a nusift::testing would shadow GoogleTest's global
// ::testing everywhere these helpers are used inside namespace nusift, and the resulting
// errors point at gtest internals rather than at the collision.
namespace nusift::synth {

// Beta-minus: Z+1 at constant A. The simplest way to build a linear chain, since successive
// members differ only in Z and cannot collide with each other.
inline constexpr double kBetaMinus = 1.0;

inline double halfLifeFor(double lambda) {
  return units::kLn2 / lambda;
}

// A linear chain Z, Z+1, ... Z+n-1 at fixed A, each decaying beta-minus into the next. The
// last member is stable, terminating the chain. Decay constants are given directly because
// every analytic solution below is written in terms of lambda.
inline StoreArrays linearChain(const std::vector<double>& lambdas, int z0 = 50, int a = 100) {
  const int n = static_cast<int>(lambdas.size()) + 1;  // + stable terminator
  StoreArrays s;
  s.provenance.version = 1;
  s.nuclideKey.reserve(n);
  s.halfLife.reserve(n);
  for (int i = 0; i < n; ++i) {
    s.nuclideKey.push_back(Zai{z0 + i, a, 0}.key());
    const bool stable = i == n - 1;
    s.halfLife.push_back(stable ? 0.0 : halfLifeFor(lambdas[static_cast<std::size_t>(i)]));
  }
  s.modeOffset.assign(static_cast<std::size_t>(n) + 1, 0);
  for (int i = 0; i < n; ++i) {
    const bool stable = i == n - 1;
    if (!stable) {
      s.modeRtyp.push_back(kBetaMinus);
      s.modeBranching.push_back(1.0);
      s.modeFinalState.push_back(0);
      s.modeIsFission.push_back(0);
    }
    s.modeOffset[static_cast<std::size_t>(i) + 1] = static_cast<int>(s.modeRtyp.size());
  }
  return s;
}

// One nuclide with two decay branches of the given branching fractions: beta-minus to Z+1
// and electron capture to Z-1. Used to check that branching splits production correctly and
// conserves atoms.
inline StoreArrays branchingChain(double lambda, double branchUp, double branchDown, int z0 = 50,
                                  int a = 100) {
  StoreArrays s;
  s.provenance.version = 1;
  // Sorted ascending by key: Z-1, Z, Z+1.
  s.nuclideKey = {Zai{z0 - 1, a, 0}.key(), Zai{z0, a, 0}.key(), Zai{z0 + 1, a, 0}.key()};
  s.halfLife = {0.0, halfLifeFor(lambda), 0.0};
  s.modeOffset = {0, 0, 2, 2};
  s.modeRtyp = {kBetaMinus, 2.0};  // 2.0 = electron capture, Z-1
  s.modeBranching = {branchUp, branchDown};
  s.modeFinalState = {0, 0};
  s.modeIsFission = {0, 0};
  return s;
}

// Attach discrete photon lines to nuclide `index`. Only the CSR bookkeeping is interesting
// here; the exposure physics that consumes them is checked separately.
inline void addLines(StoreArrays& s, int index, const std::vector<double>& energiesEv,
                     const std::vector<double>& intensities) {
  const int n = s.nuclideCount();
  if (s.lineOffset.empty()) {
    s.lineOffset.assign(static_cast<std::size_t>(n) + 1, 0);
  }
  const int at = s.lineOffset[static_cast<std::size_t>(index)];
  for (std::size_t k = 0; k < energiesEv.size(); ++k) {
    s.lineEnergyEv.insert(s.lineEnergyEv.begin() + at + static_cast<int>(k), energiesEv[k]);
    s.lineIntensity.insert(s.lineIntensity.begin() + at + static_cast<int>(k), intensities[k]);
    s.lineStyp.insert(s.lineStyp.begin() + at + static_cast<int>(k), 0);
  }
  for (int i = index + 1; i <= n; ++i) {
    s.lineOffset[static_cast<std::size_t>(i)] += static_cast<int>(energiesEv.size());
  }
}

// --- closed-form solutions -------------------------------------------------
//
// Bateman for a linear chain seeded entirely in member 0. Written with distinct decay
// constants only; the degenerate case lambda_i == lambda_j has a different closed form and
// none of the tests use it.

// N_0(t) = N_0(0) e^{-lambda_0 t}
inline double batemanN0(double n0, double lambda, double t) {
  return n0 * std::exp(-lambda * t);
}

// Two-member chain, second member.
inline double batemanN1(double n0, double l0, double l1, double t) {
  return n0 * l0 / (l1 - l0) * (std::exp(-l0 * t) - std::exp(-l1 * t));
}

// Three-member chain, third member.
inline double batemanN2(double n0, double l0, double l1, double l2, double t) {
  const double a = std::exp(-l0 * t) / ((l1 - l0) * (l2 - l0));
  const double b = std::exp(-l1 * t) / ((l0 - l1) * (l2 - l1));
  const double c = std::exp(-l2 * t) / ((l0 - l2) * (l1 - l2));
  return n0 * l0 * l1 * (a + b + c);
}

// \int_0^t N_0 dtau, in atom-seconds. The integral of a decaying exponential.
inline double batemanIntegralN0(double n0, double lambda, double t) {
  return n0 * (1.0 - std::exp(-lambda * t)) / lambda;
}

// \int_0^t N_1 dtau. Term-by-term integral of batemanN1; this is what the augmented matrix's
// bottom block must reproduce, and it is the only independent check on that mechanism.
inline double batemanIntegralN1(double n0, double l0, double l1, double t) {
  const double first = (1.0 - std::exp(-l0 * t)) / l0;
  const double second = (1.0 - std::exp(-l1 * t)) / l1;
  return n0 * l0 / (l1 - l0) * (first - second);
}

}  // namespace nusift::synth
