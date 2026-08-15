#pragma once
/**
 * @file
 * @brief The decay engine's entire output: atoms and their exact time integrals, per nuclide
 *        and per time.
 * @ingroup engine
 */
//
// The engine produces exactly two matrices and nothing else:
//
//   atoms(k, i)            = n_i(t_k)
//   integratedAtoms(k, i)  = \int_0^{t_k} n_i(tau) dtau      [atom.s]
//
// Every metric NuSIFT reports is a linear functional of one of those with a fixed per-nuclide
// or per-line weight. Activity is lambda_i * atoms; exposure is lambda_i * atoms * a photon
// weight; the same weights against integratedAtoms give the time-integrated forms, in decays
// and in roentgen respectively. Keeping the engine's contract this narrow is what makes
// ranking by nuclide, mass chain, element, or gamma line a post-multiply instead of four
// separate code paths -- and it is why the engine never sums anything.
//
// Storage is dense and row-major by time. 200 times x 4000 nuclides x 2 arrays x 8 B is
// 12.8 MB; there is nothing here worth being clever about.
//
#include <cstdint>
#include <span>
#include <string>
#include <vector>

namespace nusift {

struct DecayResult {
  // Seconds since t=0, ascending and unique.
  std::vector<double> times;

  // The index space of the two matrices below. This is the PRUNED chain -- only nuclides
  // forward-reachable from the seed -- so it is generally a subset of the NuclearData index
  // space and must not be used to index into it. Lookups go through the key.
  std::vector<std::int64_t> nuclideKeys;

  std::vector<double> atoms;            // [nT * nNuc], row-major by time
  std::vector<double> integratedAtoms;  // [nT * nNuc], row-major by time  [atom.s]

  // Carried through to every report, because for a triage answer the inputs that produced it
  // are part of the answer.
  std::string seedProvenance;

  int timeCount() const { return static_cast<int>(times.size()); }
  int nuclideCount() const { return static_cast<int>(nuclideKeys.size()); }

  std::span<const double> atomsAt(int timeIndex) const {
    const std::size_t n = nuclideKeys.size();
    return std::span<const double>(atoms.data() + static_cast<std::size_t>(timeIndex) * n, n);
  }

  std::span<const double> integratedAtomsAt(int timeIndex) const {
    const std::size_t n = nuclideKeys.size();
    return std::span<const double>(integratedAtoms.data() + static_cast<std::size_t>(timeIndex) * n,
                                   n);
  }
};

}  // namespace nusift
