// Writes a small nuclear-data store from hand-entered evaluated data.
//
//   make_test_store <out.h5>
//
// Six decaying nuclides and their stable terminators, chosen to exercise the cases that
// matter rather than to be representative of an evaluation:
//
//   Cs-137 -> Ba-137m -> Ba-137    a branching decay feeding a short-lived isomer that
//                                  carries essentially all of the chain's photon output
//   Sr-90  -> Y-90    -> Zr-90     a pure beta chain with no photons at all
//   Co-60  -> Ni-60                two strong lines, the calibration case for gamma dose
//
// That combination covers a nuclide whose activity is invisible to a photon metric (Sr-90),
// one whose photons come from a daughter rather than itself (Cs-137), and one with a
// continuum fraction large enough to trip the unmodeled-photon flag (Y-90 bremsstrahlung).
//
// Half-lives and line intensities are from ENDF/B-VIII.0 decay data. They are accurate
// enough to be worth checking answers against by hand, which is the point -- but this is a
// hand-entered fixture, not an evaluation, and nothing outside the test suite should read it.
#include <algorithm>
#include <cstdio>
#include <string>
#include <vector>

#include "nusift/core/nuclide.hpp"
#include "nusift/nucdata/data_store.hpp"
#include "nusift/units.hpp"
#include "nusift/version.hpp"

namespace {

using namespace nusift;

constexpr double kYear = units::kSecondsPerYear;

struct Line {
  double energyEv;
  double intensity;
  int styp;
};

struct Mode {
  double rtyp;
  double branching;
  int finalState;
};

struct Nuclide {
  Zai zai;
  double halfLifeSeconds;  // <= 0 for stable
  double awr;
  double emEnergyEv;
  double continuumPhotonEv;
  std::vector<Mode> modes;
  std::vector<Line> lines;
};

// Atomic weight ratios are mass / neutron mass, so molar mass is AWR * 1.00866491595.
// Entered from the nuclide's atomic mass rather than from A, which is the whole reason the
// store carries them: using A instead is wrong by ~0.1 g/mol here.
std::vector<Nuclide> evaluatedData() {
  return {
      // Cs-137: 94.7% to the Ba-137m isomer, 5.3% straight to the ground state. Its own
      // photon output is negligible -- the 662 keV line everyone attributes to Cs-137 is
      // actually emitted by Ba-137m.
      {Zai{55, 137, 0},
       30.08 * kYear,
       135.7305,
       1.0e3,
       0.0,
       {{1.0, 0.9470, 1}, {1.0, 0.0530, 0}},
       {}},
      // Ba-137m: 153 s isomeric transition carrying the 661.657 keV line.
      {Zai{56, 137, 1},
       153.0,
       135.7192,
       6.0e5,
       0.0,
       {{3.0, 1.0, 0}},
       {{661657.0, 0.8994, 0}, {31817.0, 0.0199, 9}, {32194.0, 0.0366, 9}}},
      {Zai{56, 137, 0}, 0.0, 135.7192, 0.0, 0.0, {}, {}},

      // Sr-90 -> Y-90: a pure beta chain. No photon lines at all, so it dominates activity
      // while contributing nothing to a photon-based metric -- exactly the divergence
      // between metrics that makes ranking by the right one matter.
      {Zai{38, 90, 0}, 28.79 * kYear, 89.0777, 0.0, 0.0, {{1.0, 1.0, 0}}, {}},
      // Y-90's photon output is almost entirely bremsstrahlung continuum, which NuSIFT does
      // not model; recording it is what makes the understatement visible.
      {Zai{39, 90, 0}, 64.05 * 3600.0, 89.0714, 1.7e3, 1.5e4, {{1.0, 1.0, 0}}, {}},
      {Zai{40, 90, 0}, 0.0, 89.0655, 0.0, 0.0, {}, {}},

      // Co-60: two strong lines at essentially unit intensity. The standard case for a
      // gamma-dose constant, and the easiest result to check against a published value.
      {Zai{27, 60, 0},
       5.2711 * kYear,
       59.3357,
       2.5e6,
       0.0,
       {{1.0, 1.0, 0}},
       {{1173228.0, 0.9985, 0}, {1332492.0, 0.9998, 0}}},
      {Zai{28, 60, 0}, 0.0, 59.3216, 0.0, 0.0, {}, {}},
  };
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 2) {
    std::fprintf(stderr, "usage: make_test_store <out.h5>\n");
    return 2;
  }

  std::vector<Nuclide> nuclides = evaluatedData();
  // The store's nuclide axis is sorted ascending by key, so sort here rather than relying on
  // the order they happen to be written in above.
  std::sort(nuclides.begin(), nuclides.end(),
            [](const Nuclide& a, const Nuclide& b) { return a.zai.key() < b.zai.key(); });

  StoreArrays arrays;
  arrays.provenance.library = "hand-entered test fixture (ENDF/B-VIII.0 values)";
  arrays.provenance.createdUtc = "2026-08-13T00:00:00Z";
  arrays.provenance.nusiftVersion = kVersion;
  arrays.provenance.stagedTapeCount = 0;
  arrays.provenance.decaySource = DataSource::Endf;
  arrays.provenance.linesSource = DataSource::Endf;
  arrays.provenance.yieldsSource = DataSource::None;

  arrays.modeOffset.push_back(0);
  arrays.lineOffset.push_back(0);
  for (const Nuclide& nuclide : nuclides) {
    arrays.nuclideKey.push_back(nuclide.zai.key());
    arrays.halfLife.push_back(nuclide.halfLifeSeconds);
    arrays.awr.push_back(nuclide.awr);
    arrays.emEnergyEv.push_back(nuclide.emEnergyEv);
    arrays.continuumPhotonEv.push_back(nuclide.continuumPhotonEv);

    for (const Mode& mode : nuclide.modes) {
      arrays.modeRtyp.push_back(mode.rtyp);
      arrays.modeBranching.push_back(mode.branching);
      arrays.modeFinalState.push_back(mode.finalState);
      arrays.modeIsFission.push_back(0);
    }
    arrays.modeOffset.push_back(static_cast<int>(arrays.modeRtyp.size()));

    for (const Line& line : nuclide.lines) {
      arrays.lineEnergyEv.push_back(line.energyEv);
      arrays.lineIntensity.push_back(line.intensity);
      arrays.lineStyp.push_back(line.styp);
    }
    arrays.lineOffset.push_back(static_cast<int>(arrays.lineEnergyEv.size()));
  }

  try {
    writeStore(argv[1], arrays);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "make_test_store: %s\n", e.what());
    return 1;
  }
  std::printf("wrote %s (%d nuclides, %zu photon lines)\n", argv[1], arrays.nuclideCount(),
              arrays.lineEnergyEv.size());
  return 0;
}
