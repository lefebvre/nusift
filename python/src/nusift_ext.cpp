// The nusift Python extension.
//
// The binding is thin on purpose. Everything here is a direct projection of the C++ API, with
// two adaptations that matter in Python and nowhere else:
//
//   * The big arrays -- inventories over time, response tables -- are handed out as ZERO-COPY
//     NumPy views over the C++ storage, with the owning Python object as the array's base. A
//     sixty-point run over a fission source is a 60 x 1000 matrix; copying it into a list of
//     lists would cost more than the solve did.
//
//   * Times and units are accepted as the same strings the CLI takes ("30d", "1h:100y:log:60",
//     "Sv/h"), parsed by the same functions. A notebook and a terminal should not disagree
//     about what "1.5y" means.
//
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/string_view.h>
#include <nanobind/stl/vector.h>

#include <optional>
#include <string>
#include <vector>

#include "nusift/core/error.hpp"
#include "nusift/core/nuclide_name.hpp"
#include "nusift/engine/decay_engine.hpp"
#include "nusift/engine/inventory.hpp"
#include "nusift/io/inventory_io.hpp"
#include "nusift/io/time_spec.hpp"
#include "nusift/nucdata/nuclear_data.hpp"
#include "nusift/nucdata/store_locator.hpp"
#include "nusift/seed/seed_fission.hpp"
#include "nusift/triage/forecast.hpp"
#include "nusift/triage/ranking.hpp"
#include "nusift/triage/response.hpp"
#include "nusift/version.hpp"

namespace nb = nanobind;
using namespace nb::literals;
using namespace nusift;

namespace {

// --- string forms, shared with the CLI ---------------------------------------

Metric metricFrom(const std::string& text) {
  if (text == "activity") {
    return Metric::Activity;
  }
  if (text == "exposure") {
    return Metric::Exposure;
  }
  throw InputError("metric: \"" + text + "\" is not a metric (activity or exposure)");
}

Aggregate aggregateFrom(const std::string& text) {
  if (text == "nuclide") {
    return Aggregate::Nuclide;
  }
  if (text == "mass-chain") {
    return Aggregate::MassChain;
  }
  if (text == "element") {
    return Aggregate::Element;
  }
  if (text == "line") {
    return Aggregate::GammaLine;
  }
  throw InputError("by: \"" + text + "\" is not an aggregate (nuclide, mass-chain, element, line)");
}

// A time given either as a number of seconds or as a string the CLI would accept.
double timeFrom(const nb::object& value) {
  if (nb::isinstance<nb::str>(value)) {
    return parseDuration(nb::cast<std::string>(value));
  }
  return nb::cast<double>(value);
}

// Zero-copy 2-D view over C++-owned storage, with `owner` keeping it alive. Const, because a
// consumer writing into a response table would corrupt the totals computed alongside it.
nb::ndarray<nb::numpy, const double, nb::ndim<2>> view2d(const double* data, std::size_t rows,
                                                         std::size_t cols, nb::handle owner) {
  return nb::ndarray<nb::numpy, const double, nb::ndim<2>>(data, {rows, cols}, nb::find(owner));
}

nb::ndarray<nb::numpy, const double, nb::ndim<1>> view1d(const double* data, std::size_t n,
                                                         nb::handle owner) {
  return nb::ndarray<nb::numpy, const double, nb::ndim<1>>(data, {n}, nb::find(owner));
}

// The store this wheel ships, if it ships one. Asked of the pure-Python locator rather than
// reimplemented here: the C++ search takes caller-supplied candidates precisely so a binding
// can contribute the one location only it knows about, without becoming a second search order
// that could disagree with the CLI's about which evaluation is in use.
std::vector<std::string> packagedStorePaths() {
  try {
    const nb::object located = nb::module_::import_("nusift._data").attr("default_store_path")();
    if (!located.is_none()) {
      return {nb::cast<std::string>(nb::str(located))};
    }
  } catch (const nb::python_error&) {
    // _core imported bare, without the package around it. That removes one candidate from the
    // search; it is not a reason to fail one that would have succeeded without it.
  }
  return {};
}

std::vector<std::string> nuclideNames(const std::vector<std::int64_t>& keys) {
  std::vector<std::string> names;
  names.reserve(keys.size());
  for (const std::int64_t key : keys) {
    names.push_back(formatNuclideName(Zai::fromKey(key)));
  }
  return names;
}

}  // namespace

NB_MODULE(_core, m) {
  m.doc() = "NuSIFT: which isotopes dominate activity or exposure, and when.";
  m.attr("__version__") = nusift::kVersion;

  // Every NuSIFT error becomes a Python exception, and InputError derives from NusiftError so
  // the C++ hierarchy survives the crossing. Without the explicit base they would be unrelated
  // Python types, and `except NusiftError` would silently miss every bad-input error -- which
  // is most of what a user actually hits.
  const nb::object baseError = nb::exception<NusiftError>(m, "NusiftError");
  nb::exception<InputError>(m, "InputError", baseError);

  // --- time helpers ----------------------------------------------------------
  m.def("parse_duration", &parseDuration, "text"_a,
        "Seconds from '30d', '1.5y', '90m', or a bare number of seconds.");
  m.def("parse_time_grid", &parseTimeGrid, "spec"_a,
        "Times from a grid spec such as '1h:100y:log:60'.");
  m.def(
      "logspace",
      [](const nb::object& start, const nb::object& stop, int count) {
        return logspace(timeFrom(start), timeFrom(stop), count);
      },
      "start"_a, "stop"_a, "count"_a, "Log-spaced times; endpoints exact.");
  m.def(
      "linspace",
      [](const nb::object& start, const nb::object& stop, int count) {
        return linspace(timeFrom(start), timeFrom(stop), count);
      },
      "start"_a, "stop"_a, "count"_a, "Linearly spaced times; endpoints exact.");
  m.def("format_duration", &formatDuration, "seconds"_a);

  // --- geometry --------------------------------------------------------------
  nb::class_<exposure::PointSourceGeometry>(m, "PointSource", "Unshielded point source in air.")
      .def(nb::init<>())
      .def(
          "__init__",
          [](exposure::PointSourceGeometry* self, double distance_m, double air_density,
             bool air_attenuation, double buildup) {
            new (self)
                exposure::PointSourceGeometry{distance_m, air_density, air_attenuation, buildup};
          },
          "distance_m"_a = 1.0, "air_density"_a = 1.205, "air_attenuation"_a = true,
          "buildup"_a = 1.0)
      .def_rw("distance_m", &exposure::PointSourceGeometry::distanceM)
      .def_rw("air_density", &exposure::PointSourceGeometry::airDensityKgM3)
      .def_rw("air_attenuation", &exposure::PointSourceGeometry::airAttenuation)
      .def_rw("buildup", &exposure::PointSourceGeometry::buildup);

  // --- nuclear data ----------------------------------------------------------
  m.def(
      "store_search_paths",
      [] {
        StoreSearch search;
        search.extraPaths = packagedStorePaths();
        return storeSearchPaths(search);
      },
      "Every path NuclearData.open() would try, in order. Diagnostic; touches no files.");

  nb::class_<NuclearData>(m, "NuclearData", "A staged nuclear-data store.")
      .def_static(
          "open",
          [](std::optional<std::string> path) {
            StoreSearch search;
            if (path) {
              search.explicitPath = *path;
            } else {
              // Only when no path was given: an explicit path short-circuits the search, so
              // resolving the packaged store would be work whose result is discarded.
              search.extraPaths = packagedStorePaths();
            }
            return NuclearData::open(locateStore(search));
          },
          "path"_a = nb::none(),
          "Open a store. With no path, searches $NUSIFT_DATA_STORE, the store packaged with "
          "this wheel, and the usual locations.")
      .def_prop_ro("size", &NuclearData::size)
      .def_prop_ro("staged_count", &NuclearData::stagedCount)
      .def_prop_ro("has_photon_lines", &NuclearData::hasPhotonLines)
      .def_prop_ro("has_atomic_weights", &NuclearData::hasAtomicWeights)
      .def_prop_ro("library", [](const NuclearData& d) { return d.provenance().library; })
      .def_prop_ro("staged_utc", [](const NuclearData& d) { return d.provenance().createdUtc; })
      .def_prop_ro(
          "fissionable",
          [](const NuclearData& d) {
            std::vector<std::string> names;
            for (const Zai& zai : d.fissionYields().parents()) {
              names.push_back(formatNuclideName(zai));
            }
            return names;
          },
          "Nuclides this store can seed a fission inventory from.")
      .def(
          "half_life",
          [](const NuclearData& d, const std::string& name) {
            const int i = d.indexOf(requireNuclideName(name));
            return i >= 0 ? d.halfLifeSeconds(i) : 0.0;
          },
          "nuclide"_a, "Half-life in seconds; 0 for stable or absent.")
      .def(
          "molar_mass",
          [](const NuclearData& d, const std::string& name) {
            const int i = d.indexOf(requireNuclideName(name));
            return i >= 0 ? d.molarMassGPerMol(i) : 0.0;
          },
          "nuclide"_a, "Molar mass in g/mol from the staged atomic weight; 0 if absent.")
      .def(
          "gamma_constant",
          [](const NuclearData& d, const std::string& name, double minEnergyEv) {
            const int i = d.indexOf(requireNuclideName(name));
            if (i < 0) {
              return 0.0;
            }
            if (minEnergyEv <= 0.0) {
              return exposure::gammaConstant(d.lines(i));
            }
            // Published tabulations usually state a low-energy cutoff -- 20 keV is the common
            // one -- because soft X-rays are absorbed by any real source encapsulation before
            // they reach air. Applying the same cutoff is what makes a comparison against such
            // a table a like-for-like one rather than a comparison of two conventions.
            std::vector<GammaLine> kept;
            for (const GammaLine& line : d.lines(i)) {
              if (line.energyEv >= minEnergyEv) {
                kept.push_back(line);
              }
            }
            return exposure::gammaConstant(LineSpectrum(kept.data(), kept.size()));
          },
          "nuclide"_a, "min_energy_ev"_a = 0.0,
          "Specific gamma-ray constant in R*m^2/(h*Bq), vacuum. With min_energy_ev, counts only "
          "photons at or above that energy, matching tabulations that state a cutoff.")
      .def("__repr__", [](const NuclearData& d) {
        return "<NuclearData " + d.provenance().library + ", " + std::to_string(d.stagedCount()) +
               " nuclides>";
      });

  // --- inventory -------------------------------------------------------------
  nb::class_<Inventory>(m, "Inventory", "An isotopic inventory, in atoms.")
      .def(nb::init<>())
      .def(
          "add",
          [](Inventory& inv, const std::string& name, double atoms) {
            inv.add(requireNuclideName(name), atoms);
          },
          "nuclide"_a, "atoms"_a)
      .def_prop_ro("total_atoms", &Inventory::totalAtoms)
      .def_prop_ro("provenance", &Inventory::provenance)
      .def_prop_ro("nuclides",
                   [](const Inventory& inv) {
                     std::vector<std::string> names;
                     for (const InventoryEntry& e : inv.entries()) {
                       names.push_back(formatNuclideName(Zai::fromKey(e.zaiKey)));
                     }
                     return names;
                   })
      .def_prop_ro("atoms",
                   [](const Inventory& inv) {
                     std::vector<double> values;
                     for (const InventoryEntry& e : inv.entries()) {
                       values.push_back(e.atoms);
                     }
                     return values;
                   })
      .def("__len__", &Inventory::size)
      .def("__repr__", [](const Inventory& inv) {
        return "<Inventory " + std::to_string(inv.size()) + " nuclides>";
      });

  m.def(
      "read_inventory",
      [](const std::string& path, const NuclearData& data, bool ignore_unknown) {
        InventoryReadOptions options;
        options.ignoreUnknown = ignore_unknown;
        return readInventory(path, data, options);
      },
      "path"_a, "data"_a, "ignore_unknown"_a = false, "Read an inventory CSV or JSON.");

  m.def(
      "seed_fission",
      [](const NuclearData& data, const std::string& fissile, const std::string& energy,
         std::optional<double> fissions, std::optional<double> yield_kt,
         std::optional<double> energy_j, double mev_per_fission) {
        seed::FissionSeed fissionSeed;
        fissionSeed.fissile = requireNuclideName(fissile);
        if (!parseIncidentEnergy(energy, fissionSeed.incidentEnergyEv)) {
          throw InputError("energy: \"" + energy + "\" is not an incident energy");
        }
        fissionSeed.meVPerFission = mev_per_fission;

        const int given = (fissions ? 1 : 0) + (yield_kt ? 1 : 0) + (energy_j ? 1 : 0);
        if (given != 1) {
          throw InputError(
              "give exactly one of fissions=, yield_kt=, or energy_j= to size the source");
        }
        if (fissions) {
          fissionSeed.fissions = *fissions;
        } else if (yield_kt) {
          fissionSeed.fissions = seed::fissionsFromKt(*yield_kt, mev_per_fission);
        } else {
          fissionSeed.fissions = seed::fissionsFromEnergyJ(*energy_j, mev_per_fission);
        }
        return seed::seedFromFission(data, fissionSeed);
      },
      "data"_a, "fissile"_a, "energy"_a = "thermal", "fissions"_a = nb::none(),
      "yield_kt"_a = nb::none(), "energy_j"_a = nb::none(),
      "mev_per_fission"_a = seed::kMeVPerFissionExplosiveYield,
      "Build an inventory from fission. 180 MeV per fission is the explosive-yield "
      "convention; pass 200 for total recoverable energy.");

  m.def("fissions_from_kt", &seed::fissionsFromKt, "kt"_a,
        "mev_per_fission"_a = seed::kMeVPerFissionExplosiveYield);

  // --- decay -----------------------------------------------------------------
  nb::class_<DecayResult>(m, "DecayResult", "Atoms and their exact time integrals.")
      .def_prop_ro("times",
                   [](nb::handle self) {
                     const DecayResult& r = nb::cast<const DecayResult&>(self);
                     return view1d(r.times.data(), r.times.size(), self);
                   })
      .def_prop_ro("nuclides", [](const DecayResult& r) { return nuclideNames(r.nuclideKeys); })
      .def_prop_ro(
          "atoms",
          [](nb::handle self) {
            const DecayResult& r = nb::cast<const DecayResult&>(self);
            return view2d(r.atoms.data(), static_cast<std::size_t>(r.timeCount()),
                          static_cast<std::size_t>(r.nuclideCount()), self);
          },
          "(times, nuclides) atom counts. A zero-copy view, not a copy.")
      .def_prop_ro(
          "integrated_atoms",
          [](nb::handle self) {
            const DecayResult& r = nb::cast<const DecayResult&>(self);
            return view2d(r.integratedAtoms.data(), static_cast<std::size_t>(r.timeCount()),
                          static_cast<std::size_t>(r.nuclideCount()), self);
          },
          "(times, nuclides) atom-seconds. A zero-copy view, not a copy.")
      .def("__repr__", [](const DecayResult& r) {
        return "<DecayResult " + std::to_string(r.timeCount()) + " times x " +
               std::to_string(r.nuclideCount()) + " nuclides>";
      });

  m.def(
      "decay",
      [](const NuclearData& data, const Inventory& inventory, const std::vector<double>& times,
         int threads, bool prune, int cram_order) {
        DecayOptions options;
        options.threads = threads;
        options.prune = prune;
        options.order = cram_order == 16 ? CramOrder::Order16 : CramOrder::Order48;
        return decay(data, inventory, times, options);
      },
      "data"_a, "inventory"_a, "times"_a, "threads"_a = 0, "prune"_a = true, "cram_order"_a = 48,
      "Decay an inventory to a set of times, in seconds.");

  // --- response and ranking --------------------------------------------------
  nb::class_<Contributor>(m, "Contributor")
      .def_ro("label", &Contributor::label)
      .def_ro("value", &Contributor::value)
      .def_ro("fraction", &Contributor::fraction)
      .def_ro("cumulative_fraction", &Contributor::cumulativeFraction)
      .def_ro("rank", &Contributor::rank)
      .def_ro("pinned", &Contributor::pinned)
      .def_prop_ro("key", [](const Contributor& c) { return c.id.key; })
      .def_prop_ro("line_energy_ev", [](const Contributor& c) { return c.id.lineEnergyEv; })
      .def("__repr__", [](const Contributor& c) {
        char buffer[128];
        std::snprintf(buffer, sizeof(buffer), "<Contributor %s %.4g (%.1f%%)>", c.label.c_str(),
                      c.value, c.fraction * 100.0);
        return std::string(buffer);
      });

  nb::class_<Ranking>(m, "Ranking")
      .def_ro("contributors", &Ranking::contributors)
      .def_ro("total", &Ranking::total)
      .def_ro("covered_fraction", &Ranking::coveredFraction)
      .def_ro("omitted_count", &Ranking::omittedCount)
      .def_ro("time", &Ranking::time)
      .def_ro("unmodeled_energy_fraction", &Ranking::unmodeledEnergyFraction)
      .def_prop_ro("labels",
                   [](const Ranking& r) {
                     std::vector<std::string> names;
                     for (const Contributor& c : r.contributors) {
                       names.push_back(c.label);
                     }
                     return names;
                   })
      .def("__len__", [](const Ranking& r) { return r.contributors.size(); })
      .def("__repr__", [](const Ranking& r) {
        return "<Ranking " + std::to_string(r.contributors.size()) + " of " +
               std::to_string(r.contributors.size() + r.omittedCount) + " contributors>";
      });

  nb::class_<DominanceWindow>(m, "DominanceWindow")
      .def_ro("label", &DominanceWindow::label)
      .def_ro("start_s", &DominanceWindow::startSeconds)
      .def_ro("end_s", &DominanceWindow::endSeconds)
      .def_ro("peak_fraction", &DominanceWindow::peakFraction)
      .def("__repr__", [](const DominanceWindow& w) {
        return "<DominanceWindow " + w.label + " " + formatDuration(w.startSeconds) + " to " +
               formatDuration(w.endSeconds) + ">";
      });

  nb::class_<ResponseTable>(m, "ResponseTable", "Per-contributor values over time.")
      .def_prop_ro("labels", [](const ResponseTable& t) { return t.labels; })
      .def_prop_ro("times",
                   [](nb::handle self) {
                     const ResponseTable& t = nb::cast<const ResponseTable&>(self);
                     return view1d(t.times.data(), t.times.size(), self);
                   })
      .def_prop_ro("totals",
                   [](nb::handle self) {
                     const ResponseTable& t = nb::cast<const ResponseTable&>(self);
                     return view1d(t.totals.data(), t.totals.size(), self);
                   })
      .def_prop_ro(
          "values",
          [](nb::handle self) {
            const ResponseTable& t = nb::cast<const ResponseTable&>(self);
            return view2d(t.values.data(), static_cast<std::size_t>(t.timeCount()),
                          static_cast<std::size_t>(t.contributorCount()), self);
          },
          "(times, contributors). A zero-copy view, not a copy.")
      .def_prop_ro("unit", [](const ResponseTable& t) { return std::string(unitName(t.unit)); })
      .def(
          "rank",
          [](const ResponseTable& table, const nb::object& at, int top, double coverage,
             double min_fraction, const nb::object& pin) {
            RankRequest request;
            request.topN = top;
            request.coverage = coverage;
            request.minFraction = min_fraction;

            // A bare string is one pin, not an iterable of one-character ones. Python makes
            // that mistake easy to write and impossible to notice, since "Cs-137" is a perfectly
            // good sequence -- of six spellings that name nothing.
            if (!pin.is_none()) {
              if (nb::isinstance<nb::str>(pin)) {
                request.pinned.push_back(requirePin(table, nb::cast<std::string>(pin)));
              } else {
                for (const nb::handle item : pin) {
                  request.pinned.push_back(requirePin(table, nb::cast<std::string>(item)));
                }
              }
            }

            int index = 0;
            if (!at.is_none()) {
              // Nearest grid point to the requested time, so `at="30d"` works on a log grid
              // that has no sample exactly there.
              const double wanted = timeFrom(at);
              double best = std::abs(table.times[0] - wanted);
              for (int k = 1; k < table.timeCount(); ++k) {
                const double distance = std::abs(table.times[static_cast<std::size_t>(k)] - wanted);
                if (distance < best) {
                  best = distance;
                  index = k;
                }
              }
            }
            return rank(table, index, request);
          },
          "at"_a = nb::none(), "top"_a = 10, "coverage"_a = 0.0, "min_fraction"_a = 0.0,
          "pin"_a = nb::none(),
          "Rank at the grid time nearest `at`. `pin` is a contributor name, or several, that "
          "appear whatever they rank -- appended below the ranking carrying the place they "
          "actually hold.")
      .def(
          "dominance_windows",
          [](const ResponseTable& table, int min_samples) {
            return dominanceWindows(table, min_samples);
          },
          "min_samples"_a = 2, "Who leads, and over which windows.")
      .def("__repr__", [](const ResponseTable& t) {
        return "<ResponseTable " + std::to_string(t.timeCount()) + " times x " +
               std::to_string(t.contributorCount()) + " contributors, " + unitName(t.unit) + ">";
      });

  m.def(
      "response",
      [](const NuclearData& data, const DecayResult& result, const std::string& metric,
         const std::string& by, const std::string& units,
         const exposure::PointSourceGeometry& geometry) {
        ResponseSpec spec;
        spec.metric = metricFrom(metric);
        spec.aggregate = aggregateFrom(by);
        spec.unit = requireUnit(units, spec.metric, Domain::Instant);
        spec.geometry = geometry;
        return buildResponse(data, result, spec);
      },
      "data"_a, "result"_a, "metric"_a = "activity", "by"_a = "nuclide", "units"_a = "",
      "geometry"_a = exposure::PointSourceGeometry{},
      "Turn a decay result into a table of per-contributor values.");
}
