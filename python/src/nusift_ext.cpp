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
#include <sstream>
#include <string>
#include <vector>

#include "nusift/core/error.hpp"
#include "nusift/core/nuclide_name.hpp"
#include "nusift/engine/decay_engine.hpp"
#include "nusift/engine/inventory.hpp"
#include "nusift/engine/reconcile.hpp"
#include "nusift/engine/sensitivity.hpp"
#include "nusift/exposure/dose_coefficients.hpp"
#include "nusift/io/inventory_io.hpp"
#include "nusift/io/number_format.hpp"
#include "nusift/io/source_report.hpp"
#include "nusift/io/time_spec.hpp"
#include "nusift/nucdata/coefficient_pack.hpp"
#include "nusift/nucdata/nuclear_data.hpp"
#include "nusift/nucdata/store_locator.hpp"
#include "nusift/seed/seed_fission.hpp"
#include "nusift/triage/allowable.hpp"
#include "nusift/triage/attribution.hpp"
#include "nusift/triage/events.hpp"
#include "nusift/triage/forecast.hpp"
#include "nusift/triage/intervention.hpp"
#include "nusift/triage/ranking.hpp"
#include "nusift/triage/response.hpp"
#include "nusift/triage/spectrum.hpp"
#include "nusift/triage/task_plan.hpp"
#include "nusift/triage/triage_set.hpp"
#include "nusift/triage/uncertainty.hpp"
#include "nusift/triage/yield_uncertainty.hpp"
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
  if (text == "photon") {
    return Metric::Photon;
  }
  throw InputError("metric: \"" + text + "\" is not a metric (activity, exposure, or photon)");
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
// Resolve a contributor named the way a user writes one into its column. requirePin does the
// naming half -- "Cs-137", "A=140", "Cs" -- and returns a key; a series needs the column.
int columnFor(const ResponseTable& table, const std::string& name) {
  const std::int64_t key = requirePin(table, name);
  for (int c = 0; c < table.contributorCount(); ++c) {
    if (table.contributors[static_cast<std::size_t>(c)].key == key) {
      return c;
    }
  }
  throw NusiftError("trajectory: \"" + name + "\" resolved to a key the table does not hold");
}

// Which curve an event search runs over: the total, one contributor, or a ratio of two. Shared
// by every search method so `of=` and `ratio=` cannot come to mean different things on
// different calls.
EventSeries seriesFor(const ResponseTable& table, const nb::object& of, const nb::object& ratio,
                      const ResponseEvaluator* evaluator = nullptr) {
  if (!ratio.is_none()) {
    std::vector<std::string> names;
    for (const nb::handle item : ratio) {
      names.push_back(nb::cast<std::string>(item));
    }
    if (names.size() != 2) {
      throw InputError(
          "trajectory: ratio takes exactly two contributors, e.g. (\"Zr-95\", \"Nb-95\")");
    }
    const int numerator = columnFor(table, names[0]);
    const int denominator = columnFor(table, names[1]);
    return evaluator != nullptr ? ratioSeries(table, numerator, denominator, *evaluator)
                                : ratioSeries(table, numerator, denominator);
  }
  if (!of.is_none()) {
    const int column = columnFor(table, nb::cast<std::string>(of));
    return evaluator != nullptr ? contributorSeries(table, column, *evaluator)
                                : contributorSeries(table, column);
  }
  return evaluator != nullptr ? totalSeries(table, *evaluator) : totalSeries(table);
}

// The evaluator a search was handed, or nothing. Written once because every search takes the
// same pair of optional arguments, and a `refine=` that meant something different on one of
// them would be worse than not offering it there.
const ResponseEvaluator* evaluatorFrom(const nb::object& refine) {
  if (refine.is_none()) {
    return nullptr;
  }
  // A reference cast rather than a pointer one: it fails loudly when `refine` is not an
  // evaluator, instead of handing a search something that is not a curve.
  const ResponseEvaluator& evaluator = nb::cast<const ResponseEvaluator&>(refine);
  return &evaluator;
}

EventTolerance toleranceOf(double relative) {
  EventTolerance tolerance;
  tolerance.relative = relative;
  return tolerance;
}

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

// Open the store the search resolves to, and hand anything the search had to say about the
// choice to Python's warnings module -- a notebook does not see the process's stderr, and a
// store silently picked out of several is exactly the thing a user of one needs to know.
NuclearData openLocated(StoreSearch search) {
  std::ostringstream passedOver;
  search.warnings = &passedOver;
  const std::string path = locateStore(search);
  NuclearData data = NuclearData::open(path);
  if (!passedOver.str().empty()) {
    nb::module_::import_("warnings").attr("warn")(passedOver.str());
  }
  return data;
}

// The interval counterpart of DecayResult: one row of atom-seconds over [t1, t2], in the
// pruned index space intervalIntegral() solved over. Defined here rather than in the library
// because the CLI consumes the integral and its keys directly; Python needs one object to
// hand back, so that the array can be a view with an owner to keep it alive.
struct IntervalResult {
  double t1 = 0.0;
  double t2 = 0.0;
  std::vector<std::int64_t> nuclideKeys;
  std::vector<double> integratedAtoms;  // [nNuc], atom-seconds over the window
};

const char* domainName(Domain domain) {
  return domain == Domain::Interval ? "interval" : "instant";
}

// A binning described the way a keyword call describes one. Energies are eV here rather than
// the CLI's keV: the library speaks eV, and a binding whose numbers needed a different scale
// from the API it wraps would be a conversion waiting to be forgotten.
BinningSpec binningFrom(int bins, const std::string& scale, double minEv, double maxEv,
                        const nb::object& edges) {
  BinningSpec binning;
  binning.count = bins;
  binning.minEv = minEv;
  binning.maxEv = maxEv;
  if (!edges.is_none()) {
    binning.edgesEv = nb::cast<std::vector<double>>(edges);
  }
  if (!parseBinScale(scale, binning.scale)) {
    throw InputError("spectrum: \"" + scale + "\" is not a bin scale (linear or log)");
  }
  return binning;
}

}  // namespace

NB_MODULE(_core, m) {
  m.doc() = "NuSIFT: which isotopes dominate activity, exposure, or photon output, and when.";
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
             bool air_attenuation, double buildup, const std::string& irradiation) {
            exposure::Irradiation orientation = exposure::Irradiation::AP;
            if (!exposure::parseIrradiation(irradiation, orientation)) {
              throw InputError("exposure: \"" + irradiation +
                               "\" is not an irradiation geometry (ap, pa, llat, rlat, rot, iso)");
            }
            new (self) exposure::PointSourceGeometry{distance_m, air_density, air_attenuation,
                                                     buildup, orientation};
          },
          "distance_m"_a = 1.0, "air_density"_a = 1.205, "air_attenuation"_a = true,
          "buildup"_a = 1.0, "irradiation"_a = "ap")
      .def_rw("distance_m", &exposure::PointSourceGeometry::distanceM)
      .def_rw("air_density", &exposure::PointSourceGeometry::airDensityKgM3)
      .def_rw("air_attenuation", &exposure::PointSourceGeometry::airAttenuation)
      .def_rw("buildup", &exposure::PointSourceGeometry::buildup)
      // Read and written as the spelling ICRP prints, not as an opaque enum: the string is what
      // a report shows and what the CLI takes, and a binding that spoke a different language
      // for the same fact would be a third spelling to keep in step.
      .def_prop_rw(
          "irradiation",
          [](const exposure::PointSourceGeometry& g) {
            return std::string(exposure::irradiationName(g.irradiation));
          },
          [](exposure::PointSourceGeometry& g, const std::string& text) {
            if (!exposure::parseIrradiation(text, g.irradiation)) {
              throw InputError("exposure: \"" + text +
                               "\" is not an irradiation geometry (ap, pa, llat, rlat, rot, iso)");
            }
          },
          "How the body is oriented in the field, for sievert units: AP, PA, LLAT, RLAT, ROT, "
          "ISO. Half of what an effective dose means, so it is reported beside every one.");

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
            return openLocated(std::move(search));
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
      .def(
          "effective_dose_constant",
          [](const NuclearData& d, const std::string& name, const std::string& irradiation,
             double minEnergyEv) {
            const int i = d.indexOf(requireNuclideName(name));
            if (i < 0) {
              return 0.0;
            }
            // Vacuum at one metre, so the 1/(4 pi d^2) is the whole geometry and the result is
            // distance-independent -- the same construction gammaConstant() uses, and what
            // makes this comparable with a published per-activity coefficient.
            exposure::PointSourceGeometry vacuum;
            vacuum.distanceM = 1.0;
            vacuum.airAttenuation = false;
            if (!exposure::parseIrradiation(irradiation, vacuum.irradiation)) {
              throw InputError("exposure: \"" + irradiation +
                               "\" is not an irradiation geometry (ap, pa, llat, rlat, rot, iso)");
            }
            if (minEnergyEv <= 0.0) {
              return exposure::effectiveDoseRatePerBecquerel(d.lines(i), vacuum);
            }
            std::vector<GammaLine> kept;
            for (const GammaLine& line : d.lines(i)) {
              if (line.energyEv >= minEnergyEv) {
                kept.push_back(line);
              }
            }
            return exposure::effectiveDoseRatePerBecquerel(LineSpectrum(kept.data(), kept.size()),
                                                           vacuum);
          },
          "nuclide"_a, "irradiation"_a = "ap", "min_energy_ev"_a = 0.0,
          "ICRP 116 effective dose per unit activity at 1 m in vacuum, in Sv*m^2/(h*Bq). The "
          "gamma constant's counterpart for the quantity a sievert actually names, and what a "
          "published effective-dose coefficient is comparable with.")
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
      "path"_a, "data"_a, "ignore_unknown"_a = false,
      "Read an inventory CSV or JSON. A file whose rows carry assay dates is reconciled to its "
      "latest one on the way back, so what returns always describes the material at a single "
      "instant -- see read_assays() to see the reconciliation instead of only its result.");

  // --- assays taken on different dates ---------------------------------------
  //
  // An Inventory is atoms at ONE instant, which is what makes every answer built from it well
  // posed. Sheets from different dates are not addable until they are brought to a common one,
  // and bringing them there is a decay solve rather than a bookkeeping step.
  m.def("parse_date", &parseCalendarDate, "text"_a,
        "Seconds since 1970-01-01T00:00:00Z from an ISO date: '2024-03-15', or "
        "'2024-03-15T09:30:00Z'. UTC only -- every other spelling is ambiguous somewhere.");
  m.def("format_date", &formatCalendarDate, "seconds"_a, "The inverse of parse_date.");

  nb::class_<AssayGroup>(m, "AssayGroup", "One assay: what was measured, and when.")
      .def(nb::init<>())
      .def(
          "__init__",
          [](AssayGroup* self, const nb::object& date, const Inventory& inventory,
             const std::string& label) {
            // A date given as a string is parsed the way the file reader parses one, so a
            // notebook and a spreadsheet cannot disagree about which day "2024-03-15" is.
            const double seconds = nb::isinstance<nb::str>(date)
                                       ? parseCalendarDate(nb::cast<std::string>(date))
                                       : nb::cast<double>(date);
            // A caller who passed a date has stated one, so this is always dated.
            new (self) AssayGroup{seconds, inventory, label, /*dated=*/true};
          },
          "date"_a, "inventory"_a, "label"_a = "")
      .def_prop_ro("date", [](const AssayGroup& g) { return formatCalendarDate(g.dateSeconds); })
      .def_ro("date_s", &AssayGroup::dateSeconds)
      .def_ro("inventory", &AssayGroup::inventory)
      .def_ro("label", &AssayGroup::label)
      .def_ro("dated", &AssayGroup::dated,
              "Whether the file stated a date. An undated sheet's date is a placeholder zero "
              "and cannot be carried to any epoch but zero.")
      .def("__repr__", [](const AssayGroup& g) {
        return "<AssayGroup " + formatCalendarDate(g.dateSeconds) + ", " +
               std::to_string(g.inventory.size()) + " nuclides>";
      });

  nb::class_<AssayContribution>(m, "AssayContribution", "What one assay contributed.")
      .def_prop_ro("date",
                   [](const AssayContribution& c) { return formatCalendarDate(c.dateSeconds); })
      .def_ro("date_s", &AssayContribution::dateSeconds)
      .def_ro("label", &AssayContribution::label)
      .def_ro("carried_s", &AssayContribution::carriedSeconds,
              "How far this assay was carried to reach the epoch. Zero for the one defining it.")
      .def_ro("nuclides", &AssayContribution::nuclides)
      .def_ro("atoms_at_assay", &AssayContribution::atomsAtAssay)
      .def_ro("atoms_at_epoch", &AssayContribution::atomsAtEpoch);

  nb::class_<Reconciliation>(m, "Reconciliation", "Assays carried to one epoch and merged.")
      .def_ro("inventory", &Reconciliation::inventory,
              "The merged result: an ordinary Inventory, because at the epoch there is only "
              "one date left.")
      .def_prop_ro("epoch",
                   [](const Reconciliation& r) { return formatCalendarDate(r.epochSeconds); })
      .def_ro("epoch_s", &Reconciliation::epochSeconds)
      .def_ro("contributions", &Reconciliation::contributions)
      .def_ro("span_s", &Reconciliation::spanSeconds,
              "Between the earliest and latest assay. A wide span means an old measurement was "
              "carried a long way on nothing but the decay model.")
      .def("__repr__", [](const Reconciliation& r) {
        return "<Reconciliation " + std::to_string(r.contributions.size()) + " assays to " +
               formatCalendarDate(r.epochSeconds) + ">";
      });

  m.def(
      "read_assays",
      [](const std::string& path, const NuclearData& data, bool ignore_unknown) {
        InventoryReadOptions options;
        options.ignoreUnknown = ignore_unknown;
        const DatedInventory dated = readInventoryDated(path, data, options);
        // The `dated` flag is not returned separately: a file with no dates has exactly one
        // group, and asking whether a list has more than one entry is the same question without
        // a second thing to keep in step.
        return dated.groups;
      },
      "path"_a, "data"_a, "ignore_unknown"_a = false,
      "The assays a file describes, before anything is carried anywhere. One group for a file "
      "that names no dates.");

  m.def(
      "reconcile",
      [](const NuclearData& data, const std::vector<AssayGroup>& groups, const nb::object& epoch,
         int threads, bool prune, int cram_order) {
        DecayOptions options;
        options.threads = threads;
        options.prune = prune;
        options.order = cram_order == 16 ? CramOrder::Order16 : CramOrder::Order48;
        const double at = epoch.is_none() ? latestAssayDate(groups)
                          : nb::isinstance<nb::str>(epoch)
                              ? parseCalendarDate(nb::cast<std::string>(epoch))
                              : nb::cast<double>(epoch);
        const nb::gil_scoped_release release;
        return reconcile(data, groups, at, options);
      },
      "data"_a, "groups"_a, "epoch"_a = nb::none(), "threads"_a = 0, "prune"_a = true,
      "cram_order"_a = 48,
      "Carry every assay to `epoch` and merge. Defaults to the LATEST assay date, the only "
      "choice that carries every assay forward and none backward. An earlier epoch is refused: "
      "un-growing a daughter is an inverse problem with no unique answer, not a decay solve.");

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

  // The GIL is released for the solve. Nothing below touches a Python object once the
  // arguments have been converted, and an eighty-point forecast is a second of factorization
  // that other Python threads have no reason to wait behind.
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
      nb::call_guard<nb::gil_scoped_release>(),
      "Decay an inventory to a set of times, in seconds.");

  // --- intervals -------------------------------------------------------------
  nb::class_<IntervalResult>(m, "IntervalResult",
                             "Atom-seconds over one window, per nuclide, solved in closed form.")
      .def_ro("t1", &IntervalResult::t1, "Start of the window, in seconds.")
      .def_ro("t2", &IntervalResult::t2, "End of the window, in seconds.")
      .def_prop_ro("nuclides", [](const IntervalResult& r) { return nuclideNames(r.nuclideKeys); })
      .def_prop_ro(
          "integrated_atoms",
          [](nb::handle self) {
            const IntervalResult& r = nb::cast<const IntervalResult&>(self);
            return view1d(r.integratedAtoms.data(), r.integratedAtoms.size(), self);
          },
          "Per-nuclide atom-seconds over [t1, t2]. A zero-copy view, not a copy.")
      .def("__repr__", [](const IntervalResult& r) {
        return "<IntervalResult " + formatDuration(r.t1) + " to " + formatDuration(r.t2) + ", " +
               std::to_string(r.integratedAtoms.size()) + " nuclides>";
      });

  // The guarded path. Differencing two rows of DecayResult.integrated_atoms is correct in
  // exact arithmetic and loses its digits to cancellation for a narrow window late in the
  // decay, which is exactly the window a user asking for one is looking at carefully. This is
  // the route that re-solves the window directly when that would happen; see
  // docs/interval-integration.md.
  m.def(
      "integrate",
      [](const NuclearData& data, const Inventory& inventory, const nb::object& t1,
         const nb::object& t2, int threads, bool prune, int cram_order) {
        DecayOptions options;
        options.threads = threads;
        options.prune = prune;
        options.order = cram_order == 16 ? CramOrder::Order16 : CramOrder::Order48;
        IntervalResult result;
        result.t1 = timeFrom(t1);
        result.t2 = timeFrom(t2);
        // Released only once the Python arguments have been read; the solve itself needs
        // nothing from the interpreter.
        const nb::gil_scoped_release release;
        result.integratedAtoms =
            intervalIntegral(data, inventory, result.t1, result.t2, &result.nuclideKeys, options);
        return result;
      },
      "data"_a, "inventory"_a, "t1"_a, "t2"_a, "threads"_a = 0, "prune"_a = true,
      "cram_order"_a = 48,
      "The exact per-nuclide time integral over [t1, t2], in atom-seconds. Times are seconds "
      "or the strings the CLI takes ('1h', '30d').");

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
      .def_ro("mean_optical_depth", &Ranking::meanOpticalDepth,
              "Exposure and photon fluence only: the air path in mean free paths at the "
              "energies carrying this answer. Past about 0.5 with buildup 1.0, scattered "
              "photons are a large omission.")
      .def_ro("buildup", &Ranking::buildup)
      .def_ro("pack_coverage", &Ranking::packCoverage,
              "Pack metrics only: the share of the inventory the pack carries a coefficient "
              "for. A total taken without reading this cannot be told from one computed over "
              "the whole inventory.")
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

  nb::class_<SeedShare>(m, "SeedShare")
      .def_ro("label", &SeedShare::label)
      .def_ro("key", &SeedShare::key)
      .def_ro("seed_atoms", &SeedShare::seedAtoms)
      .def_ro("importance", &SeedShare::importance,
              "dR/dn0: what one more atom of this nuclide IN THE SEED would be worth, in the "
              "ranking's unit. A potent seed and a large one are different things, and only "
              "this column separates them.")
      .def_ro("value", &SeedShare::value)
      .def_ro("fraction", &SeedShare::fraction)
      .def_ro("cumulative_fraction", &SeedShare::cumulativeFraction)
      .def_ro("rank", &SeedShare::rank)
      .def_ro("pinned", &SeedShare::pinned)
      .def("__repr__", [](const SeedShare& s) {
        char buffer[128];
        std::snprintf(buffer, sizeof(buffer), "<SeedShare %s %.4g (%.1f%%)>", s.label.c_str(),
                      s.value, s.fraction * 100.0);
        return std::string(buffer);
      });

  nb::class_<SeedAttribution>(m, "SeedAttribution",
                              "A response decomposed over the nuclides that were seeded.")
      .def_ro("shares", &SeedAttribution::shares)
      .def_ro("total", &SeedAttribution::total)
      .def_ro("covered_fraction", &SeedAttribution::coveredFraction)
      .def_ro("omitted_count", &SeedAttribution::omittedCount)
      .def_ro("time", &SeedAttribution::time)
      .def_ro("seed_provenance", &SeedAttribution::seedProvenance)
      // The same three a Ranking carries, under the same names: an attributed exposure is the
      // ranking's number seen from the other side, so it is understated by the same amount.
      .def_ro("unmodeled_energy_fraction", &SeedAttribution::unmodeledEnergyFraction)
      .def_ro("mean_optical_depth", &SeedAttribution::meanOpticalDepth,
              "Exposure and photon fluence only: the air path in mean free paths at the "
              "energies carrying this answer. Past about 0.5 with buildup 1.0, scattered "
              "photons are a large omission.")
      .def_ro("buildup", &SeedAttribution::buildup)
      .def_ro("unmodeled_continuum", &SeedAttribution::unmodeledContinuum,
              "Exposure and photon only: emitters carrying photon energy NuSIFT does not "
              "model, so their contribution to this figure is understated.")
      .def_prop_ro("labels",
                   [](const SeedAttribution& a) {
                     std::vector<std::string> names;
                     for (const SeedShare& s : a.shares) {
                       names.push_back(s.label);
                     }
                     return names;
                   })
      .def("__len__", [](const SeedAttribution& a) { return a.shares.size(); })
      .def("__repr__", [](const SeedAttribution& a) {
        return "<SeedAttribution " + std::to_string(a.shares.size()) + " of " +
               std::to_string(a.shares.size() + a.omittedCount) + " seeds>";
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
      .def_prop_ro(
          "time_ends",
          [](nb::handle self) {
            const ResponseTable& t = nb::cast<const ResponseTable&>(self);
            return view1d(t.timeEnds.data(), t.timeEnds.size(), self);
          },
          "Window ends for an interval table; empty for an instant one.")
      .def_prop_ro("domain",
                   [](const ResponseTable& t) { return std::string(domainName(t.domain)); })
      .def_prop_ro("totals",
                   [](nb::handle self) {
                     const ResponseTable& t = nb::cast<const ResponseTable&>(self);
                     return view1d(t.totals.data(), t.totals.size(), self);
                   })
      .def_prop_ro(
          "mean_optical_depth",
          [](nb::handle self) {
            const ResponseTable& t = nb::cast<const ResponseTable&>(self);
            return view1d(t.meanOpticalDepth.data(), t.meanOpticalDepth.size(), self);
          },
          "Exposure and photon fluence only, per time: the air path in mean free paths at the "
          "energies carrying the answer. Empty for activity and for photon strength, which "
          "uses no geometry.")
      .def_prop_ro(
          "values",
          [](nb::handle self) {
            const ResponseTable& t = nb::cast<const ResponseTable&>(self);
            return view2d(t.values.data(), static_cast<std::size_t>(t.timeCount()),
                          static_cast<std::size_t>(t.contributorCount()), self);
          },
          "(times, contributors). A zero-copy view, not a copy.")
      .def_prop_ro("unit",
                   [](const ResponseTable& t) {
                     return t.unitLabel.empty() ? std::string(unitName(t.unit)) : t.unitLabel;
                   })
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
          [](const ResponseTable& table, int min_samples, const nb::object& refine,
             double tolerance) {
            const ResponseEvaluator* evaluator = evaluatorFrom(refine);
            return evaluator != nullptr
                       ? dominanceWindows(table, *evaluator, min_samples, toleranceOf(tolerance))
                       : dominanceWindows(table, min_samples);
          },
          "min_samples"_a = 2, "refine"_a = nb::none(), "tolerance"_a = 1.0e-6,
          "Who leads, and over which windows. A boundary is a located event like any other: "
          "pass an `evaluator` as `refine` to place it by re-solving inside its bracket "
          "instead of interpolating across it.")
      .def(
          "crossings",
          [](const ResponseTable& table, double level, const nb::object& of,
             const nb::object& ratio, const nb::object& refine, double tolerance) {
            return crossings(seriesFor(table, of, ratio, evaluatorFrom(refine)), level,
                             toleranceOf(tolerance));
          },
          "level"_a, "of"_a = nb::none(), "ratio"_a = nb::none(), "refine"_a = nb::none(),
          "tolerance"_a = 1.0e-6,
          "Every crossing of `level` the grid observed, in time order. `of` follows one "
          "contributor instead of the total; `ratio` a pair of them. An excursion between two "
          "samples leaves no sign change and is not found -- the grid decides what is seen. "
          "`refine` takes an `evaluator`, which narrows each event by re-solving rather than "
          "interpolating; `TrajectoryEvent.refined` says which happened.")
      .def(
          "extrema",
          [](const ResponseTable& table, const nb::object& of, const nb::object& ratio,
             const nb::object& refine, double tolerance) {
            return extrema(seriesFor(table, of, ratio, evaluatorFrom(refine)),
                           toleranceOf(tolerance));
          },
          "of"_a = nb::none(), "ratio"_a = nb::none(), "refine"_a = nb::none(),
          "tolerance"_a = 1.0e-6,
          "Interior maxima and minima. A turn needs three samples to be told from a monotone "
          "run, so one in the first or last interval is not reported.")
      .def(
          "windows_above",
          [](const ResponseTable& table, double level, const nb::object& of,
             const nb::object& ratio, const nb::object& refine, double tolerance) {
            return windowsAbove(seriesFor(table, of, ratio, evaluatorFrom(refine)), level,
                                toleranceOf(tolerance));
          },
          "level"_a, "of"_a = nb::none(), "ratio"_a = nb::none(), "refine"_a = nb::none(),
          "tolerance"_a = 1.0e-6,
          "The stretches at or above `level`. An edge the grid never observed is flagged open "
          "rather than clipped to the grid's own endpoint.")
      .def(
          "windows_below",
          [](const ResponseTable& table, double level, const nb::object& of,
             const nb::object& ratio, const nb::object& refine, double tolerance) {
            return windowsBelow(seriesFor(table, of, ratio, evaluatorFrom(refine)), level,
                                toleranceOf(tolerance));
          },
          "level"_a, "of"_a = nb::none(), "ratio"_a = nb::none(), "refine"_a = nb::none(),
          "tolerance"_a = 1.0e-6,
          "The complement of windows_above, and the shape a \"when is it safe\" question takes.")
      .def_prop_ro(
          "pack_coverage",
          [](nb::handle self) {
            const ResponseTable& t = nb::cast<const ResponseTable&>(self);
            return view1d(t.packCoverage.data(), t.packCoverage.size(), self);
          },
          "Pack metrics only: the share of the inventory the pack carries a coefficient for, at "
          "each time, measured in the quantity its coefficients multiply. Empty otherwise.")
      .def("__repr__", [](const ResponseTable& t) {
        return "<ResponseTable " + std::to_string(t.timeCount()) + " times x " +
               std::to_string(t.contributorCount()) + " contributors, " + unitName(t.unit) + ">";
      });

  // A located event is the one quantity NuSIFT reports that is not exact, and the bracket is
  // why: it is the grid interval that observed the event, and `located_to_s` is how tightly the
  // instant was placed inside it. A consumer reading `time_s` alone has thrown away the error
  // bar.
  nb::class_<TrajectoryEvent>(m, "TrajectoryEvent")
      .def_prop_ro("kind",
                   [](const TrajectoryEvent& e) { return std::string(eventKindName(e.kind)); })
      .def_ro("time_s", &TrajectoryEvent::timeSeconds)
      .def_ro("value", &TrajectoryEvent::value)
      .def_ro("bracket_start_s", &TrajectoryEvent::bracketStartSeconds)
      .def_ro("bracket_end_s", &TrajectoryEvent::bracketEndSeconds)
      .def_ro("located_to_s", &TrajectoryEvent::locatedToSeconds,
              "The width the location was narrowed to: the honest error bar on time_s.")
      .def_ro("refined", &TrajectoryEvent::refined,
              "True when real evaluations placed it, False when it is an interpolation between "
              "two samples.")
      .def_ro("converged", &TrajectoryEvent::converged)
      .def("__repr__", [](const TrajectoryEvent& e) {
        return "<TrajectoryEvent " + std::string(eventKindName(e.kind)) + " at " +
               formatDuration(e.timeSeconds) + " within " + formatDuration(e.locatedToSeconds) +
               ">";
      });

  nb::class_<LevelWindow>(m, "LevelWindow")
      .def_ro("start_s", &LevelWindow::startSeconds)
      .def_ro("end_s", &LevelWindow::endSeconds)
      .def_ro("entry_observed", &LevelWindow::entryObserved,
              "False when the grid began already inside: the start is a bound, not a crossing.")
      .def_ro("exit_observed", &LevelWindow::exitObserved,
              "False when the grid ended still inside.")
      .def("__repr__", [](const LevelWindow& w) {
        return "<LevelWindow " + formatDuration(w.startSeconds) + " to " +
               formatDuration(w.endSeconds) + (w.entryObserved && w.exitObserved ? "" : " (open)") +
               ">";
      });

  // Whether an event is NARROWED or merely interpolated turns on whether the curve can be
  // asked what it does between two samples. For a response that is a single-time solve, and
  // this is the object that performs it -- passed to any search as `refine=`.
  //
  // It has to describe the same response the table does: the same store, inventory, metric,
  // aggregate, units and geometry. A mismatch is refused rather than quietly refining an event
  // on one curve with values taken from another.
  nb::class_<ResponseEvaluator>(m, "ResponseEvaluator")
      .def_prop_ro("solves", &ResponseEvaluator::solves,
                   "Single-time solves spent so far. What refinement cost is not visible in "
                   "the events it produced.")
      .def("__repr__", [](const ResponseEvaluator& e) {
        return "<ResponseEvaluator " + std::string(metricName(e.spec().metric)) + " in " +
               std::string(unitName(e.spec().unit)) + ", " + std::to_string(e.solves()) +
               " solves>";
      });

  m.def(
      "evaluator",
      [](const NuclearData& data, const Inventory& inventory, const std::string& metric,
         const std::string& by, const std::string& units,
         const exposure::PointSourceGeometry& geometry, bool prune, int cram_order) {
        ResponseSpec spec;
        spec.metric = metricFrom(metric);
        spec.aggregate = aggregateFrom(by);
        spec.unit = requireUnit(units, spec.metric, Domain::Instant);
        spec.geometry = geometry;
        DecayOptions options;
        options.prune = prune;
        options.order = cram_order == 16 ? CramOrder::Order16 : CramOrder::Order48;
        return ResponseEvaluator(data, inventory, spec, options);
      },
      "data"_a, "inventory"_a, "metric"_a = "activity", "by"_a = "nuclide", "units"_a = "",
      "geometry"_a = exposure::PointSourceGeometry{}, "prune"_a = true, "cram_order"_a = 48,
      // The evaluator reads the store and the inventory on every call, so Python must not
      // collect either while it lives.
      nb::keep_alive<0, 1>(), nb::keep_alive<0, 2>(),
      "A re-solver for the response described by these arguments, to hand to a search as "
      "`refine=`. Build it from the same arguments the table was built from.");

  // A curve to locate events on, when that curve is not a column of a response table. The
  // searches are the same functions the table methods call, so a task curve and a rate curve
  // are found the same way and reported the same way.
  nb::class_<EventSeries>(m, "EventSeries")
      .def_prop_ro("times",
                   [](nb::handle self) {
                     const EventSeries& s = nb::cast<const EventSeries&>(self);
                     return view1d(s.times.data(), s.times.size(), self);
                   })
      .def_prop_ro("values",
                   [](nb::handle self) {
                     const EventSeries& s = nb::cast<const EventSeries&>(self);
                     return view1d(s.values.data(), s.values.size(), self);
                   })
      .def_prop_ro(
          "refines", [](const EventSeries& s) { return static_cast<bool>(s.evaluate); },
          "True when the curve can be evaluated between samples, so events are "
          "narrowed rather than interpolated.")
      .def(
          "crossings",
          [](const EventSeries& series, double level, double tolerance) {
            return crossings(series, level, toleranceOf(tolerance));
          },
          "level"_a, "tolerance"_a = 1.0e-6)
      .def(
          "extrema",
          [](const EventSeries& series, double tolerance) {
            return extrema(series, toleranceOf(tolerance));
          },
          "tolerance"_a = 1.0e-6)
      .def(
          "windows_above",
          [](const EventSeries& series, double level, double tolerance) {
            return windowsAbove(series, level, toleranceOf(tolerance));
          },
          "level"_a, "tolerance"_a = 1.0e-6)
      .def(
          "windows_below",
          [](const EventSeries& series, double level, double tolerance) {
            return windowsBelow(series, level, toleranceOf(tolerance));
          },
          "level"_a, "tolerance"_a = 1.0e-6)
      .def("__repr__", [](const EventSeries& s) {
        return "<EventSeries " + std::to_string(s.times.size()) + " samples" +
               (s.evaluate ? ", refinable>" : ">");
      });

  // The fixed-duration task: what a job of a given length accrues, against WHEN it starts.
  // Every sample is an exact interval integral, so this is the expensive curve in the library
  // -- two to three solves a sample -- and the reason to pay for it is that its extrema answer
  // "when should this be done" and its crossings answer "when does it fit the budget", neither
  // of which is a question about the rate curve.
  m.def(
      "task_series",
      [](const NuclearData& data, const Inventory& inventory, const nb::object& starts,
         const nb::object& duration, const std::string& metric, const std::string& by,
         const std::string& units, const exposure::PointSourceGeometry& geometry, bool refine,
         int threads, bool prune, int cram_order) {
        std::vector<double> startTimes;
        for (const nb::handle item : starts) {
          startTimes.push_back(timeFrom(nb::borrow<nb::object>(item)));
        }
        ResponseSpec spec;
        spec.metric = metricFrom(metric);
        spec.aggregate = aggregateFrom(by);
        spec.unit = requireUnit(units, spec.metric, Domain::Interval);
        spec.geometry = geometry;
        DecayOptions options;
        options.threads = threads;
        options.prune = prune;
        options.order = cram_order == 16 ? CramOrder::Order16 : CramOrder::Order48;
        return taskSeries(data, inventory, spec, startTimes, timeFrom(duration), options, refine);
      },
      "data"_a, "inventory"_a, "starts"_a, "duration"_a, "metric"_a = "activity",
      "by"_a = "nuclide", "units"_a = "", "geometry"_a = exposure::PointSourceGeometry{},
      "refine"_a = false, "threads"_a = 0, "prune"_a = true, "cram_order"_a = 48,
      // Like the evaluator, the refining form holds the store and the inventory and reads them
      // during the search.
      nb::keep_alive<0, 1>(), nb::keep_alive<0, 2>(),
      "What a task of `duration` accrues, sampled at each start in `starts`. The units are "
      "interval units -- roentgen, decays -- because the values are accrued totals.");

  // The converse inversion. task_series fixes the length and searches for the start; this fixes
  // the start and solves for the LENGTH, which is a root-find in one window rather than a search
  // along a curve of them, and no sampling of the task curve produces it.
  nb::class_<StayTime>(m, "StayTime", "How long a stay can run before it spends a budget.")
      .def_ro("start_s", &StayTime::startSeconds)
      .def_ro("budget", &StayTime::budget)
      .def_ro("max_duration_s", &StayTime::maxDurationSeconds,
              "The longest stay considered. Part of the question: a search with no ceiling "
              "would have to invent one.")
      .def_prop_ro(
          "duration_s",
          [](const StayTime& s) -> nb::object {
            // None rather than a number when the budget was never spent. A caller reading 0.0
            // here would get "leave immediately" for "stay as long as you like", which is the
            // worst available way to be wrong about this question, so the type refuses it.
            return s.bounded ? nb::cast(s.durationSeconds) : nb::none();
          },
          "How long the budget lasts, in seconds, or None when it is not spent inside "
          "max_duration_s.")
      .def_ro("bounded", &StayTime::bounded)
      .def_ro("accrued_at_max", &StayTime::accruedAtMax,
              "What a stay of the full max_duration_s accrues. The useful number when not "
              "bounded: not 'never spent' but 'a stay this long costs this much of it'.")
      .def_ro("located_to_s", &StayTime::locatedToSeconds,
              "The width the duration was narrowed to -- the honest error bar on it.")
      .def_ro("converged", &StayTime::converged)
      .def_ro("samples", &StayTime::samples, "Interval integrals spent finding it.")
      .def("__repr__", [](const StayTime& s) {
        return "<StayTime from " + formatDuration(s.startSeconds) + ": " +
               (s.bounded ? formatDuration(s.durationSeconds) : std::string("budget not spent")) +
               ">";
      });

  m.def(
      "stay_time",
      [](const NuclearData& data, const Inventory& inventory, const nb::object& at, double budget,
         const nb::object& max_stay, const std::string& metric, const std::string& units,
         const exposure::PointSourceGeometry& geometry, double tolerance, int threads, bool prune,
         int cram_order) {
        ResponseSpec spec;
        spec.metric = metricFrom(metric);
        // The total is what a budget constrains, so the aggregate is fixed rather than offered:
        // which nuclide spends it does not change how long it lasts, and `rank` is where that
        // question lives.
        spec.aggregate = Aggregate::Nuclide;
        spec.unit = requireUnit(units, spec.metric, Domain::Interval);
        spec.geometry = geometry;
        DecayOptions options;
        options.threads = threads;
        options.prune = prune;
        options.order = cram_order == 16 ? CramOrder::Order16 : CramOrder::Order48;
        const nb::gil_scoped_release release;
        return stayTime(data, inventory, spec, timeFrom(at), budget, timeFrom(max_stay), options,
                        toleranceOf(tolerance));
      },
      "data"_a, "inventory"_a, "at"_a, "budget"_a, "max_stay"_a = "24h", "metric"_a = "activity",
      "units"_a = "", "geometry"_a = exposure::PointSourceGeometry{}, "tolerance"_a = 1.0e-6,
      "threads"_a = 0, "prune"_a = true, "cram_order"_a = 48,
      "How long a stay beginning at `at` can run before it accrues `budget`. The units are "
      "interval units -- Sv, R, decays -- because a budget is an accrued total and not a rate. "
      "Costs tens of solves: each step of the root-find is an exact interval integral.");

  // A published coefficient table, read at runtime rather than compiled in. What makes it a
  // pack rather than a dictionary is everything in its header: the quantity, the unit, what the
  // coefficient multiplies, the scenario it was tabulated under, and the version -- all of which
  // travel into the answer, because a transport index from the 2012 edition of SSR-6 is not the
  // same number as one from a later edition.
  nb::class_<CoefficientPack>(m, "CoefficientPack",
                              "A versioned table of per-nuclide coefficients.")
      .def_prop_ro("name", [](const CoefficientPack& p) { return p.provenance().name; })
      .def_prop_ro("version", [](const CoefficientPack& p) { return p.provenance().version; })
      .def_prop_ro("quantity", [](const CoefficientPack& p) { return p.provenance().quantity; })
      .def_prop_ro("unit", [](const CoefficientPack& p) { return p.provenance().unit; })
      .def_prop_ro("scenario", [](const CoefficientPack& p) { return p.provenance().scenario; })
      .def_prop_ro("source", [](const CoefficientPack& p) { return p.provenance().source; })
      .def_prop_ro(
          "basis",
          [](const CoefficientPack& p) { return std::string(packBasisName(p.provenance().basis)); })
      .def_prop_ro("folds_progeny",
                   [](const CoefficientPack& p) { return p.provenance().foldsProgeny; })
      .def_prop_ro(
          "per", [](const CoefficientPack& p) { return p.provenance().per; },
          "Concentration packs only: the denominator its coefficients are per -- m2, "
          "m3 or kg. Empty for every other basis.")
      .def_prop_ro("size", &CoefficientPack::size)
      .def(
          "covers",
          [](const CoefficientPack& p, const std::string& nuclide) {
            return p.covers(requireNuclideName(nuclide).key());
          },
          "nuclide"_a,
          "Whether the pack carries a row for this nuclide -- which is NOT the same question "
          "as whether its coefficient is zero. A zero written in the file means a limit that "
          "does not bind; a nuclide with no row means the pack cannot speak for it at all.")
      .def(
          "coefficient",
          [](const CoefficientPack& p, const std::string& nuclide) {
            return p.coefficient(requireNuclideName(nuclide).key());
          },
          "nuclide"_a, "The coefficient as the file states it, or zero when there is no row.")
      .def(
          "folded_into",
          [](const CoefficientPack& p, const std::string& nuclide) {
            return nuclideNames(p.foldedInto(requireNuclideName(nuclide).key()));
          },
          "nuclide"_a,
          "The parents whose coefficients already account for this nuclide. Several, because "
          "decay chains nest: SSR-6 folds Tl-208 into Bi-212, Pb-212, Ra-224 and Th-228.")
      .def("__repr__", [](const CoefficientPack& p) {
        return "<CoefficientPack " + p.provenance().name + " " + p.provenance().version + ", " +
               std::to_string(p.size()) + " nuclides>";
      });

  m.def(
      "load_pack", [](const std::string& path) { return CoefficientPack::open(path); }, "path"_a,
      "Read a coefficient pack: CSV with a header carrying its quantity, unit, basis, scenario "
      "and version.");

  // A pack resolved against a seed. Which coefficient applies to a folded daughter is a
  // question about the inventory -- yttrium-90 alone is limited by its own value, and yttrium-90
  // with strontium-90 is inside its parent's -- so the resolution happens once, here, and the
  // response layer post-multiplies a fixed vector as it does for every other metric.
  nb::class_<ResolvedPack>(m, "ResolvedPack", "A pack resolved against a seed inventory.")
      .def("__repr__", [](const ResolvedPack& r) {
        return "<ResolvedPack " + r.pack->provenance().name + " over " +
               std::to_string(r.weights.size()) + " nuclides>";
      });

  m.def(
      "resolve_pack",
      [](const CoefficientPack& pack, const NuclearData& data, const Inventory& seed,
         double extent) {
        PackExtent spread;
        if (extent > 0.0) {
          spread.value = extent;
          // The unit is the pack's own, so a caller cannot disagree with it about whether the
          // number is an area or a volume.
          spread.unit = pack.provenance().per;
        }
        return resolvePack(pack, data, seed, spread);
      },
      "pack"_a, "data"_a, "seed"_a, "extent"_a = 0.0, nb::keep_alive<0, 1>(),
      "Resolve a pack against the inventory it will be used with, which is what decides whether "
      "a folded daughter takes its parent's coefficient or its own. `extent` is the volume, "
      "area or mass the inventory is spread through, required by a concentration pack and "
      "refused by every other.");

  m.def(
      "response",
      [](const NuclearData& data, const DecayResult& result, const std::string& metric,
         const std::string& by, const std::string& units,
         const exposure::PointSourceGeometry& geometry, const ResolvedPack* pack) {
        ResponseSpec spec;
        // A pack IS the metric, so naming both would be two answers to one question.
        spec.metric = pack != nullptr ? Metric::Pack : metricFrom(metric);
        spec.aggregate = aggregateFrom(by);
        spec.unit = requireUnit(units, spec.metric, Domain::Instant);
        spec.geometry = geometry;
        spec.pack = pack;
        return buildResponse(data, result, spec);
      },
      "data"_a, "result"_a, "metric"_a = "activity", "by"_a = "nuclide", "units"_a = "",
      "geometry"_a = exposure::PointSourceGeometry{}, "pack"_a = nb::none(),
      "Turn a decay result into a table of per-contributor values. `pack` takes a resolved "
      "coefficient pack and becomes the metric.");

  // The same call over an interval result gives the interval domain: one row of totals --
  // decays, or roentgen accrued -- with the units gated accordingly, so `units="Bq"` is
  // refused here for the same reason the CLI refuses it on `integrate`.
  m.def(
      "response",
      [](const NuclearData& data, const IntervalResult& result, const std::string& metric,
         const std::string& by, const std::string& units,
         const exposure::PointSourceGeometry& geometry, const ResolvedPack* pack) {
        ResponseSpec spec;
        spec.metric = pack != nullptr ? Metric::Pack : metricFrom(metric);
        spec.aggregate = aggregateFrom(by);
        spec.unit = requireUnit(units, spec.metric, Domain::Interval);
        spec.geometry = geometry;
        spec.pack = pack;
        return buildIntervalResponse(data, result.nuclideKeys, result.integratedAtoms, result.t1,
                                     result.t2, spec);
      },
      "data"_a, "result"_a, "metric"_a = "activity", "by"_a = "nuclide", "units"_a = "",
      "geometry"_a = exposure::PointSourceGeometry{}, "pack"_a = nb::none(),
      "Turn an interval result into a one-row table of per-contributor totals over the window.");

  // --- what the EVALUATED data does to the answer ----------------------------
  //
  // The other parameter class. `uncertainty` propagates what the assay said at fixed nuclear
  // data; this asks how much the answer rests on the half-lives the evaluation supplied, at a
  // fixed inventory. A complete error budget wants both and they do not overlap.
  nb::class_<DecaySensitivity>(m, "DecaySensitivity")
      .def_ro("label", &DecaySensitivity::label)
      .def_ro("half_life_s", &DecaySensitivity::halfLifeSeconds)
      .def_ro("implicit", &DecaySensitivity::implicit,
              "The response moving because the decay MATRIX moved. cram's adjoint returns this "
              "term alone.")
      .def_ro("explicit_weight", &DecaySensitivity::explicitWeight,
              "... because the WEIGHT moved, since every built-in weight is lambda times "
              "something lambda-independent.")
      .def_ro("basis", &DecaySensitivity::basis,
              "... because the SEED moved, for whatever part of the row was given as an "
              "activity: n0 = A0/lambda.")
      .def_ro("total", &DecaySensitivity::total)
      .def_ro("elasticity", &DecaySensitivity::elasticity,
              "lambda * (dR/dlambda) / R. THE number to read: dimensionless, and it stays "
              "interpretable where the terms cancel and the total does not.")
      .def_prop_ro(
          "relative_uncertainty",
          [](const DecaySensitivity& d) -> nb::object {
            return d.relativeUncertainty > 0.0 ? nb::cast(d.relativeUncertainty) : nb::none();
          },
          "Evaluated 1-sigma on the half-life as a fraction, or None when none is staged.")
      .def_ro("sigma_contribution", &DecaySensitivity::sigmaContribution)
      .def("__repr__", [](const DecaySensitivity& d) {
        return "<DecaySensitivity " + d.label + " elasticity " + shortestRoundTrip(d.elasticity) +
               ">";
      });

  nb::class_<BranchingSensitivity>(m, "BranchingSensitivity")
      .def_ro("parent", &BranchingSensitivity::parent)
      .def_ro("daughter", &BranchingSensitivity::daughter)
      .def_ro("mode", &BranchingSensitivity::mode)
      .def_ro("branching", &BranchingSensitivity::branching)
      .def_ro("sigma", &BranchingSensitivity::sigma)
      .def_ro("total", &BranchingSensitivity::total)
      .def_ro("elasticity", &BranchingSensitivity::elasticity)
      .def("__repr__", [](const BranchingSensitivity& b) {
        return "<BranchingSensitivity " + b.parent + " " + b.mode + " -> " + b.daughter + ">";
      });

  nb::class_<BranchingBlock>(m, "BranchingBlock",
                             "One nuclide's branchings, whose modes sum to one.")
      .def_ro("parent", &BranchingBlock::parent)
      .def_ro("modes", &BranchingBlock::modes)
      .def_ro("diagonal_variance", &BranchingBlock::diagonalVariance,
              "What a diagonal treatment would report. Kept so the constraint's effect is "
              "visible rather than asserted.")
      .def_ro("constrained_variance", &BranchingBlock::constrainedVariance,
              "e^T C e with C the covariance the sum-to-one constraint forces. The honest one.")
      .def("__repr__", [](const BranchingBlock& b) {
        return "<BranchingBlock " + b.parent + " " + std::to_string(b.modes) + " modes>";
      });

  nb::class_<DecaySensitivities>(m, "DecaySensitivities",
                                 "How much a response rests on the evaluated half-lives.")
      .def_ro("time", &DecaySensitivities::time)
      .def_ro("response", &DecaySensitivities::response)
      .def_ro("nuclides", &DecaySensitivities::nuclides)
      .def_ro("relative_norm", &DecaySensitivities::relativeNorm,
              "Root-sum-square of the contributions, as a fraction of R. A sensitivity NORM and "
              "not an error budget: it takes Sigma diagonal, and evaluated half-lives are not "
              "independent of the branchings and yields fitted alongside them.")
      .def_ro("covered_fraction", &DecaySensitivities::coveredFraction)
      .def_ro("branchings", &DecaySensitivities::branchings)
      .def_ro("branching_blocks", &DecaySensitivities::branchingBlocks)
      .def_ro("branching_norm", &DecaySensitivities::branchingNorm,
              "Branching contribution as a fraction of R, with the sum-to-one constraint "
              "imposed. NOT a diagonal figure: the correlation is derived from the constraint "
              "rather than imported, because ENDF carries no covariance for decay data.")
      .def_ro("branching_norm_diagonal", &DecaySensitivities::branchingNormDiagonal,
              "The same taking the branchings as independent. The DIFFERENCE from the above is "
              "the result: it says what ignoring the constraint would have cost.")
      .def_ro("with_uncertainty", &DecaySensitivities::withUncertainty)
      .def_ro("without_uncertainty", &DecaySensitivities::withoutUncertainty)
      .def_ro("end_refinements", &DecaySensitivities::endRefinements)
      .def_ro("refinement_capped", &DecaySensitivities::refinementCapped,
              "True when the cap bound before the smallest quadrature piece reached the shortest "
              "removal time. The answer is then UNDER-REFINED by cram's own criterion and may be "
              "wrong by tens of percent, with nothing in the numbers to show it.")
      .def_ro("solves", &DecaySensitivities::solves)
      .def_ro("shortest_removal_s", &DecaySensitivities::shortestRemovalSeconds)
      .def("__len__", [](const DecaySensitivities& d) { return d.nuclides.size(); })
      .def("__repr__", [](const DecaySensitivities& d) {
        return "<DecaySensitivities " + std::to_string(d.nuclides.size()) + " nuclides, norm " +
               shortestRoundTrip(d.relativeNorm) + ">";
      });

  m.def(
      "decay_sensitivity",
      [](const NuclearData& data, const Inventory& inventory, const nb::object& at,
         const std::string& metric, const std::string& units,
         const exposure::PointSourceGeometry& geometry, const ResolvedPack* pack, int intervals,
         int max_refinements, int threads, bool prune, int cram_order) {
        ResponseSpec spec;
        spec.metric = pack != nullptr ? Metric::Pack : metricFrom(metric);
        spec.aggregate = Aggregate::Nuclide;
        spec.unit = requireUnit(units, spec.metric, Domain::Instant);
        spec.geometry = geometry;
        spec.pack = pack;
        SensitivityOptions sensitivity;
        sensitivity.scheduleIntervals = intervals;
        sensitivity.maxEndRefinements = max_refinements;
        DecayOptions options;
        options.threads = threads;
        options.prune = prune;
        options.order = cram_order == 16 ? CramOrder::Order16 : CramOrder::Order48;
        const double time = timeFrom(at);
        const nb::gil_scoped_release release;
        return decaySensitivities(data, inventory, time, spec, sensitivity, options);
      },
      "data"_a, "inventory"_a, "at"_a, "metric"_a = "activity", "units"_a = "",
      "geometry"_a = exposure::PointSourceGeometry{}, "pack"_a = nb::none(), "intervals"_a = 1,
      "max_refinements"_a = 30, "threads"_a = 0, "prune"_a = true, "cram_order"_a = 48,
      "dR/dlambda for every nuclide the seed reaches. Costs quadrature over the forward and "
      "adjoint trajectories -- hundreds of solves where a ranking costs one -- and the "
      "refinement is chosen automatically from the shortest removal time, because cram's "
      "defaults are silently wrong by tens of percent on a decay chain.");

  // --- the error bar the assay puts on the answer ----------------------------
  //
  // The one uncertainty question that needs no evaluated data. R is LINEAR in the seed, so
  // sigma_R^2 = g^T Sigma g is exact rather than first-order -- no expansion, no estimated
  // derivative. What is assumed is the shape of Sigma, and the report says so.
  nb::class_<SeedUncertainty>(m, "SeedUncertainty")
      .def_ro("label", &SeedUncertainty::label)
      .def_ro("assay", &SeedUncertainty::assay)
      .def_ro("carried_s", &SeedUncertainty::carriedSeconds)
      .def_ro("seed_atoms", &SeedUncertainty::seedAtoms)
      .def_prop_ro(
          "sigma_atoms",
          [](const SeedUncertainty& s) -> nb::object {
            // None rather than 0.0: "stated no uncertainty" and "stated an uncertainty of
            // zero" are different claims, and only one is ever true of a measurement.
            return s.sigmaAtoms > 0.0 ? nb::cast(s.sigmaAtoms) : nb::none();
          },
          "1-sigma on this row in atoms, or None when the row stated none.")
      .def_ro("importance", &SeedUncertainty::importance, "dR/dn0 for this row, at T + carried.")
      .def_ro("share", &SeedUncertainty::share, "This row's exact share of R.")
      .def_ro("sigma_contribution", &SeedUncertainty::sigmaContribution,
              "The standard deviation this row alone puts on R. Not additive -- the variance "
              "fractions are what sum to one.")
      .def_ro("variance_fraction", &SeedUncertainty::varianceFraction,
              "Share of the VARIANCE. The assay-planning number: it says which measurement to "
              "improve, and a row with a negligible share of R can dominate it.")
      .def("__repr__", [](const SeedUncertainty& s) {
        return "<SeedUncertainty " + s.label + " " + shortestRoundTrip(s.varianceFraction) +
               " of variance>";
      });

  nb::class_<ResponseUncertainty>(m, "ResponseUncertainty",
                                  "An error bar on a response, and what it rests on.")
      .def_prop_ro("unit",
                   [](const ResponseUncertainty& u) { return std::string(unitName(u.unit)); })
      .def_ro("time", &ResponseUncertainty::time)
      .def_ro("response", &ResponseUncertainty::response)
      .def_ro("sigma", &ResponseUncertainty::sigma)
      .def_ro("relative", &ResponseUncertainty::relative)
      .def_ro("covered_fraction", &ResponseUncertainty::coveredFraction,
              "Share of the response coming from rows that stated an uncertainty. An error bar "
              "propagated from rows holding half the answer is not an error bar on the answer.")
      .def_ro("rows_with_sigma", &ResponseUncertainty::rowsWithSigma)
      .def_ro("rows_without_sigma", &ResponseUncertainty::rowsWithoutSigma)
      .def_ro("seeds", &ResponseUncertainty::seeds, "Ranked by variance fraction, largest first.")
      .def("__repr__", [](const ResponseUncertainty& u) {
        return "<ResponseUncertainty " + shortestRoundTrip(u.response) + " +/- " +
               shortestRoundTrip(u.sigma) + " " + unitName(u.unit) + ">";
      });

  m.def(
      "uncertainty",
      [](const NuclearData& data, const std::vector<AssayGroup>& assays, const nb::object& at,
         const nb::object& epoch, const std::string& metric, const std::string& units,
         const exposure::PointSourceGeometry& geometry, const ResolvedPack* pack, int threads,
         bool prune, int cram_order) {
        ResponseSpec spec;
        spec.metric = pack != nullptr ? Metric::Pack : metricFrom(metric);
        spec.aggregate = Aggregate::Nuclide;
        spec.unit = requireUnit(units, spec.metric, Domain::Instant);
        spec.geometry = geometry;
        spec.pack = pack;
        DecayOptions options;
        options.threads = threads;
        options.prune = prune;
        options.order = cram_order == 16 ? CramOrder::Order16 : CramOrder::Order48;
        const double time = timeFrom(at);
        const double when = epoch.is_none() ? latestAssayDate(assays)
                            : nb::isinstance<nb::str>(epoch)
                                ? parseCalendarDate(nb::cast<std::string>(epoch))
                                : nb::cast<double>(epoch);
        const nb::gil_scoped_release release;
        return responseUncertainty(data, assays, when, time, spec, options);
      },
      "data"_a, "assays"_a, "at"_a, "epoch"_a = nb::none(), "metric"_a = "activity", "units"_a = "",
      "geometry"_a = exposure::PointSourceGeometry{}, "pack"_a = nb::none(), "threads"_a = 0,
      "prune"_a = true, "cram_order"_a = 48,
      "The error bar the assay uncertainties put on the response at `at`, and which row "
      "dominates it. Takes ASSAYS rather than an inventory -- read them with read_assays() -- "
      "because a sigma cannot ride on a reconciled inventory: a diagonal covariance at assay is "
      "not diagonal at the epoch. Costs one adjoint per assay.");

  // --- the error bar the EVALUATED yields put on a fission source ------------
  //
  // The other half of the same question, and the only one that needs data no evaluation
  // publishes. A fission seed is n0 = N_f Y, so the propagation is the one above with a
  // different Sigma -- and the off-diagonal of that Sigma is imported rather than assumed.
  nb::class_<YieldCorrelation>(m, "YieldCorrelation",
                               "A correlation between independent fission yields, read from a "
                               "published FYCoM matrix.")
      .def_static("read", &YieldCorrelation::read, "path"_a, "library"_a = "ENDF/B-VIII.0",
                  "Read a FYCoM CORRELATION csv (a *_corr.csv). Only the correlation is taken; "
                  "the variances come from the store's own staged sigma_Y, which is what the "
                  "paper's own 'normalized' covariance does. `library` is the evaluation the "
                  "matrix was BUILT FOR, carried into the answer so the pairing with the "
                  "store's edition is declared rather than assumed.")
      .def_prop_ro("size", &YieldCorrelation::size)
      .def_prop_ro("system", [](const YieldCorrelation& c) { return c.provenance().system; })
      .def_prop_ro("library", [](const YieldCorrelation& c) { return c.provenance().library; })
      .def_prop_ro("citation", [](const YieldCorrelation& c) { return c.provenance().citation; })
      .def("__len__", &YieldCorrelation::size)
      .def("__repr__", [](const YieldCorrelation& c) {
        return "<YieldCorrelation " + c.provenance().system + " " + std::to_string(c.size()) +
               " products from " + c.provenance().library + ">";
      });

  nb::class_<YieldContribution>(m, "YieldContribution",
                                "One fission product's yield, and what its uncertainty does.")
      .def_ro("label", &YieldContribution::label)
      .def_ro("yield_per_fission", &YieldContribution::yield)
      .def_prop_ro(
          "sigma_yield",
          [](const YieldContribution& p) -> nb::object {
            // None rather than 0.0, as the assay side does it: an evaluation stating no
            // uncertainty and one stating zero are different claims.
            return p.sigmaYield > 0.0 ? nb::cast(p.sigmaYield) : nb::none();
          },
          "Absolute 1-sigma on the yield, or None when the evaluation states none.")
      .def_ro("seed_atoms", &YieldContribution::seedAtoms)
      .def_ro("importance", &YieldContribution::importance)
      .def_ro("share", &YieldContribution::share)
      .def_ro("sigma_contribution", &YieldContribution::sigmaContribution)
      .def_ro("variance_fraction", &YieldContribution::varianceFraction,
              "Share of the DIAGONAL variance -- the only per-product share there is, since "
              "half of every off-diagonal term belongs to each of two products.")
      .def_ro("correlated", &YieldContribution::correlated,
              "Whether the imported matrix carries this product.")
      .def("__repr__", [](const YieldContribution& p) {
        return "<YieldContribution " + p.label + " " + shortestRoundTrip(p.varianceFraction) +
               " of variance>";
      });

  nb::class_<YieldUncertainty>(m, "YieldUncertainty",
                               "The error bar the evaluated yields put on a response.")
      .def_prop_ro("unit", [](const YieldUncertainty& u) { return std::string(unitName(u.unit)); })
      .def_ro("time", &YieldUncertainty::time)
      .def_ro("fissions", &YieldUncertainty::fissions)
      .def_ro("response", &YieldUncertainty::response)
      .def_ro("sigma_diagonal", &YieldUncertainty::sigmaDiagonal,
              "Over every seeded product stating a sigma: what a diagonal treatment reports.")
      .def_ro("relative_diagonal", &YieldUncertainty::relativeDiagonal)
      .def_ro("sigma_diagonal_matched", &YieldUncertainty::sigmaDiagonalMatched,
              "The same sum restricted to products the matrix carries. THIS is what the "
              "correlated figure should be compared against: the two run over one set of "
              "products, so their whole difference is the off-diagonal.")
      .def_ro("relative_diagonal_matched", &YieldUncertainty::relativeDiagonalMatched)
      .def_prop_ro(
          "sigma_correlated",
          [](const YieldUncertainty& u) -> nb::object {
            // None when the contraction came out negative. The published matrices are not
            // positive semi-definite, so that is a real outcome and a square root would be an
            // invented number.
            return u.varianceNegative ? nb::none() : nb::cast(u.sigmaCorrelated);
          },
          "1-sigma with the imported correlation, or None when the contraction was negative.")
      .def_prop_ro("relative_correlated",
                   [](const YieldUncertainty& u) -> nb::object {
                     return u.varianceNegative ? nb::none() : nb::cast(u.relativeCorrelated);
                   })
      .def_ro("variance_correlated", &YieldUncertainty::varianceCorrelated,
              "The raw quadratic form, SIGNED. Negative is possible and is not an error.")
      .def_ro("variance_negative", &YieldUncertainty::varianceNegative)
      .def_ro("variance_ratio", &YieldUncertainty::varianceRatio,
              "Correlated variance over the matched diagonal one: what the correlation costs.")
      .def_ro("smallest_eigenvalue", &YieldUncertainty::smallestEigenvalue,
              "Of the correlation block actually contracted. Negative means the published "
              "matrix is indefinite on these products, which is expected and is not repaired.")
      .def_ro("covered_fraction", &YieldUncertainty::coveredFraction,
              "Share of the response sitting on products the matrix carries.")
      .def_ro("sigma_covered_fraction", &YieldUncertainty::sigmaCoveredFraction,
              "Share of the response sitting on products whose yield states an uncertainty.")
      .def_ro("products_matched", &YieldUncertainty::productsMatched)
      .def_ro("products_unmatched", &YieldUncertainty::productsUnmatched)
      .def_ro("products_with_sigma", &YieldUncertainty::productsWithSigma)
      .def_ro("products_without_sigma", &YieldUncertainty::productsWithoutSigma)
      .def_ro("products_unused_in_matrix", &YieldUncertainty::productsUnusedInMatrix)
      .def_ro("store_library", &YieldUncertainty::storeLibrary)
      .def_prop_ro("correlation_library",
                   [](const YieldUncertainty& u) { return u.provenance.library; })
      .def_prop_ro("correlation_system",
                   [](const YieldUncertainty& u) { return u.provenance.system; })
      .def_ro("products", &YieldUncertainty::products,
              "Ranked by diagonal variance fraction, largest first.")
      .def("__repr__", [](const YieldUncertainty& u) {
        return "<YieldUncertainty " + shortestRoundTrip(u.response) + " " +
               std::string(unitName(u.unit)) + ", diagonal " +
               shortestRoundTrip(u.relativeDiagonalMatched) + " correlated " +
               (u.varianceNegative ? std::string("indefinite")
                                   : shortestRoundTrip(u.relativeCorrelated)) +
               ">";
      });

  m.def(
      "yield_uncertainty",
      [](const NuclearData& data, const YieldCorrelation& correlation, const std::string& fissile,
         const nb::object& at, const std::string& energy, std::optional<double> fissions,
         std::optional<double> yield_kt, std::optional<double> energy_j, double mev_per_fission,
         const std::string& metric, const std::string& units,
         const exposure::PointSourceGeometry& geometry, const ResolvedPack* pack, int threads,
         bool prune, int cram_order) {
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

        ResponseSpec spec;
        spec.metric = pack != nullptr ? Metric::Pack : metricFrom(metric);
        spec.aggregate = Aggregate::Nuclide;
        spec.unit = requireUnit(units, spec.metric, Domain::Instant);
        spec.geometry = geometry;
        spec.pack = pack;
        DecayOptions options;
        options.threads = threads;
        options.prune = prune;
        options.order = cram_order == 16 ? CramOrder::Order16 : CramOrder::Order48;
        const double time = timeFrom(at);
        const nb::gil_scoped_release release;
        return yieldUncertainty(data, fissionSeed, correlation, time, spec, options);
      },
      "data"_a, "correlation"_a, "fissile"_a, "at"_a, "energy"_a = "thermal",
      "fissions"_a = nb::none(), "yield_kt"_a = nb::none(), "energy_j"_a = nb::none(),
      "mev_per_fission"_a = seed::kMeVPerFissionExplosiveYield, "metric"_a = "activity",
      "units"_a = "", "geometry"_a = exposure::PointSourceGeometry{}, "pack"_a = nb::none(),
      "threads"_a = 0, "prune"_a = true, "cram_order"_a = 48,
      "The error bar the EVALUATED fission yields put on the response at `at`, with and without "
      "the imported correlation. Costs one adjoint solve plus one dense contraction. The number "
      "of fissions is taken as exact; an uncertainty on the source term scales the response and "
      "this sigma together and is the caller's to carry.");

  // --- a job with a shape ----------------------------------------------------
  //
  // task_series treats a job as one window at one distance. A real job is legs at different
  // distances, and the dose is not distributed the way the durations are -- which is the
  // question a plan answers and a single interval cannot.
  nb::class_<PlanLeg>(m, "PlanLeg", "One stretch of a job.")
      .def(
          "__init__",
          [](PlanLeg* self, const std::string& name, const nb::object& duration, double distance_m,
             double occupancy, const exposure::PointSourceGeometry& geometry) {
            PlanLeg leg;
            leg.name = name;
            leg.durationSeconds = timeFrom(duration);
            leg.occupancy = occupancy;
            leg.geometry = geometry;
            // distance_m overrides whatever the geometry carried, because the distance is the
            // thing a plan varies and passing a whole geometry to change one number is friction.
            if (distance_m > 0.0) {
              leg.geometry.distanceM = distance_m;
            }
            new (self) PlanLeg(std::move(leg));
          },
          "name"_a, "duration"_a, "distance_m"_a = 0.0, "occupancy"_a = 1.0,
          "geometry"_a = exposure::PointSourceGeometry{})
      .def_ro("name", &PlanLeg::name)
      .def_ro("duration_s", &PlanLeg::durationSeconds)
      .def_ro("occupancy", &PlanLeg::occupancy,
              "Fraction of the leg actually spent in the field. Zero makes it a BREAK: the clock "
              "runs, nothing accrues, and no solve is spent.")
      .def("__repr__", [](const PlanLeg& leg) {
        return "<PlanLeg " + leg.name + " " + formatDuration(leg.durationSeconds) + ">";
      });

  nb::class_<LegResult>(m, "LegResult")
      .def_ro("name", &LegResult::name)
      .def_ro("start_s", &LegResult::startSeconds)
      .def_ro("end_s", &LegResult::endSeconds)
      .def_ro("occupancy", &LegResult::occupancy)
      .def_ro("distance_m", &LegResult::distanceM)
      .def_ro("accrued", &LegResult::accrued, "In the plan's unit, already scaled by occupancy.")
      .def_ro("fraction", &LegResult::fraction)
      .def_ro("cumulative", &LegResult::cumulative)
      .def_ro("cumulative_fraction", &LegResult::cumulativeFraction)
      .def_ro("mean_rate", &LegResult::meanRate,
              "Accrued per second actually spent there. What separates a leg that is expensive "
              "because it is long from one that is expensive because it is close.")
      .def_ro("is_break", &LegResult::isBreak)
      .def("__repr__", [](const LegResult& leg) {
        return "<LegResult " + leg.name + " " + shortestRoundTrip(leg.accrued) + ">";
      });

  nb::class_<TaskPlan>(m, "TaskPlan", "A job as a sequence of legs, each with its own geometry.")
      .def_prop_ro("unit", [](const TaskPlan& p) { return std::string(unitName(p.unit)); })
      .def_ro("start_s", &TaskPlan::startSeconds)
      .def_ro("end_s", &TaskPlan::endSeconds)
      .def_ro("total", &TaskPlan::total)
      .def_ro("elapsed_s", &TaskPlan::elapsedSeconds, "Wall clock, breaks included.")
      .def_ro("exposed_s", &TaskPlan::exposedSeconds,
              "Sum of duration times occupancy -- the time that earns the dose.")
      .def_ro("legs", &TaskPlan::legs)
      .def_ro("budget", &TaskPlan::budget)
      .def_ro("budget_spent", &TaskPlan::budgetSpent)
      .def_prop_ro(
          "spent_in_leg",
          [](const TaskPlan& p) -> nb::object {
            // None rather than -1: a parser reading -1 as an index would be reading the last leg.
            return p.budgetSpent ? nb::cast(p.spentInLeg) : nb::none();
          },
          "Index of the leg the budget runs out in, or None.")
      .def_prop_ro(
          "spent_at_s",
          [](const TaskPlan& p) -> nb::object {
            return p.budgetSpent ? nb::cast(p.spentAtSeconds) : nb::none();
          },
          "When it runs out, located by the same root-find stay_time uses.")
      .def("__len__", [](const TaskPlan& p) { return p.legs.size(); })
      .def("__repr__", [](const TaskPlan& p) {
        return "<TaskPlan " + std::to_string(p.legs.size()) + " legs, " +
               shortestRoundTrip(p.total) + " " + unitName(p.unit) + ">";
      });

  m.def(
      "task_plan",
      [](const NuclearData& data, const Inventory& inventory, const nb::object& at,
         const std::vector<PlanLeg>& legs, const std::string& metric, const std::string& units,
         double budget, int threads, bool prune, int cram_order) {
        ResponseSpec spec;
        spec.metric = metricFrom(metric);
        spec.aggregate = Aggregate::Nuclide;
        spec.unit = requireUnit(units, spec.metric, Domain::Interval);
        DecayOptions options;
        options.threads = threads;
        options.prune = prune;
        options.order = cram_order == 16 ? CramOrder::Order16 : CramOrder::Order48;
        const double start = timeFrom(at);
        const nb::gil_scoped_release release;
        return runTaskPlan(data, inventory, spec, start, legs, budget, options);
      },
      "data"_a, "inventory"_a, "at"_a, "legs"_a, "metric"_a = "activity", "units"_a = "",
      "budget"_a = 0.0, "threads"_a = 0, "prune"_a = true, "cram_order"_a = 48,
      "What a job costs leg by leg, starting at `at`. Legs follow one another with no gap, so a "
      "wait between them is a break leg rather than an implied silence. The units are interval "
      "units, because a leg accrues a total rather than holding a rate. Pass `budget` to have "
      "the plan say where it runs out. No leg is shielded.");

  // --- the smallest list that works everywhere -------------------------------
  //
  // `rank(coverage=0.95)` is this question at one instant for one metric. This is the same
  // question asked of every metric at every time simultaneously, which is a covering problem
  // rather than a longer sort -- and the union of the per-time answers is not it.
  nb::class_<SetMember>(m, "SetMember")
      .def_ro("label", &SetMember::label)
      .def_ro("order", &SetMember::order,
              "When greedy chose it, from 1. NOT a ranking: the second member is whichever most "
              "improved the constraints still unmet given the first.")
      .def_ro("peak_fraction", &SetMember::peakFraction,
              "The largest share of any requirement's total this member ever holds. A member at "
              "0.1% is there to hold up one particular instant, not because it is ever large.")
      .def_ro("closed_shortfall", &SetMember::closedShortfall)
      .def("__repr__", [](const SetMember& s) {
        return "<SetMember " + std::to_string(s.order) + " " + s.label + ">";
      });

  nb::class_<CoveragePoint>(m, "CoveragePoint")
      .def_ro("requirement", &CoveragePoint::requirement)
      .def_ro("time_s", &CoveragePoint::timeSeconds)
      .def_ro("achieved", &CoveragePoint::achieved)
      .def_ro("required", &CoveragePoint::required)
      .def("__repr__", [](const CoveragePoint& p) {
        return "<CoveragePoint " + p.requirement + " at " + formatDuration(p.timeSeconds) + ": " +
               shortestRoundTrip(p.achieved) + " of " + shortestRoundTrip(p.required) + ">";
      });

  nb::class_<TriageSet>(m, "TriageSet", "The smallest set holding a floor everywhere at once.")
      .def_ro("members", &TriageSet::members)
      .def_ro("binding", &TriageSet::binding,
              "Where the set is closest to failing. The number to read first: a list of twenty "
              "says nothing about whether it is comfortable or exactly on the edge.")
      .def_ro("shortfalls", &TriageSet::shortfalls,
              "Constraints no set these tables can form would meet. Empty in the ordinary case.")
      .def_ro("candidates", &TriageSet::candidateCount)
      .def_ro("constraints", &TriageSet::constraintCount)
      .def("__len__", [](const TriageSet& s) { return s.members.size(); })
      .def("__repr__", [](const TriageSet& s) {
        return "<TriageSet " + std::to_string(s.members.size()) + " of " +
               std::to_string(s.candidateCount) + " contributors>";
      });

  m.def(
      "shortlist",
      [](const nb::dict& requirements, double coverage) {
        std::vector<CoverageRequirement> specs;
        // The tables are borrowed, not copied: the dict passed in owns them and is alive for
        // the whole call, which is the only lifetime a raw pointer here needs.
        for (const auto [key, value] : requirements) {
          CoverageRequirement requirement;
          requirement.label = nb::cast<std::string>(key);
          requirement.fraction = coverage;
          if (nb::isinstance<nb::tuple>(value)) {
            const nb::tuple pair = nb::cast<nb::tuple>(value);
            if (pair.size() != 2) {
              throw InputError("shortlist: a requirement is a table, or a (table, fraction) pair");
            }
            requirement.table = &nb::cast<const ResponseTable&>(pair[0]);
            requirement.fraction = nb::cast<double>(pair[1]);
          } else {
            requirement.table = &nb::cast<const ResponseTable&>(value);
          }
          specs.push_back(requirement);
        }
        return robustTriageSet(specs);
      },
      "requirements"_a, "coverage"_a = 0.95,
      "The smallest set of contributors holding `coverage` of EVERY requirement at EVERY time. "
      "`requirements` maps a name to a response table, or to a (table, fraction) pair when the "
      "floors differ. Set cover is NP-hard and this is a deterministic greedy: a small set that "
      "meets the floor, not a proof that none smaller exists.");

  // --- the source term -------------------------------------------------------
  //
  // The photon metric collapsed onto energy instead of onto nuclides. `response(..., by="line")`
  // ranks the discrete lines, which is what a person reads; this is what another code reads,
  // and the difference is not resolution but audience -- see triage/spectrum.hpp.
  nb::class_<BinnedSpectrum>(m, "BinnedSpectrum",
                             "Photon emission binned onto an energy grid: a source term.")
      .def_prop_ro("unit", [](const BinnedSpectrum& s) { return std::string(unitName(s.unit)); })
      .def_prop_ro("domain", [](const BinnedSpectrum& s) { return domainName(s.domain); })
      .def_ro("t1", &BinnedSpectrum::timeSeconds, "The instant, or the window's start.")
      .def_ro("t2", &BinnedSpectrum::timeEndSeconds,
              "The window's end; equal to t1 for an instantaneous spectrum.")
      .def_prop_ro(
          "edges_ev",
          [](nb::handle self) {
            const BinnedSpectrum& s = nb::cast<const BinnedSpectrum&>(self);
            return view1d(s.edgesEv.data(), s.edgesEv.size(), self);
          },
          "Bin boundaries in eV, ascending. One longer than `emission`.")
      .def_prop_ro(
          "emission",
          [](nb::handle self) {
            const BinnedSpectrum& s = nb::cast<const BinnedSpectrum&>(self);
            return view1d(s.values.data(), s.values.size(), self);
          },
          "Emission in each bin, in `unit`. A zero-copy view, not a copy.")
      .def_ro("total", &BinnedSpectrum::total,
              "Every evaluated line, binned or not. Equals emission.sum() + below_range + "
              "above_range, which is the identity that makes a truncated grid visible.")
      .def_ro("below_range", &BinnedSpectrum::belowRange,
              "Emission softer than the first edge, and so absent from a deck built on it.")
      .def_ro("above_range", &BinnedSpectrum::aboveRange, "Emission harder than the last edge.")
      .def_ro("line_count", &BinnedSpectrum::lineCount)
      .def_ro("emitter_count", &BinnedSpectrum::emitterCount)
      .def_ro("unmodeled_energy_fraction", &BinnedSpectrum::unmodeledEnergyFraction,
              "Share of emitted photon ENERGY in continua NuSIFT does not model. A source "
              "understated by this much understates whatever is run against it.")
      .def_ro("unmodeled_continuum", &BinnedSpectrum::unmodeledContinuum,
              "The emitters carrying that continuum, by name.")
      .def("__repr__", [](const BinnedSpectrum& s) {
        return "<BinnedSpectrum " + std::to_string(s.binCount()) + " bins, " +
               shortestRoundTrip(s.total) + " " + unitName(s.unit) + ">";
      });

  m.def(
      "binned_spectrum",
      [](const NuclearData& data, const DecayResult& result, int time_index, int bins,
         const std::string& scale, double min_ev, double max_ev, const nb::object& edges_ev) {
        return binnedSpectrum(data, result, time_index,
                              binningFrom(bins, scale, min_ev, max_ev, edges_ev));
      },
      "data"_a, "result"_a, "time_index"_a = 0, "bins"_a = 100, "scale"_a = "linear",
      "min_ev"_a = 0.0, "max_ev"_a = 0.0, "edges_ev"_a = nb::none(),
      "The photon emission spectrum at one time of a solve, in photons/s. An INDEX into the "
      "solve rather than a time, because a spectrum is built per instant and a grid point is "
      "what a caller iterating one already holds. A range of zero is taken from the lines "
      "present, which is the only default that cannot silently drop photons.");

  m.def(
      "binned_spectrum",
      [](const NuclearData& data, const IntervalResult& result, int bins, const std::string& scale,
         double min_ev, double max_ev, const nb::object& edges_ev) {
        return binnedIntervalSpectrum(data, result.nuclideKeys, result.integratedAtoms, result.t1,
                                      result.t2,
                                      binningFrom(bins, scale, min_ev, max_ev, edges_ev));
      },
      "data"_a, "result"_a, "bins"_a = 100, "scale"_a = "linear", "min_ev"_a = 0.0,
      "max_ev"_a = 0.0, "edges_ev"_a = nb::none(),
      "The same over an interval: a COUNT of photons emitted across the window, integrated "
      "exactly, which is the source term for a job that runs while the inventory decays.");

  // Renders the deck rather than handing back the pieces, because the pieces are not the
  // deliverable: the comment cards saying what has to be multiplied by what, and what the
  // source is missing, are the part a transport code will not supply for itself. A notebook
  // that assembled its own SDEF from `edges_ev` and `emission` would lose exactly those.
  m.def(
      "source_deck",
      [](const NuclearData& data, const Inventory& inventory, const BinnedSpectrum& spectrum,
         const std::string& format) {
        SourceFormat chosen = SourceFormat::McnpSdef;
        if (!parseSourceFormat(format, chosen)) {
          throw InputError("source: \"" + format +
                           "\" is not a format (text, csv, json, mcnp, openmc)");
        }
        ReportContext context;
        context.storeLibrary = data.provenance().library;
        context.storeCreatedUtc = data.provenance().createdUtc;
        context.storeNuclideCount = data.stagedCount();
        context.seedProvenance = inventory.provenance();
        std::ostringstream out;
        writeSourceSpectrum(out, spectrum, context, chosen);
        return out.str();
      },
      "data"_a, "inventory"_a, "spectrum"_a, "format"_a = "mcnp",
      "The spectrum as a transport code's source definition: 'mcnp' for an SDEF card, 'openmc' "
      "for a Python snippet, or 'text', 'csv', 'json'. Returns the deck as a string, caveats "
      "included -- they are the half a transport code cannot infer.");

  // The other attribution of the same number. `rank` says what is producing the response now;
  // this says which seeded nuclide it came from, and the two totals agree because they are
  // partitions of one quantity rather than two calculations of it.
  m.def(
      "attribute",
      [](const NuclearData& data, const Inventory& inventory, const nb::object& at,
         const std::string& metric, const std::string& by, const std::string& units,
         const exposure::PointSourceGeometry& geometry, int top, double coverage,
         double min_fraction, const nb::object& pin, bool prune, int cram_order) {
        DecayOptions options;
        options.prune = prune;
        options.order = cram_order == 16 ? CramOrder::Order16 : CramOrder::Order48;
        ResponseSpec spec;
        spec.metric = metricFrom(metric);
        // Passed through, not pinned to Nuclide: attributeToSeed refuses the others with a
        // message saying why, which is a better answer than silently attributing by nuclide
        // to someone who asked for elements.
        spec.aggregate = aggregateFrom(by);
        spec.unit = requireUnit(units, spec.metric, Domain::Instant);
        spec.geometry = geometry;

        RankRequest request;
        request.topN = top;
        request.coverage = coverage;
        request.minFraction = min_fraction;
        // A bare string is one pin, not an iterable of one-character ones -- the same trap
        // ResponseTable.rank guards, and for the same reason: "Cs-137" is a perfectly good
        // sequence of six spellings that name nothing.
        if (!pin.is_none()) {
          if (nb::isinstance<nb::str>(pin)) {
            request.pinned.push_back(requireSeedPin(data, inventory, nb::cast<std::string>(pin)));
          } else {
            for (const nb::handle item : pin) {
              request.pinned.push_back(
                  requireSeedPin(data, inventory, nb::cast<std::string>(item)));
            }
          }
        }

        // Released only once every Python argument has been read, as `integrate` does. The
        // adjoint solve -- and, for exposure, the forward one beside it -- needs nothing from
        // the interpreter, and holding the GIL across it blocks every other thread.
        const double time = timeFrom(at);
        const nb::gil_scoped_release release;
        return attributeToSeed(data, inventory, time, spec, request, options);
      },
      "data"_a, "inventory"_a, "at"_a, "metric"_a = "activity", "by"_a = "nuclide", "units"_a = "",
      "geometry"_a = exposure::PointSourceGeometry{}, "top"_a = 10, "coverage"_a = 0.0,
      "min_fraction"_a = 0.0, "pin"_a = nb::none(), "prune"_a = true, "cram_order"_a = 48,
      "Attribute the response at `at` to the nuclides in `inventory`. `at` is seconds or a "
      "duration string; `pin` is a seed name, or several.");

  // --- how much is allowed ----------------------------------------------------

  // A criterion is a quantity and a limit on it. Whether the limit is legitimate, genuinely
  // linear, or composed correctly with others is a rule layer above this one -- naming it is
  // the caller's assertion that dividing by it means something.
  nb::class_<Criterion>(m, "Criterion", "A quantity, and how much of it is allowed.")
      .def(
          "__init__",
          [](Criterion* self, const std::string& name, double limit, const std::string& metric,
             const std::string& by, const std::string& units,
             const exposure::PointSourceGeometry& geometry) {
            new (self) Criterion();
            self->name = name;
            self->limit = limit;
            self->spec.metric = metricFrom(metric);
            self->spec.aggregate = aggregateFrom(by);
            // Instant, always: an interval unit is an accrued total over a window, and a window
            // is not what a possession or transport limit constrains.
            self->spec.unit = requireUnit(units, self->spec.metric, Domain::Instant);
            self->spec.geometry = geometry;
          },
          "name"_a, "limit"_a, "metric"_a = "activity", "by"_a = "nuclide", "units"_a = "",
          "geometry"_a = exposure::PointSourceGeometry{})
      .def_ro("name", &Criterion::name)
      .def_ro("limit", &Criterion::limit)
      .def_prop_ro("units", [](const Criterion& c) { return std::string(unitName(c.spec.unit)); })
      .def_prop_ro("metric",
                   [](const Criterion& c) { return std::string(metricName(c.spec.metric)); })
      .def("__repr__", [](const Criterion& c) {
        return "<Criterion " + c.name + " <= " + shortestRoundTrip(c.limit) + " " +
               unitName(c.spec.unit) + ">";
      });

  nb::class_<LimitingContributor>(m, "LimitingContributor")
      .def_ro("label", &LimitingContributor::label)
      .def_ro("fraction", &LimitingContributor::fraction)
      .def("__repr__",
           [](const LimitingContributor& c) { return "<LimitingContributor " + c.label + ">"; });

  nb::class_<CriterionHeadroom>(m, "CriterionHeadroom")
      .def_ro("name", &CriterionHeadroom::name)
      .def_ro("response", &CriterionHeadroom::response)
      .def_ro("limit", &CriterionHeadroom::limit)
      .def_ro("fraction", &CriterionHeadroom::fraction, "R/L, the sum-of-fractions idiom.")
      .def_ro("scale", &CriterionHeadroom::scale, "L/R, the factor this criterion alone allows.")
      .def_ro("binding", &CriterionHeadroom::binding)
      .def_ro("unbounded", &CriterionHeadroom::unbounded,
              "The response is zero, so this criterion constrains nothing. Not a very large "
              "allowance -- no allowance at issue.")
      .def("__repr__", [](const CriterionHeadroom& h) {
        return "<CriterionHeadroom " + h.name + (h.unbounded ? " unbounded>" : ">");
      });

  nb::class_<AllowableScale>(m, "AllowableScale")
      .def_ro("time_s", &AllowableScale::timeSeconds)
      .def_prop_ro(
          "scale",
          [](const AllowableScale& at) {
            // None rather than a number when nothing binds: a zero here would read as the exact
            // opposite of what it means, and Python has a value for "there isn't one".
            return at.bounded ? nb::cast(at.scale) : nb::none();
          },
          "The largest factor the inventory can be multiplied by, or None when nothing "
          "constrains it at this time.")
      .def_ro("bounded", &AllowableScale::bounded)
      .def_prop_ro(
          "binding",
          [](const AllowableScale& at) {
            return at.bounded
                       ? nb::cast(at.criteria[static_cast<std::size_t>(at.bindingIndex)].name)
                       : nb::none();
          })
      .def_ro("criteria", &AllowableScale::criteria)
      .def_ro("limiting", &AllowableScale::limiting,
              "The contributors driving the binding criterion, most first.")
      .def("__repr__", [](const AllowableScale& at) {
        return "<AllowableScale at " + formatDuration(at.timeSeconds) +
               (at.bounded ? " x" + shortestRoundTrip(at.scale) : " unbounded") + ">";
      });

  m.def(
      "allowable_scale",
      [](const NuclearData& data, const DecayResult& result, const std::vector<Criterion>& criteria,
         int limiting) { return allowableScale(data, result, criteria, limiting); },
      "data"_a, "result"_a, "criteria"_a, "limiting"_a = 3,
      "By what factor the inventory can be multiplied before the first limit binds, at every "
      "time on the grid. Exact: every response is linear in the inventory, so the scale at "
      "which a criterion binds is a division rather than a search. Criteria over DIFFERENT "
      "metrics are where this earns its keep, since only then can the binding one change.");

  // --- counterfactual interventions ---------------------------------------------

  nb::class_<Removal>(m, "Removal", "Something taken out, and how much of it.")
      .def(
          "__init__",
          [](Removal* self, const std::string& selector, double fraction) {
            new (self) Removal();
            self->selector = selector;
            self->fraction = fraction;
          },
          "selector"_a, "fraction"_a = 1.0,
          "`selector` is a nuclide (\"Cs-137\"), an element (\"Cs\"), Z=55, or A=137. The "
          "element form is the one a chemical separation can actually perform.")
      .def_ro("selector", &Removal::selector)
      .def_ro("fraction", &Removal::fraction)
      .def("__repr__", [](const Removal& r) {
        return "<Removal " + r.selector + " x" + shortestRoundTrip(r.fraction) + ">";
      });

  nb::class_<Intervention>(m, "Intervention", "A set of removals applied at one instant.")
      .def(
          "__init__",
          [](Intervention* self, const std::string& name, const std::vector<Removal>& removals) {
            new (self) Intervention();
            self->name = name;
            self->removals = removals;
          },
          "name"_a, "removals"_a)
      .def_ro("name", &Intervention::name)
      .def_ro("removals", &Intervention::removals)
      .def("__repr__", [](const Intervention& i) { return "<Intervention " + i.name + ">"; });

  nb::class_<RemovedContributor>(m, "RemovedContributor")
      .def_ro("label", &RemovedContributor::label)
      .def_ro("key", &RemovedContributor::key)
      .def_ro("atoms_removed", &RemovedContributor::atomsRemoved)
      .def_ro("value", &RemovedContributor::value,
              "The response these removed atoms were carrying.")
      .def_ro("fraction", &RemovedContributor::fraction)
      .def("__repr__",
           [](const RemovedContributor& c) { return "<RemovedContributor " + c.label + ">"; });

  nb::class_<InterventionEffect>(m, "InterventionEffect")
      .def_ro("name", &InterventionEffect::name)
      .def_ro("response", &InterventionEffect::response, "R(T) once this was applied.")
      .def_ro("removed", &InterventionEffect::removed, "What it bought: baseline - response.")
      .def_ro("removed_fraction", &InterventionEffect::removedFraction)
      .def_ro("contributors", &InterventionEffect::contributors,
              "The nuclides the benefit came from, most first.")
      .def("__repr__", [](const InterventionEffect& e) {
        return "<InterventionEffect " + e.name + " -" + shortestRoundTrip(e.removedFraction) + ">";
      });

  nb::class_<InterventionStudy>(m, "InterventionStudy")
      .def_ro("baseline", &InterventionStudy::baseline, "R(T) with nothing removed.")
      .def_ro("intervention_time_s", &InterventionStudy::interventionTimeSeconds)
      .def_ro("response_time_s", &InterventionStudy::responseTimeSeconds)
      .def_ro("effects", &InterventionStudy::effects)
      .def_prop_ro("unit", [](const InterventionStudy& s) { return std::string(unitName(s.unit)); })
      .def_prop_ro("metric",
                   [](const InterventionStudy& s) { return std::string(metricName(s.metric)); })
      .def("__len__", [](const InterventionStudy& s) { return s.effects.size(); })
      .def("__repr__", [](const InterventionStudy& s) {
        return "<InterventionStudy " + std::to_string(s.effects.size()) + " alternatives>";
      });

  m.def(
      "compare_interventions",
      [](const NuclearData& data, const Inventory& inventory, const nb::object& remove_at,
         const nb::object& at, const std::vector<Intervention>& interventions,
         const std::string& metric, const std::string& by, const std::string& units,
         const exposure::PointSourceGeometry& geometry, bool prune, int cram_order) {
        ResponseSpec spec;
        spec.metric = metricFrom(metric);
        spec.aggregate = aggregateFrom(by);
        spec.unit = requireUnit(units, spec.metric, Domain::Instant);
        spec.geometry = geometry;

        DecayOptions options;
        options.prune = prune;
        options.order = cram_order == 16 ? CramOrder::Order16 : CramOrder::Order48;

        // Every Python argument read before the GIL goes, as attribute() does: the forward and
        // adjoint solves need nothing from the interpreter, and holding it across them blocks
        // every other thread.
        const double t0 = timeFrom(remove_at);
        const double target = timeFrom(at);
        const nb::gil_scoped_release release;
        return compareInterventions(data, inventory, t0, target, spec, interventions, options);
      },
      "data"_a, "inventory"_a, "remove_at"_a, "at"_a, "interventions"_a, "metric"_a = "activity",
      "by"_a = "nuclide", "units"_a = "", "geometry"_a = exposure::PointSourceGeometry{},
      "prune"_a = true, "cram_order"_a = 48,
      "What each intervention is worth to R at `at`, from ONE adjoint solve over "
      "[remove_at, at]. The interventions are alternatives against a common baseline, not a "
      "sequence; a schedule of removals on different dates is a call per date. Exact: R is "
      "linear in the inventory, so each answer is a dot product rather than a re-solve.");
}
