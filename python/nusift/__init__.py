"""NuSIFT: which isotopes dominate activity, exposure, or photon output, and when.

Given an isotopic inventory -- born from burnup, activation, or fission -- decay it forward
and rank what contributes, by nuclide, mass chain, element, or individual photon line.

    import nusift

    nd  = nusift.NuclearData.open()
    inv = nusift.seed_fission(nd, "U-235", energy="thermal", yield_kt=20)
    res = nusift.decay(nd, inv, nusift.logspace("1h", "100y", 60))

    tab = nusift.response(nd, res, metric="exposure", by="nuclide", units="Sv/h")
    for c in tab.rank(at="30d", top=5).contributors:
        print(c.label, c.value, c.fraction)

    for w in tab.dominance_windows():
        print(w.label, "leads", w.start_s, "to", w.end_s)

    # Totals over a window come from the same engine, solved in closed form:
    window = nusift.integrate(nd, inv, "1d", "30d")
    for c in nusift.response(nd, window, metric="activity").rank(top=5).contributors:
        print(c.label, c.value, "decays")

    # When a curve reaches a value. Every located event carries the grid bracket that found
    # it, because a crossing is a root of something SAMPLED and is never exact:
    for e in tab.crossings(1e-3):
        print(e.kind, "at", e.time_s, "known to", e.located_to_s)

    # ...narrowed by re-solving inside that bracket, rather than interpolated across it:
    ev = nusift.evaluator(nd, inv, metric="exposure", units="Sv/h")
    for e in tab.crossings(1e-3, refine=ev):
        print(e.time_s, e.refined, e.located_to_s)

    # When to do a fixed-length job, and when it fits a budget. Every sample is an exact
    # interval integral, so this curve costs solves the rate curve does not:
    task = nusift.task_series(nd, inv, nusift.logspace("1h", "1y", 40), "1h",
                              metric="exposure", units="R")
    for e in task.extrema():
        print("best start" if e.kind == "minimum" else "worst start", e.time_s)
    for w in task.windows_below(0.5):
        print("may start between", w.start_s, "and", w.end_s)

    # And how much is allowed, against limits you supply:
    limits = [nusift.Criterion("A2 transport", 3.7e13, metric="activity", units="Bq")]
    for at in nusift.allowable_scale(nd, res, limits):
        print(at.time_s, at.scale, at.binding)   # scale is None where nothing binds

    # And what taking something out would buy, from one adjoint solve over the window:
    plan = [nusift.Intervention("Cs separation", [nusift.Removal("Cs")])]
    study = nusift.compare_interventions(nd, inv, remove_at="30d", at="30y",
                                         interventions=plan, metric="exposure", units="Sv/h")
    for e in study.effects:
        print(e.name, "removes", e.removed_fraction, "of", study.baseline)

The large arrays -- ``DecayResult.atoms``, ``ResponseTable.values`` -- are zero-copy NumPy
views over the C++ storage rather than copies, so they are cheap to take and must not
outlive the object they came from. NumPy's own base-object tracking enforces that.
"""

from ._core import (  # noqa: F401
    AllowableScale,
    Contributor,
    Criterion,
    CriterionHeadroom,
    SeedAttribution,
    SeedShare,
    DecayResult,
    DominanceWindow,
    EventSeries,
    InputError,
    IntervalResult,
    Intervention,
    InterventionEffect,
    InterventionStudy,
    Inventory,
    LevelWindow,
    LimitingContributor,
    NuclearData,
    NusiftError,
    PointSource,
    Ranking,
    RemovedContributor,
    Removal,
    ResponseEvaluator,
    ResponseTable,
    TrajectoryEvent,
    __version__,
    allowable_scale,
    attribute,
    compare_interventions,
    decay,
    evaluator,
    fissions_from_kt,
    format_duration,
    integrate,
    linspace,
    logspace,
    parse_duration,
    parse_time_grid,
    read_inventory,
    response,
    seed_fission,
    task_series,
)
from ._data import default_store_path  # noqa: F401

__all__ = [
    "AllowableScale",
    "Contributor",
    "Criterion",
    "CriterionHeadroom",
    "SeedAttribution",
    "SeedShare",
    "DecayResult",
    "DominanceWindow",
    "EventSeries",
    "InputError",
    "IntervalResult",
    "Intervention",
    "InterventionEffect",
    "InterventionStudy",
    "Inventory",
    "LevelWindow",
    "LimitingContributor",
    "NuclearData",
    "NusiftError",
    "PointSource",
    "Ranking",
    "Removal",
    "RemovedContributor",
    "ResponseEvaluator",
    "ResponseTable",
    "TrajectoryEvent",
    "__version__",
    "allowable_scale",
    "attribute",
    "compare_interventions",
    "decay",
    "default_store_path",
    "evaluator",
    "fissions_from_kt",
    "format_duration",
    "integrate",
    "linspace",
    "logspace",
    "parse_duration",
    "parse_time_grid",
    "read_inventory",
    "response",
    "seed_fission",
    "task_series",
]
