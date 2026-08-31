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

    # And how much is allowed, against limits you supply:
    limits = [nusift.Criterion("A2 transport", 3.7e13, metric="activity", units="Bq")]
    for at in nusift.allowable_scale(nd, res, limits):
        print(at.time_s, at.scale, at.binding)   # scale is None where nothing binds

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
    InputError,
    IntervalResult,
    Inventory,
    LevelWindow,
    LimitingContributor,
    NuclearData,
    NusiftError,
    PointSource,
    Ranking,
    ResponseTable,
    TrajectoryEvent,
    __version__,
    allowable_scale,
    attribute,
    decay,
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
    "InputError",
    "IntervalResult",
    "Inventory",
    "LevelWindow",
    "LimitingContributor",
    "NuclearData",
    "NusiftError",
    "PointSource",
    "Ranking",
    "ResponseTable",
    "TrajectoryEvent",
    "__version__",
    "allowable_scale",
    "attribute",
    "decay",
    "default_store_path",
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
]
