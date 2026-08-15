"""NuSIFT: which isotopes dominate activity or exposure, and when.

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

The large arrays -- ``DecayResult.atoms``, ``ResponseTable.values`` -- are zero-copy NumPy
views over the C++ storage rather than copies, so they are cheap to take and must not
outlive the object they came from. NumPy's own base-object tracking enforces that.
"""

from ._core import (  # noqa: F401
    Contributor,
    DecayResult,
    DominanceWindow,
    InputError,
    Inventory,
    NuclearData,
    NusiftError,
    PointSource,
    Ranking,
    ResponseTable,
    __version__,
    decay,
    fissions_from_kt,
    format_duration,
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
    "Contributor",
    "DecayResult",
    "DominanceWindow",
    "InputError",
    "Inventory",
    "NuclearData",
    "NusiftError",
    "PointSource",
    "Ranking",
    "ResponseTable",
    "__version__",
    "decay",
    "default_store_path",
    "fissions_from_kt",
    "format_duration",
    "linspace",
    "logspace",
    "parse_duration",
    "parse_time_grid",
    "read_inventory",
    "response",
    "seed_fission",
]
