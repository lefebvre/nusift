"""Locating the nuclear-data store that ships inside the wheel.

The C++ store locator owns the search order -- explicit path, then ``$NUSIFT_DATA_STORE``,
then caller-supplied paths, then the install prefix, then ``./data``. This module contributes
the one location only Python knows about: the store packaged alongside the module itself.

Deliberately a contribution to that search rather than a second search of its own. Two front
ends with their own resolution order is how a CLI and a binding end up quietly reading
different evaluations and reporting different answers to the same question.
"""

from __future__ import annotations

from pathlib import Path


def default_store_path() -> Path | None:
    """The store shipped in this package, or None if the wheel carries none.

    A wheel built without a bundled evaluation is a legitimate configuration -- the store is a
    megabyte and a half and some deployments would rather point at a shared one -- so this
    answers with None instead of raising, and ``NuclearData.open()`` falls through to the rest
    of the search order.
    """
    packaged = Path(__file__).parent / "data"
    if not packaged.is_dir():
        return None
    stores = sorted(packaged.glob("*.h5"))
    return stores[0] if stores else None


def store_search_paths() -> list[str]:
    """Every place ``NuclearData.open()`` would look, in order. Diagnostic only.

    Delegates to the C++ locator, with this module's packaged store already contributed to
    it, so the list is what the search actually does rather than a second description of it
    that could drift from the first.
    """
    # Imported here rather than at module scope: the extension calls back into this module to
    # resolve the packaged store, and a top-level import would close that loop at load time.
    from ._core import store_search_paths as _paths

    return list(_paths())
