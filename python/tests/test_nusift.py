"""Tests for the Python bindings.

These check the BINDING, not the physics -- the C++ suite already covers the physics against
analytic solutions and published constants, and repeating it here would only test that
nanobind can pass a double. What can go wrong at this layer is different: an array copied when
it should be a view, an exception swallowed, a string form the CLI accepts that the binding
does not, or a number quietly transposed on the way across.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

import nusift

STORE = Path(__file__).resolve().parents[2] / "data" / "nusift_b8.1.h5"
needs_store = pytest.mark.skipif(not STORE.is_file(), reason="no staged data store")


@pytest.fixture(scope="module")
def data():
    return nusift.NuclearData.open(str(STORE))


@pytest.fixture(scope="module")
def result(data):
    inv = nusift.seed_fission(data, "U-235", energy="thermal", fissions=1e20)
    return inv, nusift.decay(data, inv, nusift.logspace("1h", "10y", 24))


def test_version_is_exposed():
    assert nusift.__version__
    assert nusift.__version__.count(".") == 2


# --- time forms, shared with the CLI -----------------------------------------


def test_durations_parse_like_the_cli():
    assert nusift.parse_duration("30d") == 30 * 86400
    assert nusift.parse_duration("2h") == 7200
    assert nusift.parse_duration("3600") == 3600
    # A year is the Julian year, and the binding must not quietly use a different one.
    assert nusift.parse_duration("1y") / nusift.parse_duration("1d") == pytest.approx(365.25)


def test_grids_hit_their_endpoints_exactly():
    times = nusift.logspace("1h", "100y", 40)
    assert len(times) == 40
    assert times[0] == pytest.approx(3600.0)
    assert times[-1] == pytest.approx(100 * 365.25 * 86400)

    spec = nusift.parse_time_grid("1h:1y:log:10")
    assert len(spec) == 10


def test_bad_time_raises_a_python_exception():
    with pytest.raises(nusift.InputError):
        nusift.parse_duration("next tuesday")


# Python hands out float("inf") far more casually than C++ does -- it is what a division by
# zero or an overflowing product produces upstream in a notebook -- so the binding is where a
# non-finite time is most likely to arrive, and it must not pass through to the solver.
def test_non_finite_times_are_refused():
    with pytest.raises(nusift.InputError):
        nusift.logspace(1.0, math.inf, 2)
    with pytest.raises(nusift.InputError):
        nusift.logspace(math.inf, 100.0, 10)
    with pytest.raises(nusift.InputError):
        nusift.linspace(0.0, math.nan, 5)
    with pytest.raises(nusift.InputError):
        nusift.parse_duration("infs")


def test_inventory_refuses_counts_that_are_not_atom_counts():
    inv = nusift.Inventory()
    with pytest.raises(nusift.InputError):
        inv.add("Cs-137", -1.0)
    with pytest.raises(nusift.InputError):
        inv.add("Cs-137", math.inf)
    with pytest.raises(nusift.InputError):
        inv.add("Cs-137", math.nan)
    assert inv.total_atoms == 0.0


# --- zero copy ----------------------------------------------------------------


@needs_store
def test_decay_arrays_are_views_not_copies(result):
    _, res = result
    atoms = res.atoms
    assert atoms.dtype == np.float64
    assert atoms.shape == (len(res.times), len(res.nuclides))
    # A NumPy array with a base is a view over someone else's memory. Without this the binding
    # would copy a (times x nuclides) matrix on every attribute access, which for a fission
    # source is megabytes per touch.
    assert atoms.base is not None
    assert res.integrated_atoms.base is not None


@needs_store
def test_response_values_are_views_not_copies(data, result):
    _, res = result
    table = nusift.response(data, res, metric="activity", by="nuclide")
    assert table.values.base is not None
    assert table.values.shape == (len(table.times), len(table.labels))


@needs_store
def test_views_are_read_only(data, result):
    # Writing into a response table would corrupt the totals computed alongside it, so the
    # views are const on the C++ side and NumPy must see that.
    _, res = result
    table = nusift.response(data, res, metric="activity")
    with pytest.raises(ValueError):
        table.values[0, 0] = 1.0


# --- the numbers agree with the C++ ------------------------------------------


@needs_store
def test_gamma_constant_matches_the_published_value(data):
    # 12.91 against Ninkovic & Adrovic's 13.05 R*cm^2/(h*mCi). The residual is the air-table
    # evaluation and the roentgen convention, not scatter -- published constants are vacuum
    # quantities by definition. The classic 13.2 is the same physics in the pre-1979 roentgen.
    #
    # This stays here as a binding check: that the number survives the trip across nanobind
    # intact. The authoritative published-value pass is the sweep in test_validation.py, which
    # compares thirty of these against their reference table.
    published = data.gamma_constant("Co-60") * 1e4 * 3.7e7
    assert published == pytest.approx(13.05, rel=0.03)


@needs_store
def test_kiloton_conversion_is_glasstones(data):
    assert nusift.fissions_from_kt(1.0) == pytest.approx(1.45e23, rel=0.01)
    # The reactor convention counts delayed energy an explosive yield does not.
    assert nusift.fissions_from_kt(1.0, 200.0) == pytest.approx(1.31e23, rel=0.01)


@needs_store
def test_activity_equals_lambda_times_atoms(data, result):
    _, res = result
    table = nusift.response(data, res, metric="activity", by="nuclide", units="Bq")
    # Total activity is sum(lambda_i n_i), and every nuclide's lambda is available, so the
    # table can be checked against the atoms it was built from.
    atoms = np.asarray(res.atoms)
    lambdas = np.array(
        [math.log(2.0) / h if h > 0 else 0.0 for h in (data.half_life(n) for n in res.nuclides)]
    )
    expected = atoms @ lambdas
    assert np.allclose(np.asarray(table.totals), expected, rtol=1e-9)


@needs_store
def test_seeding_conserves_the_yield_sum(data):
    fissions = 1e20
    inv = nusift.seed_fission(data, "U-235", energy="thermal", fissions=fissions)
    # Independent yields sum to about 2.0, so the atom count is about twice the fissions.
    assert inv.total_atoms == pytest.approx(2.0 * fissions, rel=1e-6)
    assert "sum Y_indep = 2" in inv.provenance


# --- ranking and forecasting --------------------------------------------------


@needs_store
def test_ranking_is_ordered_and_reports_its_coverage(data, result):
    _, res = result
    table = nusift.response(data, res, metric="activity")
    ranking = table.rank(at="30d", top=5)

    assert len(ranking.contributors) == 5
    values = [c.value for c in ranking.contributors]
    assert values == sorted(values, reverse=True)
    assert ranking.contributors[0].rank == 1
    # A truncated ranking must say what it left out.
    assert 0.0 < ranking.covered_fraction <= 1.0
    assert ranking.omitted_count > 0
    assert ranking.total > 0


@needs_store
def test_rank_at_accepts_the_same_strings_as_the_cli(data, result):
    _, res = result
    table = nusift.response(data, res, metric="activity")
    # Snaps to the nearest grid point, so a time with no exact sample still works.
    assert table.rank(at="30d", top=1).time == pytest.approx(nusift.parse_duration("30d"), rel=0.2)
    assert table.rank(at=nusift.parse_duration("30d"), top=1).contributors


@needs_store
def test_a_pin_reaches_past_the_cut_without_moving_it(data, result):
    _, res = result
    table = nusift.response(data, res, metric="activity")
    plain = table.rank(at="10y", top=3)
    pinned = table.rank(at="10y", top=3, pin="Sm-151")

    # The ranking is untouched; the pin is a row after it, carrying where it really stands.
    assert [c.label for c in pinned.contributors[:3]] == [c.label for c in plain.contributors]
    assert len(pinned.contributors) == 4
    tail = pinned.contributors[3]
    assert tail.label == "Sm-151"
    assert tail.pinned
    assert not any(c.pinned for c in pinned.contributors[:3])
    assert tail.rank > 3
    assert pinned.covered_fraction > plain.covered_fraction


@needs_store
def test_several_pins_are_accepted_and_a_bare_string_is_one_pin(data, result):
    _, res = result
    table = nusift.response(data, res, metric="activity", by="mass-chain")

    # "A=147" is one pin, not five one-character ones -- a string is a perfectly good sequence
    # in Python, which is what makes that mistake worth a test.
    assert len(table.rank(at="10y", top=1, pin="A=147").contributors) == 2
    assert len(table.rank(at="10y", top=1, pin=["A=147", "A=85"]).contributors) == 3
    # A pin that already ranked is not a second row.
    assert len(table.rank(at="10y", top=1, pin="A=137").contributors) == 1

    with pytest.raises(Exception, match="pin"):
        table.rank(at="10y", top=1, pin="not-a-chain")


@needs_store
def test_dominance_windows_partition_the_grid(data, result):
    _, res = result
    table = nusift.response(data, res, metric="activity")
    windows = table.dominance_windows()

    assert windows
    assert windows[0].start_s == pytest.approx(np.asarray(table.times)[0])
    assert windows[-1].end_s == pytest.approx(np.asarray(table.times)[-1])
    for earlier, later in zip(windows, windows[1:]):
        assert earlier.end_s == pytest.approx(later.start_s)
        assert earlier.start_s < earlier.end_s


@needs_store
def test_exposure_and_activity_rank_differently(data, result):
    _, res = result
    by_activity = nusift.response(data, res, metric="activity").rank(at="1h", top=5)
    by_exposure = nusift.response(data, res, metric="exposure", units="R/h").rank(at="1h", top=5)
    # Not merely reordered -- pure beta emitters lead activity and contribute no exposure.
    assert {c.label for c in by_activity.contributors} != {
        c.label for c in by_exposure.contributors
    }


@needs_store
def test_geometry_scales_as_inverse_square(data, result):
    _, res = result
    near = nusift.response(
        data, res, metric="exposure", units="R/h",
        geometry=nusift.PointSource(distance_m=1.0, air_attenuation=False),
    )
    far = nusift.response(
        data, res, metric="exposure", units="R/h",
        geometry=nusift.PointSource(distance_m=2.0, air_attenuation=False),
    )
    assert np.asarray(far.totals)[0] == pytest.approx(np.asarray(near.totals)[0] / 4.0, rel=1e-9)


@needs_store
def test_per_line_columns_carry_their_energy(data, result):
    _, res = result
    lines = nusift.response(data, res, metric="exposure", by="line", units="R/h")
    top = lines.rank(at="1h", top=3).contributors
    assert top
    for contributor in top:
        assert contributor.line_energy_ev > 0
        # The label names the emitter and the energy in keV.
        assert "keV" in contributor.label


# --- errors -------------------------------------------------------------------


@needs_store
def test_errors_surface_as_python_exceptions(data, result):
    _, res = result

    with pytest.raises(nusift.InputError):
        nusift.response(data, res, metric="not-a-metric")
    with pytest.raises(nusift.InputError):
        nusift.response(data, res, by="not-an-aggregate")
    with pytest.raises(nusift.InputError):
        # Becquerel does not measure exposure; a category error, refused at the boundary.
        nusift.response(data, res, metric="exposure", units="Bq")
    with pytest.raises(nusift.InputError):
        # Cs-137 does not fission, and the message lists what the store does carry.
        nusift.seed_fission(data, "Cs-137", fissions=1e20)


@needs_store
def test_fission_source_size_must_be_given_exactly_once(data):
    with pytest.raises(nusift.InputError):
        nusift.seed_fission(data, "U-235")
    with pytest.raises(nusift.InputError):
        nusift.seed_fission(data, "U-235", fissions=1e20, yield_kt=1.0)


def test_missing_store_raises_rather_than_returning_none():
    with pytest.raises(nusift.NusiftError):
        nusift.NuclearData.open("definitely_not_a_store_12345.h5")


# --- store discovery ----------------------------------------------------------


@needs_store
def test_packaged_store_is_contributed_to_the_search(monkeypatch, tmp_path):
    """``NuclearData.open()`` with no argument must find the store the wheel ships.

    The packaged location is the one thing only Python knows, and the extension picks it up by
    calling back into ``nusift._data`` -- so this stands a real store in for the packaged one
    and checks it is found from a directory with no ``./data`` and no environment variable.
    Without that wiring the no-argument workflow in the README works only where it happens to
    be run from.
    """
    monkeypatch.delenv("NUSIFT_DATA_STORE", raising=False)
    monkeypatch.setattr(nusift._data, "default_store_path", lambda: STORE)
    monkeypatch.chdir(tmp_path)

    assert str(STORE) in nusift._data.store_search_paths()
    assert nusift.NuclearData.open().size > 0


def test_search_paths_are_the_ones_the_search_actually_uses(monkeypatch, tmp_path):
    # Reported by the C++ locator rather than described a second time in Python, so the
    # diagnostic cannot drift from the search it is meant to explain.
    monkeypatch.setenv("NUSIFT_DATA_STORE", str(tmp_path / "from_env.h5"))
    paths = nusift._data.store_search_paths()
    assert paths[0] == str(tmp_path / "from_env.h5")


def test_input_error_is_a_nusift_error():
    assert issubclass(nusift.InputError, nusift.NusiftError)
    assert issubclass(nusift.NusiftError, Exception)
