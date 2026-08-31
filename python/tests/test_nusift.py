"""Tests for the Python bindings.

These check the BINDING, not the physics -- the C++ suite already covers the physics against
analytic solutions and published constants, and repeating it here would only test that
nanobind can pass a double. What can go wrong at this layer is different: an array copied when
it should be a view, an exception swallowed, a string form the CLI accepts that the binding
does not, or a number quietly transposed on the way across.
"""

from __future__ import annotations

import math
import shutil
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


# --- intervals ----------------------------------------------------------------


@needs_store
def test_a_window_from_zero_is_the_cumulative_row(data):
    inv = nusift.Inventory()
    inv.add("Cs-137", 1e20)
    thirty_years = nusift.parse_duration("30y")
    window = nusift.integrate(data, inv, 0.0, "30y")
    cumulative = nusift.decay(data, inv, [thirty_years])

    assert window.t1 == 0.0
    assert window.t2 == thirty_years
    assert list(window.nuclides) == list(cumulative.nuclides)
    # From zero the window IS the cumulative integral, from the same solve: identical bits,
    # not merely close -- and a view, like every other array the binding hands out.
    np.testing.assert_array_equal(
        np.asarray(window.integrated_atoms), np.asarray(cumulative.integrated_atoms)[0]
    )
    assert window.integrated_atoms.base is not None


@needs_store
def test_integrate_reaches_the_guarded_path_for_a_narrow_late_window(data):
    # One second at thirty years. Differencing the two cumulative rows loses about ten of the
    # sixteen digits to cancellation, which is the reason the binding exists; it has to reach
    # intervalIntegral()'s re-solve rather than reproduce that subtraction.
    inv = nusift.Inventory()
    inv.add("Cs-137", 1e20)
    t1 = nusift.parse_duration("30y")
    window = nusift.integrate(data, inv, t1, t1 + 1.0)
    i = list(window.nuclides).index("Cs-137")

    lam = math.log(2.0) / data.half_life("Cs-137")
    expected = 1e20 / lam * math.exp(-lam * t1) * -math.expm1(-lam * 1.0)
    got = float(np.asarray(window.integrated_atoms)[i])
    assert got == pytest.approx(expected, rel=1e-9)

    rows = np.asarray(nusift.decay(data, inv, [t1, t1 + 1.0]).integrated_atoms)
    naive = float(rows[1, i] - rows[0, i])
    assert abs(got - expected) < abs(naive - expected)


@needs_store
def test_an_interval_response_is_in_the_interval_domain(data):
    inv = nusift.Inventory()
    inv.add("Cs-137", 1e20)
    window = nusift.integrate(data, inv, "1d", "30d")
    table = nusift.response(data, window, metric="activity")

    assert table.domain == "interval"
    assert table.unit == "decays"
    assert np.asarray(table.times)[0] == pytest.approx(nusift.parse_duration("1d"))
    assert np.asarray(table.time_ends)[0] == pytest.approx(nusift.parse_duration("30d"))
    # Ba-137m in equilibrium decays very nearly once per Cs-137 decay, so both lead.
    leaders = {c.label for c in table.rank(top=2).contributors}
    assert leaders == {"Cs-137", "Ba-137m"}

    with pytest.raises(nusift.InputError):
        # A rate cannot express a total; refused here as the CLI refuses it on `integrate`.
        nusift.response(data, window, metric="activity", units="Bq")


# --- the air path -------------------------------------------------------------


@needs_store
def test_exposure_tables_carry_their_air_path(data, result):
    _, res = result
    far = nusift.response(
        data, res, metric="exposure", units="R/h", geometry=nusift.PointSource(distance_m=100.0)
    )
    depth = np.asarray(far.mean_optical_depth)
    assert depth.shape == (len(far.times),)
    # A hundred metres of air is most of a mean free path at fission-product energies.
    assert (depth > 0.1).all()
    assert far.rank(at="1h", top=1).mean_optical_depth == pytest.approx(depth[0])

    vacuum = nusift.response(
        data, res, metric="exposure", units="R/h",
        geometry=nusift.PointSource(distance_m=100.0, air_attenuation=False),
    )
    assert not np.asarray(vacuum.mean_optical_depth).any()
    # Activity has no air path, and the table does not carry one that reads as a measurement.
    assert len(nusift.response(data, res, metric="activity").mean_optical_depth) == 0


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


@needs_store
def test_a_directory_holding_several_stores_is_warned_about(monkeypatch, tmp_path):
    """The first store by name is opened, and the ones passed over are named.

    Through the warnings module rather than the process's stderr, because a notebook never
    shows the latter -- and a store silently picked out of several is exactly what a user of
    one needs to hear about.
    """
    monkeypatch.delenv("NUSIFT_DATA_STORE", raising=False)
    monkeypatch.setattr(nusift._data, "default_store_path", lambda: None)
    (tmp_path / "data").mkdir()
    shutil.copy(STORE, tmp_path / "data" / "a_fixture.h5")
    shutil.copy(STORE, tmp_path / "data" / "b_evaluation.h5")
    monkeypatch.chdir(tmp_path)

    with pytest.warns(UserWarning, match="2 stores") as record:
        assert nusift.NuclearData.open().size > 0
    assert "b_evaluation.h5" in str(record[0].message)


def test_search_paths_are_the_ones_the_search_actually_uses(monkeypatch, tmp_path):
    # Reported by the C++ locator rather than described a second time in Python, so the
    # diagnostic cannot drift from the search it is meant to explain.
    monkeypatch.setenv("NUSIFT_DATA_STORE", str(tmp_path / "from_env.h5"))
    paths = nusift._data.store_search_paths()
    assert paths[0] == str(tmp_path / "from_env.h5")


def test_input_error_is_a_nusift_error():
    assert issubclass(nusift.InputError, nusift.NusiftError)
    assert issubclass(nusift.NusiftError, Exception)


# --- seed attribution ------------------------------------------------------
#
# The binding, not the physics: that the shares cross intact, that the pin path resolves
# against the SEED rather than a response table, and that the total is the same number `rank`
# reports rather than a second calculation of it.


@needs_store
def test_attribute_shares_partition_the_total(data):
    inv = nusift.seed_fission(data, "U-235", energy="thermal", yield_kt=20)
    at = nusift.parse_duration("30d")

    full = nusift.attribute(data, inv, at=at, top=0)
    assert full.omitted_count == 0
    assert math.isclose(sum(s.value for s in full.shares), full.total, rel_tol=1e-10)
    assert math.isclose(full.covered_fraction, 1.0, rel_tol=1e-10)


@needs_store
def test_attribute_total_matches_the_forward_ranking(data):
    inv = nusift.seed_fission(data, "U-235", energy="thermal", yield_kt=20)
    at = nusift.parse_duration("30d")

    attributed = nusift.attribute(data, inv, at=at, top=5)
    forward = nusift.response(data, nusift.decay(data, inv, [at])).rank(at=at, top=5)

    assert math.isclose(attributed.total, forward.total, rel_tol=1e-9)
    # Same number, different partition -- the two lists need not agree at all.
    assert attributed.labels != forward.labels


@needs_store
def test_attribute_truncation_reports_what_it_omitted(data):
    inv = nusift.seed_fission(data, "U-235", energy="thermal", yield_kt=20)
    at = nusift.parse_duration("30d")

    top5 = nusift.attribute(data, inv, at=at, top=5)
    assert len(top5) == 5
    assert top5.omitted_count > 0
    assert 0.0 < top5.covered_fraction < 1.0


@needs_store
def test_attribute_pins_reach_past_the_cut(data):
    inv = nusift.seed_fission(data, "U-235", energy="thermal", yield_kt=20)
    at = nusift.parse_duration("30d")

    pinned = nusift.attribute(data, inv, at=at, top=3, pin=["Cs-137"])
    assert len(pinned) == 4
    assert pinned.shares[-1].label == "Cs-137"
    assert pinned.shares[-1].pinned
    assert pinned.shares[-1].rank > 3


@needs_store
def test_attribute_refuses_a_pin_that_was_never_seeded(data):
    inv = nusift.Inventory()
    inv.add("Cs-137", 1.0e18)
    with pytest.raises(nusift.InputError):
        nusift.attribute(data, inv, at=0.0, pin=["Co-60"])


@needs_store
def test_attribute_importance_is_shared_within_a_decay_chain(data):
    """Xe-140 (13.6 s) decays into Cs-140 (63.7 s) long before 30 days, so at that time an atom
    seeded as either is worth very nearly the same. Not exactly: Xe-140 has a small
    delayed-neutron branch that leaves the A=140 chain, and the ~1e-5 gap is that branch rather
    than solver noise -- which is why the tolerance is loose enough to admit it and tight enough
    that a genuinely transposed importance would still fail."""
    inv = nusift.seed_fission(data, "U-235", energy="thermal", yield_kt=20)
    at = nusift.parse_duration("30d")

    shares = {s.label: s.importance for s in nusift.attribute(data, inv, at=at, top=0).shares}
    assert math.isclose(shares["Xe-140"], shares["Cs-140"], rel_tol=1e-3)


# The input forms the other time and pin APIs already take. `attribute` took a bare float and a
# list of strings, so `at="30d"` and `pin="Cs-137"` -- both of which work everywhere else --
# raised TypeError instead of answering.


@needs_store
def test_attribute_takes_a_duration_string_for_at(data):
    inv = nusift.seed_fission(data, "U-235", energy="thermal", yield_kt=20)

    by_string = nusift.attribute(data, inv, at="30d", top=5)
    by_seconds = nusift.attribute(data, inv, at=nusift.parse_duration("30d"), top=5)

    assert by_string.time == by_seconds.time
    assert by_string.labels == by_seconds.labels


@needs_store
def test_attribute_takes_a_bare_string_as_one_pin(data):
    """Not as an iterable of one-character pins. "Cs-137" is a perfectly good sequence of six
    spellings that name nothing, which is why this has to be handled rather than iterated."""
    inv = nusift.seed_fission(data, "U-235", energy="thermal", yield_kt=20)
    at = nusift.parse_duration("30d")

    one = nusift.attribute(data, inv, at=at, top=3, pin="Cs-137")
    listed = nusift.attribute(data, inv, at=at, top=3, pin=["Cs-137"])

    assert one.labels == listed.labels
    assert one.shares[-1].label == "Cs-137"
    assert one.shares[-1].pinned


@needs_store
def test_attribute_refuses_an_aggregate_it_cannot_answer(data):
    """An inventory row names a nuclide, so that is what a share can name. Answering `by=
    "element"` with a nuclide attribution would return a different table than the one asked
    for, silently."""
    inv = nusift.seed_fission(data, "U-235", energy="thermal", yield_kt=20)
    for by in ("element", "mass-chain", "line"):
        with pytest.raises(nusift.InputError):
            nusift.attribute(data, inv, at="30d", by=by)


@needs_store
def test_attribute_exposure_carries_the_same_caveats_the_ranking_does(data):
    """The same figure reached two ways cannot be better characterised one way than the other.
    `attribute` used to print the exposure `rank` warns about, with no warning."""
    inv = nusift.seed_fission(data, "U-235", energy="thermal", yield_kt=20)
    at = nusift.parse_duration("30d")
    geometry = nusift.PointSource(distance_m=100.0)

    attributed = nusift.attribute(
        data, inv, at=at, metric="exposure", units="R/h", geometry=geometry, top=5
    )
    forward = nusift.response(
        data, nusift.decay(data, inv, [at]), metric="exposure", units="R/h", geometry=geometry
    ).rank(at=at, top=5)

    assert math.isclose(attributed.total, forward.total, rel_tol=1e-9)
    assert math.isclose(
        attributed.unmodeled_energy_fraction, forward.unmodeled_energy_fraction, rel_tol=1e-12
    )
    assert math.isclose(attributed.mean_optical_depth, forward.mean_optical_depth, rel_tol=1e-12)
    assert attributed.buildup == forward.buildup
    # A caveat of zero would satisfy the equalities above and prove nothing.
    assert attributed.mean_optical_depth > 0.0
    assert attributed.unmodeled_continuum


@needs_store
def test_attribute_activity_carries_no_exposure_caveats(data):
    inv = nusift.seed_fission(data, "U-235", energy="thermal", yield_kt=20)
    a = nusift.attribute(data, inv, at="30d", top=5)
    assert a.unmodeled_energy_fraction == 0.0
    assert a.mean_optical_depth == 0.0
    assert a.unmodeled_continuum == []


@needs_store
def test_attribute_does_not_count_an_inert_seed_as_omitted(data):
    """A stable seed places real atoms worth zero becquerel forever. Ranked, it spends a `top`
    slot; counted as omitted, it reports a gap that showing it would not close."""
    inv = nusift.Inventory()
    inv.add("Cs-137", 1.0e18)
    inv.add("Cs-133", 5.0e18)  # stable

    a = nusift.attribute(data, inv, at="30d", top=1)
    assert a.labels == ["Cs-137"]
    assert a.omitted_count == 0
    assert math.isclose(a.covered_fraction, 1.0, rel_tol=1e-10)

    # Asked about directly, it answers -- rankless, which is what "contributes nothing here"
    # looks like in a table of ranks.
    pinned = nusift.attribute(data, inv, at="30d", top=1, pin="Cs-133")
    assert pinned.labels == ["Cs-137", "Cs-133"]
    assert pinned.shares[-1].rank == 0
    assert pinned.shares[-1].value == 0.0
    assert pinned.shares[-1].seed_atoms > 0.0
    assert pinned.omitted_count == 0


# --- located events ----------------------------------------------------------


@needs_store
def test_crossings_carry_the_bracket_that_found_them(data, result):
    """A crossing is a root of a SAMPLED curve, so `time_s` alone is not the answer. The binding
    has to hand over the grid interval and the width the instant was placed to, or a consumer
    cannot tell a located event from a grid artefact."""
    _, res = result
    table = nusift.response(data, res, metric="activity")

    # A level the total certainly passes on the way down.
    level = float(table.totals[0]) / 10.0
    events = table.crossings(level)
    assert events, "the total falls by more than a decade over this grid"

    first = events[0]
    assert first.kind == "falling"
    assert first.bracket_start_s <= first.time_s <= first.bracket_end_s
    assert first.located_to_s > 0.0
    # A table carries samples and no evaluator, so every event over it is interpolated.
    assert first.refined is False
    assert first.converged is True
    assert "falling" in repr(first)


@needs_store
def test_crossings_follow_a_contributor_or_a_ratio(data, result):
    _, res = result
    table = nusift.response(data, res, metric="activity")
    labels = table.labels

    # `of` narrows to one column; the total and one contributor are different curves and may
    # cross a level at different times, or a different number of times.
    name = labels[0]
    of_events = table.crossings(1.0, of=name)
    assert isinstance(of_events, list)

    ratio_events = table.crossings(1.0, ratio=(labels[0], labels[1]))
    assert isinstance(ratio_events, list)


@needs_store
def test_a_ratio_needs_exactly_two_contributors(data, result):
    _, res = result
    table = nusift.response(data, res, metric="activity")
    with pytest.raises(nusift.InputError):
        table.crossings(1.0, ratio=(table.labels[0],))


@needs_store
def test_windows_flag_the_edges_the_grid_never_observed(data, result):
    """An edge outside the grid is a bound, not a crossing. Clipping it to the grid's own
    endpoint would report an instant that was never seen."""
    _, res = result
    table = nusift.response(data, res, metric="activity")

    # Far below anything on the curve, so the whole grid is inside and neither edge is real.
    windows = table.windows_above(1.0)
    assert len(windows) == 1
    assert windows[0].entry_observed is False
    assert windows[0].exit_observed is False
    assert "open" in repr(windows[0])

    # Above and below partition the grid, so a level nothing reaches gives one window below.
    below = table.windows_below(1.0)
    assert below == [] or below[0].start_s >= 0.0


@needs_store
def test_an_evaluator_narrows_events_a_table_alone_interpolates(data, result):
    """The table carries samples; the evaluator carries a way to ask what happens between two of
    them. Same grid, same crossing, placed by solving instead of by interpolating -- and
    `refined` says which happened rather than leaving it to be inferred."""
    inv, res = result
    table = nusift.response(data, res, metric="activity")
    level = float(table.totals[0]) / 10.0

    ev = nusift.evaluator(data, inv, metric="activity")
    wide = table.crossings(level)[0]
    narrow = table.crossings(level, refine=ev)[0]

    assert wide.refined is False
    assert narrow.refined is True
    assert narrow.located_to_s < wide.located_to_s / 100.0
    assert wide.bracket_start_s <= narrow.time_s <= wide.bracket_end_s
    assert ev.solves > 0

    # A boundary is a located event too, and takes refinement on the same terms.
    refined_windows = table.dominance_windows(refine=ev)
    assert len(refined_windows) == len(table.dominance_windows())


@needs_store
def test_an_evaluator_for_another_curve_is_refused(data, result):
    inv, res = result
    table = nusift.response(data, res, metric="activity")
    with pytest.raises(nusift.InputError):
        table.crossings(1.0, refine=nusift.evaluator(data, inv, metric="exposure", units="R/h"))


@needs_store
def test_a_task_curve_answers_when_to_start_and_when_it_fits(data, result):
    """What a fixed-length job accrues, against when it starts. Not a slice of the rate curve:
    every sample is an exact interval integral over its own window."""
    inv, _ = result
    starts = nusift.logspace("1h", "10y", 20)
    task = nusift.task_series(data, inv, starts, "1h", metric="activity", units="decays")

    assert len(task.values) == len(starts)
    assert task.refines is False
    assert all(v > 0.0 for v in task.values)

    # Waiting makes this job cheaper, so a budget between the first and last sample is met from
    # some start onwards and never stops being met.
    budget = (float(task.values[0]) + float(task.values[-1])) / 2.0
    allowed = task.windows_below(budget)
    assert allowed
    assert allowed[-1].exit_observed is False

    events = task.crossings(budget)
    assert events and events[0].kind == "falling"
    assert events[0].time_s == pytest.approx(allowed[-1].start_s)


@needs_store
def test_extrema_find_an_ingrowth_peak(data):
    """Y-90 grows into equilibrium with Sr-90 and then follows its parent down, so the curve
    turns. A monotone one does not, and reporting a turn on it would be an artefact."""
    inv = nusift.Inventory()
    inv.add("Sr-90", 1.0e20)
    res = nusift.decay(data, inv, nusift.logspace("1h", "100y", 80))
    table = nusift.response(data, res, metric="activity")

    turns = table.extrema(of="Y-90")
    assert turns, "Y-90 ingrowth turns over once it reaches equilibrium"
    assert turns[0].kind == "maximum"
    assert turns[0].bracket_start_s <= turns[0].time_s <= turns[0].bracket_end_s

    # The parent only decays, so it has no interior turn at all.
    assert table.extrema(of="Sr-90") == []


# --- how much is allowed ------------------------------------------------------


@needs_store
def test_allowable_scale_is_exact_against_the_limit(data, result):
    """The scale is a division, not a search: R is linear in the inventory, so R times the
    scale is the limit. Checked here at the binding layer because a transposed or rescaled
    number would still look plausible."""
    _, res = result
    limit = 3.7e13
    criteria = [nusift.Criterion("A2 transport", limit, metric="activity", units="Bq")]
    scaled = nusift.allowable_scale(data, res, criteria)

    assert len(scaled) == len(res.times)
    for at in scaled:
        if not at.bounded:
            continue
        assert at.binding == "A2 transport"
        headroom = at.criteria[0]
        assert math.isclose(headroom.response * at.scale, limit, rel_tol=1e-9)
        # The two idioms are the same number inverted.
        assert math.isclose(headroom.fraction, 1.0 / at.scale, rel_tol=1e-9)


@needs_store
def test_an_unconstrained_time_has_no_scale_rather_than_a_zero(data):
    """`scale` is None where nothing binds. A zero would read as the exact opposite of what it
    means, which is the one mistake this value must not invite."""
    inv = nusift.Inventory()
    inv.add("Cs-133", 1.0e20)  # stable: no activity at any time
    res = nusift.decay(data, inv, nusift.logspace("1h", "10y", 8))

    criteria = [nusift.Criterion("possession", 1.0e10, metric="activity", units="Bq")]
    scaled = nusift.allowable_scale(data, res, criteria)

    assert scaled
    for at in scaled:
        assert at.bounded is False
        assert at.scale is None
        assert at.binding is None
        assert at.criteria[0].unbounded is True
        assert at.limiting == []


@needs_store
def test_allowable_scale_names_what_drives_the_binding_criterion(data, result):
    _, res = result
    criteria = [nusift.Criterion("possession", 1.0e15, metric="activity", units="Bq")]
    scaled = nusift.allowable_scale(data, res, criteria, limiting=2)

    bounded = [at for at in scaled if at.bounded]
    assert bounded
    assert len(bounded[0].limiting) <= 2
    assert bounded[0].limiting[0].label
    assert bounded[0].limiting[0].fraction > 0.0
    # Ordered by share of the binding criterion's total.
    fractions = [c.fraction for c in bounded[0].limiting]
    assert fractions == sorted(fractions, reverse=True)


@needs_store
def test_a_criterion_is_refused_when_it_cannot_mean_anything(data, result):
    _, res = result
    # An interval unit is an accrued total over a window, which is not what a possession limit
    # constrains. Refused where the C++ refuses it -- when the criteria are used, not when the
    # spelling is parsed -- so the binding and the library agree about which call fails.
    accrued = nusift.Criterion("decays", 1.0e15, metric="activity", units="decays")
    with pytest.raises(nusift.InputError, match="total accrued over a window"):
        nusift.allowable_scale(data, res, [accrued])

    with pytest.raises(nusift.InputError):
        nusift.allowable_scale(data, res, [])

    with pytest.raises(nusift.InputError):
        nusift.allowable_scale(
            data, res, [nusift.Criterion("", 1.0e15, metric="activity", units="Bq")]
        )


# --- counterfactual interventions ---------------------------------------------


@needs_store
def test_an_intervention_is_exact_against_a_forward_solve(data):
    """The benefit is a dot product against the adjoint, so it has to equal what a full forward
    solve of the counterfactual gives. Checked at the binding layer because a slip in the time
    arguments would still produce a plausible-looking number."""
    inv = nusift.Inventory()
    inv.add("Cs-137", 1.0e20)
    inv.add("Sr-90", 5.0e19)

    plan = [nusift.Intervention("strip Cs", [nusift.Removal("Cs")])]
    study = nusift.compare_interventions(
        data, inv, remove_at="30d", at="30y", interventions=plan
    )

    # Rebuild the counterfactual by hand: decay to t0, drop every caesium atom, decay on.
    res = nusift.decay(data, inv, [nusift.parse_duration("30d")])
    after = nusift.Inventory()
    for label, atoms in zip(res.nuclides, res.atoms[0]):
        if not label.startswith("Cs-"):
            after.add(label, max(0.0, float(atoms)))
    window = nusift.parse_duration("30y") - nusift.parse_duration("30d")
    independently = nusift.response(data, nusift.decay(data, after, [window])).rank(top=0).total

    assert study.effects[0].response == pytest.approx(independently, rel=1e-9)


@needs_store
def test_the_benefit_is_booked_against_what_was_removed(data):
    """Removing Cs-137 takes with it the Ba-137m it would have fed -- but the benefit belongs to
    the caesium that was removed, not to the barium that would have emitted."""
    inv = nusift.Inventory()
    inv.add("Cs-137", 1.0e20)

    plan = [nusift.Intervention("strip Cs", [nusift.Removal("Cs")])]
    study = nusift.compare_interventions(
        data, inv, remove_at="0", at="30y", interventions=plan, metric="exposure", units="Sv/h"
    )

    effect = study.effects[0]
    assert effect.removed_fraction == pytest.approx(1.0, rel=1e-6)
    labels = [c.label for c in effect.contributors]
    assert "Cs-137" in labels
    assert "Ba-137m" not in labels


@needs_store
def test_what_a_separation_buys_depends_on_the_metric(data):
    """Sr-90 and its Y-90 daughter are beta emitters. They carry real activity and almost no
    photon exposure, so the same separation is worth wholly different amounts depending on which
    question is being asked -- which is the reason both metrics exist."""
    inv = nusift.Inventory()
    inv.add("Cs-137", 1.0e20)
    inv.add("Sr-90", 1.0e20)

    plan = [nusift.Intervention("strip Sr", [nusift.Removal("Sr")])]
    by_activity = nusift.compare_interventions(
        data, inv, remove_at="30d", at="30y", interventions=plan
    )
    by_exposure = nusift.compare_interventions(
        data, inv, remove_at="30d", at="30y", interventions=plan,
        metric="exposure", units="Sv/h",
    )

    assert by_activity.effects[0].removed_fraction > 0.2
    assert by_exposure.effects[0].removed_fraction < 0.01


@needs_store
def test_alternatives_share_one_baseline(data):
    inv = nusift.Inventory()
    inv.add("Cs-137", 1.0e20)
    inv.add("Sr-90", 1.0e20)

    plan = [
        nusift.Intervention("Cs", [nusift.Removal("Cs")]),
        nusift.Intervention("half the Cs", [nusift.Removal("Cs", 0.5)]),
        nusift.Intervention("both", [nusift.Removal("Cs"), nusift.Removal("Sr")]),
    ]
    study = nusift.compare_interventions(data, inv, remove_at="30d", at="30y", interventions=plan)

    assert len(study) == 3
    full, half, both = study.effects
    # Linear in the inventory: removing half the caesium is worth exactly half of removing it all.
    assert half.removed == pytest.approx(full.removed / 2.0, rel=1e-9)
    assert both.removed > full.removed
    for effect in study.effects:
        assert effect.response == pytest.approx(study.baseline - effect.removed, rel=1e-9)


@needs_store
def test_an_intervention_that_cannot_mean_anything_is_refused(data):
    inv = nusift.Inventory()
    inv.add("Cs-137", 1.0e20)

    # A response before the removal: taking something out later cannot change an earlier number.
    with pytest.raises(nusift.InputError):
        nusift.compare_interventions(
            data, inv, remove_at="30y", at="30d",
            interventions=[nusift.Intervention("late", [nusift.Removal("Cs")])],
        )

    # A selector naming nothing the chain reaches reads as "worth nothing" when in fact the
    # question never arrived.
    with pytest.raises(nusift.InputError):
        nusift.compare_interventions(
            data, inv, remove_at="0", at="30y",
            interventions=[nusift.Intervention("absent", [nusift.Removal("Pu")])],
        )

    with pytest.raises(nusift.InputError):
        nusift.compare_interventions(
            data, inv, remove_at="0", at="30y",
            interventions=[nusift.Intervention("over", [nusift.Removal("Cs", 1.5)])],
        )

    # The same nuclide taken out twice is not a stated quantity.
    with pytest.raises(nusift.InputError):
        nusift.compare_interventions(
            data, inv, remove_at="0", at="30y",
            interventions=[
                nusift.Intervention("twice", [nusift.Removal("Cs", 0.5),
                                              nusift.Removal("Cs-137", 0.5)])
            ],
        )
