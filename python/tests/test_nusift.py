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


PACKS = Path(__file__).resolve().parents[2] / "data" / "packs"
PACK = PACKS / "iaea-ssr6-a2.csv"
needs_pack = pytest.mark.skipif(not PACK.is_file(), reason="no shipped pack")


@needs_store
@needs_pack
def test_a_pack_is_a_metric_carrying_its_own_provenance(data):
    """A pack answer is only interpretable with its version and scenario, so the binding has to
    hand those over rather than just the numbers."""
    pack = nusift.load_pack(str(PACK))
    assert pack.name == "iaea-ssr6-a2"
    assert pack.version == "2012 edition"
    assert pack.unit == "1"
    assert pack.basis == "activity"
    assert pack.folds_progeny is True
    assert pack.size > 300
    assert "SSR-6" in pack.source

    # 1/A2 for Cs-137, whose A2 is 6e-1 TBq.
    assert pack.coefficient("Cs-137") == pytest.approx(1.0 / 6.0e11, rel=1e-6)
    # SSR-6 has no Ba-137m row: its contribution lives inside its parent's value.
    assert pack.covers("Ba-137m") is False
    assert pack.folded_into("Ba-137m") == ["Cs-137"]
    # Chains nest, so one daughter can sit under several parents.
    assert len(pack.folded_into("Tl-208")) > 1


@needs_store
@needs_pack
def test_the_intake_packs_carry_the_published_coefficients(data):
    """Values published widely enough to be recognised on sight. They also guard the extraction:
    every exponent's minus sign is a glyph the PDF's text layer drops, so a transcription that
    lost it would put these out by twenty-two orders of magnitude rather than a little."""
    ingestion = nusift.load_pack(str(PACKS / "icrp119-ingestion-worker.csv"))
    assert ingestion.unit == "Sv"
    assert ingestion.basis == "activity"
    assert ingestion.folds_progeny is False
    assert "worker" in ingestion.scenario

    assert ingestion.coefficient("Cs-137") == pytest.approx(1.3e-8, rel=1e-6)
    assert ingestion.coefficient("I-131") == pytest.approx(2.2e-8, rel=1e-6)
    assert ingestion.coefficient("Sr-90") == pytest.approx(2.8e-8, rel=1e-6)

    inhalation = nusift.load_pack(str(PACKS / "icrp119-inhalation-worker-5um.csv"))
    assert inhalation.coefficient("Pu-239") == pytest.approx(3.2e-5, rel=1e-6)
    assert inhalation.coefficient("Am-241") == pytest.approx(2.7e-5, rel=1e-6)
    # ICRP gives these no particulate coefficient at all -- they are gases, in another table --
    # so they are absent rather than approximated, and coverage will say so.
    assert inhalation.covers("H-3") is False
    assert ingestion.covers("H-3") is True


@needs_store
@needs_pack
def test_h10_is_the_conservative_estimator_of_effective_dose(data):
    """Two kernels transcribed from two publications by two different routes -- ICRP 116 out of a
    text layer, ICRP 74 off an image scan -- reproducing the relationship between the quantities
    they describe. H*(10) is designed to be a conservative estimate of effective dose for photons
    up to about 10 MeV, and on a caesium field it exceeds it by the expected fifth."""
    pack = nusift.load_pack(str(PACKS / "icrp74-ambient-dose-h10.csv"))
    assert pack.unit == "Sv/s"
    assert pack.size == 0, "a kernel carries no nuclides; it applies to whatever lines exist"

    inv = nusift.Inventory()
    inv.add("Cs-137", 1.0e15)
    res = nusift.decay(data, inv, [86400.0])
    geometry = nusift.PointSource(distance_m=1.0)

    ambient = nusift.response(data, res, geometry=geometry,
                              pack=nusift.resolve_pack(pack, data, inv))
    effective = nusift.response(data, res, metric="exposure", units="Sv/h", geometry=geometry)

    per_hour = float(ambient.totals[0]) * 3600.0
    ratio = per_hour / float(effective.totals[0])
    assert 1.1 < ratio < 1.3, f"H*(10) should exceed effective dose by about a fifth here: {ratio}"

    # Almost all of the answer comes from lines inside the curve's range, but not quite all:
    # caesium's L X-rays sit below ICRP 74's 10 keV floor and take a clamped value. That is the
    # coverage figure doing its job -- it is a share of the ANSWER, so a few soft lines carrying
    # a ten-thousandth of the dose show up as a ten-thousandth rather than as a warning.
    coverage = float(ambient.pack_coverage[0])
    assert 0.999 < coverage < 1.0


@needs_store
@needs_pack
def test_a_concentration_pack_needs_the_extent_it_is_spread_through(data):
    """An inventory is atoms and a concentration is atoms over an extent, so the extent is an
    input rather than a property of the material -- and it is half of what the answer means."""
    pack = nusift.load_pack(str(PACKS / "fgr15-air-submersion-adult.csv"))
    assert pack.basis == "concentration"
    assert pack.per == "m3"
    assert pack.unit == "Sv/s"
    assert pack.folds_progeny is False

    inv = nusift.Inventory()
    inv.add("Cs-137", 1.0e14)
    res = nusift.decay(data, inv, [0.0, 86400.0])

    with pytest.raises(nusift.InputError):
        nusift.resolve_pack(pack, data, inv)

    small = nusift.response(data, res, pack=nusift.resolve_pack(pack, data, inv, extent=1.0e6))
    large = nusift.response(data, res, pack=nusift.resolve_pack(pack, data, inv, extent=2.0e6))
    # Twice the volume is half the concentration and half the dose rate.
    assert float(large.totals[1]) == pytest.approx(float(small.totals[1]) / 2.0)

    # Caesium's own coefficient is negligible; the dose comes from the barium the chain grows in,
    # which is exactly the combination FGR-15 tells its readers to make for themselves.
    labels = list(small.labels)
    assert small.rank(at=86400.0, top=1).contributors[0].label == "Ba-137m"
    assert "Cs-137" in labels


@needs_store
@needs_pack
def test_intake_and_external_hazard_rank_differently(data):
    """The whole argument for ranking by the metric you care about, a fourth time. Caesium leads
    the ingestion hazard and strontium the inhalation hazard, on one inventory at one instant."""
    inv = nusift.Inventory()
    inv.add("Cs-137", 1.0e18)
    inv.add("Sr-90", 1.0e17)
    res = nusift.decay(data, inv, [0.0])

    def leader(name):
        pack = nusift.load_pack(str(PACKS / name))
        table = nusift.response(data, res, pack=nusift.resolve_pack(pack, data, inv))
        return table.rank(top=1).contributors[0].label

    assert leader("icrp119-ingestion-worker.csv") == "Cs-137"
    assert leader("icrp119-inhalation-worker-5um.csv") == "Sr-90"


@needs_store
@needs_pack
def test_a_folded_daughter_is_not_counted_twice(data):
    """The regulation's own rule, which only a chain-tracking tool can apply: yttrium-90 shipped
    alone is limited by its own A2, and yttrium-90 accompanying strontium-90 is already inside
    its parent's."""
    pack = nusift.load_pack(str(PACK))

    together = nusift.Inventory()
    together.add("Sr-90", 1.0e18)
    alone = nusift.Inventory()
    alone.add("Y-90", 1.0e18)

    times = [0.0]
    with_parent = nusift.response(
        data, nusift.decay(data, together, times),
        pack=nusift.resolve_pack(pack, data, together))
    on_its_own = nusift.response(
        data, nusift.decay(data, alone, times),
        pack=nusift.resolve_pack(pack, data, alone))

    labels = list(with_parent.labels)
    yttrium = labels.index("Y-90")
    assert float(with_parent.values[0][yttrium]) == 0.0
    assert float(on_its_own.totals[0]) > 0.0

    # Covered either way, and the coverage figure says so: folded is accounted for, not missing.
    assert float(with_parent.pack_coverage[0]) == pytest.approx(1.0)


@needs_store
@needs_pack
def test_pack_coverage_says_what_the_pack_could_not_speak_for(data):
    """A sum of fractions over a quarter of the activity looks identical to one over all of it."""
    pack = nusift.load_pack(str(PACK))
    inv = nusift.seed_fission(data, "U-235", energy="thermal", fissions=1e20)
    res = nusift.decay(data, inv, [3600.0])
    table = nusift.response(data, res, pack=nusift.resolve_pack(pack, data, inv))

    assert table.unit == "1"
    coverage = float(table.pack_coverage[0])
    assert 0.0 < coverage < 0.5, "a fresh fission mix is mostly nuclides SSR-6 does not list"
    assert table.rank(top=3).pack_coverage == pytest.approx(coverage)


@needs_store
def test_the_sievert_column_is_effective_dose_and_names_its_geometry(data, result):
    """Sv is ICRP 116 effective dose, Gy is air kerma, and they are different quantities rather
    than two spellings of one. On Am-241 -- soft, and the nuclide the repair was for -- they are
    a factor of four apart."""
    _, res = result
    ap = nusift.PointSource(distance_m=1.0)
    assert ap.irradiation == "AP"

    sievert = nusift.response(data, res, metric="exposure", units="Sv/h", geometry=ap)
    gray = nusift.response(data, res, metric="exposure", units="Gy/h", geometry=ap)
    assert float(sievert.totals[0]) != float(gray.totals[0])

    # The orientation reaches the numbers, and is carried on the table that was built with it.
    rot = nusift.PointSource(distance_m=1.0, irradiation="rot")
    assert rot.irradiation == "ROT"
    rotated = nusift.response(data, res, metric="exposure", units="Sv/h", geometry=rot)
    assert float(rotated.totals[0]) < float(sievert.totals[0])

    with pytest.raises(nusift.InputError):
        nusift.PointSource(irradiation="sideways")


@needs_store
def test_effective_dose_constant_matches_the_published_coefficient(data):
    """The gamma constant's counterpart for the quantity a sievert names. Checked against ICRP
    116 as tabulated by Peplow (2020), in mSv/(h.MBq) at 1 m: the same three nuclides the
    validation suite gates, one of which used to be off by a factor of five."""
    for nuclide, published in (("Co-60", 3.062e-4), ("Ba-137m", 8.228e-5), ("Am-241", 5.413e-6)):
        computed = data.effective_dose_constant(nuclide) * 1.0e9
        assert computed == pytest.approx(published, rel=0.03), nuclide

    # Air kerma has not gone anywhere; it is simply no longer what a sievert means. For Am-241
    # it is still five times the effective dose, which is the size of what was repaired.
    kerma = data.gamma_constant("Am-241") * 0.00876 * 1.0e9
    assert kerma / (data.effective_dose_constant("Am-241") * 1.0e9) > 4.0

    assert data.effective_dose_constant("Am-241", irradiation="iso") < data.effective_dose_constant(
        "Am-241", irradiation="ap"
    )
    with pytest.raises(nusift.InputError):
        data.effective_dose_constant("Am-241", irradiation="sideways")


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


# --- the binned source term --------------------------------------------------


@needs_store
def test_binned_spectrum_arrays_are_views_and_edges_bound_the_bins(data, result):
    _, res = result
    spectrum = nusift.binned_spectrum(data, res, time_index=0, bins=32)

    assert spectrum.edges_ev.base is not None
    assert spectrum.emission.base is not None
    # One more boundary than bin, in every direction the binding is read from. A snippet built
    # on the other convention is accepted by OpenMC and silently one bin wrong.
    assert len(spectrum.edges_ev) == len(spectrum.emission) + 1
    assert np.all(np.diff(spectrum.edges_ev) > 0)
    assert spectrum.unit == "photons/s"
    assert spectrum.domain == "instant"


@needs_store
def test_binned_spectrum_conserves_every_photon(data, result):
    _, res = result
    # A grid narrow enough that emission falls off both ends, which is the case the identity
    # exists to make visible.
    spectrum = nusift.binned_spectrum(data, res, bins=16, min_ev=1.0e5, max_ev=4.0e5)
    accounted = spectrum.emission.sum() + spectrum.below_range + spectrum.above_range
    assert accounted == pytest.approx(spectrum.total, rel=1e-12)
    assert spectrum.above_range > 0.0


@needs_store
def test_the_source_total_is_the_photon_metric_total(data, result):
    _, res = result
    spectrum = nusift.binned_spectrum(data, res, time_index=3)
    table = nusift.response(data, res, metric="photon", units="photons/s")
    # Two front doors onto one quantity. If these ever diverge, the number that leaves in a
    # deck is the one nobody re-derives.
    assert spectrum.total == pytest.approx(table.totals[3], rel=1e-12)


@needs_store
def test_an_interval_source_is_a_count_of_photons(data):
    inv = nusift.Inventory()
    inv.add("Cs-137", 1.0e20)
    window = nusift.integrate(data, inv, "1h", "30d")
    spectrum = nusift.binned_spectrum(data, window, bins=8)

    assert spectrum.unit == "photons"
    assert spectrum.domain == "interval"
    assert spectrum.t1 == nusift.parse_duration("1h")
    assert spectrum.t2 == nusift.parse_duration("30d")


@needs_store
def test_explicit_edges_win_and_cannot_be_combined_with_a_range(data, result):
    _, res = result
    edges = [1.0e4, 1.0e5, 1.0e6, 3.0e6]
    spectrum = nusift.binned_spectrum(data, res, edges_ev=edges)
    assert list(spectrum.edges_ev) == edges

    with pytest.raises(nusift.InputError):
        nusift.binned_spectrum(data, res, edges_ev=edges, min_ev=500.0)
    with pytest.raises(nusift.InputError):
        nusift.binned_spectrum(data, res, edges_ev=[1.0e5])
    with pytest.raises(nusift.InputError):
        nusift.binned_spectrum(data, res, scale="quadratic")


@needs_store
def test_source_decks_carry_the_caveats_a_transport_code_cannot_infer(data, result):
    inv, res = result
    spectrum = nusift.binned_spectrum(data, res, bins=12)

    sdef = nusift.source_deck(data, inv, spectrum, format="mcnp")
    assert "SDEF PAR=P ERG=D1" in sdef
    assert "SHAPE, NOT A STRENGTH" in sdef
    assert sdef.count("\nSI1 H") == 1

    snippet = nusift.source_deck(data, inv, spectrum, format="openmc")
    assert 'interpolation="histogram"' in snippet
    assert "strength=" in snippet
    # The snippet has to be Python, not merely Python-shaped: a deck that does not compile is
    # found by whoever pastes it, at the worst moment.
    compile(snippet, "<openmc deck>", "exec")

    with pytest.raises(nusift.InputError):
        nusift.source_deck(data, inv, spectrum, format="serpent")


# --- the stay time -----------------------------------------------------------


@needs_store
def test_stay_time_inverts_the_task_curve_it_is_the_converse_of(data):
    inv = nusift.Inventory()
    inv.add("Cs-137", 1.0e18)
    geometry = nusift.PointSource(distance_m=2.0)

    stay = nusift.stay_time(
        data, inv, at="30d", budget=1.0e-4, metric="exposure", units="Sv", geometry=geometry,
    )
    assert stay.bounded
    assert stay.duration_s > 0.0
    assert stay.converged
    assert 0.0 < stay.located_to_s < stay.duration_s

    # A task of exactly that length, starting at the same instant, must cost exactly the budget.
    # The two inversions are one integral read two ways.
    series = nusift.task_series(
        data, inv, starts=["30d", "31d"], duration=stay.duration_s,
        metric="exposure", units="Sv", geometry=geometry,
    )
    assert series.values[0] == pytest.approx(1.0e-4, rel=1e-5)


@needs_store
def test_an_unspent_budget_is_none_and_not_zero(data):
    inv = nusift.Inventory()
    inv.add("Cs-137", 1.0e18)

    stay = nusift.stay_time(
        data, inv, at="30d", budget=1.0e9, metric="exposure", units="Sv",
        geometry=nusift.PointSource(distance_m=2.0),
    )
    assert not stay.bounded
    # The whole point of the property. A 0.0 here would read as "leave immediately", which is
    # the exact inverse of what an unspent budget means.
    assert stay.duration_s is None
    assert 0.0 < stay.accrued_at_max < stay.budget


@needs_store
def test_waiting_buys_a_longer_stay(data):
    inv = nusift.Inventory()
    inv.add("Co-60", 1.0e18)
    geometry = nusift.PointSource(distance_m=1.0)

    def stay(at):
        return nusift.stay_time(
            data, inv, at=at, budget=5.0e-4, max_stay="30d",
            metric="exposure", units="Sv", geometry=geometry,
        ).duration_s

    # Co-60 has no ingrowth feeding it, so the source only gets weaker and the same budget can
    # only buy longer. Monotone in the direction that makes the root unique.
    assert stay("1d") < stay("5y") < stay("20y")


@needs_store
def test_a_stay_budget_must_be_an_accrued_total(data):
    inv = nusift.Inventory()
    inv.add("Cs-137", 1.0e18)

    # A rate is the wrong dimension for a budget, refused where every other unit mismatch is.
    with pytest.raises(nusift.InputError):
        nusift.stay_time(data, inv, at="30d", budget=0.02, metric="exposure", units="Sv/h")
    with pytest.raises(nusift.InputError):
        nusift.stay_time(data, inv, at="30d", budget=0.0, metric="exposure", units="Sv")
    with pytest.raises(nusift.InputError):
        nusift.stay_time(
            data, inv, at="30d", budget=0.02, max_stay="0s", metric="exposure", units="Sv",
        )


# --- assays on different dates -----------------------------------------------


def test_dates_parse_and_round_trip():
    assert nusift.parse_date("1970-01-01") == 0.0
    assert nusift.format_date(nusift.parse_date("2024-03-15")) == "2024-03-15"
    # A leap day that exists and one that does not.
    assert nusift.parse_date("2024-02-29")
    with pytest.raises(nusift.InputError):
        nusift.parse_date("2023-02-29")
    # Ambiguous anywhere is refused everywhere: 15/03 is two different days by country.
    with pytest.raises(nusift.InputError):
        nusift.parse_date("15/03/2024")


@needs_store
def test_reading_a_dated_inventory_reconciles_it(data, tmp_path):
    path = tmp_path / "assays.csv"
    path.write_text(
        "nuclide,quantity,unit,assayed\n"
        "Cs-137,1.0e14,Bq,2024-03-15\n"
        "Sr-90,5.0e13,Bq,2023-01-10\n"
    )

    groups = nusift.read_assays(str(path), data)
    assert [g.date for g in groups] == ["2023-01-10", "2024-03-15"]

    reconciled = nusift.reconcile(data, groups)
    assert reconciled.epoch == "2024-03-15"
    assert reconciled.span_s > 0.0
    assert [c.carried_s > 0.0 for c in reconciled.contributions] == [True, False]

    # read_inventory does the same thing on its own, because a date the tool ignored would be
    # worse than one it refused.
    assert len(nusift.read_inventory(str(path), data)) == len(reconciled.inventory)
    # Y-90 grew in during the carry and is in the seed; it was in neither sheet.
    assert "Y-90" in reconciled.inventory.nuclides


@needs_store
def test_an_epoch_earlier_than_an_assay_is_refused(data, tmp_path):
    path = tmp_path / "assays.csv"
    path.write_text(
        "Cs-137,1.0e14,Bq,2024-03-15\n"
        "Sr-90,5.0e13,Bq,2023-01-10\n"
    )
    groups = nusift.read_assays(str(path), data)

    # Un-growing a daughter has no unique answer, so this is refused rather than caveated.
    with pytest.raises(nusift.InputError):
        nusift.reconcile(data, groups, epoch="2023-06-01")
    # Later than every assay is fine: that is just carrying everything further forward.
    assert nusift.reconcile(data, groups, epoch="2030-01-01").epoch == "2030-01-01"


@needs_store
def test_a_file_dates_every_row_or_none(data, tmp_path):
    mixed = tmp_path / "mixed.csv"
    mixed.write_text("Cs-137,1.0e14,Bq,2024-03-15\nSr-90,5.0e13,Bq\n")
    with pytest.raises(nusift.InputError):
        nusift.read_assays(str(mixed), data)

    undated = tmp_path / "undated.csv"
    undated.write_text("Cs-137,1.0e14,Bq\nSr-90,5.0e13,Bq\n")
    groups = nusift.read_assays(str(undated), data)
    assert len(groups) == 1


# --- the robust triage set ---------------------------------------------------


@needs_store
def test_a_shortlist_holds_its_floor_at_every_time_and_metric(data, result):
    _, res = result
    tables = {
        "activity": nusift.response(data, res, metric="activity", units="Bq"),
        "exposure": nusift.response(data, res, metric="exposure", units="R/h"),
    }
    chosen = nusift.shortlist(tables, coverage=0.95)

    assert len(chosen) == len(chosen.members)
    assert chosen.constraints == 2 * len(tables["activity"].times)
    assert not chosen.shortfalls

    # Checked against the tables rather than against what the search reported: at every time of
    # every requirement, the chosen columns hold the floor.
    keep = {m.label for m in chosen.members}
    for table in tables.values():
        columns = [i for i, name in enumerate(table.labels) if name in keep]
        for k, total in enumerate(table.totals):
            if total > 0.0:
                assert table.values[k][columns].sum() / total >= 0.95 - 1e-12

    # The binding point is the honest headline: it says whether the set clears its floor
    # comfortably or by a thousandth at one instant.
    assert chosen.binding.achieved >= chosen.binding.required - 1e-12
    assert chosen.binding.required == 0.95


@needs_store
def test_a_shortlist_beats_the_union_of_the_per_time_answers(data, result):
    _, res = result
    tables = {
        "activity": nusift.response(data, res, metric="activity", units="Bq"),
        "exposure": nusift.response(data, res, metric="exposure", units="R/h"),
    }
    chosen = nusift.shortlist(tables, coverage=0.95)

    # What `rank --coverage 0.95` repeated at every time and every metric would give. It is not
    # wrong -- it over-delivers -- but it is not minimal, and nothing about it was chosen.
    union = set()
    for table in tables.values():
        for k, total in enumerate(table.totals):
            if total <= 0.0:
                continue
            row = table.values[k]
            order = np.argsort(-row)
            cumulative = np.cumsum(row[order]) / total
            union.update(order[: int(np.searchsorted(cumulative, 0.95)) + 1].tolist())

    assert len(chosen.members) < len(union)


@needs_store
def test_a_shortlist_is_reproducible_and_refuses_what_it_cannot_answer(data, result):
    _, res = result
    table = nusift.response(data, res, metric="activity", units="Bq")

    first = nusift.shortlist({"activity": table}, coverage=0.9)
    second = nusift.shortlist({"activity": table}, coverage=0.9)
    # A monitoring list that changed between runs would be worse than no list.
    assert [m.label for m in first.members] == [m.label for m in second.members]

    with pytest.raises(nusift.InputError):
        nusift.shortlist({}, coverage=0.9)
    with pytest.raises(nusift.InputError):
        nusift.shortlist({"activity": table}, coverage=1.5)
    # A nuclide column and a mass-chain column are not the same kind of thing.
    with pytest.raises(nusift.InputError):
        nusift.shortlist({
            "by nuclide": table,
            "by chain": nusift.response(data, res, metric="activity", units="Bq", by="mass-chain"),
        })


# --- a job with a shape ------------------------------------------------------


@needs_store
def test_a_plan_partitions_its_total_across_the_legs(data):
    inv = nusift.Inventory()
    inv.add("Cs-137", 1.0e16)
    legs = [
        nusift.PlanLeg("approach", "3m", distance_m=4.0),
        nusift.PlanLeg("valve work", "20m", distance_m=0.8),
        nusift.PlanLeg("break", "10m", occupancy=0.0),
        nusift.PlanLeg("retreat", "3m", distance_m=4.0),
    ]
    plan = nusift.task_plan(data, inv, at="30d", legs=legs, metric="exposure", units="Sv")

    assert len(plan) == 4
    assert sum(leg.accrued for leg in plan.legs) == pytest.approx(plan.total, rel=1e-12)
    assert sum(leg.fraction for leg in plan.legs) == pytest.approx(1.0, rel=1e-12)
    # The clock runs through the break; the time that earns the dose does not.
    assert plan.elapsed_s == 36 * 60
    assert plan.exposed_s == 26 * 60
    assert plan.legs[2].is_break
    assert plan.legs[2].accrued == 0.0
    # The legs abut: no gap the plan did not name.
    for before, after in zip(plan.legs, plan.legs[1:]):
        assert after.start_s == before.end_s


@needs_store
def test_the_close_leg_is_the_expensive_one_even_when_it_is_short(data):
    inv = nusift.Inventory()
    inv.add("Co-60", 1.0e15)
    legs = [
        nusift.PlanLeg("far and long", "60m", distance_m=5.0),
        nusift.PlanLeg("near and short", "5m", distance_m=0.5),
    ]
    plan = nusift.task_plan(data, inv, at="1y", legs=legs, metric="exposure", units="Sv")

    long_leg, short_leg = plan.legs
    # Twelve times the time at ten times the distance: inverse square wins, and the mean-rate
    # column is what says so. This is the whole reason to break a job into legs.
    assert short_leg.accrued > long_leg.accrued
    assert short_leg.mean_rate > long_leg.mean_rate * 50


@needs_store
def test_a_plan_says_where_a_budget_runs_out(data):
    inv = nusift.Inventory()
    inv.add("Cs-137", 1.0e16)
    legs = [
        nusift.PlanLeg("a", "10m", distance_m=1.0),
        nusift.PlanLeg("b", "10m", distance_m=1.0),
        nusift.PlanLeg("c", "10m", distance_m=1.0),
    ]
    whole = nusift.task_plan(data, inv, at="30d", legs=legs, metric="exposure", units="Sv")
    budget = whole.legs[0].accrued + 0.5 * whole.legs[1].accrued

    plan = nusift.task_plan(
        data, inv, at="30d", legs=legs, metric="exposure", units="Sv", budget=budget,
    )
    assert plan.budget_spent
    assert plan.spent_in_leg == 1
    # Cs-137 barely decays over ten minutes, so half the leg's dose is half the leg's time.
    assert plan.spent_at_s == pytest.approx(plan.legs[1].start_s + 300.0, abs=5.0)

    fits = nusift.task_plan(
        data, inv, at="30d", legs=legs, metric="exposure", units="Sv", budget=whole.total * 10,
    )
    assert not fits.budget_spent
    # None rather than -1: a parser reading -1 as an index would be reading the last leg.
    assert fits.spent_in_leg is None
    assert fits.spent_at_s is None


@needs_store
def test_a_plan_refuses_a_leg_that_is_not_one(data):
    inv = nusift.Inventory()
    inv.add("Cs-137", 1.0e16)
    good = nusift.PlanLeg("a", "10m", distance_m=1.0)

    with pytest.raises(nusift.InputError):
        nusift.task_plan(data, inv, at="30d", legs=[], metric="exposure", units="Sv")
    with pytest.raises(nusift.InputError):
        nusift.task_plan(data, inv, at="30d", legs=[nusift.PlanLeg("", "10m", distance_m=1.0)],
                         metric="exposure", units="Sv")
    with pytest.raises(nusift.InputError):
        nusift.task_plan(data, inv, at="30d",
                         legs=[nusift.PlanLeg("a", "10m", distance_m=1.0, occupancy=1.5)],
                         metric="exposure", units="Sv")
    # A leg accrues a total, so a rate unit is the wrong dimension for it.
    with pytest.raises(nusift.InputError):
        nusift.task_plan(data, inv, at="30d", legs=[good], metric="exposure", units="Sv/h")


# --- the error bar the assay puts on the answer ------------------------------


@needs_store
def test_uncertainty_is_exact_and_ranks_by_variance(data, tmp_path):
    path = tmp_path / "assay.csv"
    path.write_text(
        "nuclide,quantity,unit,uncertainty\n"
        "Cs-137,1.0e14,Bq,3%\n"
        "Sr-90,5.0e13,Bq,25%\n"
    )
    assays = nusift.read_assays(str(path), data)
    u = nusift.uncertainty(data, assays, at="30d")

    assert u.sigma > 0.0
    assert u.covered_fraction == pytest.approx(1.0)
    assert u.rows_without_sigma == 0

    # sigma_R is the quadrature sum of the per-row contributions, and the variance fractions
    # partition. Shares of R do not: they are a different decomposition of a different quantity.
    assert sum(s.variance_fraction for s in u.seeds) == pytest.approx(1.0, rel=1e-12)
    assert math.hypot(*[s.sigma_contribution for s in u.seeds]) == pytest.approx(u.sigma, rel=1e-12)

    # The point of the ordering: the loosely measured row leads the variance even though the
    # well-measured one is the larger share of the answer.
    assert u.seeds[0].label == "Sr-90"
    assert u.seeds[0].share < u.seeds[1].share


@needs_store
def test_the_propagated_response_is_the_ordinary_total(data, tmp_path):
    # The identity that says the per-assay adjoints ran at the right times. For a dated file
    # these are genuinely different solves reaching the same number.
    path = tmp_path / "dated.csv"
    path.write_text(
        "nuclide,quantity,unit,assayed,uncertainty\n"
        "Cs-137,1.0e14,Bq,2024-03-15,3%\n"
        "Sr-90,5.0e13,Bq,2023-01-10,25%\n"
    )
    assays = nusift.read_assays(str(path), data)
    u = nusift.uncertainty(data, assays, at="30d")

    merged = nusift.read_inventory(str(path), data)
    table = nusift.response(
        data, nusift.decay(data, merged, [nusift.parse_duration("30d")]),
        metric="activity", units="Bq",
    )
    assert u.response == pytest.approx(table.totals[0], rel=1e-9)
    # The older sheet was carried; the one defining the epoch was not.
    carried = {s.label: s.carried_s for s in u.seeds}
    assert carried["Sr-90"] > 0.0
    assert carried["Cs-137"] == 0.0


@needs_store
def test_importance_reorders_the_error_bar_by_metric(data, tmp_path):
    path = tmp_path / "assay.csv"
    path.write_text(
        "nuclide,quantity,unit,uncertainty\n"
        "Co-60,2.0e13,Bq,12%\n"
        "Sr-90,5.0e13,Bq,25%\n"
    )
    assays = nusift.read_assays(str(path), data)

    by_activity = nusift.uncertainty(data, assays, at="30d")
    by_exposure = nusift.uncertainty(data, assays, at="30d", metric="exposure", units="Sv/h")

    # Sr-90 is a pure beta emitter, so it dominates the activity error bar and contributes
    # essentially nothing to the photon-dose one. Which measurement to improve depends on the
    # question being asked, which is the whole reason importance is in the product.
    assert by_activity.seeds[0].label == "Sr-90"
    assert by_exposure.seeds[0].label == "Co-60"
    assert by_exposure.seeds[-1].variance_fraction < 1e-6


@needs_store
def test_rows_stating_no_uncertainty_are_reported_not_guessed(data, tmp_path):
    path = tmp_path / "partial.csv"
    path.write_text(
        "nuclide,quantity,unit,uncertainty\n"
        "Cs-137,1.0e14,Bq,3%\n"
        "Sr-90,5.0e13,Bq,\n"
    )
    assays = nusift.read_assays(str(path), data)
    u = nusift.uncertainty(data, assays, at="30d")

    assert u.rows_with_sigma == 1
    assert u.rows_without_sigma == 1
    # An error bar propagated from rows holding part of the answer is not an error bar on the
    # answer, and the covered fraction is what lets a reader see that.
    assert 0.0 < u.covered_fraction < 1.0
    silent = next(s for s in u.seeds if s.label == "Sr-90")
    # None rather than 0.0, and no contribution invented for it.
    assert silent.sigma_atoms is None
    assert silent.sigma_contribution == 0.0


@needs_store
def test_uncertainty_refuses_a_question_it_cannot_answer(data, tmp_path):
    path = tmp_path / "assay.csv"
    path.write_text("nuclide,quantity,unit,uncertainty\nCs-137,1.0e14,Bq,3%\n")
    assays = nusift.read_assays(str(path), data)

    # The shares are instantaneous; an accrued total needs the integrated adjoint.
    with pytest.raises(nusift.InputError):
        nusift.uncertainty(data, assays, at="30d", units="decays")
    with pytest.raises(nusift.InputError):
        nusift.uncertainty(data, assays, at=-1.0)
    # An undated sheet has no date to carry FROM, so it cannot be carried to any epoch. The
    # placeholder zero taken at face value would age it by fifty-four years, silently.
    with pytest.raises(nusift.InputError):
        nusift.uncertainty(data, assays, at="30d", epoch="2024-01-01")

    # And a dated one cannot be carried BACKWARD, which is the inverse problem reconcile()
    # refuses on the same terms.
    dated_path = tmp_path / "dated.csv"
    dated_path.write_text(
        "nuclide,quantity,unit,assayed,uncertainty\nCs-137,1.0e14,Bq,2024-03-15,3%\n"
    )
    dated = nusift.read_assays(str(dated_path), data)
    with pytest.raises(nusift.InputError):
        nusift.uncertainty(data, dated, at="30d", epoch="2000-01-01")
