"""Validation of the shipped evaluation against published values and an independent code.

These are physics tests, unlike test_nusift.py next door, which deliberately tests only the
binding. The rows they assert on are computed by validation/checks.py -- the same module the
report generator uses -- so a number that appears in docs/validation.md and a number CI gates
on can never be two different numbers.

What each row is compared against, and why its band is what it is, lives with the reference
data in validation/references/README.md rather than here.

    pytest python/tests -m validation

Rows marked `report` in the reference tables are computed and printed but not asserted. That is
how a documented convention difference -- a published constant that applies an energy cutoff
NuSIFT does not, say -- stays visible instead of being buried under a widened tolerance.
"""

from __future__ import annotations

import pytest

from validation import checks

needs_store = pytest.mark.skipif(not checks.STORE.is_file(), reason="no staged data store")

pytestmark = [pytest.mark.validation, needs_store]


@pytest.fixture(scope="module")
def data():
    return checks.open_store()


def _assert_gated(rows):
    """Every gated row inside its band, reported one failure per row rather than one per run."""
    failures = [
        f"{row['key']}: computed {row['computed']:.6g} against {row['published']:.6g} "
        f"{row['unit']} -- {row['residual'] * 100:+.2f}%, band {row['tolerance'] * 100:.1f}% "
        f"[{row['source']}]"
        for row in rows
        if row["gate"] and not row["within"]
    ]
    assert not failures, "\n".join(failures)


def _assert_has_gated_rows(rows, minimum):
    """A reference table emptied by a bad edit would otherwise make every check above vacuous."""
    gated = sum(1 for row in rows if row["gate"])
    assert gated >= minimum, f"only {gated} gated rows; the reference table looks truncated"


# --- published values --------------------------------------------------------


def test_gamma_constants_match_published_values(data):
    rows = checks.gamma_constant_rows(data)
    _assert_has_gated_rows(rows, 15)
    _assert_gated(rows)


def test_half_lives_match_the_evaluated_values(data):
    rows = checks.half_life_rows(data)
    _assert_has_gated_rows(rows, 15)
    _assert_gated(rows)


def test_molar_masses_match_ame2020(data):
    rows = checks.molar_mass_rows(data)
    _assert_has_gated_rows(rows, 10)
    _assert_gated(rows)


def test_chain_yields_match_the_evaluated_cumulative_yields(data):
    rows = checks.chain_yield_rows(data)
    _assert_has_gated_rows(rows, 8)
    _assert_gated(rows)


# The sievert column is air kerma with a photon weighting factor of 1, not effective dose to a
# person. Co-60 and Ba-137m agree with ICRP 116 by coincidence of energy; Am-241 does not, and
# asserting that it still diverges by about five is what stops the caveat quietly expiring.
def test_the_sievert_column_tracks_and_departs_from_effective_dose(data):
    rows = checks.icrp116_rows(data)
    _assert_has_gated_rows(rows, 3)
    _assert_gated(rows)

    americium = next(row for row in rows if row["key"] == "Am-241")
    assert americium["ratio"] > 3.0, (
        "Am-241 air kerma should overstate ICRP 116 effective dose several-fold; "
        f"ratio is {americium['ratio']:.2f}"
    )


# --- an empirical law nothing here was fitted to ------------------------------


def test_fission_product_activity_follows_the_way_wigner_slope(data):
    result = checks.way_wigner(data)
    low, high = result["band"]
    assert result["within"], (
        f"log-log slope {result['slope']:.3f} outside [{low}, {high}] over "
        f"{checks.WAY_WIGNER_START} to {checks.WAY_WIGNER_END}"
    )


# --- an independent implementation --------------------------------------------


def test_activities_agree_with_radioactivedecay(data):
    pytest.importorskip(
        "radioactivedecay",
        reason=f"pip install radioactivedecay=={checks.RADIOACTIVEDECAY_PIN}",
    )
    _, summaries = checks.cross_code_rows(data)
    assert summaries, "no cross-code cases ran"

    failures = [
        f"{s['key']}: worst disagreement {s['max_residual'] * 100:.2f}% over "
        f"{s['comparisons']} comparisons, band {s['tolerance'] * 100:.0f}%"
        for s in summaries
        if not s["within"]
    ]
    assert not failures, "\n".join(failures)

    assert sum(s["comparisons"] for s in summaries) >= 50, "too few nuclides actually compared"
