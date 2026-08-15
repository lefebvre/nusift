"""The validation checks, computed once and consumed twice.

Both python/tests/test_validation.py and validation/make_report.py import this module. The
tests assert on the rows; the report prints them. Computing them in one place is what keeps a
green suite and a published table from ever disagreeing about the same number.

Every function returns plain dicts with a fixed set of keys, so the report needs no knowledge
of any particular check:

    key            what is being checked -- a nuclide name, a mass number
    published      the reference value, as its source prints it
    computed       what NuSIFT gives
    unit           the unit both are in
    residual       (computed - published) / published, a fraction
    tolerance      the band this row is accepted within
    gate           True if a test asserts on it, False if it is reported only
    source         short key resolved in validation/references/README.md
    note           why this row deviates, or why it is not gated
    within         whether |residual| <= tolerance

Only numpy and the nusift extension are required. radioactivedecay is imported lazily by the
cross-code check alone, so everything else runs without it.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np

import nusift

REFERENCES = Path(__file__).parent / "references"
STORE = Path(__file__).resolve().parent.parent / "data" / "nusift_b8.1.h5"

# The cross-code comparison is pinned to one release of radioactivedecay. Its ICRP-107 decay
# data is frozen, but a package upgrade could still move a branching or add a nuclide, and an
# unpinned comparison would turn that into a mysterious CI failure on an unrelated commit.
RADIOACTIVEDECAY_PIN = "0.6.1"

# R*m^2/(h*Bq) -> R*cm^2/(h*mCi), the units the rest of this repository quotes gamma constants
# in -- the CLI's `data nuclide`, the unit suite, and docs/exposure.md.
TO_PUBLISHED_GAMMA = 1.0e4 * 3.7e7

# R*m^2/(h*Bq) -> uGy*m^2/(GBq*h), the units the air-kerma reference table is printed in. The
# comparison is done in the source's units and NuSIFT's number is the one converted, so the
# reference table holds exactly what the paper prints and nothing is transcribed twice.
# 0.00876 Gy/R, then micro-gray, then per GBq. For orientation: 309.0 in these units is the
# 13.05 R*cm^2/(h*mCi) quoted elsewhere in this repository.
TO_AIR_KERMA = 0.00876 * 1.0e6 * 1.0e9

# Each reference table is compared in its own units, converting NuSIFT's number rather than the
# published one, so a reference file holds exactly what its source prints.
GAMMA_UNIT_SCALE = {
    "uGy.m2/(GBq.h)": TO_AIR_KERMA,
    "R.cm2/(h.mCi)": TO_PUBLISHED_GAMMA,
}

# Every tabulation of these constants states a low-energy cutoff -- 20 keV for the air-kerma
# table, 15 keV for the exposure one -- below which photons are not counted, because any real
# source encapsulation absorbs them before they reach air. Summing the whole spectrum against
# such a table compares two different quantities: for the X-ray emitters the disagreement
# reaches a factor of ten, all of it convention. The cutoff is carried per row.

# Above this, a constant is carried by photons the air table handles well. Below it the
# tabulation is sparse, the energy-absorption coefficient is turning over at its Compton
# minimum, and X-ray intensity evaluations differ between compilations. A row drawing more than
# SOFT_FRACTION_LIMIT of its constant from below this energy is reported rather than gated --
# stated as a property of the spectrum, decided before the residual is looked at.
SOFT_PHOTON_EV = 100.0e3
SOFT_FRACTION_LIMIT = 0.30

# R*m^2/(h*Bq) -> mSv/(h*MBq) at 1 m, for comparison against tabulated effective-dose
# constants: 0.00876 Gy/R with a photon radiation weighting factor of 1, then per MBq, then
# milli. See docs/exposure.md section 6 -- this column is air kerma wearing a sievert label.
TO_MSV_PER_H_MBQ = 0.00876 * 1.0e6 * 1.0e3

SECONDS_PER = {"s": 1.0, "min": 60.0, "h": 3600.0, "d": 86400.0, "y": 365.25 * 86400.0}

AVOGADRO = 6.02214076e23
BQ_PER_CI = 3.7e10


def open_store():
    return nusift.NuclearData.open(str(STORE))


def mass_number_of(name):
    """Mass number out of a canonical nuclide name: Cs-137 and Ba-137m are both 137."""
    digits = ""
    for character in name.split("-")[1]:
        if not character.isdigit():
            break
        digits += character
    return int(digits)


def load_reference(name):
    """Rows of a reference CSV, in file order. Order is never sorted or grouped downstream, so
    the committed file controls how the report reads.

    Leading `#` lines are comments, as they are in an inventory CSV, and they carry the reason
    each table's band is what it is -- which belongs next to the numbers rather than only in a
    README nobody opens while reading a diff.
    """
    path = REFERENCES / f"{name}.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        lines = [line for line in handle if not line.lstrip().startswith("#")]
    rows = [r for r in csv.DictReader(lines) if r.get("key", "").strip()]
    if not rows:
        raise SystemExit(f"reference table {path.name} has no rows")
    return rows


def _row(key, published, computed, unit, tolerance, gate, source, note):
    residual = (computed - published) / published if published else math.nan
    return {
        "key": key,
        "published": published,
        "computed": computed,
        "unit": unit,
        "residual": residual,
        "tolerance": tolerance,
        "gate": gate,
        "source": source,
        "note": note,
        "within": abs(residual) <= tolerance if published else False,
    }


def _common(reference):
    return (
        float(reference["value"]),
        reference["unit"],
        float(reference["tolerance_rel"]),
        reference["gate"].strip() == "gate",
        reference["source"],
        reference.get("note", "").strip(),
    )


# --- the sweeps -------------------------------------------------------------


def gamma_constant_rows(data):
    """Specific gamma-ray constants, computed from the shipped photon spectra.

    This is the broadest published-value comparison in the suite and the one that exercises the
    most machinery at once: the staged line energies and intensities, the NIST air table, the
    log-log interpolation, and the roentgen conversion, all against numbers NuSIFT was not
    fitted to.
    """
    rows = []
    for reference in load_reference("gamma_constants"):
        published, unit, tolerance, gate, source, note = _common(reference)
        key = reference["key"]
        cutoff = float(reference["cutoff_ev"])
        scale = GAMMA_UNIT_SCALE[unit]

        computed = data.gamma_constant(key, cutoff)
        whole = data.gamma_constant(key)
        # A few entries are tabulated for a parent in equilibrium with its daughter rather than
        # for either nuclide alone -- Cs-137 is the familiar one, where the photons are the
        # daughter's. Reconstructing the pair is the only way to compare against such a row
        # without folding the branch into a constant, which is the practice docs/exposure.md
        # section 5 explains is a live bug rather than a tidiness question.
        partner = reference.get("pair_with", "").strip()
        if partner:
            branch = float(reference["pair_branch"])
            computed += branch * data.gamma_constant(partner, cutoff)
            whole += branch * data.gamma_constant(partner)

        hard = data.gamma_constant(key, SOFT_PHOTON_EV)
        if partner:
            hard += float(reference["pair_branch"]) * data.gamma_constant(partner, SOFT_PHOTON_EV)
        soft_fraction = 1.0 - hard / computed if computed > 0.0 else 0.0

        row = _row(key, published, computed * scale, unit, tolerance, gate, source, note)
        row["soft_fraction"] = soft_fraction
        row["cutoff_ev"] = cutoff
        # What the same spectrum gives with no cutoff, which is what NuSIFT reports by default.
        # Printing it next to the cutoff value is how the report shows the size of the
        # convention rather than merely asserting that one exists.
        row["whole_spectrum"] = whole * scale
        rows.append(row)
    return rows


def molar_mass_rows(data):
    """Staged atomic weights against AME2020.

    The store carries ENDF's atomic weight ratio, and everything expressed per gram goes
    through it: an inventory given in grams, and every specific activity. AME2020 is an
    independent evaluation of the same masses, so this is the one check that pins that
    conversion against something outside the ENDF pipeline.

    It also stands in for a specific-activity sweep, which would be ln(2) N_A / (T_half M) over
    two quantities that already have their own tables here. Published specific-activity tables
    disagree with each other by more than the staging error being looked for -- Am-241 is
    tabulated anywhere between 3.2 and 3.5 Ci/g depending on which half-life the compiler used
    -- so checking the two inputs separately says more than checking their quotient.
    """
    rows = []
    for reference in load_reference("molar_masses"):
        published, unit, tolerance, gate, source, note = _common(reference)
        computed = data.molar_mass(reference["key"])
        rows.append(_row(reference["key"], published, computed, unit, tolerance, gate, source, note))
    return rows


def half_life_rows(data):
    """Staged half-lives against the evaluated compilations.

    Every decay constant, and so every activity and every exposure, is ln(2) over one of these.
    """
    rows = []
    for reference in load_reference("half_lives"):
        published, unit, tolerance, gate, source, note = _common(reference)
        computed = data.half_life(reference["key"]) / SECONDS_PER[unit]
        rows.append(_row(reference["key"], published, computed, unit, tolerance, gate, source, note))
    return rows


def specific_activity(data, nuclide):
    """ln(2) N_A / (T_half M), in Ci/g, from the store's own half-life and atomic weight.

    There is no sweep against a published specific-activity table, because such a table is
    itself this derivation and a sweep would mostly re-test the half-lives that already have
    their own table. The keystone check lives in tests/validation, where it is the only thing
    in the suite that exercises the staged atomic weight at all; see
    validation/references/README.md.
    """
    return math.log(2.0) * AVOGADRO / (data.half_life(nuclide) * data.molar_mass(nuclide)) \
        / BQ_PER_CI


def chain_yield_rows(data, fissions=1.0e20):
    """Cumulative fission-product chain yields for U-235 thermal fission.

    Beta decay preserves mass number, so the cumulative yield of a mass chain is the sum of the
    independent yields of every nuclide on it -- available at t = 0, with no decay solve needed
    and no dependence on the time grid. What this checks is that seeding applies the yield set
    correctly and that the products land on the mass chains they belong to.

    Delayed-neutron emission moves a nucleus off its chain and is the one term this identity
    misses; it is below a tenth of a percent for the chains here and is noted per row where it
    matters.
    """
    inventory = nusift.seed_fission(data, "U-235", energy="thermal", fissions=fissions)
    per_chain = {}
    for name, atoms in zip(inventory.nuclides, inventory.atoms):
        per_chain[mass_number_of(name)] = per_chain.get(mass_number_of(name), 0.0) + atoms

    rows = []
    for reference in load_reference("cumulative_yields"):
        published, unit, tolerance, gate, source, note = _common(reference)
        computed = per_chain.get(int(reference["key"]), 0.0) / fissions * 100.0
        rows.append(_row(reference["key"], published, computed, unit, tolerance, gate, source, note))
    return rows


def icrp116_rows(data):
    """NuSIFT's sievert column against tabulated ICRP 116 effective-dose constants.

    This is the honesty check rather than an agreement check. NuSIFT reports air kerma with a
    photon radiation weighting factor of 1, which is not effective dose to a person; near 1 MeV
    the two happen to coincide, and going soft they do not. Am-241 is gated on DIVERGING by
    about a factor of five, so the caveat can never quietly stop being true.
    """
    rows = []
    for reference in load_reference("icrp116_ratios"):
        published, unit, tolerance, gate, source, note = _common(reference)
        expected_ratio = float(reference["expected_ratio"])
        computed = data.gamma_constant(reference["key"]) * TO_MSV_PER_H_MBQ
        ratio = computed / published
        row = _row(reference["key"], published, computed, unit, tolerance, gate, source, note)
        row["ratio"] = ratio
        row["expected_ratio"] = expected_ratio
        row["residual"] = (ratio - expected_ratio) / expected_ratio
        row["within"] = abs(row["residual"]) <= tolerance
        rows.append(row)
    return rows


# --- the empirical law ------------------------------------------------------

WAY_WIGNER_START = "1h"
WAY_WIGNER_END = "30d"
WAY_WIGNER_POINTS = 25
# Way and Wigner's rule is a fit to gross fission-product behaviour, quoted as good to roughly
# 25% over its validity window rather than as an exact exponent, and the local slope genuinely
# moves across that window -- measured here between -1.07 and -1.23 decade by decade. The band
# is set to admit that spread while still rejecting anything qualitatively wrong: a single
# exponential decaying away, or a chain that never turns over.
WAY_WIGNER_BAND = (-1.35, -1.05)


def way_wigner(data, fissions=1.0e20):
    """The t^-1.2 decay of mixed fission-product activity.

    Way and Wigner's rule is empirical: gross fission-product activity falls off as roughly
    t^-1.2 from minutes to months, and no part of NuSIFT was built to reproduce it. Getting the
    exponent right therefore exercises the yield set, the chain topology, the solve, and the
    activity weights together against something none of them knows about.

    The window is 1 hour to 30 days, inside the rule's classical validity range. Local slopes
    per decade are returned as well, so a future drift shows where the curve bent rather than
    only that it did.
    """
    inventory = nusift.seed_fission(data, "U-235", energy="thermal", fissions=fissions)
    times = nusift.logspace(WAY_WIGNER_START, WAY_WIGNER_END, WAY_WIGNER_POINTS)
    result = nusift.decay(data, inventory, times)
    table = nusift.response(data, result, metric="activity", units="Bq")

    t = np.asarray(times, dtype=float)
    activity = np.asarray(table.totals, dtype=float)
    log_t, log_a = np.log10(t), np.log10(activity)
    slope = float(np.polyfit(log_t, log_a, 1)[0])

    segments = []
    for lo, hi, label in ((3600.0, 36000.0, "1 h - 10 h"), (36000.0, 360000.0, "10 h - 100 h"),
                          (360000.0, 30.0 * 86400.0, "100 h - 30 d")):
        mask = (t >= lo) & (t <= hi)
        if mask.sum() >= 2:
            segments.append((label, float(np.polyfit(log_t[mask], log_a[mask], 1)[0])))

    low, high = WAY_WIGNER_BAND
    return {
        "times": t,
        "activity": activity,
        "slope": slope,
        "segments": segments,
        "band": WAY_WIGNER_BAND,
        "within": low <= slope <= high,
        "source": "glasstone1977",
    }


def equilibrium_curves(data):
    """Two real decay chains solved by the engine, against their closed forms.

    The unit suite proves the solver against Bateman on chains it invents; this runs the same
    comparison on chains built out of the shipped store, so the decay constants and the
    branching are the evaluation's rather than chosen. Sr-90/Y-90 reaches secular equilibrium
    because the parent outlives the daughter four-thousand-fold; Mo-99/Tc-99m reaches transient
    equilibrium and only matches if the 87.6% branch to the isomer was staged correctly.

    The C++ suite gates these; here they are drawn.
    """
    curves = []
    for parent, daughter, span, label in (
        ("Sr-90", "Y-90", ("1h", "60d"), "secular"),
        ("Mo-99", "Tc-99m", ("1h", "14d"), "transient"),
    ):
        seed = 1.0e20
        times = np.asarray(nusift.logspace(span[0], span[1], 60), dtype=float)
        inventory = nusift.Inventory()
        inventory.add(parent, seed)
        result = nusift.decay(data, inventory, list(times))

        names = list(result.nuclides)
        atoms = np.asarray(result.atoms)
        lambda_p = math.log(2.0) / data.half_life(parent)
        lambda_d = math.log(2.0) / data.half_life(daughter)

        engine_parent = atoms[:, names.index(parent)] * lambda_p
        engine_daughter = atoms[:, names.index(daughter)] * lambda_d

        # The two-member Bateman solution, times the branch to this daughter. The branch is
        # recovered from the solve rather than asserted, so the curve says which of the two is
        # wrong if they part company.
        unbranched = seed * lambda_p / (lambda_d - lambda_p) * (
            np.exp(-lambda_p * times) - np.exp(-lambda_d * times))
        branch = float(engine_daughter[0] / (lambda_d * unbranched[0]))
        closed_parent = seed * np.exp(-lambda_p * times) * lambda_p
        closed_daughter = branch * lambda_d * unbranched

        curves.append({
            "parent": parent,
            "daughter": daughter,
            "kind": label,
            "branch": branch,
            "times": times,
            "engine_parent": engine_parent,
            "engine_daughter": engine_daughter,
            "closed_parent": closed_parent,
            "closed_daughter": closed_daughter,
            "worst_residual": float(np.max(np.abs(engine_daughter - closed_daughter)
                                           / np.maximum(closed_daughter, 1e-300))),
        })
    return curves


# --- the independent implementation ------------------------------------------
#
# Five single-parent cases plus a mixed one. The parents are chosen for having well-established
# decay data in both evaluations; nuclides where ICRP-107 and ENDF/B-VIII.1 are known to
# disagree substantially are deliberately absent, because this check is meant to catch an error
# in NuSIFT rather than to survey differences between data libraries.

CROSS_CODE_CASES = [
    {"name": "Co-60", "seed": {"Co-60": 1.0e20}, "parent": "Co-60", "tolerance": 0.02},
    {"name": "Cs-137", "seed": {"Cs-137": 1.0e20}, "parent": "Cs-137", "tolerance": 0.02},
    {"name": "Sr-90", "seed": {"Sr-90": 1.0e20}, "parent": "Sr-90", "tolerance": 0.02},
    {"name": "I-131", "seed": {"I-131": 1.0e20}, "parent": "I-131", "tolerance": 0.02},
    {"name": "Mo-99", "seed": {"Mo-99": 1.0e20}, "parent": "Mo-99", "tolerance": 0.02},
    {
        "name": "mixed",
        "seed": {"Cs-137": 1.0e20, "Sr-90": 5.0e19, "Co-60": 2.0e19, "I-131": 1.0e18,
                 "Mo-99": 1.0e18},
        "parent": None,
        "tolerance": 0.03,
    },
]

# Nuclides below this share of a case's total activity are dropped: a species at 1e-12 of the
# total is numerically meaningless in both codes and its ratio is noise, not disagreement.
CROSS_CODE_FLOOR = 1.0e-6

# Where the two evaluations themselves disagree, and the disagreement is documented rather than
# smoothed away with a wider band. These rows are still computed and still printed with their
# residuals; they are held out of the gate because gating them would be gating ICRP-107 against
# ENDF/B-VIII.1, which is not a statement about NuSIFT.
#
# Each entry names the measured cause. Adding one without measuring the cause first would turn
# this into a place for inconvenient failures to go, which is exactly what it must not be.
KNOWN_EVALUATION_DIFFERENCES = {
    "Xe-131m": (
        "I-131 feeds the metastable state with branching 0.0108477 in ENDF/B-VIII.1 against "
        "0.011759 in ICRP-107, a ratio of 0.9225. The residual is that ratio, is the same at "
        "every time, and is therefore the branching rather than the solve."
    ),
}


def _case_times(data, case):
    if case["parent"] is None:
        return list(nusift.logspace("1h", "30y", 8))
    half_life = data.half_life(case["parent"])
    cap = 100.0 * SECONDS_PER["y"]
    return [min(f * half_life, cap) for f in (0.1, 0.5, 1.0, 2.0, 5.0)]


def cross_code_curve(data, parent, span=("1h", "14d"), points=40):
    """One chain on a dense grid, decayed by both codes, for plotting.

    cross_code_rows compares at a handful of times chosen as multiples of the parent's
    half-life, which is the right sampling for a gate and too sparse to draw. This runs the same
    two codes over a fine grid so the agreement can be seen as well as measured.
    """
    import radioactivedecay as rd

    seed = {parent: 1.0e20}
    times = np.asarray(nusift.logspace(span[0], span[1], points), dtype=float)

    inventory = nusift.Inventory()
    inventory.add(parent, seed[parent])
    result = nusift.decay(data, inventory, list(times))
    names = list(result.nuclides)
    atoms = np.asarray(result.atoms)

    reference = rd.Inventory(seed, "num")
    series = {}
    for name in names:
        half_life = data.half_life(name)
        if half_life <= 0.0:
            continue
        ours = math.log(2.0) / half_life * atoms[:, names.index(name)]
        if ours.max() < 1.0e-3 * (math.log(2.0) / data.half_life(parent) * seed[parent]):
            continue
        theirs = np.array([float(reference.decay(float(t), "s").activities("Bq").get(name, 0.0))
                           for t in times])
        series[name] = (ours, theirs)
    return {"parent": parent, "times": times, "series": series}


def cross_code_rows(data):
    """Per-nuclide activities against radioactivedecay, an independent implementation.

    Different solver (matrix exponential via sympy/numpy rather than CRAM), different decay data
    (ICRP-107 rather than ENDF/B-VIII.1), different author. Agreement therefore says something
    neither an analytic check nor a published constant can: that two unrelated implementations
    of the same physics land in the same place.

    The residual is not expected to be zero and should not be. ICRP-107 and ENDF/B-VIII.1 differ
    on half-lives at the tenth-of-a-percent level and on some branchings by a few tenths -- the
    Cs-137 half-life alone differs by 0.29% -- so a few percent is data, and anything larger is
    NuSIFT.

    Returns (fine_rows, case_summaries). Raises ImportError if radioactivedecay is absent; the
    caller decides whether that is a skip or a failure.
    """
    import radioactivedecay as rd

    fine, summaries = [], []
    for case in CROSS_CODE_CASES:
        times = _case_times(data, case)

        inventory = nusift.Inventory()
        for name, atoms in case["seed"].items():
            inventory.add(name, atoms)
        result = nusift.decay(data, inventory, times)
        lambdas = {n: (math.log(2.0) / data.half_life(n) if data.half_life(n) > 0 else 0.0)
                   for n in result.nuclides}
        atoms = np.asarray(result.atoms)

        reference_inventory = rd.Inventory(case["seed"], "num")
        worst = 0.0
        residuals = []
        for k, t in enumerate(times):
            ours = {n: lambdas[n] * atoms[k][i] for i, n in enumerate(result.nuclides)}
            total = sum(ours.values())
            theirs = reference_inventory.decay(float(t), "s").activities("Bq")

            for name, mine in sorted(ours.items()):
                if total <= 0.0 or mine < CROSS_CODE_FLOOR * total:
                    continue
                other = float(theirs.get(name, 0.0))
                if other <= 0.0:
                    continue
                residual = (mine - other) / other
                excluded = name in KNOWN_EVALUATION_DIFFERENCES
                if not excluded:
                    residuals.append(abs(residual))
                    worst = max(worst, abs(residual))
                fine.append({
                    "case": case["name"],
                    "time_s": float(t),
                    "key": name,
                    "published": other,
                    "computed": mine,
                    "unit": "Bq",
                    "residual": residual,
                    "tolerance": case["tolerance"],
                    "gate": False,
                    "excluded": excluded,
                    "source": "radioactivedecay",
                    "note": KNOWN_EVALUATION_DIFFERENCES.get(name, ""),
                    "within": excluded or abs(residual) <= case["tolerance"],
                })

        summaries.append({
            "key": case["name"],
            "max_residual": worst,
            "median_residual": float(np.median(residuals)) if residuals else math.nan,
            "comparisons": len(residuals),
            "tolerance": case["tolerance"],
            "gate": True,
            "within": worst <= case["tolerance"],
            "source": "radioactivedecay",
        })

    return fine, summaries
