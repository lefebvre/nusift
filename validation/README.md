# Validation

Checks of the shipped evaluation against values NuSIFT was not fitted to, and the generator that
turns them into [`docs/validation.md`](../docs/validation.md).

The unit suite proves the code does what the physics says on chains it invents. This asks a
different question — whether the evaluation NuSIFT actually ships reproduces published physics —
and it can only be asked against a specific data file.

## Layout

| | |
| --- | --- |
| `checks.py` | every check, computed once. Imported by both the tests and the report generator, which is what stops a published number and a gated number drifting apart. |
| `make_report.py` | writes `docs/validation.md` and its figures |
| `references/*.csv` | the published values, one file per table |
| `references/README.md` | full citations, and why every acceptance band is what it is |

The gates live outside this directory, next to the suites they belong to:
[`tests/validation/`](../tests/validation/) for the C++ side and
[`python/tests/test_validation.py`](../python/tests/test_validation.py) for the sweeps.

## Running it

```bash
# the C++ gate: census, published constants through the store, real chains vs closed forms
ctest --test-dir build -L validation --output-on-failure

# the sweeps and the cross-code comparison
PYTHONPATH=python python -m pytest python/tests -m validation

# regenerate the report and its figures, then confirm nothing moved
PYTHONPATH=python python validation/make_report.py
git diff --exit-code -- docs/validation.md docs/figures/
```

The cross-code comparison needs `radioactivedecay`, pinned in `checks.py`:

```bash
pip install radioactivedecay==0.6.1
```

Without it, pytest skips that one test and `make_report.py` fails outright — the committed
report must always contain the section.

## When the store is restaged

`tests/validation/test_store_census.cpp` pins what the store contains: how many nuclides are
staged, how many are unstable, how many emit photons with no evaluated spectrum, how many lines
fall outside the tabulated air-coefficient range. Those numbers are quoted in
[`docs/nuclear-data.md`](../docs/nuclear-data.md) and printed by `nusift data info`, so they are
published facts rather than internal details.

They are **expected to change** when the store is restaged, and updating them is a step of that
work rather than a failure:

1. Restage, and replace `data/nusift_b8.1.h5`.
2. `ctest -L validation` — the census test fails, listing what moved.
3. Update the pinned counts, and check the new ones against `nusift data info`. A count that
   moved by thousands when the evaluation moved by a point release is worth understanding
   before it is written down.
4. `PYTHONPATH=python python validation/make_report.py`, and read the diff. The residuals are
   the interesting part: a band that was 20% used and is now 90% used means the new evaluation
   disagrees with a published value more than the old one did.
5. Update `docs/nuclear-data.md` where it quotes the census.

The pins are tripwires on a deliberately-updated file, not golden files for the physics. Nothing
about *what the numbers should be* is asserted from a stored baseline anywhere in this suite —
every value is checked against a published source, an analytic solution, or an independent code.

## Adding a check

Add the reference value to the CSV for its table, with a band and a reason, and cite the source
in `references/README.md`. Two rules keep the suite honest:

**Choose the band before looking at the residual.** A band chosen afterwards records what the
code does, which is what a golden file does, and proves nothing. Every band in
`references/README.md` is justified from the physics that limits the comparison.

**Mark a row `report` rather than widening a band to swallow it.** A published value and a
computed one are sometimes not the same quantity — a different energy cutoff, a constant
tabulated for a parent in equilibrium with its daughter, an evaluation that simply disagrees.
Those are computed, printed with their residuals, and left ungated with the cause named. The
same applies to `KNOWN_EVALUATION_DIFFERENCES` in `checks.py`: nothing goes in it without the
cause being measured first, or it becomes a place for inconvenient failures to go.
