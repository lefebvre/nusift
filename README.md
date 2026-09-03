# NuSIFT — Nuclear Source-term Isotope Forecasting and Triage

Given an isotopic inventory — born from burnup, activation, or fission — NuSIFT answers a
single question: **which isotopes or mass chains are the top contributors to activity or
exposure, at a specific time or over a time interval?**

It decays the inventory forward with [CRAM](https://github.com/lefebvre/cram-depletion) and
then *ranks* the result, rather than collapsing it to a scalar. That distinction is the whole
design: the engine produces per-nuclide inventories and their exact time integrals, and every
metric is a linear functional of those with a fixed per-nuclide or per-line weight. Ranking by
nuclide, mass chain, element, or individual gamma line — instantaneously or integrated, in Bq
or R/h or Sv/h — is therefore a post-multiply, not a separate code path.

**Status: early.** Ranking, forecasting, and seed attribution work end to end, for activity and
exposure, from an inventory file or from fission, from the CLI or from Python. Activation seeding
is still to come.

## Using it

```bash
# What dominates 30 days after shutdown?
nusift rank -i inventory.csv --at 30d

# Which mass chains are 95% of the activity at 30 years, in curies?
nusift rank -i inventory.csv --at 30y --by mass-chain --units Ci --top 0 --coverage 0.95

# Top 10, and wherever Cs-137 and the A=90 chain happen to fall
nusift rank -i inventory.csv --at 30d --pin Cs-137 --pin A=90

# How many decays occur in the first year?
nusift integrate -i inventory.csv --interval 0,1y

# What dominates the effective dose rate at 2 m, in Sv/h? (ICRP 116, AP by default)
nusift rank -i inventory.csv --at 30d --metric exposure --distance 2 --units Sv/h

# The same field as air kerma instead, which for a soft emitter is a different number
nusift rank -i inventory.csv --at 30d --metric exposure --distance 2 --units Gy/h

# Build the inventory from fission instead of reading one
nusift rank --seed-fission U-235 --energy thermal --yield-kt 20 --at 1h --metric exposure

# Which of the nuclides I seeded is the answer riding on?
nusift attribute --seed-fission U-235 --yield-kt 20 --at 30d --top 8

# Who dominates, and when does it change?
nusift forecast -i inventory.csv --times 1d:300y:log:70 --metric exposure

# Which individual photon lines drive the dose?
nusift spectrum --seed-fission U-235 --yield-kt 20 --at 1h --top 10

# Hand the same photons to a transport code, binned, with the caveats it cannot
# infer written into the deck: what to multiply the tally by, and what is missing
nusift source -i inventory.csv --at 30d --bins 100 --format mcnp > source.i
nusift source -i inventory.csv --interval 0,8h --format openmc > source.py

# When does the total fall below a level? Every answer carries the grid bracket
# that found it -- an event between two samples is not found, and it says so.
nusift when -i inventory.csv --times 1h:300y:log:80 --level 5e13

# When does an ingrowth-fed nuclide stop getting worse?
nusift when -i inventory.csv --times 1s:100y:log:80 --of Y-90

# When should a one-hour job be done, and from when does it fit a budget?
# (each sample is an exact interval integral, so this costs solves the rate curve does not)
nusift when -i inventory.csv --times 1h:50y:log:60 --task 1h --metric exposure --level 0.5

# Decay heat: the power the inventory releases, and the energy over a window.
# Geometry-free -- a property of the material, not of where it is standing.
nusift rank --seed-fission U-235 --yield-kt 20 --at 30d --metric heat

# And the other parameter class: how much does the answer rest on the evaluated
# half-lives? (elasticities, with the refinement chosen from the chain itself)
nusift sensitivity -i inventory.csv --at 30d --metric exposure --units Sv/h

# Given what the sheet says about each row, what is the error bar -- and which
# measurement is it resting on? (exact: R is linear in the seed)
nusift uncertainty -i assay.csv --at 30d --metric exposure --units Sv/h

# And for a fission source, the same question about the EVALUATED yields -- with the
# one correlation no evaluation publishes, imported (data/fycom/fetch_fycom.sh)
nusift uncertainty --seed-fission U-235 --energy thermal --fissions 1e20 --at 30d \
        --yield-covariance data/fycom/ENDF/U235T_corr.csv

# What does the job actually cost, leg by leg -- and where does the budget run out?
nusift plan -i inventory.csv --at 30d --metric exposure --units Sv --budget 2e-3 \
        --leg "approach,3m,4" --leg "valve work,20m,0.8" --leg "break:10m" \
        --leg "reassemble,15m,1.2,0.7" --leg "retreat,3m,4"

# The smallest monitoring list holding 95% of BOTH metrics at every time on the grid
nusift shortlist --seed-fission U-235 --yield-kt 20 --times 1h:100y:log:40 \
        --metrics activity,exposure --coverage 0.95

# Three assay sheets from different dates, brought to one and merged
nusift reconcile -i assays.csv --write merged.csv --write-units Bq

# And the converse: going in at 30 d, how long may someone stay on a dose budget?
nusift stay -i inventory.csv --at 30d --at 1y --at 30y --metric exposure \
        --units Sv --budget 0.02 --distance 2

# Place an event by re-solving inside its bracket instead of interpolating across it
nusift when -i inventory.csv --times 1h:300y:log:80 --level 5e13 --refine

# By what factor can this be scaled before a limit binds, and when does it ship?
nusift allowable -i inventory.csv --times 1h:300y:log:60 --limit 3.7e13 \
        --limit-name "A2 transport"

# What would a Cs or Sr separation before storage actually buy at 30 years?
nusift intervene -i inventory.csv --remove-at 30d --at 30y --remove Cs --remove Sr \
        --metric exposure --units Sv/h

# Rank by a published coefficient table instead -- here the IAEA transport A2 values,
# so the total IS the sum of fractions and 1.0 is where a Type A package binds
nusift rank -i inventory.csv --pack data/packs/iaea-ssr6-a2.csv --at 30d

# What does the data store actually cover?
nusift data info

# What does the store know about one nuclide?
nusift data nuclide Cs-137 Ba-137m
```

An inventory is a three-column CSV; comments, blank lines, a BOM, and an optional header row
are all accepted, and each row carries its own unit:

```
# quantities in whatever unit each source came in
nuclide, quantity, unit
Cs-137,  1.2e14,   Bq
Sr-90,   3.5,      g
Co-60,   0.8,      Ci
```

Every report states the total over *all* contributors and the fraction the shown rows cover,
so a top-10 worth 40% and one worth 99% can never look alike.

A ranking says what is producing the response *now*, which is why a Cs-137 source's exposure
lands on its Ba-137m daughter. `nusift attribute` answers the complementary question — which of
the nuclides you *seeded* the answer is riding on — by running the same weights backwards
through one adjoint solve. For a fission seed the two lists barely overlap: at 30 days the
emitters are La-140 and Pr-143, while the seeds carrying them are Xe-140 and Ba-143, both long
gone. Because decay is linear the shares are an exact partition of the same total, not an
estimate, so the attribution carries a coverage figure like any other ranking.

Every other way of shortening a ranking truncates it; `--pin` is the one that reaches past the
cut. A pinned nuclide, mass chain, or element appears below the ranking whatever it ranks,
carrying the place it actually holds — `27  Cs-137  0.063%  99.6%` — so following one specific
isotope never means printing the whole chain or guessing a `--top` large enough to reach it.
Pinning changes nothing about the ranking above it, and a pin that resolves to nothing is
refused rather than answered with a row of zeros.

Activity and exposure routinely give different answers, which is the point of ranking by the
one you care about. A pure beta emitter can dominate activity and contribute no exposure at
all; a nuclide can dominate exposure through a daughter that emits the photons rather than
itself. NuSIFT attaches photon lines to the nuclide that actually emits them, so a Cs-137
source's exposure is correctly attributed to its Ba-137m daughter and the published gamma
constant falls out of the equilibrium ratio rather than being folded into a table.

Exposure is modelled as an unshielded point source in air: inverse-square spreading, air
attenuation applied per photon line, and air kerma converted to roentgen. Because attenuation
is energy-dependent it sits *inside* the sum over lines, which is why the data store keeps
whole spectra rather than one constant per nuclide -- no single constant is right at more than
one distance. Scatter buildup, source self-absorption, bremsstrahlung, and beta/neutron dose
are not modelled; what a nuclide emits as continuum is recorded and reported, so an
understated row says so rather than looking merely small.

Solves are parallel across time points and deterministic: `--threads` defaults to every core,
and the result is bit-for-bit what the serial path gives. An 80-point forecast over a full
ENDF/B-VIII.1 chain takes well under a second.

## From Python

```python
import nusift

nd  = nusift.NuclearData.open()
inv = nusift.seed_fission(nd, "U-235", energy="thermal", yield_kt=20)
res = nusift.decay(nd, inv, nusift.logspace("1h", "100y", 60))

tab = nusift.response(nd, res, metric="exposure", units="Sv/h",
                      geometry=nusift.PointSource(distance_m=2.0))

for c in tab.rank(at="30d", top=5).contributors:
    print(f"{c.label:10s} {c.value:.3e} Sv/h  {c.fraction:.1%}")

for w in tab.dominance_windows():
    print(f"{w.label} leads {nusift.format_duration(w.start_s)} "
          f"to {nusift.format_duration(w.end_s)}")

# Which seeded nuclide is the 30-day answer riding on? (an exact partition, not an estimate)
att = nusift.attribute(nd, inv, at=nusift.parse_duration("30d"), top=5)
for s in att.shares:
    print(f"{s.label:8s} {s.fraction:6.1%}  one more atom is worth {s.importance:.3e} Bq")

# Totals over a window -- decays, or roentgen accrued -- in closed form, guard included
window = nusift.integrate(nd, inv, "1d", "30d")
for c in nusift.response(nd, window, metric="activity").rank(top=5).contributors:
    print(f"{c.label:10s} {c.value:.3e} decays")
```

`res.atoms` and `tab.values` are zero-copy NumPy views over the C++ storage rather than
copies, so they are cheap to take and compose directly with NumPy. Times and units accept the
same strings the CLI does, parsed by the same code — a notebook and a terminal never disagree
about what `1.5y` means.

Build the extension with `-DNUSIFT_BUILD_PYTHON=ON` (needs `pip install nanobind`), or
`pip install .` to go through scikit-build-core.

A wheel carries the staged store inside the package, which is what lets `NuclearData.open()`
take no argument from anywhere. An explicit path or `$NUSIFT_DATA_STORE` still wins over it,
so a shared evaluation can be used instead of the packaged one.

## Documentation

[**docs/**](docs/README.md) documents the methodology stage by stage — what is evaluated
exactly, what is approximated and by how much, and what is not modelled at all.

| | |
| --- | --- |
| [Nuclear data](docs/nuclear-data.md) | Staging ENDF into a store, chain closure, and why whole photon spectra are kept |
| [Inventory and seeding](docs/inventory.md) | Unit conversions to atoms, and sizing a fission source |
| [The decay solve](docs/decay-solve.md) | CRAM, exact pruning, time grids, threading, determinism |
| [Interval integration](docs/interval-integration.md) | Time-integrated answers in closed form, and the cancellation guard |
| [Exposure](docs/exposure.md) | The point-source photon model, and what it excludes |
| [Ranking and forecasting](docs/ranking.md) | Weights, aggregation, coverage, and dominance windows |
| [Seed attribution](docs/attribution.md) | The other attribution of the same number, from one adjoint solve |
| [**Validation**](docs/validation.md) | Computed against published constants, evaluated data, an empirical law, and an independent code — regenerated and diff-checked in CI |

## Building

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build -L unit --output-on-failure
```

Formatting is a CI gate, and the same check runs locally as a target:

```bash
cmake --build build --target format-check   # what CI runs, same files and same flags
cmake --build build --target format         # the same, rewriting the files in place
```

Both take their file list from `git ls-files` at the moment they run, so a source added since
the last configure is still checked — but one never `git add`ed is invisible to them, as it is
to CI.

On Windows, point CMake at a vcpkg toolchain for HDF5 and Eigen:

```bash
cmake -B build -DCMAKE_TOOLCHAIN_FILE=$VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake
```

### Build options

All options are `NUSIFT_`-prefixed so they cannot collide with the identically-named
`CRAM_*` options when cram is pulled in as a sub-build.

| Option | Default | Purpose |
| --- | --- | --- |
| `NUSIFT_ENABLE_TESTS` | `ON` | Build the GoogleTest suite |
| `NUSIFT_BUILD_CLI` | `ON` | Build the `nusift` command-line tool |
| `NUSIFT_WITH_STAGING` | `OFF` | Build `nusift_stage_data` (requires ENDFtk via cram) |
| `NUSIFT_BUILD_PYTHON` | `OFF` | Build the nanobind extension |
| `NUSIFT_INSTALL` | *auto* | Generate install/export rules — see below |
| `NUSIFT_ENABLE_COVERAGE` | `OFF` | Instrument for gcov |
| `NUSIFT_ENABLE_SANITIZERS` | `OFF` | ASan + UBSan |
| `NUSIFT_CLANG_TIDY` | `OFF` | Run clang-tidy during the build |

### Two build configurations

NuSIFT has two configurations, and they are mutually exclusive by construction:

| | Eigen | ENDFtk | cram | Installable |
| --- | --- | --- | --- | --- |
| **Runtime** (default) | system | off | `find_package(cram)` | **yes** |
| **Staging** (`NUSIFT_WITH_STAGING=ON`) | either | fetched | fetched, `CRAM_WITH_ENDFTK=ON` | no |

The reason is structural, not incidental: a *fetched* dependency is a target in this build
tree that is never installed, so exporting a target whose interface names it is rejected
outright by `install(EXPORT)`. cram documents the same constraint for the same reason, and
NuSIFT mirrors its policy. `NUSIFT_INSTALL` therefore defaults `ON` only for a top-level
build that found both cram and Eigen and is not staging.

This costs nothing in practice: staging is a one-time offline step that reads ENDF tapes and
writes an HDF5 store, which is then committed and shipped. Production runs never link ENDFtk.

**Staging on MSVC needs cram 2.0.0 or later**, which is what `cmake/dep_versions.cmake` pins.
Earlier cram releases built ENDFtk with `SPDLOG_USE_STD_FORMAT`, and njoy's `tools::Log`
forwarded its arguments by value into spdlog, so the format string reached
`std::format_string` as a runtime value rather than a compile-time constant; MSVC rejects that
(`error C7595`) where libstdc++ accepted it. cram 2.0.0 takes an explicit `fmt::format_string`
through njoy/tools 0.4.4, and the staging tool has been built and run against cram's test
tapes on MSVC 19.44. CI still exercises staging on Linux only, because that job compiles
ENDFtk, spdlog, and range-v3 from source, and doubling that cost would buy a toolchain
difference nothing in NuSIFT's own staging code is sensitive to.

## Consuming

```cmake
find_package(nusift REQUIRED)
target_link_libraries(my_target PRIVATE nusift::nusift)
```

Public headers are Eigen-free: `cram::DepletionChain` and Eigen are PIMPL'd out of every
installed header, so a consumer calling the NuSIFT API never compiles an Eigen template.

## Scope

Implemented:

- Inventory input in atoms, moles, mass, or activity units
- Assays dated per row, reconciled to a common epoch and merged -- forward only, by refusal
- Error bars from the assay's own uncertainties, exact rather than first order, ranked by which
  measurement dominates the variance
- Decay-constant sensitivities: which evaluated half-lives the answer rests on, as elasticities
- Evaluated 1-sigma on half-lives, branchings and fission yields, staged from the ENDF tapes
- Decay to a set of cooling times, with exact time integrals over an interval
- Ranking by nuclide, mass chain, element, or gamma line, with any of them pinnable
- Point-source gamma exposure and air kerma in R/h or Gy/h, at any distance
- ICRP 116 effective dose in Sv/h, in any of the six irradiation geometries
- Coefficient packs: a published per-nuclide table as a metric, with its version, scenario and
  coverage carried into the answer
- Staging a data store from ENDF decay and fission-yield tapes
- Seeding an inventory from fission, by fission count, kilotons, or joules
- Dominance forecasting: who leads, and when that changes
- Located events: crossings, turns, level windows, and when a fixed-length task fits a budget
- Stay times: how long a stay beginning at a given time can run before it spends a dose budget
- Robust triage sets: the smallest list holding a coverage floor across every metric and time
- Piecewise task plans: a job as legs with their own distance, occupancy, and breaks
- Seed attribution: which seeded nuclide a response is riding on, as an exact partition
- A binned photon emission spectrum, written as an MCNP `SDEF` card or an OpenMC source
- Python bindings, with zero-copy NumPy views

Planned:

- Seeding an inventory from neutron activation
- Shielding, and decay heat
- Wheels for the three platforms

## License

BSD-3-Clause. See `LICENSE.md`.
