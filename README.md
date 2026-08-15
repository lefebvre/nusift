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

**Status: early.** Ranking and forecasting work end to end, for activity and exposure, from an
inventory file or from fission, from the CLI or from Python. Activation seeding is still to
come.

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

# What dominates the exposure rate at 2 m, in Sv/h?
nusift rank -i inventory.csv --at 30d --metric exposure --distance 2 --units Sv/h

# Build the inventory from fission instead of reading one
nusift rank --seed-fission U-235 --energy thermal --yield-kt 20 --at 1h --metric exposure

# Who dominates, and when does it change?
nusift forecast -i inventory.csv --times 1d:300y:log:70 --metric exposure

# Which individual photon lines drive the dose?
nusift spectrum --seed-fission U-235 --yield-kt 20 --at 1h --top 10

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

## Building

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build -L unit --output-on-failure
```

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

**The staging configuration does not build with MSVC.** It is a GCC/Clang path only. cram
builds ENDFtk with `SPDLOG_USE_STD_FORMAT`, and njoy's `tools::Log` forwards its arguments by
value into spdlog, so the format string reaches `std::format_string` as a runtime value rather
than a compile-time constant. MSVC rejects that (`error C7595`); libstdc++ accepts it. The
failure is in cram's own translation units, before anything of NuSIFT's is reached, so there
is nothing to work around on this side. Stage on Linux, WSL, or macOS. The runtime
configuration — the library, the CLI, and the tests — builds fine on MSVC.

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
- Decay to a set of cooling times, with exact time integrals over an interval
- Ranking by nuclide, mass chain, element, or gamma line, with any of them pinnable
- Point-source gamma exposure in R/h, Gy/h, or Sv/h, at any distance
- Staging a data store from ENDF decay and fission-yield tapes
- Seeding an inventory from fission, by fission count, kilotons, or joules
- Dominance forecasting: who leads, and when that changes
- Python bindings, with zero-copy NumPy views

Planned:

- A binned photon spectrum for detector-response work
- Seeding an inventory from neutron activation
- Shielding, and decay heat
- Wheels for the three platforms

## License

BSD-3-Clause. See `LICENSE.md`.
