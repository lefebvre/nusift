# Figure data

Committed so `make_figures.py` regenerates the figures without a built CLI or a staged store.
Every file is verbatim CLI output from ENDF/B-VIII.1 (`data/nusift_b8.1.h5`, staged
2026-08-13), with `S="--store data/nusift_b8.1.h5"` and `F="--seed-fission U-235 --yield-kt 20"`:

| File | Command |
| --- | --- |
| `fission-exposure-timeline.csv` | `nusift rank $F $S --times 1m:100y:log:70 --metric exposure --top 8 --format csv` |
| `fission-dominance-windows.csv` | `nusift forecast $F $S --times 1m:100y:log:70 --metric exposure --format csv` |
| `fission-1d-activity.csv` | `nusift rank $F $S --at 1d --metric activity --top 8 --format csv` |
| `fission-1d-exposure.csv` | `nusift rank $F $S --at 1d --metric exposure --top 8 --format csv` |

The air-coefficient figures need no data file: the script parses the NIST table straight out of
`nusift/exposure/air_coefficients.cpp`, so those figures cannot drift from the values the code
actually interpolates.
