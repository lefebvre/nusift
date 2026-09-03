#!/usr/bin/env bash
#
# Regenerate the committed CLI output the scenario figures are drawn from.
#
#   ./docs/figures/data/make_data.sh [path/to/nusift] [path/to/store.h5]
#
# Every scenario figure in docs/scenarios.md is built from one of these files rather than from
# numbers typed into a plotting script, so a figure cannot claim something the tool does not
# print. The files are committed for the same reason the validation report is: `make_figures.py`
# then needs neither a built CLI nor a staged store, and a figure that changed is a diff on the
# DATA rather than an unexplained change in a picture.
#
# Run this after anything that could move the numbers -- a restage, a new metric, a changed
# weight -- and commit whatever moved.
set -euo pipefail

NUSIFT="${1:-./build/dev-linux/nusift_apps/nusift}"
STORE="${2:-data/nusift_b8.1.h5}"
HERE="$(cd "$(dirname "$0")" && pwd)"
# Commands run from the repository ROOT with RELATIVE paths, because several of them echo the
# input path into their own provenance line. An absolute path there would bake this checkout's
# location into committed output and make the regenerate-and-diff check fail on every other
# machine -- which is the same reason the validation report publishes bounds rather than measured
# round-off.
ROOT="$(cd "$HERE/../../.." && pwd)"
cd "$ROOT"
EXAMPLES="examples"
DATA="docs/figures/data"

command -v "$NUSIFT" >/dev/null 2>&1 || [ -x "$NUSIFT" ] || {
  echo "error: no nusift binary at $NUSIFT" >&2; exit 1; }

S=(--store "$STORE")
# Three objects, one per track, and each is the one the tool's own model actually fits.
#
#   F  a 20 kt fission source -- the SOURCE TERM questions, asked in geometry-free quantities
#      (activity, photon strength, decay heat) because debris on the ground is a distributed
#      source and the point kernel is not a model of it;
#   C  a contaminated valve body, drum-scale and compact, where a point source at a few metres
#      IS the model, so the dose-at-a-distance questions belong here;
#   I  a tonne of irradiated fuel, for the waste and disposal questions.
F=(--seed-fission U-235 --energy thermal --yield-kt 20)
C=(-i "$EXAMPLES/contaminated_component.csv")
I=(-i "$EXAMPLES/spent_fuel.csv")

run() { # run <output file> <args...>
  local out="$1"; shift
  "$NUSIFT" "$@" --format csv -o "$DATA/$out"
  printf '  %-38s %s\n' "$out" "$(wc -l < "$DATA/$out") rows"
}

text() { # text <output file> <args...>
  local out="$1"; shift
  "$NUSIFT" "$@" > "$DATA/$out"
  printf '  %-38s %s\n' "$out" "$(wc -l < "$DATA/$out") lines"
}

echo "regenerating scenario data with $NUSIFT against $STORE"

# --- track 0: knowing what you have ------------------------------------------
text scenario-store-info.txt       data info "${S[@]}"
run scenario-decay.csv             decay "${S[@]}" "${F[@]}" --times 1h:100y:log:60 --top 8
# Deliberately misspelled on the way in, to show what the canonicalizer accepts and what it
# refuses to blur: an isomer is a different nuclide, however casually it was typed.
text scenario-nuclide.txt          nuclide Cs-137 ba137m "TC-99M" u235 Co60 Ba-137
# The same object in the unit each row was reported in, and again in atoms. The conversion is
# per nuclide because a becquerel is a RATE and an atom count is an amount.
"$NUSIFT" inventory convert "${S[@]}" "${C[@]}" --units atoms > "$DATA/scenario-inventory-atoms.csv"
printf '  %-38s %s\n' scenario-inventory-atoms.csv "$(wc -l < "$DATA/scenario-inventory-atoms.csv") rows"

# --- track 1: the source term, in quantities that need no geometry -----------
run scenario-rank-activity.csv     rank "${S[@]}" "${F[@]}" --at 1d --metric activity --top 10
run scenario-rank-heat.csv         rank "${S[@]}" "${F[@]}" --at 1d --metric heat --units W --top 10
run scenario-spectrum.csv          spectrum "${S[@]}" "${F[@]}" --at 1d --by line --top 14 \
                                     --metric photon --units photons/s
run scenario-forecast.csv          forecast "${S[@]}" "${F[@]}" --times 1m:100y:log:70 --metric activity
# Decay heat crossing a handling threshold: geometry-free, so the date is a property of the
# material rather than of where anyone is standing.
run scenario-when.csv              when "${S[@]}" "${F[@]}" --times 1h:100y:log:60 \
                                     --metric heat --units W --level 1000 --refine --peaks
run scenario-shortlist.csv         shortlist "${S[@]}" "${F[@]}" --times 1h:100y:log:40 \
                                     --metrics activity,heat --coverage 0.95
run scenario-source.csv            source "${S[@]}" "${F[@]}" --at 1d --bins 40
run scenario-attribute.csv         attribute "${S[@]}" "${F[@]}" --at 30d --metric activity --top 10

# --- track 2: working around a compact source --------------------------------
run scenario-rank-component.csv    rank "${S[@]}" "${C[@]}" --at 5y --metric exposure \
                                     --units Sv/h --distance 2 --top 8
run scenario-stay.csv              stay "${S[@]}" "${C[@]}" --metric exposure --units Sv \
                                     --distance 2 --budget 0.02 --max-stay 8h \
                                     --at 1y --at 5y --at 10y --at 20y --at 30y --at 50y --at 75y
run scenario-plan.csv              plan "${S[@]}" "${C[@]}" --at 5y --metric exposure --units Sv \
                                     --budget 0.02 \
                                     --leg "approach,3m,4" --leg "survey,15m,2" \
                                     --leg "break:10m" --leg "cut and cap,45m,1.2,0.8" \
                                     --leg "retreat,3m,4"
run scenario-intervene.csv         intervene "${S[@]}" "${C[@]}" --at 30y --metric exposure \
                                     --units Sv/h --distance 2 --remove Cs --remove Co --remove Eu
run scenario-sensitivity.csv       sensitivity "${S[@]}" "${C[@]}" --at 30y --metric exposure \
                                     --units Sv/h --distance 2 --top 12

# --- track 3: waste, transport and disposal ----------------------------------
run scenario-allowable.csv         allowable "${S[@]}" "${I[@]}" --pack data/packs/iaea-ssr6-a2.csv \
                                     --limit 1.0 --times 1y:1000y:log:50
run scenario-uncertainty.csv       uncertainty "${S[@]}" -i "$DATA/scenario-assay.csv" --at 5y \
                                     --metric exposure --units Sv/h --distance 2
run scenario-reconcile.csv         reconcile "${S[@]}" -i "$DATA/scenario-assay.csv"
# The rate at the end of a window against the total accrued over it. Same source, same metric,
# same instant of asking -- and disjoint answers, which is what the interval path is for.
#
# In DECAY HEAT rather than exposure, and that is not a stylistic choice: an exposure answer here
# would be a point-source kernel pointed at dispersed debris, which is the pairing this guide
# tells the reader not to make. Heat needs no geometry, so the contrast survives the correction.
run scenario-integrate-rate.csv    rank "${S[@]}" "${F[@]}" --at 10y --metric heat --units W --top 8
run scenario-integrate-total.csv   integrate "${S[@]}" "${F[@]}" --interval 0,10y --metric heat \
                                     --units J --top 8
# How far inside the window the endpoint contrast holds, so the claim about it can be checked
# rather than generalized from one instant.
run scenario-integrate-1h.csv      rank "${S[@]}" "${F[@]}" --at 1h --metric heat --units W --top 8
run scenario-integrate-1d.csv      rank "${S[@]}" "${F[@]}" --at 1d --metric heat --units W --top 8

echo "done"
