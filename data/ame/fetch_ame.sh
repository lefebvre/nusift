#!/usr/bin/env bash
#
# Fetch the AME2020 atomic mass table NuSIFT stages masses from, into $1 (default ./AME2020).
#
#   ./fetch_ame.sh [destination]
#
# One file, and deliberately only one:
#
#   mass_1.mas20.txt  the unrounded atomic mass table of
#                     W.J. Huang, M. Wang, F.G. Kondev, G. Audi, S. Naimi, "The AME2020 atomic
#                     mass evaluation (I)", Chinese Physics C 45 (2021) 030002, and
#                     M. Wang et al., "(II)", Chinese Physics C 45 (2021) 030003.
#
# WHY THIS IS HERE AT ALL. A molar mass is the one thing NuSIFT needs that the ENDF decay
# sublibrary does not always carry. ENDF states an atomic weight ratio in the head record of
# every MF8/MT457 section, so every nuclide with a decay evaluation is covered -- but a chain
# also contains nuclides registered by closure or named by a fission-yield set that have no
# decay evaluation at all, and those have no AWR anywhere in the three sublibraries staged.
# Without a mass they cannot be given or reported in grams. AME2020 is the evaluation that
# covers them: it is the source ENDF's own AWR values derive from, which is why the two agree
# to within rounding wherever both exist (staging cross-checks exactly that and reports drift).
#
# The rounded table (massround.mas20.txt) is NOT used: it is the published-precision copy, and
# there is no reason to stage fewer digits than the evaluation has.
#
# ESTIMATED VALUES ARE TAKEN. AME marks extrapolated masses with '#' in place of the decimal
# point, and about a thousand of the table's entries are extrapolations rather than
# measurements. They are staged, and flagged as estimated so the store can say so, because the
# quantity being derived is a molar mass: even a 2 MeV extrapolation uncertainty is a part in
# 1e5 of the mass of a fission product, which is far below anything a gram conversion resolves.
#
# The unpacked file is gitignored and regenerable; only the staged .h5 store is committed.
set -euo pipefail

BASE="https://www-nds.iaea.org/amdc/ame2020"
TABLE="mass_1.mas20.txt"
DEST="${1:-$(dirname "$0")/AME2020}"

command -v curl >/dev/null || { echo "error: curl is required" >&2; exit 1; }

mkdir -p "$DEST"
if [ -s "$DEST/$TABLE" ]; then
  echo "already present: $DEST/$TABLE"
else
  tmp="$(mktemp)"
  trap 'rm -f "$tmp"' EXIT
  echo "fetching $BASE/$TABLE"
  curl -fsS -o "$tmp" "$BASE/$TABLE"
  # A truncated download is a silently short mass table, which would stage as a coverage gap
  # rather than an error. The table is ~470 kB and ~3550 nuclides; anything tiny is a failure.
  if [ "$(wc -c < "$tmp")" -lt 400000 ]; then
    echo "error: $TABLE came back short -- refusing to stage a truncated mass table" >&2
    exit 1
  fi
  mv "$tmp" "$DEST/$TABLE"
  trap - EXIT
fi

cat <<EOF

Stage from it with:

  nusift_stage_data \\
    --decay-dir "data/endf/B-VIII.1/decay" \\
    --nfy-dir   "data/endf/B-VIII.1/nfpy" \\
    --sfy-dir   "data/endf/B-VIII.1/sfpy" \\
    --ame       "$DEST/$TABLE" \\
    -o data/nusift_b8.1.h5 --library "ENDF/B-VIII.1"
EOF
