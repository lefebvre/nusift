#!/usr/bin/env bash
#
# Fetch the ENDF/B-VIII.1 sublibraries NuSIFT stages from, into $1 (default ./B-VIII.1).
#
#   ./fetch_endf.sh [destination]
#
# Three sublibraries, and deliberately only three:
#
#   decay  MF8/MT457 -- half-lives, decay modes, branching ratios, average decay energies,
#                       and the discrete photon spectra that exposure is computed from
#   nfpy   MF8/MT454 -- independent neutron-induced fission yields
#   sfpy   MF8/MT454 -- independent spontaneous fission yields
#
# The neutron transport sublibrary is NOT fetched. It is multi-gigabyte, and NuSIFT has no use
# for it: cross sections would only matter for activation, and a defensible one-group value
# cannot be produced from MF3 alone anyway (the resonance region lives in MF2 and needs
# reconstruction, which is NJOY's job).
#
# Idempotent: a tape already unpacked at the destination is skipped, so an interrupted run can
# simply be repeated. Downloads run 16-way parallel because the tapes are small and round-trip
# latency dominates.
#
# The unpacked tree is gitignored and regenerable; only the staged .h5 store is committed.
set -euo pipefail

BASE="https://www-nds.iaea.org/public/download-endf/ENDF-B-VIII.1"
DEST="${1:-$(dirname "$0")/B-VIII.1}"
SUBLIBRARIES=(decay nfpy sfpy)
PARALLEL=16

for tool in curl unzip xargs; do
  command -v "$tool" >/dev/null || { echo "error: $tool is required" >&2; exit 1; }
done

fetch_one() {
  local sublibrary="$1" archive="$2"
  local unpacked="$DEST/$sublibrary/${archive%.zip}.dat"
  [ -s "$unpacked" ] && return 0
  local tmp
  tmp="$(mktemp)"
  if ! curl -sf --max-time 60 -o "$tmp" "$BASE/$sublibrary/$archive"; then
    echo "  warn: could not fetch $archive" >&2
    rm -f "$tmp"
    return 0   # one missing tape costs that nuclide, not the whole run
  fi
  unzip -oqj "$tmp" -d "$DEST/$sublibrary" 2>/dev/null || echo "  warn: bad archive $archive" >&2
  rm -f "$tmp"
}
export -f fetch_one
export BASE DEST

for sublibrary in "${SUBLIBRARIES[@]}"; do
  mkdir -p "$DEST/$sublibrary"
  echo "listing $sublibrary ..."
  listing="$(mktemp)"
  curl -sf --max-time 120 "$BASE/$sublibrary/" \
    | grep -oE 'href="[^"]*\.zip"' | sed 's/href="//;s/"//' | sort -u > "$listing"
  echo "  $(wc -l < "$listing") tapes"

  export sublibrary
  xargs -a "$listing" -P "$PARALLEL" -I{} bash -c 'fetch_one "$sublibrary" "{}"'
  rm -f "$listing"
  echo "  $sublibrary: $(find "$DEST/$sublibrary" -type f | wc -l) files"
done

cat <<EOF

Fetched to $DEST

Stage a store with (a staging build -- Linux, WSL, or macOS; see the README):

  nusift_stage_data \\
    --decay-dir "$DEST/decay" \\
    --nfy-dir   "$DEST/nfpy" \\
    --sfy-dir   "$DEST/sfpy" \\
    -o data/nusift_b8.1.h5 --library "ENDF/B-VIII.1"
EOF
