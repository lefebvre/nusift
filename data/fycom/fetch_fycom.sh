#!/usr/bin/env bash
#
# Fetch the published fission-yield CORRELATION matrices NuSIFT imports, into $1 (default ./ENDF).
#
#   ./fetch_fycom.sh [destination] [system ...]
#
# The matrices are from:
#
#   E.F. Matthews, L.A. Bernstein, W. Younes, "Stochastically estimated covariance matrices for
#   independent and cumulative fission yields in the ENDF/B-VIII.0 and JEFF-3.3 evaluations",
#   Atomic Data and Nuclear Data Tables 140 (2021) 101441, doi:10.1016/j.adt.2021.101441
#
# published as FYCoM (doi:10.5281/zenodo.5985936). The author releases the DATA for public use on
# the condition that it is cited, and reserves all rights to the code that generated it. So this
# script fetches rather than vendors, and `nusift uncertainty` prints the citation with every
# answer built on a matrix. Neither is optional.
#
# THREE THINGS ARE DELIBERATELY NARROW ABOUT WHAT THIS TAKES.
#
#   INDEPENDENT, never cumulative. NuSIFT decays the chain explicitly, so a cumulative yield --
#   which already includes everything a precursor decays into -- would count every precursor
#   decay twice. Same reason the ENDF staging reads MT454 and offers no option to do otherwise.
#
#   CORRELATION, never covariance. Only the off-diagonal is imported; the variances come from the
#   sigma_Y the store stages from ENDF/B-VIII.1. That is the paper's own "normalized covariance"
#   construction -- correlation times evaluated variance -- one edition later, and it makes the
#   diagonal of the assembled matrix identical to the figure a diagonal treatment reports, so the
#   difference between the two answers is the correlation and nothing else. `*_normed_corr.csv`
#   is byte-identical to `*_corr.csv` for this reason: renormalizing variances cannot change a
#   correlation.
#
#   ENDF, not JEFF, by default. The store is an ENDF/B-VIII.1 artifact and the nearer pairing is
#   the defensible one. JEFF matrices exist at the same paths under matrices/JEFF and can be
#   fetched by editing EVALUATION below -- the pairing is declared to the tool with
#   --covariance-library either way, so nothing downstream has to guess.
#
# The systems default to the six NuSIFT ships examples for. Every system in the publication is
# available; the codes are <nuclide><spelling of the energy>, where T is thermal, F is fast, H is
# 14 MeV (DT), DD is 2 MeV, and SF is spontaneous. See yields/systems.txt in the repository for
# the full list of 58.
#
# Each file is 15-30 MB, the download is not resumable, and the destination is gitignored: these
# are a fetched input like the ENDF tapes, not a checked-in artifact.
set -euo pipefail

BASE="https://raw.githubusercontent.com/efmatthews/FYCoM/master/matrices"
EVALUATION="ENDF"
DEST="${1:-$(dirname "$0")/ENDF}"
shift || true

# U-235 and Pu-239 thermal are the reactor and criticality workhorses; U-235 fast and 14 MeV and
# U-238 fast and 14 MeV are the energies a weapons-relevant seed reaches. U-238 has no thermal
# entry, because U-238 does not fission thermally -- the evaluation tabulates what exists.
SYSTEMS=("$@")
if [ ${#SYSTEMS[@]} -eq 0 ]; then
  SYSTEMS=(U235T U235F U235H Pu239T U238F U238H)
fi

command -v curl >/dev/null || { echo "error: curl is required" >&2; exit 1; }

mkdir -p "$DEST"
echo "fetching ${#SYSTEMS[@]} $EVALUATION correlation matrices into $DEST"

failed=0
for system in "${SYSTEMS[@]}"; do
  target="$DEST/${system}_corr.csv"
  # Idempotent, like the ENDF fetch: a file already here is left alone, so an interrupted run is
  # simply repeated. Delete one to re-fetch it.
  if [ -s "$target" ]; then
    echo "  $system: present, skipped"
    continue
  fi
  tmp="$(mktemp)"
  if curl -sfL --max-time 600 -o "$tmp" \
      "$BASE/$EVALUATION/independent/${system}_corr.csv"; then
    mv "$tmp" "$target"
    echo "  $system: $(du -h "$target" | cut -f1)"
  else
    rm -f "$tmp"
    # One missing system costs that system, not the whole run -- and a mistyped code is the
    # likeliest cause, so it is named rather than swallowed.
    echo "  warn: could not fetch $system (is that a system code? see yields/systems.txt)" >&2
    failed=$((failed + 1))
  fi
done

cat <<EOF

Fetched to $DEST

Use one with a fission seed:

  nusift uncertainty --seed-fission U-235 --energy thermal --fissions 1e20 --at 30d \\
    --yield-covariance "$DEST/U235T_corr.csv"

The matrix must be the one for the fissioning system being seeded: it is keyed by fission
product, so another system's file shares products with this one and would contract against
correlations that were never estimated for it. NuSIFT refuses a file sharing no product at all,
which catches the gross case and not the subtle one.
EOF

exit $((failed > 0 ? 1 : 0))
