#!/usr/bin/env bash
# Rerun the 26Q3 addendum. Does not download: the DepMap portal is behind human verification.
# Place CRISPRGeneEffect.csv and Model.csv in ../data/26Q3/ (or set DEPMAP_26Q3_DIR to a folder that contains them).
set -euo pipefail
cd "$(dirname "$0")"
DEST="$(cd .. && pwd)/data/26Q3"
mkdir -p "$DEST" results
if [ -n "${DEPMAP_26Q3_DIR:-}" ]; then
  for f in CRISPRGeneEffect.csv Model.csv; do
    [ -f "$DEST/$f" ] || ln -sfn "$(cd "$DEPMAP_26Q3_DIR" && pwd)/$f" "$DEST/$f"
  done
fi
[ -f "$DEST/CRISPRGeneEffect.csv" ] && [ -f "$DEST/Model.csv" ] || {
  echo "Missing 26Q3 files in $DEST (CRISPRGeneEffect.csv, Model.csv). Set DEPMAP_26Q3_DIR or copy them there." >&2
  exit 1
}
{ sha256sum -c SHA256SUMS.pinned; python3 analyze_26Q3.py > /dev/null; python3 -m json.tool results/comparison_22Q2_vs_26Q3.json; } 2>&1 | tee results/run.log
