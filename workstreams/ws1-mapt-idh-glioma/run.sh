#!/usr/bin/env bash
# One-command rerun: build snapshot if missing, verify checksum, run analysis, log.
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p results
[ -f data/glioma_mapt_idh_os.tsv ] || python3 fetch.py
{ sha256sum -c SHA256SUMS.pinned; python3 analyze.py; } 2>&1 | tee results/run.log
