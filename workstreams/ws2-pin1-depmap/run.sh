#!/usr/bin/env bash
# One-command rerun: fetch pinned files if missing, verify checksums, analyze, log.
set -euo pipefail
cd "$(dirname "$0")"; mkdir -p results
./fetch.sh
{ sha256sum -c SHA256SUMS.pinned; python3 analyze.py; } 2>&1 | tee results/run.log
