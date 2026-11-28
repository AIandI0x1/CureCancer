#!/usr/bin/env bash
# One-command rerun: search GEO/Synapse (dated, results differ over time), fetch pinned GSE319253 counts, verify, test H2.
set -euo pipefail
cd "$(dirname "$0")"; mkdir -p data results
python3 search.py
[ -f data/GSE319253_featureCounts_count_matrix_human.txt.gz ] || curl -sSfL -o data/GSE319253_featureCounts_count_matrix_human.txt.gz \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE319nnn/GSE319253/suppl/GSE319253_featureCounts_count_matrix_human.txt.gz
{ sha256sum -c SHA256SUMS.pinned; python3 test_h2.py; } 2>&1 | tee results/run.log
