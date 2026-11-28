#!/usr/bin/env bash
# One-command reproduction. Needs ~25 GB scratch, ~15 GB RAM, ~1 h (download 6 GB + 20 GB unpacked, streaming extraction).
# Large data stay in data/ (git-ignored). Final output: results/results.json + results/run.log
set -euo pipefail
cd "$(dirname "$0")"; mkdir -p data/raw155114 data/raw261807 data/c155114 results
G=https://ftp.ncbi.nlm.nih.gov/geo/series
# GSE155114: 8 CTE vs 8 control donors, inDrops dense counts
[ -f data/raw155114/GSE155114_RAW.tar ] || curl -sSfL -o data/raw155114/GSE155114_RAW.tar $G/GSE155nnn/GSE155114/suppl/GSE155114_RAW.tar
(cd data/raw155114 && tar -xf GSE155114_RAW.tar && for f in GSM*.tar.gz; do tar -xzf "$f" -C ../c155114; done)
# GSE261807: stream the RAW tar, keep only the CellRanger *filtered* samples (skip *_raw_* matrices)
[ -f data/raw261807/.done ] || { curl -sSfL $G/GSE261nnn/GSE261807/suppl/GSE261807_RAW.tar | tar -x -C data/raw261807 --exclude='*_raw_*'; touch data/raw261807/.done; }
# verify pinned inputs (paths in SHA256SUMS.pinned are relative to data/)
(cd data && sha256sum -c ../SHA256SUMS.pinned --quiet)
./run_extract_155114.sh data/c155114
./run_extract_261807.sh data/raw261807
{ python3 analyze.py > /dev/null; python3 summarize.py; } 2>&1 | tee results/run.log
