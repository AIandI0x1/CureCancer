#!/usr/bin/env bash
# usage: run_extract_155114.sh <dir with GSE155114 *.counts.tsv> ; ~45 s per file, 6 in parallel
set -euo pipefail; cd "$(dirname "$0")"; mkdir -p data/n155
ls "$1"/*.counts.tsv | xargs -P 6 -I{} sh -c 'o=data/n155/$(basename {} .counts.tsv).tsv.gz; [ -s $o ] || python3 extract_gse155114.py {} $o'
