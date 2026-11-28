#!/usr/bin/env bash
# usage: run_extract_261807.sh <dir with filtered *_matrix.mtx.gz, *_features.tsv.gz, *_barcodes.tsv.gz>
set -euo pipefail; cd "$(dirname "$0")"; mkdir -p data/n261
ls "$1"/*_matrix.mtx.gz | sed 's#.*/##; s#_matrix.mtx.gz##' | xargs -P 2 -I{} sh -c 'o=data/n261/{}.tsv.gz; [ -s $o ] || python3 extract_gse261807.py '"$1"' {} $o'
