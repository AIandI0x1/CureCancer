#!/usr/bin/env bash
# Download the pinned DepMap 22Q2 files from the Bioconductor ExperimentHub mirror (depmap package, EH7554/EH7558).
set -euo pipefail
cd "$(dirname "$0")"; mkdir -p data
B=https://mghp.osn.xsede.org/bir190004-bucket01/ExperimentHub/depmap
for f in crispr_22Q2.rda metadata_22Q2.rda; do [ -f data/$f ] || curl -sSfL -o data/$f $B/$f; done
