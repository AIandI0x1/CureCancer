# WS2 — PIN1 CRISPR dependency by lineage (DepMap)

Issue: #3 · Background: `research/notes/cte-and-cancer.md` (Pin1 theme). Not medical advice.

## Data (pinned) — and why not 26Q3
- The issue asked for the newest public release (26Q3). From this analysis machine, the DepMap portal download API sits behind a human-verification page and the figshare API returned HTTP 403, so no 26Q3 file could be fetched without bypassing a bot check (not done).
- Used instead: **DepMap 22Q2** CRISPR gene effect (17,386 genes × 1,086 cell lines) and cell-line metadata, as redistributed by the Bioconductor `depmap` package on ExperimentHub (EH7554 `crispr_22Q2.rda`, EH7558 `metadata_22Q2.rda`). sha256 pinned in `SHA256SUMS.pinned`.
- Rerunning on 26Q3 is the obvious next step once a sanctioned download path works.

## Reproduce
```
pip install -r requirements.txt
./run.sh        # download (~139 MB) → verify sha256 → analyze (~8 GB RAM) → results/run.log
```

## Results (22Q2)
| Gene | Lines | Median effect | IQR | Lines < −0.5 | Genome-wide rank of median (most essential first) |
|---|---|---|---|---|---|
| PIN1 | 1,086 | −0.046 | −0.108 to 0.008 | 0 (0.0%) | 42.9th percentile |
| PIN4 | 1,081 | −0.407 | −0.508 to −0.315 | 26.6% | 10.9th |
| CDK1 | 1,086 | −2.042 | — | 99.9% | 0.7th |
| RPL3 | 1,086 | −2.044 | — | 100% | 0.7th |
| PCNA | 1,086 | −2.414 | — | 100% | 0.1th |

- Strongest single-line PIN1 effect in the release: −0.41 (see `results/results.json`).
- By lineage (24 lineages with ≥10 lines): medians range −0.082 (breast) to −0.001 (eye); central nervous system −0.055 (n 79). No line in any lineage reaches −0.5. Kruskal-Wallis H 46.4, p 0.0026 — statistically detectable, biologically negligible differences.
- PPIA and POLR2A are absent from this release's gene set.

## Reading
- In CRISPR knockout screens, cancer cell lines do not depend on PIN1 for growth, in any lineage. Literature reports PIN1 over-expression in many cancers (Driver 2015, PMID 25583562); over-expression is not dependency.
- This weakens any reading of the CTE–cancer Pin1 link that requires Pin1 to be a tumour-cell survival gene. It does not test Pin1's role in tumour initiation, in vivo growth, drug response, or the tumour microenvironment.

## What this does NOT show
- 2D cell-line proliferation over ~3 weeks; no in vivo, no 3D, no stress conditions.
- 22Q2, not the newest release; scoring changed across releases.
- Would falsify "no dependency": a newer release or orthogonal screen (RNAi, degron) showing PIN1 effects ≤ −0.5 in a defined lineage.
