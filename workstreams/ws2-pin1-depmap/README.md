# WS2 — PIN1 CRISPR dependency by lineage (DepMap)

Issue: #3 · Background: `research/notes/cte-and-cancer.md` (Pin1 theme). Not medical advice.

## Data (pinned)
- Original analysis: **DepMap 22Q2** CRISPR gene effect (17,386 genes × 1,086 cell lines) and cell-line metadata, as redistributed by the Bioconductor `depmap` package on ExperimentHub (EH7554 `crispr_22Q2.rda`, EH7558 `metadata_22Q2.rda`). sha256 pinned in `SHA256SUMS.pinned`.
- **26Q3 addendum** (same metrics, side-by-side comparison): see `addendum_26Q3/`. The portal files were obtained 2026-11-28 after human verification; they are not committed. First-pass attempts from this machine hit the verification page and a figshare HTTP 403, which is why 22Q2 shipped first.

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

## 26Q3 addendum (short)
In 26Q3 (1,361 models) PIN1 median is −0.007 (42.9th → 58.3rd percentile, slightly *less* essential than in 22Q2). 14 lines are below −0.5 (0 in 22Q2); 10 of those 14 are organoids (post-hoc split). Kruskal-Wallis across 24 Oncotree lineages: H 44.23, p 0.0049; strongest lineage median is Thyroid −0.069. Full table, the 14 lines, and caveats: `addendum_26Q3/README.md`.

## Reading
- In the 22Q2 CRISPR knockout screen, cancer cell lines do not depend on PIN1 for growth in any lineage (median −0.046, 0 lines < −0.5). 26Q3 keeps the same median-level conclusion but is no longer a literal zero tail — see the addendum. Literature reports PIN1 over-expression in many cancers (Driver 2015, PMID 25583562); over-expression is not dependency.
- This weakens any reading of the CTE–cancer Pin1 link that requires Pin1 to be a tumour-cell survival gene. It does not test Pin1's role in tumour initiation, in vivo growth, drug response, or the tumour microenvironment.

## What this does NOT show
- 2D cell-line proliferation over ~3 weeks; no in vivo, no 3D, no stress conditions.
- The 22Q2 section above is the original release. A 26Q3 rerun is in `addendum_26Q3/` (PIN1 median −0.007, rank 58.3rd, 14/1,361 lines < −0.5, mostly organoids; no lineage median near −0.5). Scoring changed across releases; read both.
- Would more strongly falsify "no dependency": a lineage whose median is ≤ −0.5, or an orthogonal screen (RNAi, degron) that replicates the same tail.
