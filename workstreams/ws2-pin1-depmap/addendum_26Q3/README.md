# Addendum — PIN1 CRISPR dependency, DepMap 26Q3 vs 22Q2

This sits next to the original 22Q2 analysis (issue #3, PR #6). Same metrics, same reference genes, same Kruskal-Wallis rule (lineages with ≥10 lines). Not medical advice.

## Data
- **26Q3** portal files `CRISPRGeneEffect.csv` (1,361 models × 18,434 genes) and `Model.csv`, downloaded 2026-11-28 after human verification of the DepMap portal. sha256 in `SHA256SUMS.pinned`. Not committed (large; also still behind a verification wall for a fresh fetch).
- Lineage labels are `Model.OncotreeLineage`. These are not the same strings as the 22Q2 `lineage` field, so the Kruskal-Wallis tests are parallel, not paired.
- Chronos scores are not on an identical scale across releases (common essential genes are more negative in 26Q3). Absolute values should be compared with that in mind.

## Reproduce
```
# after placing the two portal files in ../data/26Q3/, or:
DEPMAP_26Q3_DIR=/path/to/folder ./run.sh
```
Needs pandas, numpy, scipy (same pins as the parent workstream). ~2 GB RAM. Writes `results/results.json`, `results/comparison_22Q2_vs_26Q3.json`, `results/pin1_by_lineage.csv`, `results/run.log`.

## Side-by-side (PIN1 and the same reference genes)

| Gene | 22Q2 median (n) | 22Q2 < −0.5 | 22Q2 rank | 26Q3 median (n) | 26Q3 < −0.5 | 26Q3 min | 26Q3 rank |
|---|---|---|---|---|---|---|---|
| PIN1 | −0.046 (1,086) | 0 (0.0%) | 42.9th | −0.007 (1,361) | 14 (1.0%) | −1.312 | 58.3th |
| PIN4 | −0.407 (1,081) | 26.6% | 10.9th | −0.214 (1,361) | 5.9% | −1.103 | 16.4th |
| PPIA | absent | — | — | −0.707 (1,243) | 80.9% | −1.769 | 6.9th |
| CDK1 | −2.042 (1,086) | 99.9% | 0.7th | −2.625 (1,361) | 100% | −4.354 | 0.26th |
| RPL3 | −2.044 (1,086) | 100% | 0.7th | −2.568 (1,361) | 98.6% | −3.915 | 0.33rd |
| PCNA | −2.414 (1,086) | 100% | 0.1th | −3.025 (1,361) | 99.6% | −4.263 | 0.08th |
| POLR2A | absent | — | — | absent | — | — | — |

Rank = percentile of that gene's median among all genes, most essential first. 22Q2 PIN1 min (−0.41) is taken from the original top-10 list (the overall min was not stored separately).

- Lineage Kruskal-Wallis (n ≥ 10): 22Q2, 24 lineages, H 46.44, p 0.00264; 26Q3, 24 lineages, H 44.23, p 0.00494. In 26Q3, lineage medians run from Thyroid −0.069 (n 11) to Eye +0.061 (n 15). CNS/Brain is −0.035 (n 133); no lineage median is close to −0.5.
- The 14 models with PIN1 < −0.5 are listed in `results/results.json`. They are not concentrated in one lineage (Bowel 4, Esophagus/Stomach 3, Bone 2, Breast 2, then one each in Kidney, Pancreas, Prostate).
- **Post-hoc (not pre-specified):** after seeing those 14, we split by `ModelType`. 10 are organoids (10/108, 9.3%; median −0.085; strongest −1.312, two colorectal organoids). 4 are classic cell lines (4/1,253, 0.3%; LoVo, HA1E, OS457, TC-106; strongest −0.653). HA1E is labelled Non-Cancerous.

## Reading
- The 22Q2 headline still holds for the typical line: PIN1's median effect is near zero, and in 26Q3 it is even slightly closer to zero (rank 42.9th → 58.3rd). No lineage looks dependent.
- The original 22Q2 claim of *zero* lines below −0.5 does not hold in 26Q3: 14 models do, including one organoid below −1. That is a real tail, mostly in organoids, not a lineage-wide dependency. Scoring also changed (CDK1/PCNA/RPL3 medians are more extreme in 26Q3), so a few cell lines crossing −0.5 can be a mix of new models, scale shift, and biology.
- This still does not make PIN1 a general tumour-cell survival gene. It also does not test initiation, in vivo growth, drug response, or the microenvironment.

## What this does NOT show
- 2D (or organoid) proliferation screens only.
- No covariate adjustment; organoid vs cell-line is a post-hoc split.
- Would more strongly falsify "no dependency": a lineage whose *median* is ≤ −0.5, or an orthogonal screen (RNAi, degron) with a similar tail that replicates in the same models.
