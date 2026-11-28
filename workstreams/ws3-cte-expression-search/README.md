# WS3 — Public CTE brain expression data search, and a first H2 test

Issue: #4 · Background: hypothesis H2 in `research/notes/cte-and-cancer.md`. Not medical advice.

## Search (`search.py`, dated output in `results/hits.json`)
NCBI GEO (E-utilities, db=gds) and Synapse full-text search, run 2026-11-28.
- GEO, `"chronic traumatic encephalopathy"`: 8 records. Human brain expression series found:
  - **GSE319253** — bulk RNA-seq, superior frontal cortex, 6 neuropathologically confirmed CTE vs 6 age-matched controls, 2 technical replicates each, featureCounts matrix public (released 2026-02-11).
  - GSE155114 — single-nucleus RNA-seq, CTE vs matched controls (PMID 34019156); raw files only.
  - GSE261807 — young athletes: controls, head-impact-exposed without CTE, and CTE (PMID 40963024); raw files only.
  - GSE319802 — PBMC (blood), not brain.
- Synapse keyword search returned hundreds to thousands of loosely matching entities (695 / 5,339 / 10,000 for three queries); none in the first 50 hits per query was a CTE brain expression dataset. Synapse was not searched further.
- **Outcome: a usable open dataset exists (GSE319253), so H2 was tested rather than recorded as blocked.**

## H2 test (`test_h2.py`, gene sets fixed in `genesets.py` before group results were computed)
Prediction (inverse-comorbidity direction): PIN1, cell-cycle and DNA-repair expression lower in CTE.
Method: log2 CPM, technical replicates averaged per donor (n 6 vs 6), Welch t and Mann-Whitney; set score = mean per-gene z.

| Readout | CTE vs control | p |
|---|---|---|
| PIN1 (primary) | log2FC −0.66 | Welch 0.0027; MWU 0.0087 |
| PIN1 adjusted for neuronal marker score | −0.43 (95% CI −0.70 to −0.16) | 0.0062 |
| Cell-cycle set (12 genes) | mean z −0.14 vs +0.14; 9/12 genes lower | MWU 0.24 |
| DNA-repair set (10 genes) | mean z −0.10 vs +0.10; 6/10 lower | MWU 0.70 |
| Composition (descriptive) | neuronal score lower (MWU 0.065); astrocyte AQP4 up (log2FC +1.15) | — |

## Reading
- PIN1 is lower in CTE cortex, in the direction H2 predicts, and the difference survives a crude adjustment for neuronal content.
- Cell-cycle and DNA-repair sets do not move detectably. H2 is supported only for PIN1, not for the broader programme.

## What this does NOT show
- 6 vs 6 donors from one study; no replication yet. CPM, not TMM/DESeq2 normalisation. 30 genes reported without multiplicity correction (PIN1 was the single pre-specified primary readout).
- Bulk tissue: cell-type composition differs (astrocyte marker up, neuronal markers down); the adjustment uses three marker genes only.
- Expression is not Pin1 activity or cis/trans P-tau.
- Would falsify the PIN1 finding: snRNA-seq (GSE155114 / GSE261807) showing no within-neuron PIN1 decrease in CTE, or an independent bulk cohort with log2FC ≈ 0.
