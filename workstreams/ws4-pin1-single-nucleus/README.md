# WS4 — Does the CTE PIN1 decrease hold inside cell types? Single-nucleus test (GSE155114, GSE261807)

Background: WS3 (PR #7) found PIN1 lower in CTE bulk cortex (GSE319253: log2FC −0.66; −0.43 after adjusting for a neuronal marker score) and named this test as its falsifier. Not medical advice.

## Method (`analyze.py`; genes in `genes.py`; thresholds fixed in code)
- Per-nucleus counts for 42 pre-specified genes only were streamed out of the GEO matrices (`extract_*.py`): PIN1, the WS3 cell-cycle set (12 genes), the WS3 DNA-repair set (10 genes) and 4 markers per cell type.
- Cell type per nucleus by marker score (mean log-normalised markers → z across nuclei; call = top z ≥ 1 and ≥ 0.5 above runner-up; else "Unassigned"). No clustering, no integration, no batch correction.
- Pseudobulk = summed counts per donor × cell type (≥ 50 nuclei; ≥ 3 donors/group). Library-size offset.
- Per gene: quasi-Poisson GLM (log2 rate ratio, 95% CI, p) and Mann-Whitney on donor log2 CPM. Gene sets: mean per-gene z of log2 CPM, Mann-Whitney. BH across the 5 cell types for PIN1.
- GSE261807 (10x, CellRanger-filtered): CTE vs control primary; RHI (head-impact exposure, no CTE) vs control secondary. GSE155114 (inDrops): CTE vs control. Datasets analysed separately; the meta-estimate in `results.json` is descriptive only.

## Results — PIN1, CTE vs control (log2 rate ratio, quasi-Poisson; MWU on donor CPM)
| Cell type | GSE261807 donors CTE/ctl | log2RR (95% CI) | p qPois / MWU | BH | GSE155114 donors | log2RR (95% CI) | p qPois |
|---|---|---|---|---|---|---|---|
| Excitatory neurons | 11/8 | −0.36 (−0.69, −0.03) | 0.031 / 0.091 | 0.039 | 6/8 | +0.35 (−0.35, 1.05) | 0.32 |
| Inhibitory neurons | 11/8 | −0.46 (−0.73, −0.19) | 0.0008 / 0.051 | 0.002 | 7/8 | −0.03 (−0.69, 0.63) | 0.93 |
| Astrocytes | 11/8 | −0.39 (−0.81, 0.03) | 0.067 / 0.109 | 0.067 | 7/6 | −0.35 (−1.39, 0.69) | 0.52 |
| Microglia | 9/8 | −0.64 (−0.98, −0.30) | 0.0002 / 0.021 | 0.001 | 5/7 | −0.13 (−0.90, 0.64) | 0.74 |
| Oligodendrocytes | 10/8 | −0.54 (−0.92, −0.17) | 0.005 / 0.012 | 0.008 | 5/6 | +0.12 (−0.68, 0.91) | 0.77 |

Other readouts (details in `results/results.json`):
- GSE261807 cell-cycle set: lower in CTE in every cell type with enough genes (7–8 of 9–10 genes down in neurons), MWU p 0.09 (excitatory), 0.09 (inhibitory), 0.0079 (microglia, only 4 genes), 0.002 (oligodendrocytes); not astrocytes (0.90). DNA-repair set: mostly lower, but p ≥ 0.15 except microglia (0.036). Neither set was adjusted for multiplicity.
- RHI vs control (GSE261807): PIN1 also lower but weaker (excitatory −0.25, p 0.12; inhibitory −0.32, p 0.043 / BH 0.11; oligodendrocyte −0.58, BH 0.038). No clear CTE-over-RHI step.
- GSE155114 shows no PIN1 difference in any cell type; its CIs (≈ ±0.7 to ±1 log2) are wide enough to include both zero and the GSE261807 effects.

## Reading
- In GSE261807, PIN1 is lower in CTE in **both neuron classes** (−0.36 and −0.46 log2, same size as the neuron-adjusted bulk estimate of −0.43), so the bulk drop is not explained by neuron loss alone.
- But the drop is not neuron-specific: it appears in all five cell types, the cell-cycle and DNA-repair sets also lean down, and RHI donors lean the same way. A global shift (depth, ambient RNA, batch, tissue quality) fits this as well as a biological PIN1 change.
- GSE155114 does not replicate it; it is too sparse to refute it (median 114–600 nuclei per donor × type, PIN1 total counts 58–478 per cell type).

## What this does NOT show
- No causality or Pin1 *activity*; RNA only. Only CTE vs control; no staging, tau burden, age, sex, PMI or APOE covariates are modelled (not provided in GEO metadata used here). Group differences in age and tissue quality are unadjusted.
- Few donors (≤ 11/8), and nuclei per donor differ between groups (CTE donors have fewer assigned nuclei in GSE261807: e.g. astrocytes median 380 vs 1,108). Quasi-Poisson p-values are more favourable than the rank tests (MWU p 0.05–0.12 for neurons); read the neuron result as suggestive, not established.
- Marker-based typing is crude; 31% (GSE261807) and 55% (GSE155114) of nuclei stayed "Unassigned"; low-depth inDrops typing is noisy.
- **Analytic change disclosed:** my first pass used QC ≥ 1000 UMI and ≥ 400 genes for both datasets. That kept only 10,643 GSE155114 nuclei (the GEO matrices include empty droplets), so GSE155114's cut was lowered to 500 / 300 *after* seeing that pass (39,403 nuclei); the strict version is kept in `results.json` (`..._strictQC_sensitivity`) and gives the same conclusion (underpowered, no signal). GSE261807 settings were not changed.
- 5 cell types × 2 datasets for PIN1; BH is applied within dataset across cell types only. The gene-set p-values are uncorrected.
- Would strengthen / refute: an independent snRNA cohort with age-matched controls and a depth-matched negative-control gene set; an independent bulk cohort with log2FC ≈ 0 would refute the WS3 finding.

## Reproduce
`./run.sh` (≈ 6 GB download, ≈ 25 GB scratch, ~15 GB RAM, ~1 h; inputs checked against `SHA256SUMS.pinned`). `results/run.log` was produced from the same scripts. Data files are not committed. GSE155114 files are processed counts from GEO; no patient identifiers are used.
