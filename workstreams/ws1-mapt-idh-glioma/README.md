# WS1 — Tau (MAPT) expression vs IDH status and survival in TCGA diffuse glioma

Issue: #2 · Background: `research/notes/cte-and-cancer.md` (tau theme) · Replication target: Gargini et al. 2020, Sci Transl Med, PMID 31969485.
Not medical advice. Open-tier data only.

## Data (pinned)
- Source: public cBioPortal API, TCGA PanCancer Atlas studies `lgg_tcga_pan_can_atlas_2018` (import 2026-10-30) and `gbm_tcga_pan_can_atlas_2018` (import 2026-10-29).
- `fetch.py` builds one derived table (primary tumours: MAPT RSEM, IDH1 R132 / IDH2 R140/R172 hotspot status among sequenced samples, OS, age). Its sha256 is pinned in `SHA256SUMS.pinned`; `run.sh` refuses to run if the snapshot differs.
- The table is not committed (patient-level rows, even though de-identified and public); it is rebuilt from the API.

## Reproduce
```
pip install -r requirements.txt
./run.sh        # fetch if missing → verify sha256 → analyze → results/run.log
```

## Results (snapshot sha256 0d1abda6…063a; 667 primary tumours: LGG 415 IDH-mut / 99 IDH-wt; GBM 7 / 145; 1 GBM not sequenced)
| Analysis | Result |
|---|---|
| MAPT, IDH-mut vs IDH-wt, all gliomas | median log2 13.54 vs 12.01 (Δ 1.54); AUC 0.933; Mann-Whitney p 1.8e-77 |
| Same, LGG only | 13.55 vs 12.46 (Δ 1.10); AUC 0.865; p 1.7e-29 |
| Cox OS, MAPT per SD, unadjusted (n 662, 245 deaths) | HR 0.54 (0.50–0.58), p 8.9e-58 |
| Cox, adjusted for IDH, age, GBM vs LGG (n 615, 199 deaths) | HR 0.83 (0.67–1.03), p 0.096 |
| Cox within IDH-mut, adjusted for age (n 417, 75 deaths) | HR 0.47 (0.26–0.85), p 0.012 |
| Cox within IDH-wt, adjusted for age + cohort (n 198, 124 deaths) | HR 0.95 (0.75–1.21), p 0.68 |
| Log-rank, IDH-mut, MAPT median split | p 0.052 |
| MAPT vs age (Spearman, n 616) | rho −0.40, p 1.1e-24 |

## Reading
- Replicates the direction reported by Gargini 2020: MAPT is much higher in IDH-mutant gliomas.
- Most of MAPT's crude survival association is explained by IDH status, age and tumour grade/cohort (HR moves from 0.54 to 0.83, CI crosses 1).
- Within IDH-mutant tumours, higher MAPT is associated with longer survival after age adjustment (HR 0.47 per SD) — but only 75 deaths, and the median-split log-rank is borderline (p 0.052).
- No association within IDH-wild-type tumours.

## What this does NOT show
- Association, not mechanism; bulk RNA cannot say which cells express MAPT (neurons infiltrated by tumour can contribute).
- No adjustment for 1p/19q codeletion, MGMT, treatment, or tumour purity. Per-SD HRs use the SD of the whole cohort.
- EGFR status (part of Gargini's triad) not analysed.
- Would falsify the IDH-mut finding: an independent cohort (e.g., CGGA) showing HR ≈ 1 after adjustment for 1p/19q and grade.
