# Research Landscape — open cancer data, tooling, and first workstreams

Surveyed 2026-11-27. Every link below was verified live during the survey; each entry states its access model honestly, including where approval gates exist.

Terminal reminder: **this goal completes only when cancer is actually cured.** This file is scaffolding — ranked pointers, not progress.

## Ranked open datasets

Ranked by openness × reproducibility × access cost. 5 = fully open / trivially reproducible.

| # | Source | What it is | Access | Openness | Reproducibility | Cost |
|---|--------|-----------|--------|----------|-----------------|------|
| 1 | [GDC / TCGA](https://portal.gdc.cancer.gov/) | NCI Genomic Data Commons: harmonized TCGA + TARGET + CPTAC genomic/clinical data, documented REST API ([docs](https://docs.gdc.cancer.gov/API/Users_Guide/Search_and_Retrieval/)) | Open-access tiers need no approval; controlled tiers need dbGaP DAC approval | 5 | 5 | Free |
| 2 | [cBioPortal](https://www.cbioportal.org/) | MSK/DFCI portal over TCGA + published studies; visual + programmatic access ([API docs](https://docs.cbioportal.org/web-api-and-clients/), [source](https://github.com/cBioPortal/cbioportal/)) | Fully open; open-source stack | 5 | 5 | Free |
| 3 | [DepMap](https://depmap.org/portal/) | Broad Cancer Dependency Map: CRISPR/RNAi screens + omics across 1000+ cell lines; quarterly releases (26Q3 current); bulk downloads, Bioconductor `depmap` package | Fully open downloads | 5 | 5 | Free |
| 4 | [ICGC ARGO](https://platform.icgc-argo.org/) | International Cancer Genome Consortium: clinical + molecular data, 7,056 donors / 24 programmes ([Release 14](https://www.icgc-argo.org/news/42/release-14)); code at [icgc-argo](https://github.com/icgc-argo) | Open tools; controlled data needs [DACO approval](https://daco.icgc-argo.org/) | 4 | 5 | Free, approval-gated |
| 5 | [SEER](https://seer.cancer.gov/data-software/) | NCI surveillance: US population incidence/survival 1975–2023; [REST API](https://api.seer.cancer.gov/) (free key), [SEER*Stat DBs](https://seer.cancer.gov/data-software/documentation/seerstat/nov2025/) (Nov 2025 submission) | Public API + downloads; limited-use files via request | 4 | 5 | Free |

## Ranked open tooling

| # | Tool | What it is | Access | Openness | Reproducibility | Cost |
|---|------|-----------|--------|----------|-----------------|------|
| 1 | [Bioconductor](https://www.bioconductor.org/) | Versioned R packages for genomics incl. `depmap`, `cBioPortalData`, `TCGAretriever`; release-versioned, Docker images | Fully open | 5 | 5 | Free |
| 2 | [nf-core/sarek](https://nf-co.re/sarek/latest/docs/usage/) + [nf-core/oncoanalyser](https://nf-co.re/oncoanalyser/latest/) | Nextflow pipelines for variant calling (sarek) and Hartwig-style oncogenomic analysis (oncoanalyser); containers, CI-tested, DOIs per release | Fully open | 5 | 5 | Free (compute only) |
| 3 | [Galaxy](https://usegalaxy.org/) (+ [Cancer hub](https://cancer.usegalaxy.org/)) | Web + API workflow execution; published histories re-run end-to-end; ToolShed versioned | Free public server (quotas apply); self-hostable | 5 | 4 | Free |

## Candidate first workstreams

Each has a definition of done — a stranger must be able to re-run it with one command.

### WS1 — One-command TCGA oncoprint reproducer (recommended first)

- **Question:** can a stranger reproduce a published TCGA alteration summary for one cancer type with one command?
- **Data:** GDC open-access MAF snapshot (pinned release + checksum) or cBioPortal bulk download.
- **Method:** pinned script (Bioconductor `maftools` or Python) emitting an oncoprint + alteration frequencies.
- **Done when:** `reproduce.sh` in the repo re-creates the figure byte-identical from the pinned snapshot; run log committed beside it.
- **Why first:** smallest credible reproducibility win; exercises GDC + cBioPortal open tiers; no approval needed.

### WS2 — DepMap selective-dependency snapshot

- **Question:** which genes are selectively essential in one cancer lineage (e.g. KRAS-mutant lung)?
- **Data:** pinned DepMap quarterly release (26Q3 at survey time) + lineage metadata.
- **Method:** pinned notebook/script computing top-10 selective dependencies by lineage enrichment.
- **Done when:** one command re-runs the ranking from the pinned release; run log + ranked table committed.
- **Why:** directly therapeutic-adjacent (what kills which cancer lines); quarterly releases give a natural refresh cadence.

### WS3 — SEER survival-trend tracker

- **Question:** are 5-year survival trends improving for 5 major cancers over 20 years?
- **Data:** SEER API pull or SEER*Stat export (Nov 2025 submission), pinned snapshot.
- **Method:** pinned notebook producing trend curves + annual percent change.
- **Done when:** versioned snapshot + notebook re-run with one command; figure + table committed.
- **Why:** population-level signal (are we winning?); free API; smallest data footprint.

## Pointers only

No patient-level data is ever committed to this repo. Snapshots are aggregated/open tiers; anything controlled stays behind its DAC with a pointer + access instructions.
