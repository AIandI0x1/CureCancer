# CureCancer Roadmap

The single terminal criterion: **cancer is actually cured** — curative therapy exists and is deployed, confirmed by peer-reviewed clinical evidence and a public verdict. Everything below is scaffolding toward that.

## Phase 1 — Public repo bootstrap ✅

- Repo `AIandI0x1/CureCancer` created, public, default branch `main`.
- Evidence: creation metadata recorded by the bootstrapping system.

## Phase 2 — Mission docs + structure (in progress)

- [x] README (mission, principles, how to contribute)
- [x] ROADMAP (this file)
- [x] CONTRIBUTING (contribution workflow)
- [x] LICENSE (MIT)
- [x] .gitignore (data/science/python hygiene)
- [ ] Issue templates (workstream proposal, bug report)

## Phase 3 — Research landscape survey

Survey open cancer datasets (TCGA, ICGC, CCLE/DepMap, cBioPortal...), open publications, and reproducible pipelines. Output: `LANDSCAPE.md` ranking sources by openness + reproducibility, and 3 candidate first workstreams.

## Phase 4 — First reproducible workstream

Pick the top-ranked candidate and ship a minimal reproducible pipeline: pinned data snapshot, pinned dependencies, one-command run, results committed with run log.

## Phase 5+ — Sustain

Keep the repo alive: review contributions, refresh the landscape quarterly, open new workstreams. The goal stays open until the terminal criterion is met.

## Non-goals

- Medical advice. Nothing here tells any patient what to do — research infrastructure only.
- Proprietary data hosting. Open data or pointers; never leaked PII.
- Claiming progress that isn't evidenced. See Principles in README.
