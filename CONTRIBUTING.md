# Contributing to CureCancer

## Propose work

1. Open an **issue** describing the workstream: question, data source, method sketch, definition of done.
2. Wait for maintainer sign-off (or self-assign small docs fixes directly).
3. One workstream = one PR. One PR = one reviewable unit.

## Land work

Every PR must include:

- **Data**: pinned snapshot, checksum, and source link (or pointer if the data can't live here).
- **Code**: pinned dependencies (`requirements.txt` / lockfile / container digest).
- **Reproduction**: one command that re-runs the result, plus the run log committed beside it.
- **Honesty**: what was NOT shown, what could confound it, what would falsify it.

## Standards

- No medical advice, no patient data, no PII — ever.
- No paywalled-only sources without an open alternative noted.
- Cite primary sources. Link DOIs, not headlines.
- Small, reviewable diffs. Rebase onto `main` before requesting review.

## Review

A maintainer checks: reproducibility (does the command run?), evidence (data + log present?), honesty (limits stated?), scope (one unit, not a bundle?). Two approvals for anything touching methods; one for docs.

Thank you for building in the open.
