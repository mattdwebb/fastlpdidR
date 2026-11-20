# Release candidate — 20 November 2026

- Add a complete school-policy walkthrough with bundled data, executable code,
  verified numerical output and interpretation in outcome units.
- Check the installed-package README example against its documented output.
- Preserve the existing local surviving-aggregate jackknife corrections while
  integrating the current upstream changes; all required deletions remain required.
- Run Windows build/check steps under an explicit shell so PowerShell does not
  interpret R as its history alias.

# Changelog

## 0.1.0 Development (2026-09-29)

Native absorbing-panel pooled and horizon LP-DiD in Python/R/Stata.
Explicit variance/equal-treated weighting, covariate regression adjustment,
cluster inference, full-recomputation JK and RI, deterministic supplied assignments.
Official lpdid 1.1.0 comparisons and independent lm regression tests.
Unsupported features fail explicitly. This is a validated reference subset,
not a full-feature parity or full-certification release.
