# Earlier development notes

These are the previous landing-page notes, preserved on 19 November 2026 for
provenance. They may describe superseded restrictions. Use the [current README](../README.md),
[inference guide](../INFERENCE.md), and [change history](../CHANGELOG.md) for current behavior.

# fastlpdidR

The standalone jackknife refits the pooled effect on surviving cohort and calendar-time support. Every input inference cluster contributes a deletion. Output reports inactive clusters and the native estimation-sample cluster count separately.

Native R reference DGJT LP-DiD. MIT original code, base R only; no Python/Stata
runtime dependency. Tested development subset, not full-feature certification.

```sh
R CMD INSTALL .
Rscript tests/reference.R
```

```r
library(fastlpdidR)
d <- read.csv("tests/panel.csv")
r <- fit(d, "y", "t", unit="i", gvar="g", event_window=c(-2,2),
         inference=c("native","jackknife"))
print(c(estimate=r$estimate, se=r$native_se))
```