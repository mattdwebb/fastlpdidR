# Local-projection DiD in R

Estimate policy effects over a chosen horizon using the Dube–Girardi–Jordà–Taylor
local-projection DiD approach. The package compares outcome changes for newly
treated units with units that remain untreated through the relevant horizon.
It reports horizon effects and a pooled effect, with native clustered and
delete-cluster jackknife inference.

## Install

Use R 4.1 or later. Estimation uses base R; it does not call Python or Stata.
Install this source version with:

```r
install.packages("remotes")  # once, if needed
remotes::install_github("mattdwebb/fastlpdidR")
```

Alternatively, download the repository, open a terminal in its directory, and run
`R CMD INSTALL .`. For reproducibility, record the installed version
with `packageVersion("fastlpdidR")` and the source commit.

## A school-policy example

The included `school_panel.csv` contains school-year means for the same 12
schools over 2010–2017: 96 observations. Four schools adopt in 2013, four in 2015,
and four never adopt. `test_score` is the mean score and `baseline_score` is a
score covariate. The panel unit and inference cluster are both the school.

These are fictional teaching data. The generator adds a two-point policy effect
and varying disturbances; the example is not evidence about a real education
policy. [Data definitions and the generator](inst/extdata/DATA_GUIDE.md) are
included, and the same input is used in the Python and Stata examples.

```r
library(fastlpdidR)
school_data <- read.csv(system.file("extdata", "school_panel.csv",
                                    package = "fastlpdidR"))
fit <- fit(school_data, "test_score", "year", unit = "school",
           gvar = "first_treat", cluster = "school", event_window = c(0, 1),
           covariates = "baseline_score", inference = c("native", "jackknife"))
example_result <- data.frame(estimate = fit$estimate, cluster_se = fit$native_se,
                             jackknife_se = fit$inference$jackknife$se,
                             jackknife_p = fit$inference$jackknife$p)
print(example_result, digits = 7, row.names = FALSE)
```

Actual output, rounded to six decimal places:

| Estimate | Clustered SE | Jackknife SE | Jackknife p-value |
| --- | --- | --- | --- |
| 1.601279 | 0.331234 | 0.365699 | 0.001102 |

The estimated adjusted effect is **1.60 score points**. Its jackknife
standard error is **0.37 points**. The p-value tests a zero effect using a
Student t distribution with 11 degrees of freedom, one less than the 12 schools.
The 96 school-year means do not create 96 independent policy assignments.
The two standard errors describe uncertainty about the same reported estimate.
The clustered calculation uses this estimator's native variance formula. The
jackknife recomputes the estimate after deleting each school, including nuisance
models and aggregation weights.

## Adapt the example

Omit `covariates` to estimate without additional controls. To change the
covariate specification, pass the names of numeric columns, for example
`covariates = c("baseline_score", "class_size")` after adding `class_size` to
your own data. Encode categorical controls before passing them. Covariates should
be chosen to suit the identification assumptions; a variable affected by the
policy is not made appropriate by including it in a regression.

`first_treat` repeats each school's first policy year on every row. Zero means
never treated. Treatment must remain on after adoption. Keep all necessary
pre-policy observations and check that the required control comparisons exist.
Missing or unsupported inputs produce explicit errors instead of silently
changing the specification.

Save the displayed results with `write.csv(example_result, "school_results.csv",
row.names = FALSE)`. The full `fit` object retains the estimate, inference results,
and the estimator-specific diagnostics described in `?fit`.

## Inference and supported scope

LP-DiD requires a genuine panel: the same unit must be observed over time.
The example uses schools as panel units; it never links unrelated pupils across
years. `c(0, 1)` pools the adoption year and the following year. Reweighted
aggregation is the default; `weighting = "variance"` selects variance weighting.
Covariates default to their lag-one values. This implementation does not support
repeated cross sections, observation weights, nonabsorbing treatment, or PMD.

The jackknife uses all input inference clusters. If deleting a school removes
one adoption cohort but leaves an estimable aggregate, the supported effects and
their weights are recomputed. If a required deletion has no estimable aggregate,
the jackknife fails explicitly; failed deletions are never discarded.

## Help, validation and citation

Method background: Dube, Girardi, Jordà and Taylor, [A Local Projections Approach to Difference-in-Differences](https://doi.org/10.1002/jae.70000).

Use `help("fit", package = "fastlpdidR")` for arguments and return values.
The [method specification](METHOD_SPECIFICATION.md), [validation record](VALIDATION.md),
and [developer reference](docs/DEVELOPER_REFERENCE.md) give additional detail.
These are supported subsets, not promises of every feature in the original
estimators' software. Cite the relevant method and use `citation("fastlpdidR")`
for the software version you used. Original code is distributed under the
[package license](LICENSE); third-party notices retain their own terms.

From a source checkout, run `R CMD INSTALL .` followed by
`Rscript tests/readme.R` to reproduce and check this example. The executable
source is [school_policy.R](inst/examples/school_policy.R).

Current candidate checks: [release validation](RELEASE_VALIDATION.md).
