\name{fit}
\alias{fit}
\title{Reference Local Projections Difference in Differences}
\usage{
fit(data, y, time, unit = NULL, gvar = NULL, treatment = NULL,
    cluster = NULL, covariates = NULL, event_window = NULL,
    data_mode = "auto", weighting = "reweighted", covariate_lag = 1L,
    inference = "native", assignments = NULL, reps = 99L, seed = 0L,
    null = 0, ...)
}
\arguments{
\item{data}{Panel data.frame.}
\item{y, time, unit}{Outcome, integer calendar time and panel ID column names.}
\item{gvar, treatment}{Cohort column (zero/missing means never), or absorbing binary treatment.}
\item{cluster}{Cluster column, defaults to unit.}
\item{covariates}{Covariate column names.}
\item{event_window}{Required ordered integer endpoints, including post-treatment periods.}
\item{data_mode}{auto or panel; rcs unsupported.}
\item{weighting}{reweighted or variance.}
\item{covariate_lag}{One by default; zero explicitly uses supplied entry-row values.}
\item{inference}{Character vector of native, jackknife, ri_beta, ri_t_native, , or point.}
\item{assignments}{Alternative cohort assignments: rows by sorted cluster columns.}
\item{reps, seed}{Unique randomization draw count and seed.}
\item{null}{Only zero supported.}
\item{...}{Unknown options rejected.}
}
\value{Named list of pooled estimate, target metadata, native inference, horizon results and resampling inference.}
\description{Native regression-based DGJT LP-DiD on genuine absorbing-treatment panels.}
