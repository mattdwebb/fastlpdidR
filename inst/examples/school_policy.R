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
