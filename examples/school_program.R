# Run from the repository's top-level folder:
# Rscript examples/school_program.R
# These are fictional teaching data, not evidence about a real program.

source("R/summclust.R")
d <- read.csv("examples/school_program.csv")

fit <- summclust(
  d,
  outcome = "test_score",
  focal = "program",
  controls = "baseline_score",
  cluster = "school",
  fevar = "year"
)

options(digits = 5)
print(fit)

diagnostics <- data.frame(
  school = names(fit$ng),
  N = as.numeric(fit$ng),
  partial_leverage = as.numeric(fit$partlev),
  estimate_without_school = as.numeric(fit$betanog)
)
diagnostics$change <- diagnostics$estimate_without_school - fit$beta_full
top <- diagnostics[order(diagnostics$partial_leverage, decreasing = TRUE), ]
print(head(top, 5), row.names = FALSE, digits = 4)

# Save a table of all schools if you want to inspect it in a spreadsheet:
# write.csv(diagnostics, "school_diagnostics.csv", row.names = FALSE)
