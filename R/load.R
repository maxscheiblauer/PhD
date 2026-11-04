# Source all project functions. Run scripts from the project root.
for (f in c("R/matnorm.R", "R/tff.R", "R/metrics.R", "R/data.R", "R/cellmmcd.R", "R/archive/r1.R")) source(f)
dir.create("output", showWarnings = FALSE)
has_robustmatrix <- requireNamespace("robustmatrix", quietly = TRUE)
