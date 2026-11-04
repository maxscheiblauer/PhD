# Source all project functions. Run scripts from the project root.
for (f in c("R/matnorm.R", "R/tff.R", "R/metrics.R", "R/data.R", "R/cellmmcd.R", "R/archive/r1.R")) source(f)
dir.create("output", showWarnings = FALSE)
has_robustmatrix <- requireNamespace("robustmatrix", quietly = TRUE)
# C++ cell step (block 4): compiled once per machine (cache in output/cpp_cache); falls back to the R loop.
has_cpp <- tryCatch({ dir.create("output/cpp_cache", showWarnings = FALSE)
  Rcpp::sourceCpp("src/cm_wstep.cpp", cacheDir = "output/cpp_cache", showOutput = FALSE); TRUE },
  error = function(e) { message("cm_wstep.cpp not compiled (", conditionMessage(e), "); using the R cell step"); FALSE })
options(cellmmcd.engine = if (has_cpp) "cpp" else "R")
