# Block 4, step 11: operating range in n. p = q = 4, clean data, n in {8, 10, 12, 16, 24, 40}, 20 reps, seed 4100 + rep.
# degen = fraction of default fits ending with min eigenvalue ratio < 1e-6 or F < F_0 - 50 n (indicator of cm_step5),
# compared with the MLE's eigenvalue ratio (a legitimate small-n fit also has small eigenvalues). Median KL of both.
# Output: output/cm_step11_smalln.log
source("R/load.R"); library(parallel)
p <- 4; q <- 4; S <- true_sigma(p, q, "full"); ev <- function(A) { e <- eigen(A, TRUE, TRUE)$values; min(e) / max(e) }
one <- function(seed, n) { set.seed(seed); X <- rmatnorm(n, Sr = S$Sr, Sc = S$Sc)
  f <- tryCatch(cellmmcd(X), error = function(e) NULL); m <- flipflop(X, center = TRUE)
  if (is.null(f)) return(c(err = 1, degen = NA, kl = NA, evr = NA, evm = ev(m$Sr) * ev(m$Sc), klm = kl_mn(m$Sr, m$Sc, S$Sr, S$Sc), tr = NA))
  c(err = 0, degen = as.numeric(min(ev(f$Sr), ev(f$Sc)) < 1e-6 || f$raw$obj[length(f$raw$obj)] < f$raw$obj[1] - 50 * n),
    kl = kl_mn(f$Sr, f$Sc, S$Sr, S$Sc), evr = ev(f$Sr) * ev(f$Sc), evm = ev(m$Sr) * ev(m$Sc), klm = kl_mn(m$Sr, m$Sc, S$Sr, S$Sc), tr = f$trimmed) }
sink("output/cm_step11_smalln.log"); cat("p = q = 4, clean, 20 reps. err / degen fractions; median KL cellMMCD vs MLE; median eigenvalue-ratio product (cellMMCD, MLE); trimmed.\n")
for (n in c(8, 10, 12, 16, 24, 40)) { r <- do.call(rbind, mclapply(1:20, function(s) one(4100 + s, n), mc.cores = 4))
  cat(sprintf("n = %2d: err %.2f degen %.2f | KL %.3f vs MLE %.3f | evr %.1e vs MLE %.1e | trimmed %.2f\n", n, mean(r[, "err"]), mean(r[, "degen"], na.rm = TRUE),
              median(r[, "kl"], na.rm = TRUE), median(r[, "klm"]), median(r[, "evr"], na.rm = TRUE), median(r[, "evm"]), median(r[, "tr"], na.rm = TRUE))) }
sink(); cat("done\n")
