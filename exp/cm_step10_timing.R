# Block 4, step 10 (rerun in block 6 with the final default, plus the reweighted hard quota; the block-4 log is cm_step10_timing_block4.log): run time against dimension. p = q in {5, 10, 15, 20}, n = 100, 5 % cells at 6 sd, seed 40 + p.
# Default cellmmcd (raw .99 + shadow + reweighting at .999); the raw fit's iterations; flip-flop MLE and
# robustmatrix::mmcd for reference. Output: output/cm_step10_timing.log (one line per size, written as it completes).
source("R/load.R"); suppressMessages(library(robustmatrix))
sink("output/cm_step10_timing.log"); cat("n = 100, 5 % cells at 6 sd. Seconds (elapsed, 1 core) / KL to truth / raw iterations.\n"); sink()
for (p in c(5, 10, 15, 20)) { q <- p; n <- 100; set.seed(40 + p); S <- true_sigma(p, q, "full"); sdc <- sqrt(outer(diag(S$Sr), diag(S$Sc)))
  X <- rmatnorm(n, Sr = S$Sr, Sc = S$Sc); bad <- array(runif(p * q * n) < 0.05, c(p, q, n))
  X[bad] <- X[bad] + 6 * sdc[array(1:(p * q), c(p, q, n))[bad]] * sign(rnorm(sum(bad)))
  t0 <- proc.time()[3]; f <- cellmmcd(X); tc <- proc.time()[3] - t0
  t0 <- proc.time()[3]; fh <- cellmmcd(X, quant_obs = NULL); th <- proc.time()[3] - t0
  t0 <- proc.time()[3]; m <- tryCatch(mmcd(X, alpha = 0.75), error = function(e) NULL); tm <- proc.time()[3] - t0
  klm <- if (is.null(m)) NA else kl_mn(m$cov_row, m$cov_col, S$Sr, S$Sc)
  t0 <- proc.time()[3]; ml <- flipflop(X, center = TRUE); tl <- proc.time()[3] - t0
  sink("output/cm_step10_timing.log", append = TRUE)
  cat(sprintf("p = q = %2d: cellMMCD %7.1fs KL %.3f (raw it %d, rw it %d, recall %.2f, fp %.3f) | hard rw %7.1fs KL %.3f | MMCD %6.1fs KL %.3f | MLE %5.2fs KL %.3f\n",
              p, tc, kl_mn(f$Sr, f$Sc, S$Sr, S$Sc), f$raw$iter, f$iter, mean(!f$W[bad]), mean(!f$W[!bad]), th, kl_mn(fh$Sr, fh$Sc, S$Sr, S$Sc), tm, klm, tl, kl_mn(ml$Sr, ml$Sc, S$Sr, S$Sc)))
  sink() }
cat("done\n")
