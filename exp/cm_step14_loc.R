# cellMMCD step 14 (block 5): joint location (loc = TRUE, EM step for M on the kept cells) vs the cellwise median plug-in. Seed 14.
# p = q = 5, n = 100, 10 reps, true M = 1 in every cell. Settings: clean; 5 % / 10 % cells at 6 / 4 sd; 20 % matrices shifted by 5 sd; mixed.
# Reported: median KL, median squared error of M (sum over cells), trimmed fraction, false flags. Output: output/cm_step14_loc.log
source("R/load.R"); library(parallel); set.seed(14)
p <- q <- 5; n <- 100; S <- true_sigma(p, q, "full"); sd_cell <- sqrt(outer(diag(S$Sr), diag(S$Sc)))
gen <- function(seed, type) { set.seed(seed); X <- rmatnorm(n, Sr = S$Sr, Sc = S$Sc) + 1; bad <- array(FALSE, c(p, q, n))
  if (type %in% c("cells6", "cells4", "mixed")) { fr <- if (type == "cells4") 0.10 else 0.05; sh <- if (type == "cells4") 4 else 6
    bad <- array(runif(p * q * n) < fr, c(p, q, n)); X[bad] <- X[bad] + sh * sd_cell[array(1:(p * q), c(p, q, n))[bad]] * sign(rnorm(sum(bad))) }
  if (type %in% c("shiftmat", "mixed")) { ib <- 1:(if (type == "shiftmat") 20 else 10); X[, , ib] <- X[, , ib] + 5 * sd_cell[rep(TRUE, p * q)]; bad[, , ib] <- TRUE }
  list(X = X, bad = bad) }
one <- function(seed, type) { d <- gen(seed, type); out <- c()
  for (v in c("median", "joint")) { f <- cellmmcd(d$X, loc = (v == "joint")); fl <- !f$W | array(rep(!f$v, each = p * q), c(p, q, n))
    out[paste0(v, "_kl")] <- kl_mn(f$Sr, f$Sc, S$Sr, S$Sc); out[paste0(v, "_Merr")] <- sum((f$M - 1)^2); out[paste0(v, "_tr")] <- mean(!f$v); out[paste0(v, "_fp")] <- mean(fl[!d$bad]) }
  ff <- flipflop(d$X); out["mle_Merr"] <- sum((apply(d$X, c(1, 2), mean) - 1)^2); out }
sink("output/cm_step14_loc.log"); cat("p = q = 5, n = 100, 10 reps, M = 1. Median KL / squared error of M / trimmed / false flags; MLE mean's M error for reference.\n")
for (type in c("clean", "cells6", "cells4", "shiftmat", "mixed")) {
  r <- do.call(rbind, mclapply(1:10, function(s) one(1400 + s, type), mc.cores = 4)); v <- apply(r, 2, median)
  cat(sprintf("%-8s: median KL %.3f Merr %.3f tr %.2f fp %.3f | joint KL %.3f Merr %.3f tr %.2f fp %.3f | sample mean Merr %.3f\n", type,
              v["median_kl"], v["median_Merr"], v["median_tr"], v["median_fp"], v["joint_kl"], v["joint_Merr"], v["joint_tr"], v["joint_fp"], v["mle_Merr"])) }
sink(); cat("done\n")
