# Block 4, step 9: interplay of the cell step and the observation step. Pitfall from cm_step8: with 10 % of
# cells at 3-4 sd in every matrix the observation step trims its full 25 % of matrices, because at the first
# observation step (W = Wfix) the still-unflagged cells push d_i beyond mu_i, and a trimmed matrix's cells were
# never re-decided (lock-in). Variants: block-3 default (shadow_it = 0), shadow cell step for trimmed matrices
# (shadow_it = 10, new default), no observation trimming (alpha_obs = 1), quant_obs = 0.999.
# p = q = 5, n = 100, 10 reps, seeds 3900 + rep. Median KL, recall, false flags, trimmed fraction, seconds.
# Output: output/cm_step9_interplay.log
source("R/load.R"); library(parallel); set.seed(39)
sink("output/cm_step9_interplay.log")
p <- 5; q <- 5; n <- 100; S <- true_sigma(p, q, "full"); sd_cell <- sqrt(outer(diag(S$Sr), diag(S$Sc)))
vars <- list(noshadow = list(shadow_it = 0), shadow = list(shadow_it = 10), aobs1 = list(alpha_obs = 1), q999 = list(quant_obs = 0.999))
gen <- function(seed, type, shift) { set.seed(seed); X <- rmatnorm(n, Sr = S$Sr, Sc = S$Sc); bad <- array(FALSE, c(p, q, n))
  if (type %in% c("cells", "mixed")) { fr <- if (type == "cells") 0.10 else 0.05; sh <- if (type == "cells") shift else 6
    bad <- array(runif(p * q * n) < fr, c(p, q, n)); X[bad] <- X[bad] + sh * sd_cell[array(1:(p * q), c(p, q, n))[bad]] * sign(rnorm(sum(bad))) }
  if (type %in% c("shiftmat", "mixed")) { ib <- 1:(if (type == "shiftmat") 20 else 10); X[, , ib] <- X[, , ib] + shift * sd_cell[rep(TRUE, p * q)]; bad[, , ib] <- TRUE }
  list(X = X, bad = bad) }
one <- function(seed, type, shift) { d <- gen(seed, type, shift); out <- c()
  for (v in names(vars)) { t0 <- proc.time()[3]; f <- do.call(cellmmcd, c(list(X = d$X), vars[[v]]))
    fl <- !f$W | array(rep(!f$v, each = p * q), c(p, q, n))
    out[paste0(v, "_kl")] <- kl_mn(f$Sr, f$Sc, S$Sr, S$Sc); out[paste0(v, "_rec")] <- if (any(d$bad)) mean(fl[d$bad]) else NA
    out[paste0(v, "_fp")] <- mean(fl[!d$bad]); out[paste0(v, "_tr")] <- mean(!f$v); out[paste0(v, "_t")] <- proc.time()[3] - t0 }
  out }
cat("p = q = 5, n = 100, 10 reps. Median KL / recall / false flag rate (cells flagged or matrix trimmed) / trimmed matrices / seconds.\n")
for (sc in list(c("none", 0), c("cells", 3), c("cells", 4), c("cells", 6), c("shiftmat", 5), c("mixed", 5))) {
  r <- do.call(rbind, mclapply(1:10, function(s) one(3900 + s, sc[1], as.numeric(sc[2])), mc.cores = 4)); v <- apply(r, 2, median, na.rm = TRUE)
  cat(sprintf("%-9s shift %s: ", sc[1], sc[2])); for (nm in names(vars)) cat(sprintf("%s KL %.3f rec %.2f fp %.3f tr %.2f %.1fs | ", nm, v[paste0(nm, "_kl")], v[paste0(nm, "_rec")], v[paste0(nm, "_fp")], v[paste0(nm, "_tr")], v[paste0(nm, "_t")])); cat("\n") }
sink(); cat("done\n")
