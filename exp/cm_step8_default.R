# Block 3: which default? quant .99 (biased at the Gaussian, cm_step6), quant .999, or raw .99 + reweighting at
# .999. Robustness against cellwise shifts of 3, 4, 6 sd (10 % of cells, random sign) and against the
# 20 % shifted-matrix casewise setting; p = q = 5, n = 100, 10 reps. Median KL, recall, false flags, seconds.
# Output: output/cm_step8_default.log
source("R/load.R"); library(parallel); set.seed(38)
sink("output/cm_step8_default.log")
p <- 5; q <- 5; n <- 100; S <- true_sigma(p, q, "full"); sd_cell <- sqrt(outer(diag(S$Sr), diag(S$Sc)))
vars <- list(q99 = list(), q999 = list(quant = 0.999), rw = list(reweight = TRUE), rw995 = list(reweight = TRUE, quant_rw = 0.995))
gen <- function(seed, type, shift) { set.seed(seed); X <- rmatnorm(n, Sr = S$Sr, Sc = S$Sc); bad <- array(FALSE, c(p, q, n))
  if (type == "cells") { bad <- array(runif(p * q * n) < 0.10, c(p, q, n)); X[bad] <- X[bad] + shift * sd_cell[array(1:(p * q), c(p, q, n))[bad]] * sign(rnorm(sum(bad))) }
  if (type == "shiftmat") { ib <- 1:20; X[, , ib] <- X[, , ib] + shift * sd_cell[rep(TRUE, p * q)]; bad[, , ib] <- TRUE }
  list(X = X, bad = bad) }
one <- function(seed, type, shift) { d <- gen(seed, type, shift); out <- c()
  for (v in names(vars)) { t0 <- proc.time()[3]; f <- do.call(cellmmcd, c(list(X = d$X), vars[[v]]))
    fl <- !f$W | array(rep(!f$v, each = p * q), c(p, q, n))
    out[paste0(v, "_kl")] <- kl_mn(f$Sr, f$Sc, S$Sr, S$Sc); out[paste0(v, "_rec")] <- if (any(d$bad)) mean(fl[d$bad]) else NA
    out[paste0(v, "_fp")] <- mean(fl[!d$bad]); out[paste0(v, "_t")] <- proc.time()[3] - t0 }
  out }
cat("p = q = 5, n = 100, 10 reps. Median KL / recall / false flag rate (cells flagged or matrix trimmed) / seconds.\n")
for (sc in list(c("none", 0), c("cells", 3), c("cells", 4), c("cells", 6), c("shiftmat", 2), c("shiftmat", 5))) {
  r <- do.call(rbind, mclapply(1:10, function(s) one(3800 + s, sc[1], as.numeric(sc[2])), mc.cores = 4)); v <- apply(r, 2, median, na.rm = TRUE)
  cat(sprintf("%-9s shift %s: ", sc[1], sc[2])); for (nm in names(vars)) cat(sprintf("%s KL %.3f rec %.2f fp %.3f %.1fs | ", nm, v[paste0(nm, "_kl")], v[paste0(nm, "_rec")], v[paste0(nm, "_fp")], v[paste0(nm, "_t")])); cat("\n") }
sink()
