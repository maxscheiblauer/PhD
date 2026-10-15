# Step 8 (R1 convergence, counterexample): find clean Gaussian data on which R1, started at the
# plain flip-flop MLE, ends in an exact limit cycle far from any fixed point.
# Run from project root: Rscript exp/archive/step8_r1_cycles.R      (~4 min)
# Search: p = q = 2 (so each factor's shape is one correlation and the whole state is a point in
# the plane (corr Sr, corr Sc)), Sr = AR(.5), Sc = AR(.7), n in 6..40, beta in {.5, .75, .9}.
# Amplitude of a cycle = max shape-KL between any two products visited on it; compared to the
# KL_shape of the MLE to the truth (the size of the ordinary estimation error).
source("R/load.R")
set.seed(808)
S <- true_sigma(2, 2, "full")
corr <- function(M) M[1, 2] / sqrt(M[1, 1] * M[2, 2])

# States visited on the cycle (det Sc = 1), by re-running one period from the final state.
cycle_states <- function(X, beta, f) {
  st <- list(list(Sr = f$Sr, Sc = f$Sc))
  for (k in seq_len(f$period - 1)) {
    g <- r1_fit(X, beta, st[[k]]$Sr, st[[k]]$Sc, maxit = 1, diag = FALSE, tol = -1)
    st[[k + 1]] <- g[c("Sr", "Sc")]
  }
  st
}
amplitude <- function(st) max(outer(seq_along(st), seq_along(st), Vectorize(function(i, j)
  kl_shape(st[[i]]$Sr, st[[i]]$Sc, st[[j]]$Sr, st[[j]]$Sc))))

t0 <- proc.time()[3]; found <- list(); tried <- data.frame()
while (proc.time()[3] - t0 < 200) {
  n <- sample(6:40, 1); beta <- sample(c(0.5, 0.75, 0.9), 1)
  X <- rmatnorm(n, Sr = S$Sr, Sc = S$Sc)
  ff <- flipflop(X, scale_id = "det")
  f <- r1_fit(X, beta, ff$Sr, ff$Sc, maxit = 2000, diag = FALSE)
  tried <- rbind(tried, data.frame(n = n, beta = beta, status = f$status))
  if (f$status != "cycle") next
  st <- cycle_states(X, beta, f)
  found[[length(found) + 1]] <- list(X = X, beta = beta, n = n, period = f$period,
    amp = amplitude(st), err_mle = kl_shape(ff$Sr, ff$Sc, S$Sr, S$Sc), states = st)
}
cat(sprintf("%d clean Gaussian data sets (p = q = 2, ff start), %.0f s\n", nrow(tried), proc.time()[3] - t0))
tab <- with(tried, table(n_band = cut(n, c(5, 10, 20, 40)), beta, status))
print(ftable(tab, row.vars = c("n_band", "beta")))
cyc <- do.call(rbind, lapply(found, function(z) data.frame(n = z$n, beta = z$beta, period = z$period,
  amp = z$amp, err_mle = z$err_mle)))
print(cyc[order(-cyc$amp), ][seq_len(min(15, nrow(cyc))), ], row.names = FALSE)

## the showcase: largest amplitude relative to the MLE's own error
ex <- found[[which.max(cyc$amp / cyc$err_mle)]]
X <- ex$X; beta <- ex$beta; n <- ex$n
ff <- flipflop(X, scale_id = "det")
run <- r1_fit(X, beta, ff$Sr, ff$Sc, maxit = 2000, diag = TRUE)
own_r <- -diff(run$crit[, "Lr"])[seq(1, 2 * run$iter, 2)]
own_c <- -diff(run$crit[, "Lc"])[seq(2, 2 * run$iter, 2)]
on_cycle <- tail(seq_len(run$iter), ex$period)
# independent check with tff (its own cycle detection)
tf <- tff(X, "R1", beta = beta, init = ff[c("Sr", "Sc")], kappa = FALSE, scale_id = "det",
          maxit = 2000, trace = FALSE)
# all other fixed points / cycles: many elemental and random starts
others <- lapply(1:60, function(s) {
  s0 <- if (s <= 30) init_elemental(X) else list(Sr = rWishart(1, 3, diag(2))[, , 1], Sc = rWishart(1, 3, diag(2))[, , 1])
  r1_fit(X, beta, s0$Sr, s0$Sc, maxit = 2000, diag = FALSE)
})
oth <- data.frame(status = sapply(others, `[[`, "status"),
                  J = sapply(others, function(f) sum(r1_crit(X, f$Sr, f$Sc, f$hr, f$hc))),
                  cr = sapply(others, function(f) corr(f$Sr)), cc = sapply(others, function(f) corr(f$Sc)))
Jcyc <- sapply(ex$states, function(s) sum(r1_crit(X, s$Sr, s$Sc, run$hr, run$hc)))

cat(sprintf("\nShowcase: n = %d, beta = %.2f, period %d; amplitude (max shape-KL on the cycle) %.3f,",
            n, beta, ex$period, ex$amp))
cat(sprintf(" KL_shape(MLE, truth) %.3f\n", ex$err_mle))
cat(sprintf("tff (independent implementation) from the same start: %s, period %s\n", tf$status, tf$period))
cat(sprintf("On the cycle every half-step strictly improves its own criterion: min gain Lr %.3g, Lc %.3g\n",
            min(own_r[on_cycle]), min(own_c[on_cycle])))
cat(sprintf("Residual gap on the cycle: %.3g .. %.3g (0 at a fixed point)\n",
            min(run$gap[on_cycle]), max(run$gap[on_cycle])))
cat("Other starts (60):\n"); print(table(oth$status))
cat(sprintf("J~ = Lr~ + Lc~ on the cycle: %.3f .. %.3f; best fixed point J~ %.3f\n",
            min(Jcyc), max(Jcyc), min(oth$J[oth$status == "converged"], Inf)))
saveRDS(list(found = found, tried = tried, cyc = cyc, example = ex, run = run, others = oth, tff = tf[c("status", "period")]),
        "output/step8.rds")

## plots: (1) the state plane, (2) criteria and residual along the run
png("output/step8_cycle.png", 1500, 520, res = 110)
par(mfrow = c(1, 3), mar = c(4.2, 4.2, 2.5, 0.8))
tr <- list(list(Sr = ff$Sr, Sc = ff$Sc)); s0 <- tr[[1]]
for (k in seq_len(run$iter)) { s0 <- r1_fit(X, beta, s0$Sr, s0$Sc, maxit = 1, diag = FALSE, tol = -1)[c("Sr", "Sc")]; tr[[k + 1]] <- s0 }
xr <- sapply(tr, function(s) corr(s$Sr)); xc <- sapply(tr, function(s) corr(s$Sc))
plot(xr, xc, type = "n", xlab = "corr of Sigma^r", ylab = "corr of Sigma^c",
     xlim = range(xr, oth$cr), ylim = range(xc, oth$cc),
     main = sprintf("R1 state, n = %d, beta = %.2f", n, beta))
points(oth$cr[oth$status == "converged"], oth$cc[oth$status == "converged"], pch = 4, col = "grey40", cex = 1.3)
lines(xr, xc, col = "grey70"); points(xr[1], xc[1], pch = 15, col = "#0072B2", cex = 1.5)
cyc_i <- run$iter + 1 - rev(seq_len(ex$period)) + 1
polygon(xr[cyc_i], xc[cyc_i], border = "#D55E00", lwd = 2); points(xr[cyc_i], xc[cyc_i], pch = 19, col = "#D55E00")
points(corr(S$Sr), corr(S$Sc), pch = 8, col = "#009E73", cex = 1.6)
legend("topleft", c("start (MLE)", "limit cycle", "fixed points (other starts)", "truth"),
       pch = c(15, 19, 4, 8), col = c("#0072B2", "#D55E00", "grey40", "#009E73"), bty = "n", cex = 0.85)
hs <- seq(0, run$iter, by = 0.5)
matplot(hs, run$crit, type = "l", lty = 1, col = c("#0072B2", "#D55E00"), lwd = 1.5,
        xlab = "iteration (half-steps)", ylab = "profiled criterion", main = "Lr~ and Lc~ along the run")
legend("topright", c("Lr~ (row step improves it)", "Lc~ (column step improves it)"),
       col = c("#0072B2", "#D55E00"), lwd = 1.5, bty = "n", cex = 0.85)
plot(seq_len(run$iter), pmax(run$gap, 1e-16), log = "y", type = "b", pch = 20, cex = 0.6,
     xlab = "iteration", ylab = "residual gap (log)", main = "one-step residual never reaches 0")
dev.off()
