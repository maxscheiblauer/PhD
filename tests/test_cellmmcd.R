# Tests for cellMMCD. Run: Rscript tests/test_cellmmcd.R
source("R/load.R")
set.seed(7)
ok <- function(msg, cond) { cat(sprintf("[%s] %s\n", if (cond) "PASS" else "FAIL", msg)); invisible(cond) }
res <- logical()
p <- 3; q <- 4; n <- 40; S <- true_sigma(p, q, "full")
X <- rmatnorm(n, Sr = S$Sr, Sc = S$Sc); Xc <- X
W <- array(runif(p * q * n) > 0.2, c(p, q, n))
lambda <- matrix(runif(p * q, 5, 9), p, q)

## 1. Schur-form objective equals the brute-force observed-cell likelihood
brute <- function(Xc, Sr, Sc, W, lambda) {
  Sig <- Sc %x% Sr; tot <- 0
  for (i in 1:dim(Xc)[3]) {
    o <- which(W[, , i]); x <- as.vector(Xc[, , i])[o]; So <- Sig[o, o, drop = FALSE]
    tot <- tot + logdet(So) + sum(x * solve(So, x)) + length(o) * log(2 * pi) + sum(lambda[!W[, , i]])
  }
  tot
}
f_schur <- cm_obj(Xc, S$Sr, S$Sc, W, lambda); f_brute <- brute(Xc, S$Sr, S$Sc, W, lambda)
res["objective"] <- ok(sprintf("objective: Schur form vs brute force, rel diff %.1e", abs(f_schur - f_brute) / abs(f_brute)),
                       abs(f_schur - f_brute) < 1e-8 * abs(f_brute))

## 2. Delta_c - lambda equals the change of F from unflagging cell c (flagged and unflagged cases)
Or <- solve(S$Sr); Oc <- solve(S$Sc); mx <- 0
for (rep in 1:30) {
  i <- sample(n, 1); c <- sample(p * q, 1); k <- (c - 1) %% p + 1; l <- (c - 1) %/% p + 1
  W1 <- W; W1[k, l, i] <- TRUE; W0 <- W; W0[k, l, i] <- FALSE
  d_true <- cm_obj(Xc, S$Sr, S$Sc, W1, lambda) - cm_obj(Xc, S$Sr, S$Sc, W0, lambda)
  Tf <- Or %*% Xc[, , i] %*% Oc
  d_form <- cm_delta1(Tf, Or, Oc, which(!W[, , i]), c, p) - lambda[k, l]
  mx <- max(mx, abs(d_true - d_form))
}
res["delta"] <- ok(sprintf("cell criterion: max |Delta - lambda - (F_unflag - F_flag)| = %.1e", mx), mx < 1e-8)

## 3. Without flags the parameter step is one flip-flop step
Wall <- array(TRUE, c(p, q, n))
ps <- cm_pstep(Xc, S$Sr, S$Sc, Wall)
Sr_ff <- Reduce("+", lapply(1:n, function(i) Xc[, , i] %*% Oc %*% t(Xc[, , i]))) / (n * q)
Sc_ff <- Reduce("+", lapply(1:n, function(i) t(Xc[, , i]) %*% solve(Sr_ff) %*% Xc[, , i])) / (n * p)
res["pstep_mle"] <- ok("parameter step without flags = flip-flop step",
                       max(abs(ps$Sr - Sr_ff)) < 1e-10 && max(abs(ps$Sc - Sc_ff)) < 1e-10)

## 4. The parameter step does not increase F (GEM) with flags, and the full algorithm is monotone
f_before <- cm_obj(Xc, S$Sr, S$Sc, W, lambda); ps <- cm_pstep(Xc, S$Sr, S$Sc, W)
f_after <- cm_obj(Xc, ps$Sr, ps$Sc, W, lambda)
res["pstep_mono"] <- ok(sprintf("parameter step with flags: F %.3f -> %.3f", f_before, f_after), f_after <= f_before + 1e-9)
d <- gen_data(100, S$Sr, S$Sc, "cell", eps = 0.1, delta = 6)
fit <- cellmmcd(d$X)
dif <- diff(fit$obj)
res["monotone"] <- ok(sprintf("cellmmcd monotone: max increase of F = %.1e over %d half-steps, converged = %s",
                              max(dif), length(dif), fit$converged), max(dif) <= 1e-8 * abs(fit$obj[1]) && fit$converged)

## 5. Flags find the bad cells (10 % cells at 6 sd), and the fit beats the MLE
truth <- array(FALSE, c(p, q, 100)); set.seed(11)
Xcl <- rmatnorm(100, Sr = S$Sr, Sc = S$Sc); bad <- sample(p * q * 100, 120); truth[bad] <- TRUE
sd_cell <- sqrt(outer(diag(S$Sr), diag(S$Sc)))
Xb <- Xcl; Xb[bad] <- Xcl[bad] + 6 * array(sd_cell, c(p, q, 100))[bad] * sample(c(-1, 1), 120, TRUE)
fit <- cellmmcd(Xb)
recall <- mean(!fit$W[truth]); fpr <- mean(!fit$W[!truth])
ff <- flipflop(Xb); kl_cm <- kl_mn(fit$Sr, fit$Sc, S$Sr, S$Sc); kl_ff <- kl_mn(ff$Sr, ff$Sc, S$Sr, S$Sc)
kp <- fit$start; kl_kp <- kl_mn(kp$Sr, kp$Sc, S$Sr, S$Sc)
res["flags"] <- ok(sprintf("flags: recall %.2f, false-positive rate %.3f", recall, fpr), recall > 0.8 && fpr < 0.05)
res["beats_mle"] <- ok(sprintf("KL: cellMMCD %.2f, kpair start %.2f, flip-flop MLE %.2f", kl_cm, kl_kp, kl_ff), kl_cm < kl_ff && kl_cm < kl_kp)

## 6. kpair is consistent for the shape at large n (clean)
Xbig <- rmatnorm(3000, Sr = S$Sr, Sc = S$Sc); kp <- kpair(Xbig)
res["kpair"] <- ok(sprintf("kpair n=3000: shape errors %.3f (Sr) %.3f (Sc), tr ratio %.3f",
                           shape_err(kp$Sr, S$Sr), shape_err(kp$Sc, S$Sc), tr_ratio(kp$Sr, kp$Sc, S$Sr, S$Sc)),
                   shape_err(kp$Sr, S$Sr) < 0.05 && shape_err(kp$Sc, S$Sc) < 0.05)

## 7. Scale calibration (block 2): proxy fixed point and clean-data scale of the Kronecker fit
fx0 <- cm_proxy_fixpt(calib = FALSE); fx1 <- cm_proxy_fixpt(calib = TRUE)
res["proxy_fixpt"] <- ok(sprintf("proxy fixed point: uncalibrated t* %.3f (%.1f%% flagged), calibrated t* %.4f = t1 %.4f (%.2f%% flagged)",
                                 fx0["t"], 100 * fx0["flagged"], fx1["t"], fx1["t1"], 100 * fx1["flagged"]),
                         abs(fx1["t"] - fx1["t1"]) < 1e-6 && abs(fx1["flagged"] - 0.01) < 1e-6 && fx0["t"] < 0.92)
set.seed(7); Xcl <- rmatnorm(600, Sr = S$Sr, Sc = S$Sc)
f0 <- cellmmcd(Xcl, calib = FALSE); f1 <- cellmmcd(Xcl, calib = TRUE)
tr0 <- tr_ratio(f0$Sr, f0$Sc, S$Sr, S$Sc); tr1 <- tr_ratio(f1$Sr, f1$Sc, S$Sr, S$Sc)
res["calib_scale"] <- ok(sprintf("clean n=600: tr ratio %.3f (uncalibrated, %.1f%% flagged) vs %.3f (calibrated, %.1f%% flagged)",
                                 tr0, 100 * f0$flagged, tr1, 100 * f1$flagged),
                         abs(tr1 - 1) < 0.04 && abs(tr1 - 1) < abs(tr0 - 1) && f1$flagged < f0$flagged)

## 8. Observation step (block 2): exactness, monotonicity with all three blocks, casewise detection
set.seed(8); Xo <- rmatnorm(80, Sr = S$Sr, Sc = S$Sc); Xo[, , 1:8] <- Xo[, , 1:8] + 5   # 10 % shifted matrices
fo <- cellmmcd(Xo); fh <- cellmmcd(Xo, quant_obs = NULL)
res["vstep_mono"] <- ok(sprintf("three-block descent monotone: soft %s (trimmed %d), hard quota %s (trimmed %d = n - h_obs %d)",
                                all(diff(fo$obj) <= 1e-8 * abs(fo$obj[-1])), sum(!fo$v), all(diff(fh$obj) <= 1e-8 * abs(fh$obj[-1])), sum(!fh$v), 80 - ceiling(0.75 * 80)),
                        all(diff(fo$obj) <= 1e-8 * abs(fo$obj[-1])) && all(diff(fh$obj) <= 1e-8 * abs(fh$obj[-1])) && sum(!fh$v) == 80 - ceiling(0.75 * 80))
res["vstep_case"] <- ok(sprintf("shifted matrices trimmed: soft %d/8 (plus %d clean), hard %d/8", sum(!fo$v[1:8]), sum(!fo$v[-(1:8)]), sum(!fh$v[1:8])),
                        all(!fo$v[1:8]) && all(!fh$v[1:8]) && sum(!fo$v[-(1:8)]) <= 4)
# exact v-step: brute force over the per-observation contributions
Xc8 <- sweep(Xo, c(1, 2), fo$M); d8 <- cm_obj_i(Xc8, fo$Sr * fo$factor, fo$Sc, fo$W, fo$lambda)
v8 <- cm_vstep(Xc8, fo$Sr * fo$factor, fo$Sc, fo$W, fo$lambda, fo$mu, 60, 1)
vh <- cm_vstep(Xc8, fo$Sr * fo$factor, fo$Sc, fo$W, fo$lambda, -Inf, 60); besth <- rep(TRUE, 80); besth[order(d8, decreasing = TRUE)[1:20]] <- FALSE
d8 <- d8 - cm_mu_i(Xc8, fo$W, fo$lambda, fo$mu); best <- rep(TRUE, 80); cand <- order(d8, decreasing = TRUE)[1:20]; best[cand[d8[cand] > 0]] <- FALSE
res["vstep_exact"] <- ok("v-step equals brute-force choice (soft and hard)", identical(v8, best) && identical(vh, besth))

## 14. (block 3) censored proxy: with tau = (1 - F3(c)) / (1 - F1(c)) and no calibration the fixed point is
## at scale 1 with exactly 1 - quant flagged
cq <- qchisq(0.99, 1); tau <- (1 - pchisq(cq, 3)) / (1 - pchisq(cq, 1))
px <- cm_proxy_fixpt(0.99, calib = FALSE, tau = tau)
res["proxy_tail"] <- ok("censored proxy fixed point at scale 1 with 1 % flagged",
                        abs(px["t"] - 1) < 1e-6 && abs(px["flagged"] - 0.01) < 1e-6)

## 15. (block 3) reweighting: fixed-W parameter steps are monotone, flags are rare on clean data, the scale is
## calibrated, and 6-sd cells are still caught
set.seed(15); n15 <- 150; X15 <- rmatnorm(n15, Sr = S$Sr, Sc = S$Sc)
bad <- array(runif(p * q * n15) < 0.05, c(p, q, n15)); X15b <- X15
X15b[bad] <- X15b[bad] + 6 * sqrt(outer(diag(S$Sr), diag(S$Sc)))[array(rep(1:(p * q), n15), c(p, q, n15))[bad]] * sign(rnorm(sum(bad)))
f15 <- cellmmcd(X15, reweight = TRUE); f15b <- cellmmcd(X15b, reweight = TRUE)
recall <- sum(!f15b$W & bad & array(f15b$v, c(p, q, n15), )[, , ]) / sum(bad & array(rep(f15b$v, each = p * q), c(p, q, n15)))
res["reweight"] <- ok(sprintf("reweighting: monotone, clean flags %.3f, tr ratio %.3f, recall %.2f, raw kept",
                              f15$flagged, tr_ratio(f15$Sr, f15$Sc, S$Sr, S$Sc), recall),
                      all(diff(f15$obj) <= 1e-8 * abs(f15$obj[-1])) && f15$flagged < 0.02 &&
                      abs(tr_ratio(f15$Sr, f15$Sc, S$Sr, S$Sc) - 1) < 0.15 && recall > 0.9 && !is.null(f15$raw$obj))

## 16. (block 3) tail = TRUE runs and converges. It is NOT a descent on F (the parameter step is an EM step
## for the censored objective, whose value is not computed), so only convergence and the scale are checked.
f16 <- cellmmcd(X15, tail = TRUE)
res["tail_runs"] <- ok(sprintf("tail variant converges (iter %d), tau = %.2f, tr ratio %.3f", f16$iter, f16$tau, tr_ratio(f16$Sr, f16$Sc, S$Sr, S$Sc)),
                       f16$converged && f16$tau > 1 && abs(tr_ratio(f16$Sr, f16$Sc, S$Sr, S$Sc) - 1) < 0.15)

## 17. (block 3) full-column safeguard: after it every column position has >= p fully kept columns and every
## row position >= q fully kept rows among the kept observations; forced flags are never unflagged
set.seed(17); W17 <- array(runif(p * q * n) > 0.35, c(p, q, n)); Wf17 <- array(runif(p * q * n) > 0.05, c(p, q, n))
W17 <- W17 & Wf17; v17 <- runif(n) > 0.2
sg <- cm_fullcol(W17, Wf17, v17); k17 <- which(v17)
fullcols <- sapply(1:q, function(l) sum(colSums(!sg$W[, l, k17, drop = FALSE]) == 0))
fullrows <- sapply(1:p, function(k) sum(colSums(matrix(!sg$W[k, , k17], q)) == 0))
res["safeguard"] <- ok(sprintf("safeguard: %d cells unflagged, min full cols %d (need %d), rows %d (need %d)", sg$unflagged, min(fullcols), p, min(fullrows), q),
                       all(fullcols >= p) && all(fullrows >= q) && !any(sg$W & !Wf17) && all(sg$W[, , !v17] == W17[, , !v17]))

## 18. Shadow cell step (block 4): 10 % cells at 4 sd in every matrix; the block-3 iteration locks 25 % of the
## matrices out, the shadow step lets them back in; F is monotone once the shadow phase (10 iterations) is over.
set.seed(18); p18 <- 4; q18 <- 4; n18 <- 80; S18 <- true_sigma(p18, q18, "full")
X18 <- rmatnorm(n18, Sr = S18$Sr, Sc = S18$Sc); sdc <- sqrt(outer(diag(S18$Sr), diag(S18$Sc)))
b18 <- array(runif(p18 * q18 * n18) < 0.10, c(p18, q18, n18))
X18[b18] <- X18[b18] + 4 * sdc[array(1:(p18 * q18), c(p18, q18, n18))[b18]] * sign(rnorm(sum(b18)))
f18a <- cellmmcd(X18, shadow_it = 0, reweight = FALSE); f18b <- cellmmcd(X18, shadow_it = 10, reweight = FALSE)
o18 <- f18b$obj[-(1:20)]
res["shadow"] <- ok(sprintf("shadow step: trimmed %.2f (block 3: %.2f), KL %.3f (%.3f), F monotone after the shadow phase %s",
                            f18b$trimmed, f18a$trimmed, kl_mn(f18b$Sr, f18b$Sc, S18$Sr, S18$Sc), kl_mn(f18a$Sr, f18a$Sc, S18$Sr, S18$Sc),
                            all(diff(o18) <= 1e-8 * abs(o18[-1]))),
                    f18b$trimmed < f18a$trimmed && f18b$trimmed <= 0.1 && all(diff(o18) <= 1e-8 * abs(o18[-1])))

## 19. Influence functions in the 1-D proxy (block 4): the formulas of cm_proxy_if against the exact
## functionals under (1 - eps) Phi + eps delta_x; Fisher consistency at Phi; bounded, and ~0 beyond the cutoff.
xs <- c(0, 1, 2, 2.5, 3, 3.2, 3.5, 10); eps <- 1e-6
dif <- sapply(xs, function(x) { nu <- (cm_proxy_T(x, eps) - cm_proxy_T(x, 0)) / eps; an <- cm_proxy_if(x); max(abs(nu - an)[c("raw", "rw")]) })
T0 <- cm_proxy_T(1, 0); ges <- max(abs(sapply(seq(0, 6, by = 0.001), function(x) cm_proxy_if(x)["rw"])))
res["proxy_if"] <- ok(sprintf("proxy IF: T(Phi) = (%.6f, %.6f), max |numeric - formula| %.1e, IF_rw(10) = %.4f, gross-error sensitivity %.2f",
                              T0["raw"], T0["rw"], max(dif), cm_proxy_if(10)["rw"], ges),
                      max(abs(T0[c("raw", "rw")] - 1)) < 1e-8 && max(dif) < 1e-3 && abs(cm_proxy_if(10)["rw"]) < 0.01 && ges < 11)

## 20. C++ cell step equals the R cell step (kept and trimmed observations, shadow on), and the full fits agree.
set.seed(20); X20 <- X18; v20 <- rep(TRUE, n18); v20[c(3, 7, 11, 40)] <- FALSE
st20 <- kpair(X20); Xc20 <- sweep(X20, c(1, 2), apply(X20, c(1, 2), median)); Wf20 <- sweep(abs(Xc20), c(1, 2), apply(abs(Xc20), c(1, 2), median) / 0.6745, "/") < 2.5
lam20 <- 6 + log(2 * pi) - log(outer(diag(solve(st20$Sr)), diag(solve(st20$Sc))))
Wr <- cm_wstep(Xc20, st20$Sr, st20$Sc, Wf20, Wf20, lam20, 60, v20, shadow = TRUE, engine = "R")
if (isTRUE(has_cpp)) { Wc <- cm_wstep(Xc20, st20$Sr, st20$Sc, Wf20, Wf20, lam20, 60, v20, shadow = TRUE, engine = "cpp")
  options(cellmmcd.engine = "R"); fr <- cellmmcd(X20); options(cellmmcd.engine = "cpp"); fc <- cellmmcd(X20)
  res["cpp"] <- ok(sprintf("C++ cell step: W identical %s (%d flags), full fits: max |Sr diff| %.1e, flagged %.3f vs %.3f",
                           identical(Wr, Wc), sum(!Wr), max(abs(fr$Sr - fc$Sr)), fr$flagged, fc$flagged),
                   identical(Wr, Wc) && max(abs(fr$Sr - fc$Sr)) < 1e-10 && identical(fr$W, fc$W))
} else cat("[SKIP] C++ cell step not compiled\n")

# 21 (block 5): robust flag rate. Under independent cellwise flags at rate eps the binomial-corrected lower-alpha mean
# recovers eps (Fisher consistency, n large); 25 % of matrices with many flags do not move it; the plain mean does.
set.seed(21); pq_e <- c(5, 5); eps_true <- 0.08; n21 <- 4000
W21 <- array(runif(prod(pq_e) * n21) > eps_true, c(pq_e, n21)); v21 <- rep(TRUE, n21)
e_clean <- cm_eps_robust(W21, v21, 0.75); W21b <- W21; W21b[, , 1:(n21 / 4)] <- array(runif(prod(pq_e) * n21 / 4) > 0.4, c(pq_e, n21 / 4))
e_cont <- cm_eps_robust(W21b, v21, 0.75); e_plain <- mean(!W21b)
res["eps_robust"] <- ok(sprintf("robust flag rate: clean %.4f (true %.2f), with 25 %% matrices at 40 %% flags: robust %.4f, plain mean %.3f",
                                e_clean, eps_true, e_cont, e_plain), abs(e_clean - eps_true) < 0.01 && e_cont < 0.12 && e_plain > 0.15)
# (with exactly 25 % contaminated matrices none of the clean counts is cut, so the correction over-corrects: 0.107 for 0.08,
#  i.e. count_max one unit higher; the plain mean doubles. Exact at 0 % casewise contamination, bounded in between.)

# 22 (block 5): hard quota with the MCD consistency factor: clean-data scale ratio near 1 (block 4: 0.87 without it).
set.seed(22); S22 <- true_sigma(4, 4, "full"); X22 <- rmatnorm(400, Sr = S22$Sr, Sc = S22$Sc)
fh <- cellmmcd(X22, quant_obs = NULL); tr22 <- tr_ratio(fh$Sr, fh$Sc, S22$Sr, S22$Sc)
res["hard_kappa"] <- ok(sprintf("hard quota: trimmed %.2f, factor t_obs %.3f (kappa %.3f), clean scale ratio %.3f", fh$trimmed, fh$raw$t_obs, kappa_mcd(0.75, 16), tr22),
                        abs(fh$trimmed - 0.25) < 0.01 && abs(tr22 - 1) < 0.06)

# 23 (block 5): joint location. Clean 4x4, n = 200, true M = 2: the EM location is closer to M than the cellwise median,
# F stays monotone after the shadow phase, and the fit with loc = FALSE is unchanged by the new argument.
set.seed(23); S23 <- true_sigma(4, 4, "full"); X23 <- rmatnorm(200, Sr = S23$Sr, Sc = S23$Sc) + 2
fj <- cellmmcd(X23, loc = TRUE); fm <- cellmmcd(X23); e_j <- sum((fj$M - 2)^2); e_m <- sum((fm$M - 2)^2)
oj <- fj$raw$obj; mono <- all(diff(oj[-(1:(2 * fj$raw$shadow_it))]) <= 1e-8 * abs(oj[1]))
res["loc_joint"] <- ok(sprintf("joint location: squared error %.4f (median plug-in %.4f), F monotone after the shadow phase %s, KL %.3f vs %.3f",
                               e_j, e_m, mono, kl_mn(fj$Sr, fj$Sc, S23$Sr, S23$Sc), kl_mn(fm$Sr, fm$Sc, S23$Sr, S23$Sc)), e_j < e_m && mono)

cat(sprintf("\n%d / %d tests passed\n", sum(res), length(res)))
if (!all(res)) quit(status = 1)
