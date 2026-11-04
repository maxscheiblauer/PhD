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

cat(sprintf("\n%d / %d tests passed\n", sum(res), length(res)))
if (!all(res)) quit(status = 1)
