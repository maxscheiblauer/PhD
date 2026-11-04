# Block 3: first-order bias formula Sigma* = Sigma - kappa diag(1/Omega_cc), kappa = F1(c) - F3(c), projected
# onto the Kronecker set by the population flip-flop, against the empirical limit of the default
# (cm_step6_bias.log: KL 0.0116 at n = 8000, 3x4 full shape). Also the fixed point of the population map
# iterated (plug the biased Sigma back into the conditional variances). Output: output/cm_step6b_biasformula.log
source("R/load.R"); sink("output/cm_step6b_biasformula.log")
kron_proj <- function(Sig, p, q, iter = 200) {     # population flip-flop: Kronecker MLE projection
  Sr <- diag(p); Sc <- diag(q); idx <- function(l) (l - 1) * p + 1:p
  for (it in 1:iter) {
    Oc <- solve(Sc); Sr <- matrix(0, p, p)
    for (l in 1:q) for (m in 1:q) Sr <- Sr + Oc[l, m] * Sig[idx(l), idx(m)]
    Sr <- Sr / q; Or <- solve(Sr)
    for (l in 1:q) for (m in 1:q) Sc[l, m] <- sum(Or * Sig[idx(l), idx(m)]) / p
  }
  list(Sr = Sr, Sc = Sc)
}
for (quant in c(0.99, 0.999)) {
  cq <- qchisq(quant, 1); kappa <- pchisq(cq, 1) - pchisq(cq, 3)
  for (pq in list(c(3, 4), c(5, 5))) { p <- pq[1]; q <- pq[2]; S <- true_sigma(p, q, "full")
    Sig <- kronecker(S$Sc, S$Sr); Om <- solve(Sig)
    pc <- -Om / sqrt(outer(diag(Om), diag(Om))); diag(pc) <- 0
    Sb <- Sig - kappa * diag(1 / diag(Om)); pr <- kron_proj(Sb, p, q)
    # iterate the population map: conditional variances at the current (biased) estimate
    Sk <- Sig; for (it in 1:50) { Ok <- solve(Sk); Sk <- kron_proj(Sig - kappa * diag(1 / diag(Ok)), p, q); Sk <- kronecker(Sk$Sc, Sk$Sr) }
    pk <- kron_proj(Sk, p, q)
    cat(sprintf("quant %.3f  %dx%d  kappa %.4f  sum pcor^2 %.2f  shapeKL(first order) %.5f  approx kappa^2/4 sum pcor^2 %.5f  shapeKL(iterated map) %.5f  tr ratio %.3f\n",
                quant, p, q, kappa, sum(pc^2), kl_shape(pr$Sr, pr$Sc, S$Sr, S$Sc), kappa^2 / 4 * sum(pc^2), kl_shape(pk$Sr, pk$Sc, S$Sr, S$Sc), tr_ratio(pk$Sr, pk$Sc, S$Sr, S$Sc)))
  }
}
sink()
