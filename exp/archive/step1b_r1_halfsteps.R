# Step 1b (Q1, R1 only): does any likelihood-type criterion improve in every half-step of R1?
# Run from project root: Rscript exp/archive/step1b_r1_halfsteps.R
# R1 with kappa off, denominator h, no rescaling, symmetric roots; data: clean / col / row / wcol (eps = .1).
# Criteria are -2 log-likelihood type (smaller = better), each minimised over the overall scale
# first, so the scale ping-pong between the two half-steps does not count as a change:
#   Lr = Sr-step criterion (kept columns), Lc = Sc-step criterion (kept rows),
#   J = Lr + Lc (spec section 6), Js = RJ trimmed likelihood evaluated at R1's kept columns.
# Lr and Lc also depend on the other factor's scale, so only J and Js are fully scale-free.
source("R/load.R"); set.seed(7)
p <- q <- 5; n <- 25; beta <- 0.75; S <- true_sigma(p, q, "full")
hr <- round(beta * n * q); hc <- round(beta * n * p)
prof <- function(H, Q) H * log(Q / H) + H          # min_t  H log t + Q / t
crit <- function(Sr, Sc, a, b) {
  D <- distances(X, Sr, Sc); sd <- sum(D$d[a]); se <- sum(D$e[b])
  c(Lr = hr * logdet(Sr) + prof(hr * p, sd),
    Lc = hc * logdet(Sc) + prof(hc * q, se),
    J  = hr * logdet(Sr) + hc * logdet(Sc) + prof(hr * p, sd + se),
    Js = hr * logdet(Sr) + hr * p / q * logdet(Sc) + prof(hr * p, sd))
}
up <- list(); tol <- 1e-9; runs <- list()
for (cd in c("clean", "col", "row", "wcol")) for (ds in 1:10) for (st in c("identity", "ff", "elemental")) {
  X <- gen_data(n, S$Sr, S$Sc, cd, if (cd == "clean") 0 else 0.1)$X
  s0 <- make_init(X, st); Sr <- s0$Sr; Sc <- s0$Sc
  D <- distances(X, Sr, Sc); a <- csel(D$d, hr); b <- csel(D$e, hc)
  old <- crit(Sr, Sc, a, b); worse <- NULL
  for (k in 1:40) {
    Z <- zmat(X, root(Sc)); a <- csel(matrix(qdist(Z, Sr), n), hr)
    Sr <- tcrossprod(Z[, as.vector(a)]) / hr
    new <- crit(Sr, Sc, a, b); w1 <- new > old + tol * abs(old); old <- new
    U <- umat(X, t(root(Sr))); b <- csel(matrix(qdist(U, Sc), n), hc)
    Sc <- tcrossprod(U[, as.vector(b)]) / hc
    new <- crit(Sr, Sc, a, b); w2 <- new > old + tol * abs(old); old <- new
    up[[length(up) + 1]] <- data.frame(cond = cd, step = c("Sr step", "Sc step"), rbind(w1, w2))
    worse <- rbind(worse, w1, w2)
  }
  runs[[length(runs) + 1]] <- data.frame(cond = cd, t(colSums(worse) > 0))
}
up <- do.call(rbind, up); runs <- do.call(rbind, runs)
saveRDS(list(halfsteps = up, runs = runs), "output/step1b.rds")
cat("Fraction of HALF-STEPS in which the criterion got worse (= likelihood went down):\n")
print(aggregate(cbind(Lr, Lc, J, Js) ~ step, up, mean), digits = 2)
cat("\nFraction of RUNS (40 iterations) with at least one worsening half-step:\n")
print(colMeans(runs[, -1]), digits = 2)
print(aggregate(cbind(Lr, Lc, J, Js) ~ cond, runs, mean), digits = 2)
