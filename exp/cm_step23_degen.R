# Block 7, step 23: numerical check of the configuration criterion for boundedness of F (JOURNAL block 7, Propositions A/B,
# Conjecture C). For a kept pattern (W, v) define c_T = #{(i, l): v_i = 1, T x {l} subset of O_i} for row sets T and
# d_U = #{(i, k): {k} x U subset of O_i} for column sets U. Prop. A: 1 <= c_T <= |T| - 1 (or the same for d_U) implies
# inf F = -inf. Conjecture C: otherwise F is bounded below. Check: (1) raw fits (no reweighting) on clean 5x5 data at
# n = 8, 9, 10, 12 (where cm_step5 found degenerate descents) and at n = 100 (clean, 20 % cells): does the criterion
# separate the degenerate runs from the others? (2) In degenerate runs, is the null direction a of the factor supported
# on a set T whose kept T-columns are orthogonal to a (the mechanism of Prop. A)? Seed 23. Output: output/cm_step23_degen.log
source("R/load.R"); library(parallel); set.seed(23)
subsets <- function(m) { s <- list(); for (k in 1:m) s <- c(s, combn(m, k, simplify = FALSE)); s }
criterion <- function(W, v) { p <- dim(W)[1]; q <- dim(W)[2]; Wv <- W[, , v, drop = FALSE]
  cT <- sapply(subsets(p), function(T) sum(apply(Wv[T, , , drop = FALSE], c(2, 3), all)))
  dU <- sapply(subsets(q), function(U) sum(apply(Wv[, U, , drop = FALSE], c(1, 3), all)))
  sT <- sapply(subsets(p), length); sU <- sapply(subsets(q), length)
  c(viol = as.numeric(any(cT >= 1 & cT <= sT - 1) || any(dU >= 1 & dU <= sU - 1)), c_full = cT[length(cT)], d_full = dU[length(dU)],
    min_slack_T = min((cT - sT)[cT >= 1]), min_slack_U = min((dU - sU)[dU >= 1])) }
nulldir <- function(f, X) { p <- dim(X)[1]; q <- dim(X)[2]; Y <- sweep(X, 1:2, f$M); out <- c(fac = NA, supp = NA, orth = NA, n_T = NA)
  er <- eigen(f$Sr, symmetric = TRUE); ec <- eigen(f$Sc, symmetric = TRUE); rr <- min(er$values) / max(er$values); rc <- min(ec$values) / max(ec$values)
  if (rr <= rc) { a <- er$vectors[, p]; T <- which(abs(a) > 1e-3 * max(abs(a))); a <- a / sqrt(sum(a^2)); cols <- list()
    for (i in which(f$v)) for (l in 1:q) if (all(f$W[T, l, i])) cols[[length(cols) + 1]] <- Y[T, l, i]
    out <- c(fac = 1, supp = length(T), orth = if (length(cols)) max(sapply(cols, function(x) abs(sum(x * a[T])) / sqrt(sum(x^2)))) else NA, n_T = length(cols)) }
  else { b <- ec$vectors[, q]; U <- which(abs(b) > 1e-3 * max(abs(b))); b <- b / sqrt(sum(b^2)); rows <- list()
    for (i in which(f$v)) for (k in 1:p) if (all(f$W[k, U, i])) rows[[length(rows) + 1]] <- Y[k, U, i]
    out <- c(fac = 2, supp = length(U), orth = if (length(rows)) max(sapply(rows, function(x) abs(sum(x * b[U])) / sqrt(sum(x^2)))) else NA, n_T = length(rows)) }
  c(out, evr = min(rr, rc)) }
one <- function(seed, n, eps) { set.seed(seed); p <- q <- 5; S <- true_sigma(p, q, "full"); X <- rmatnorm(n, Sr = S$Sr, Sc = S$Sc)
  if (eps > 0) { bad <- array(runif(p * q * n) < eps, c(p, q, n)); sdc <- sqrt(outer(diag(S$Sr), diag(S$Sc))); X[bad] <- X[bad] + 6 * sdc[array(1:(p * q), c(p, q, n))[bad]] * sign(rnorm(sum(bad))) }
  f <- tryCatch(suppressWarnings(cellmmcd(X, maxit = 300, reweight = FALSE)), error = function(e) NULL)
  if (is.null(f)) return(c(err = 1, degen = NA, viol = NA, c_full = NA, d_full = NA, min_slack_T = NA, min_slack_U = NA, fac = NA, supp = NA, orth = NA, n_T = NA, evr = NA))
  nd <- nulldir(f, X); degen <- as.numeric(nd["evr"] < 1e-6 || f$obj[length(f$obj)] < f$obj[1] - 50 * n)
  c(err = 0, degen = degen, criterion(f$W, f$v), nd) }
grid <- rbind(expand.grid(n = c(8, 9, 10, 12), eps = 0), data.frame(n = 100, eps = c(0, 0.2)))
sink("output/cm_step23_degen.log"); cat("cellMMCD step 23 (block 7): configuration criterion on the raw fit, p = q = 5, alpha = .75, 8 seeds per row.\n",
  "degen = eigenvalue ratio < 1e-6 or F dropped by 50 n; viol = criterion of Prop. A violated (some 1 <= c_T <= |T|-1 or 1 <= d_U <= |U|-1);\n",
  "c_full / d_full = fully kept columns / rows; slack = min over T with c_T >= 1 of c_T - |T| (negative = violation); orth = max |x'a|/|x| over kept T-columns, T = support of the null direction a.\n")
for (g in seq_len(nrow(grid))) { R <- if (grid$n[g] == 100) 4 else 8
  r <- do.call(rbind, mclapply(1:R, function(s) one(23000 + 10 * g + s, grid$n[g], grid$eps[g]), mc.cores = 4))
  cat(sprintf("\nn = %3d eps %.1f: errors %d\n", grid$n[g], grid$eps[g], sum(r[, "err"]))); print(round(r[, -1], 4)) }
sink(); cat("done\n")
