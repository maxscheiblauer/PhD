# True parameters and contamination schemes.

ar_cov <- function(m, rho) rho^abs(outer(seq_len(m), seq_len(m), "-"))

# Sigma^r = AR(0.5). Sigma^c either "full" = AR(0.7) or "diag" = diag(1, ..., 2).
true_sigma <- function(p, q, sc = c("full", "diag")) {
  sc <- match.arg(sc)
  list(Sr = ar_cov(p, 0.5),
       Sc = if (sc == "full") ar_cov(q, 0.7) else diag(seq(1, 2, length.out = q), q))
}

# Draw clean data and contaminate. eps is the fraction of the natural units of each type:
#   obs:       round(eps n) whole X_i, mean shift delta in every cell
#   col:       round(eps n q) original columns X_i[, j] := delta + N(0, 1)
#   row:       round(eps n p) original rows    X_i[k, ] := delta + N(0, 1)
#   cell:      round(eps n p q) cells          X_i[k, j] := delta + N(0, 1)
#   wcol:      round(eps n q) whitened columns z_ij := z_ij + delta 1 (true symmetric root)
# truth: n x q logical, TRUE where an original column (wcol: whitened column) is contaminated.
gen_data <- function(n, Sr, Sc, type = c("clean", "obs", "col", "row", "cell", "wcol"),
                     eps = 0, delta = 5) {
  type <- match.arg(type)
  p <- nrow(Sr); q <- nrow(Sc)
  truth <- matrix(FALSE, n, q)
  if (type == "wcol") {
    Z <- rmatnorm(n, Sr = Sr, Sc = diag(q))
    m <- round(eps * n * q); idx <- sample(n * q, m)
    truth[idx] <- TRUE
    for (ij in idx) {
      i <- (ij - 1) %% n + 1; j <- (ij - 1) %/% n + 1
      Z[, j, i] <- Z[, j, i] + delta
    }
    Sch <- solve(root(Sc, "sym"))           # (A^c)^{-1} = Sc^{1/2}
    X <- array(apply(Z, 3, function(z) z %*% Sch), c(p, q, n))
    return(list(X = X, truth = truth))
  }
  X <- rmatnorm(n, Sr = Sr, Sc = Sc)
  if (type == "clean" || eps == 0) return(list(X = X, truth = truth))
  if (type == "obs") {
    idx <- sample(n, round(eps * n))
    X[, , idx] <- X[, , idx] + delta
    truth[idx, ] <- TRUE
  } else if (type == "col") {
    idx <- sample(n * q, round(eps * n * q))
    for (ij in idx) {
      i <- (ij - 1) %% n + 1; j <- (ij - 1) %/% n + 1
      X[, j, i] <- delta + rnorm(p); truth[i, j] <- TRUE
    }
  } else if (type == "row") {
    idx <- sample(n * p, round(eps * n * p))
    for (ik in idx) {
      i <- (ik - 1) %% n + 1; k <- (ik - 1) %/% n + 1
      X[k, , i] <- delta + rnorm(q); truth[i, ] <- TRUE
    }
  } else if (type == "cell") {
    idx <- sample(n * p * q, round(eps * n * p * q))
    for (c0 in idx - 1) {
      k <- c0 %% p + 1; j <- (c0 %/% p) %% q + 1; i <- c0 %/% (p * q) + 1
      X[k, j, i] <- delta + rnorm(1); truth[i, j] <- TRUE
    }
  }
  list(X = X, truth = truth)
}
