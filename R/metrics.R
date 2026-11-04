# Metrics. All comparisons are on the Kronecker product or on shape-normalized factors.

# KL( N(0, Sc_hat %x% Sr_hat) || N(0, Sc %x% Sr) ), using the Kronecker structure.
kl_mn <- function(Sr_hat, Sc_hat, Sr, Sc) {
  p <- nrow(Sr); q <- nrow(Sc)
  Mr <- solve(Sr, Sr_hat); Mc <- solve(Sc, Sc_hat)
  0.5 * (sum(diag(Mr)) * sum(diag(Mc)) - p * q - q * logdet(Mr) - p * logdet(Mc))
}

# min_c KL( N(0, c Sigma_hat) || N(0, Sigma) ): shape-only KL, free of the scale/kappa choice.
kl_shape <- function(Sr_hat, Sc_hat, Sr, Sc) {
  p <- nrow(Sr); q <- nrow(Sc)
  Mr <- solve(Sr, Sr_hat); Mc <- solve(Sc, Sc_hat)
  0.5 * (p * q * log(sum(diag(Mr)) * sum(diag(Mc)) / (p * q)) - q * logdet(Mr) - p * logdet(Mc))
}

# Symmetrized KL between two fits (lists with Sr, Sc).
skl <- function(f1, f2) 0.5 * (kl_mn(f1$Sr, f1$Sc, f2$Sr, f2$Sc) + kl_mn(f2$Sr, f2$Sc, f1$Sr, f1$Sc))

shape <- function(S) S / exp(logdet(S) / nrow(S))

# Relative Frobenius error of shape-normalized matrices.
shape_err <- function(S_hat, S) norm(shape(S_hat) - shape(S), "F") / norm(shape(S), "F")

# tr(Sigma^{-1} Sigma_hat) / (pq): scale ratio of the product.
tr_ratio <- function(Sr_hat, Sc_hat, Sr, Sc)
  sum(diag(solve(Sr, Sr_hat))) * sum(diag(solve(Sc, Sc_hat))) / (nrow(Sr) * nrow(Sc))

# Flagged whitened columns (a == 0) vs. truly contaminated columns (n x q logical).
# Only meaningful when Sigma^c is diagonal (then z_ij is a rescaled original column).
prec_recall <- function(a, truth) {
  fl <- !a
  tp <- sum(fl & truth)
  c(precision = if (sum(fl)) tp / sum(fl) else NA, recall = if (sum(truth)) tp / sum(truth) else NA)
}
