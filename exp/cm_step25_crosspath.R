# Block 8, step 25: Conjecture C' (JOURNAL block 8, result 1) beyond the transversal case. p = q = 3, n = 3, pattern:
# matrix 1 keeps column 1 without (1,1) and row 1 without (1,1) (a cross with the crossing cell missing); matrix 2 keeps
# (2,1),(3,1); matrix 3 keeps (1,2),(1,3). Every c_T, d_U is 0 or >= |T| (Prop. A silent, cm_check_pattern viol 0), but
# the cross vector w = e_1 (x) a + beta b (x) e_1 with a_1 + beta b_1 = 0 is supported on O_1, and (a, b) can be chosen
# with y_1 _|_ w. Along Sigma^r = I - (1-d) aa', Sigma^c = I - (1-d) bb' C' predicts F = k log d + O(1) -> -inf. Seed 25.
# Output: output/cm_step25_crosspath.log
source("R/load.R"); set.seed(25)
Fobs <- function(X, Sr, Sc, W) { Sig <- Sc %x% Sr; tot <- 0; for (i in seq_len(dim(X)[3])) { o <- which(W[, , i]); x <- as.vector(X[, , i])[o]; So <- Sig[o, o, drop = FALSE]
  tot <- tot + logdet(So) + sum(x * solve(So, x)) + length(o) * log(2 * pi) }; tot }   # direct observed-cell form (the Schur form inverts Omega_MM, singular along the path)
p <- q <- 3; W <- array(FALSE, c(p, q, 3)); W[2:3, 1, 1] <- TRUE; W[1, 2:3, 1] <- TRUE; W[2:3, 1, 2] <- TRUE; W[1, 2:3, 3] <- TRUE
X <- array(rnorm(p * q * 3), c(p, q, 3)); y <- X[, , 1]
a <- c(1, rnorm(2)); b <- c(-1, rnorm(1), 0)                      # beta = 1 after scaling; cancellation a_1 + b_1 = 0 at (1,1)
b[3] <- -(sum(a[2:3] * y[2:3, 1]) + b[2] * y[1, 2]) / y[1, 3]      # y_1 _|_ w: a_2 y_21 + a_3 y_31 + b_2 y_12 + b_3 y_13 = 0
sa <- sqrt(sum(a^2)); sb <- sqrt(sum(b^2)); a <- a / sa; b <- b / sb  # unit vectors; w = sa e_1 (x) a + sb b (x) e_1 still supported on O_1
sink("output/cm_step25_crosspath.log"); cat("cellMMCD step 25 (block 8): cross path on a Prop.-A-silent 3x3 pattern, n = 3\n")
cat("cm_check_pattern:", unlist(cm_check_pattern(W)[1:5]), "\n")
w <- sa * (diag(3)[, 1] %x% a) + sb * (b %x% diag(3)[, 1]); cat("cross vector w (vec order), support on kept cells of matrix 1:", all(abs(w[!as.vector(W[, , 1])]) < 1e-12), "; w'y_1 =", sum(w * as.vector(y)), "\n")
for (d in 10^-(1:8)) { Sr <- diag(3) - (1 - d) * tcrossprod(a); Sc <- diag(3) - (1 - d) * tcrossprod(b)
  f <- sapply(1:3, function(i) tryCatch(Fobs(X[, , i, drop = FALSE], Sr, Sc, W[, , i, drop = FALSE]), error = function(e) NA))
  So <- (Sc %x% Sr)[as.vector(W[, , 1]), as.vector(W[, , 1])]; ev <- eigen(So, symmetric = TRUE)$values
  cat(sprintf("  d = %.0e: F = %9.3f  (matrix 1 %8.3f, 2 %7.3f, 3 %7.3f); log d = %7.3f; eigenvalues of Sigma_O1O1: %s\n", d, sum(f), f[1], f[2], f[3], log(d), paste(signif(ev, 3), collapse = " "))) }
cat("Control: same path with a generic (a, b) (no orthogonality):\n"); a2 <- rnorm(3); a2 <- a2 / sqrt(sum(a2^2)); b2 <- rnorm(3); b2 <- b2 / sqrt(sum(b2^2))
for (d in 10^-c(2, 4, 6)) { Sr <- diag(3) - (1 - d) * tcrossprod(a2); Sc <- diag(3) - (1 - d) * tcrossprod(b2)
  cat(sprintf("  d = %.0e: F = %9.3f\n", d, Fobs(X, Sr, Sc, W))) }
sink(); cat("done\n")
