## Numerical check of the Kronecker cellwise identities in notes/00_setting_matrix_cellwise.md (Sec. 4). Base R only.
set.seed(1)
p <- 4; q <- 3
rpd <- function(d){A <- matrix(rnorm(d*d), d); crossprod(A) + diag(d)}
Sr <- rpd(p); Sc <- rpd(q); Or <- solve(Sr); Oc <- solve(Sc)
M <- matrix(rnorm(p*q), p); X <- matrix(rnorm(p*q), p) * 3
S <- kronecker(Sc, Sr); Q <- solve(S); x <- c(X); mu <- c(M)
R <- Or %*% (X - M) %*% Oc            # = matrix form of Q (x - mu)
## (1) conditional mean/variance of a single cell given all other cells
j <- 2; k <- 3; a <- (k-1)*p + j
cm_vec <- mu[a] + S[a,-a] %*% solve(S[-a,-a], x[-a] - mu[-a])
cv_vec <- S[a,a] - S[a,-a] %*% solve(S[-a,-a], S[-a,a])
cat("(1) cond var  :", cv_vec, " vs 1/(or_jj*oc_kk) =", 1/(Or[j,j]*Oc[k,k]), "\n")
cat("    cond mean :", cm_vec, " vs x - R_jk/(or_jj*oc_kk) =", X[j,k] - R[j,k]/(Or[j,j]*Oc[k,k]), "\n")
## (2) Shapley = (x-m) * R, and standardized cellwise residual
phi <- (X - M) * R
r <- R / sqrt(outer(diag(Or), diag(Oc)))
cat("(2) sum(phi) =", sum(phi), " MMD^2 =", sum(diag(Oc %*% t(X-M) %*% Or %*% (X-M))), "\n")
cat("    phi_jk = (x-m)*r*sqrt(w): ", all.equal(phi, (X-M) * r * sqrt(outer(diag(Or), diag(Oc)))), "\n")
## (3) several flagged cells m: E-step via small |m|x|m| precision block
m <- c(2, 7, 11); o <- setdiff(1:(p*q), m)
cm_full <- mu[m] + S[m,o] %*% solve(S[o,o], x[o] - mu[o])
Xt <- X; Xt[m] <- M[m]                               # cells set to mean
Rt <- c(Or %*% (Xt - M) %*% Oc)
cm_fast <- mu[m] - solve(Q[m,m], Rt[m])              # Q[m,m] = kron entries, cheap
cat("(3) E-step fast == full:", all.equal(c(cm_full), c(cm_fast)), "\n")
Qmm <- outer(m, m, function(a,b){ja<-(a-1)%%p+1;ka<-(a-1)%/%p+1;jb<-(b-1)%%p+1;kb<-(b-1)%/%p+1; Or[cbind(ja,jb)]*Oc[cbind(ka,kb)]})
cat("    Q_mm from factors     :", all.equal(Qmm, Q[m,m]), "\n")
## (4) log det of observed block = q log|Sr| + p log|Sc| + log|Q_mm|
cat("(4) logdet S_oo =", determinant(S[o,o])$modulus, " vs ",
    q*determinant(Sr)$modulus + p*determinant(Sc)$modulus + determinant(Qmm)$modulus, "\n")
## (5) observed-part MD^2 = MD^2 of conditionally imputed matrix
Xi <- X; Xi[m] <- cm_fast
md_obs <- t(x[o]-mu[o]) %*% solve(S[o,o], x[o]-mu[o])
cat("(5) MD^2 observed =", md_obs, " vs MMD^2(imputed) =", sum(diag(Oc %*% t(Xi-M) %*% Or %*% (Xi-M))), "\n")
## (6) motivation: fraction of fully clean matrices under ICM with cell rate eps
for (pq in c(25, 100, 400, 10000)) cat(sprintf("(6) pq=%5d  eps=1%%: P(clean)=%.3f   eps=0.1%%: %.3f\n", pq, 0.99^pq, 0.999^pq))
