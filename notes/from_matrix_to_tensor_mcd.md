# From Matrix-Variate to Tensor-Variate Data: Mahalanobis Distance, Flip-Flop and MCD

*A study note built mainly on*
- **Obukhov (2025)**, *Robust covariance estimation for tensor-valued observations*, Master's thesis, TU Wien (`references/obukhov2025.pdf`), which introduces **TMCD**, and
- **Mayrhofer, Radojičić & Filzmoser (2025)**, *Robust covariance estimation and explainable outlier detection for matrix-valued data*, Technometrics (`references/mayrhofer2025a.pdf`), which introduces **MMCD**.

Supporting references: Ohlson, Ahmad & von Rosen (2013) on the multilinear normal (`ohlson2013.pdf`), Manceur & Dutilleul (2013) for the tensor flip-flop, Kolda & Bader (2009) for tensor algebra.

---

## 0. The one idea to take away

> **A matrix-variate observation has 2 "directions" (modes), each with its own covariance matrix. A tensor-variate observation of order $K$ has $K$ modes, each with its own covariance matrix. Every matrix formula that contains "$X$" and "$X^\top$" is really a formula about "mode 1" and "mode 2". Rewrite it in terms of modes (unfoldings and mode-wise products) and it generalises to any $K$ with essentially no new ideas.**

In the matrix case you have formulas such as

$$
\Sigma_1 \propto \sum_i (X_i - M)\,\Sigma_2^{-1}(X_i-M)^\top, \qquad
\Sigma_2 \propto \sum_i (X_i - M)^\top\Sigma_1^{-1}(X_i-M).
$$

The transpose in the second formula is not a special trick. It is the **mode-2 unfolding** of $X$, the matrix whose rows are indexed by the 2nd mode. Once you see $X = X_{(1)}$ and $X^\top = X_{(2)}$, the tensor version follows directly:

$$
\Sigma_k \propto \sum_i (X_i - M)_{(k)}\,\Big(\textstyle\bigotimes_{j\neq k}\Sigma_j\Big)^{-1}(X_i-M)_{(k)}^\top, \qquad k = 1,\dots,K .
$$

The rest of this note makes that precise and walks through each object: distribution, Mahalanobis distance, MLE, MCD.

---

## 1. What is a "tensor" here?

Forget the physics/differential-geometry meaning. In statistics and data science a tensor is simply a **multi-dimensional array**:

| order $K$ | object | one entry | example |
|---|---|---|---|
| 1 | vector $x\in\mathbb R^{p}$ | $x[i]$ | a patient's lab values |
| 2 | matrix $X\in\mathbb R^{p_1\times p_2}$ | $X[i,j]$ | grayscale image, variables × time |
| 3 | tensor $\mathcal X\in\mathbb R^{p_1\times p_2\times p_3}$ | $\mathcal X[i,j,l]$ | RGB image (row × col × channel), video (row × col × time), sensor × variable × time |
| $K$ | $\mathcal X\in\mathbb R^{p_1\times\cdots\times p_K}$ | $\mathcal X[i_1,\dots,i_K]$ | e.g. RGB video (row × col × channel × time) |

Vocabulary:

- **Mode** $k$: the $k$-th index or direction. A matrix has 2 modes (rows, columns), an RGB image has 3.
- **Fibre**: fix all indices except one. In a matrix the mode-1 fibres are the columns and the mode-2 fibres are the rows.
- **Slice**: fix one index. For an RGB image, fixing the channel gives three $p_1\times p_2$ matrices.
- Shorthand used throughout:
  $$p = \prod_{k=1}^K p_k \quad(\text{total number of cells}),\qquad p_{-k} = p/p_k \quad(\text{number of cells per mode-}k\text{ fibre's "complement"}).$$

A **sample** is $\mathcal X_1,\dots,\mathcal X_n$. Obukhov stores it as one $(K{+}1)$-way array of size $n\times p_1\times\cdots\times p_K$, which matters a lot for speed (§7).

---

## 2. Three tools that turn tensor algebra into matrix algebra

All tensor computations below reduce to ordinary matrix products through three operations.

### 2.1 Unfolding (matricisation) $\mathcal X_{(k)}$

The **mode-$k$ unfolding** $\mathcal X_{(k)}\in\mathbb R^{p_k\times p_{-k}}$ lays all mode-$k$ fibres side by side as columns. Row index = the $k$-th index. Column index = all other indices, with the **earliest remaining index running fastest** (column-major, which is R's native order; Kolda & Bader 2009, Obukhov §2.2).

Example for $p_1\times p_2\times p_3 = 3\times4\times2$:
- $\mathcal X_{(1)}$ is $3\times 8$. Its columns are $(\cdot,1,1),(\cdot,2,1),\dots,(\cdot,4,1),(\cdot,1,2),\dots$
- $\mathcal X_{(2)}$ is $4\times 6$. Its columns are $(1,\cdot,1),(2,\cdot,1),(3,\cdot,1),(1,\cdot,2),\dots$
- $\mathcal X_{(3)}$ is $2\times 12$.

> **Matrix special case:** $X_{(1)} = X$ and $X_{(2)} = X^\top$. Every transpose in the matrix-variate literature is an unfolding in disguise.

In R this is `matrix(aperm(A, c(k, setdiff(1:K, k))), dim(A)[k])`.

### 2.2 Mode-$k$ product $\mathcal X\times_k A$

For $A\in\mathbb R^{q\times p_k}$, the mode-$k$ product multiplies **every mode-$k$ fibre** by $A$:

$$
(\mathcal X\times_k A)_{(k)} = A\,\mathcal X_{(k)} .
$$

Intuition: "apply the linear map $A$ along direction $k$ and leave the other directions alone." Products along different modes commute, and you can chain them: $\mathcal X\times_1 A_1\times_2 A_2\times_3 A_3$.

> **Matrix special case:** $X\times_1 A\times_2 B = A\,X\,B^\top$. Left multiplication acts on mode 1 (rows), right multiplication acts on mode 2 (columns).

### 2.3 Vectorisation and the Kronecker product

$\operatorname{vec}(\mathcal X) = \operatorname{vec}(\mathcal X_{(1)})$ stacks all cells with $i_1$ fastest, then $i_2$, and so on. The bridge identity is

$$
\boxed{\;\operatorname{vec}\big(\mathcal X\times_1 A_1\times_2\cdots\times_K A_K\big) = (A_K\otimes\cdots\otimes A_2\otimes A_1)\,\operatorname{vec}(\mathcal X)\;}
$$

and its unfolded form (Kolda & Bader 2009, Prop. 3.7)

$$
\big(\mathcal X\times_1 A_1\cdots\times_K A_K\big)_{(k)} = A_k\,\mathcal X_{(k)}\,\Big(A_K\otimes\cdots\otimes A_{k+1}\otimes A_{k-1}\otimes\cdots\otimes A_1\Big)^\top .
$$

> **Matrix special case:** $\operatorname{vec}(AXB^\top) = (B\otimes A)\operatorname{vec}(X)$, the identity behind $\operatorname{vec}(X)\sim N(\operatorname{vec} M, \Sigma_2\otimes\Sigma_1)$.

**Why the reversed order?** Column-major vec makes index 1 run fastest. In a Kronecker product $B\otimes A$ the right factor varies fastest, so the factor of mode 1 must sit on the right.

Notation for the "everything except mode $k$" Kronecker product, used constantly below:

$$
\Sigma_{-k} := \Sigma_K\otimes\cdots\otimes\Sigma_{k+1}\otimes\Sigma_{k-1}\otimes\cdots\otimes\Sigma_1 \in\mathbb R^{p_{-k}\times p_{-k}} .
$$

Kronecker facts needed (Obukhov Table 2.2):
$(A\otimes B)^{-1} = A^{-1}\otimes B^{-1}$, $(A\otimes B)^{1/2} = A^{1/2}\otimes B^{1/2}$, $(A\otimes B)(C\otimes D) = AC\otimes BD$, and
$$
\det(\Sigma_K\otimes\cdots\otimes\Sigma_1) = \prod_{k=1}^K \det(\Sigma_k)^{p_{-k}}
\quad\Longrightarrow\quad
\log\det(\textstyle\bigotimes\Sigma) = \sum_k p_{-k}\log\det\Sigma_k .
$$
For $K=2$ this is the familiar $\det(\Sigma_2\otimes\Sigma_1) = \det(\Sigma_1)^{p_2}\det(\Sigma_2)^{p_1}$ from the MMCD objective.

---

## 3. The tensor normal distribution

### 3.1 Definition

$$
\mathcal X\sim TN_{p_1,\dots,p_K}(\mathcal M;\Sigma_1,\dots,\Sigma_K)
\iff
\operatorname{vec}(\mathcal X)\sim N_p\big(\operatorname{vec}\mathcal M,\ \Sigma_K\otimes\cdots\otimes\Sigma_1\big).
$$

Matrix case ($K=2$): $\operatorname{vec}(X)\sim N_{p_1p_2}(\operatorname{vec}M,\Sigma_2\otimes\Sigma_1)$ with $\Sigma_1=\Sigma_{\text{row}}$, $\Sigma_2=\Sigma_{\text{col}}$.

### 3.2 Generative picture (best intuition)

Start from i.i.d. $N(0,1)$ noise $\mathcal Z$ and "colour" it one mode at a time:

$$
\mathcal X = \mathcal M + \mathcal Z\times_1\Sigma_1^{1/2}\times_2\Sigma_2^{1/2}\cdots\times_K\Sigma_K^{1/2}.
$$

For a matrix this is $X = M + \Sigma_1^{1/2}Z\,\Sigma_2^{1/2}$: correlate the rows, then correlate the columns. For an RGB image, additionally correlate the three colour channels through a $3\times 3$ matrix $\Sigma_3$. Obukhov's sampler (§3.1) does exactly this, one reshape + matrix product per mode, without ever forming a Kronecker product.

### 3.3 What "separable" means cell by cell

$$
\operatorname{Cov}\big(\mathcal X[i_1,\dots,i_K],\ \mathcal X[j_1,\dots,j_K]\big) = \prod_{k=1}^K \Sigma_k[i_k,j_k].
$$

The covariance between two cells is the **product of per-mode covariances**. In an RGB video, the covariance between pixel $(i,j)$ in red at time $t$ and pixel $(i',j')$ in green at time $t'$ is (row-cov) × (col-cov) × (channel-cov) × (time-cov). That is a strong structural assumption, and it is the reason estimation is feasible at all.

### 3.4 Parameter count and sample size

| model | covariance parameters | RGB frame $128\times160\times3$ |
|---|---|---|
| unstructured on $\operatorname{vec}$ | $p(p+1)/2$ | $\approx 1.9\times10^{9}$ |
| separable (tensor normal) | $\sum_k p_k(p_k+1)/2 - (K-1)$ | $8256+12880+6-2 \approx 2.1\times10^{4}$ |

Minimum sample size for a positive-definite estimate:
- vector: $n > p$, so $n>61440$ in the example;
- matrix (Soloveychik & Trushin 2016, used in MMCD): $n > \frac{p_1}{p_2}+\frac{p_2}{p_1}+1$ (any $n\ge 4$ for square matrices);
- tensor: the literature gives several bounds (Obukhov eqs. 1.1–1.3). The general pattern is
  $$n \gtrsim \sum_k \frac{p_k}{p_{-k}} = \sum_k\frac{p_k^2}{p},$$
  and Obukhov uses the conservative elemental subset size $s = \big\lceil \sum_k p_k^2/p\big\rceil + 2$.

  Matrix check: $p_1^2/(p_1p_2)+p_2^2/(p_1p_2) = p_1/p_2+p_2/p_1$, which recovers the MMCD bound.

**Intuition:** each $\Sigma_k$ is estimated from all $n\,p_{-k}$ mode-$k$ fibres pooled together, not from $n$ observations. More modes means more fibres per observation, so you need fewer observations.

### 3.5 Identifiability

$(a_1\Sigma_1,\dots,a_K\Sigma_K)$ gives the same distribution whenever $\prod_k a_k = 1$. A matrix has 1 free scale. An order-$K$ tensor has $K-1$ free scales. Only the Kronecker product is identified, so:
- normalise, e.g. $\Sigma_k[1,1]=1$ for $k<K$ and push the leftover scale into $\Sigma_K$ (Obukhov §3.4);
- consistency and rescaling factors are applied to **one** factor only (Obukhov rescales $\Sigma_3$);
- compare estimates through scale-invariant quantities such as the Kronecker product, the KL divergence, or normalised matrices.

---

## 4. The tensor Mahalanobis distance

### 4.1 Three equivalent forms

Let $\mathcal C = \mathcal X-\mathcal M$. The squared tensor Mahalanobis distance is the ordinary Mahalanobis distance of $\operatorname{vec}\mathcal X$:

$$
\mathrm{TMD}^2(\mathcal X)
= \underbrace{\operatorname{vec}(\mathcal C)^\top(\Sigma_K\otimes\cdots\otimes\Sigma_1)^{-1}\operatorname{vec}(\mathcal C)}_{\text{(a) vectorised definition}}
= \underbrace{\big\lVert\ \mathcal C\times_1\Sigma_1^{-1/2}\times_2\cdots\times_K\Sigma_K^{-1/2}\big\rVert_F^2}_{\text{(b) whiten every mode, then Euclidean norm}}
= \underbrace{\operatorname{tr}\!\big(\Sigma_k^{-1}\,\mathcal C_{(k)}\,\Sigma_{-k}^{-1}\,\mathcal C_{(k)}^\top\big)}_{\text{(c) mode-}k\text{ trace form, any }k}.
$$

Matrix case: (b) is $\lVert\Sigma_1^{-1/2}C\,\Sigma_2^{-1/2}\rVert_F^2$, and (c) with $k=1$ is $\operatorname{tr}(\Sigma_1^{-1}C\,\Sigma_2^{-1}C^\top)$, which is MMCD's eq. (7) $\operatorname{tr}(\Omega_{\text{col}}C^\top\Omega_{\text{row}}C)$.

**Proof sketch.** (a)⇔(b): $(\bigotimes\Sigma)^{-1/2} = \Sigma_K^{-1/2}\otimes\cdots\otimes\Sigma_1^{-1/2}$, so by the bridge identity of §2.3 the whitened vector is the vec of the whitened tensor. (b)⇔(c): whiten all modes except $k$ and call the result $\mathcal W$. By the unfolded identity, $\mathcal W_{(k)} = \mathcal C_{(k)}\Sigma_{-k}^{-1/2}$ (symmetric roots). Then (b) $=\lVert\Sigma_k^{-1/2}\mathcal W_{(k)}\rVert_F^2 = \operatorname{tr}(\Sigma_k^{-1}\mathcal W_{(k)}\mathcal W_{(k)}^\top)$. ∎

**Intuition.** Form (b) is the one to picture. Whitening along each mode removes the correlation within that direction. After whitening all $K$ directions, every cell is an independent $N(0,1)$ variable under the model, and the distance is just the sum of squares. Form (c) is the one that drives the algebra of the MLE.

### 4.2 Distribution and cutoff

Under the model with known parameters, $\mathrm{TMD}^2\sim\chi^2_{p}$ with $p=\prod_k p_k$, exactly as for vectors and matrices. Outliers are flagged when $\mathrm{TMD}^2 > \chi^2_{p,0.975}$.

Caveat (Obukhov §5.2): $p$ is typically huge (61 440 for one RGB frame). The $\chi^2_p$ cutoff then sits very close to its mean, so mild non-normality alone can push many observations over the line. Obukhov's video examples show the cutoff sometimes needs manual or relative tuning.

### 4.3 Cellwise contributions: two different decompositions

There are two natural ways to split $\mathrm{TMD}^2$ into per-cell pieces.

1. **Squared whitened cells** (Obukhov §3.2): $\mathcal D = \mathcal W\odot\mathcal W$ with $\mathcal W$ from form (b). These sum to $\mathrm{TMD}^2$, but they live in the **whitened coordinate system**. Because whitening mixes cells along each mode, cell $(i,j,l)$ of $\mathcal W$ is not purely about cell $(i,j,l)$ of $\mathcal X$.
2. **Shapley values** (Mayrhofer & Filzmoser 2023; MMCD §5): for vectors $\phi_j = c_j\,(\Sigma^{-1}c)_j$. MMCD writes the matrix version as $\Phi(X) = C\odot\Sigma_1^{-1}C\,\Sigma_2^{-1}$. The direct tensor extension, obtained by applying the vector formula to $\operatorname{vec}$ with the bridge identity, is
   $$
   \Phi(\mathcal X) = \mathcal C\ \odot\ \big(\mathcal C\times_1\Sigma_1^{-1}\times_2\cdots\times_K\Sigma_K^{-1}\big).
   $$
   These also sum to $\mathrm{TMD}^2$, are tied to the **original** cells, and inherit the Shapley axioms. Row, column, channel or frame contributions are obtained by summing over the other modes.

Obukhov says he "followed the approach of [MRF25]", but the formula he implements is (1), not (2). They are different numbers; in the R check (§10) the maximum cellwise difference is about $0.025$ when the total distance is $0.76$. If you want explanations comparable to MMCD's Shapley values, use (2). That formula is my own extension and is not stated in either source, though it is verified numerically below.

---

## 5. Maximum likelihood: the flip-flop algorithm

### 5.1 Log-likelihood

$$
\ell = -\frac{np}{2}\log 2\pi \;-\;\frac n2\sum_{k=1}^K p_{-k}\log\det\Sigma_k \;-\;\frac12\sum_{i=1}^n \mathrm{TMD}^2(\mathcal X_i).
$$

### 5.2 Deriving the update

Fix all $\Sigma_j$ with $j\neq k$. With form (c) of the distance:

$$
\ell(\Sigma_k) = \text{const} - \frac{n\,p_{-k}}{2}\log\det\Sigma_k - \frac12\operatorname{tr}\!\big(\Sigma_k^{-1}S_k\big),
\qquad S_k = \sum_{i=1}^n \mathcal C_{i(k)}\,\Sigma_{-k}^{-1}\,\mathcal C_{i(k)}^\top .
$$

This has the shape of a Gaussian log-likelihood for a $p_k$-dimensional covariance with "sample size" $n p_{-k}$ and scatter $S_k$. Its maximiser is

$$
\boxed{\;\hat\Sigma_k = \frac{1}{n\,p_{-k}}\sum_{i=1}^n \hat{\mathcal C}_{i(k)}\,\hat\Sigma_{-k}^{-1}\,\hat{\mathcal C}_{i(k)}^\top\;}\qquad(\hat{\mathcal M} = \bar{\mathcal X}).
$$

- $K=2$, $k=1$: $\hat\Sigma_1 = \frac1{np_2}\sum_i C_i\hat\Sigma_2^{-1}C_i^\top$, Dutilleul's row update.
- $K=2$, $k=2$: $\hat\Sigma_2 = \frac1{np_1}\sum_i C_i^\top\hat\Sigma_1^{-1}C_i$, the column update, because $C_{(2)} = C^\top$ and $\Sigma_{-2}=\Sigma_1$.
- $K=3$: Obukhov's eqs. (2.13)–(2.15), first given by Manceur & Dutilleul (2013).

**Intuition.** To estimate the covariance along mode $k$, first whiten the data along all other modes using the current estimates. Then treat every mode-$k$ fibre of every observation as one $p_k$-dimensional vector, giving $n\,p_{-k}$ pseudo-observations, and take their ordinary sample covariance. The algorithm cycles $k=1,\dots,K$ until convergence ("flip-flop" becomes "flip-flop-flop…").

In code: $\hat\Sigma_k = \frac{1}{np_{-k}}\,\mathbf W_{(k)}\mathbf W_{(k)}^\top$, where $\mathbf W$ is the whole $n\times p_1\times\cdots\times p_K$ sample array whitened along all modes except $k$ and unfolded along $k$. That is one big matrix product with no loop over $i$.

### 5.3 A useful identity

Right after updating any $\Sigma_k$ with the formula above,
$$
\sum_{i=1}^n \mathrm{TMD}^2(\mathcal X_i) = \operatorname{tr}(\hat\Sigma_k^{-1}S_k) = \operatorname{tr}(n\,p_{-k}\,I_{p_k}) = n\,p .
$$
This is the tensor version of the "$\sum_{i\in H} d_i^2 = hpq$" fact in MMCD (proof of their eq. 23), and it is the key to the C-step argument in §6.2. It is verified numerically in §10.

### 5.4 Existence and convergence

Existence and uniqueness (up to the scale ambiguity) of the Kronecker MLE are treated by Drton, Kuriki & Hoff (2021) and Derksen, Makam & Walter (2022). Each flip-flop step maximises $\ell$ in one block, so $\ell$ is non-decreasing. In practice 2–3 sweeps usually suffice (Dutilleul 1999; Werner et al. 2008 show 2 iterations already reach asymptotic efficiency in the matrix case).

---

## 6. From MMCD to TMCD

### 6.1 The objective

MCD finds the $h$-subset whose sample covariance has the smallest determinant. Equivalently, it maximises the Gaussian likelihood over $h$-subsets (Raymaekers & Rousseeuw 2023). MMCD (Prop. 2.0.1) keeps this likelihood view and plugs in the matrix-normal likelihood. TMCD plugs in the tensor-normal likelihood:

| | objective over $H\subset\{1..n\}$, $\lvert H\rvert=h$ |
|---|---|
| MCD | $\log\det\hat\Sigma_H$ |
| MMCD | $p_2\log\det\hat\Sigma_{1,H} + p_1\log\det\hat\Sigma_{2,H} = \log\det(\hat\Sigma_{2,H}\otimes\hat\Sigma_{1,H})$ |
| TMCD | $\displaystyle\sum_{k=1}^K p_{-k}\log\det\hat\Sigma_{k,H} = \log\det(\hat\Sigma_{K,H}\otimes\cdots\otimes\hat\Sigma_{1,H})$ |

Here $(\hat{\mathcal M}_H,\hat\Sigma_{k,H})$ are the flip-flop MLEs on subset $H$. Obukhov calls the TMCD objective `ld`.

### 6.2 C-steps and why they decrease the objective

The C-step is identical in all three settings:
1. Given a subset $H_\text{old}$, compute MLEs on it (flip-flop for MMCD and TMCD).
2. Compute $\mathrm{TMD}^2$ for **all** $n$ observations using these estimates.
3. Let $H_\text{new}$ be the $h$ observations with the smallest distances.
4. Recompute MLEs on $H_\text{new}$ (flip-flop, warm-started from the old estimates).

The monotonicity argument of MMCD §4.1 carries over unchanged. By §5.3, $\sum_{i\in H_\text{old}}\mathrm{TMD}^2_i = hp$. Re-selecting the $h$ smallest distances can only lower that sum, while the determinant term stays fixed, so the likelihood does not decrease. The flip-flop on $H_\text{new}$ then increases it again. Hence

$$
\sum_k p_{-k}\log\det\hat\Sigma_{k,H_\text{new}} \le \sum_k p_{-k}\log\det\hat\Sigma_{k,H_\text{old}},
$$

and since there are finitely many subsets, the C-steps converge (to a local optimum). Obukhov (§3.5) describes the convergence as heuristic even in the matrix case. The MMCD paper does give this argument (its eq. 23), and it uses only "flip-flop is likelihood-monotone" plus the §5.3 identity, both of which hold for any $K$. The one caveat is that the flip-flop must actually be run to (or near) its fixed point for the bound to be exact.

### 6.3 The full algorithm (Fast-MCD skeleton, Obukhov §3.6)

1. **Elemental subsets.** Draw $m$ random subsets of size $s = \lceil\sum_k p_k^2/p\rceil + 2$. This is tiny: for $30\times30\times30$ data it is $\lceil 3\cdot 900/27000\rceil+2 = 3$. Run the flip-flop on each starting from identities, with only about 2 iterations.
2. **Short C-steps.** Run about 2 C-steps on each candidate.
3. **Full C-steps.** Keep the best ~10 candidates by `ld` and iterate C-steps to convergence.
4. **Raw estimator.** Take the candidate with the smallest `ld`.
5. **Consistency scaling** at the normal model, applied to one factor only (§3.5):
   $$c(\alpha) = \frac{\alpha}{F_{\chi^2_{p+2}}\!\big(\chi^2_{p,\alpha}\big)},\qquad \alpha = h/n,\quad p=\textstyle\prod_k p_k .$$
   This is the Croux–Haesbroeck factor from MCD and MMCD with $p$ replaced by the total number of cells. MMCD (Thm 3.0.2) proves that the MCD and MMCD factors must coincide; for TMCD this is assumed by analogy.
6. **Reweighting.** Flag observations with scaled $\mathrm{TMD}^2>\chi^2_{p,0.975}$. Run the flip-flop on the rest and rescale with $c(\tilde h/n)$. Obukhov also keeps the members of the raw $h$-subset in the clean set, a small deviation from MMCD's pure cutoff rule.

Choosing $h$: MMCD attains its maximal breakdown point at $h=\lfloor(n+d+2)/2\rfloor$ with $d=p_1/p_2+p_2/p_1$. The natural tensor guess is $d=\sum_k p_k^2/p$. Obukhov works with $\alpha=h/n\in[0.5,1)$, e.g. $\alpha\approx0.5$ for up to 50 % contamination.

---

## 7. Computation: never build the Kronecker product

A whitening or colouring step $\operatorname{vec}(\mathcal X)\mapsto(\bigotimes_k A_k)\operatorname{vec}(\mathcal X)$ is computed as $K$ successive mode products. Each one is a reshape plus a single matrix product.

Obukhov's R pattern for the whole sample array $\mathbf C$ of shape $(n,p_1,p_2,p_3)$ (§3.1, §3.7):

```r
dim(C) <- c(n * p1 * p2, p3); C <- C %*% A3     # mode 3: reshape is free (no copy)
dim(C) <- c(n, p1, p2, p3);  C <- aperm(C, c(1, 2, 4, 3))  # aperm = real copy
dim(C) <- c(n * p1 * p3, p2); C <- C %*% A2     # mode 2
# ... same for mode 1
```

Obukhov's "accordion" metaphor: reshaping (squeezing or expanding the accordion) is free, while permuting axes (`aperm`) copies memory and is the expensive operation to minimise.

Cost per observation of computing $\mathrm{TMD}^2$:
- vectorised: $O(p^2)$, plus storing a $p\times p$ matrix;
- mode-wise: $O\big(p\sum_k p_k\big)$.

For a $128\times160\times3$ frame: $p^2\approx3.8\times10^9$ versus $p\sum p_k\approx1.8\times10^7$, roughly 200× less work and no 30 GB matrix. The R package `Tlasso::Trnorm` forms the Kronecker product explicitly and runs out of memory at $30\times40\times50$ (Obukhov Table 4.1).

---

## 8. Matrix → tensor translation table

| concept | matrix ($K=2$) | tensor (order $K$) |
|---|---|---|
| observation | $X\in\mathbb R^{p_1\times p_2}$ | $\mathcal X\in\mathbb R^{p_1\times\cdots\times p_K}$ |
| "rows / columns" | $X$, $X^\top$ | unfoldings $\mathcal X_{(1)},\dots,\mathcal X_{(K)}$ |
| two-sided multiplication | $AXB^\top$ | $\mathcal X\times_1A_1\cdots\times_KA_K$ |
| vec identity | $(B\otimes A)\operatorname{vec}X$ | $(A_K\otimes\cdots\otimes A_1)\operatorname{vec}\mathcal X$ |
| covariance of vec | $\Sigma_2\otimes\Sigma_1$ | $\Sigma_K\otimes\cdots\otimes\Sigma_1$ |
| "other-mode" covariance | $\Sigma_{-1}=\Sigma_2$, $\Sigma_{-2}=\Sigma_1$ | $\Sigma_{-k}=\bigotimes_{j\ne k}\Sigma_j$ (reverse order) |
| cell covariance | $\Sigma_1[i,i']\,\Sigma_2[j,j']$ | $\prod_k\Sigma_k[i_k,i'_k]$ |
| log-determinant | $p_2\log\lvert\Sigma_1\rvert+p_1\log\lvert\Sigma_2\rvert$ | $\sum_k p_{-k}\log\lvert\Sigma_k\rvert$ |
| scale ambiguity | 1 degree of freedom | $K-1$ degrees of freedom |
| squared MD | $\operatorname{tr}(\Sigma_1^{-1}C\Sigma_2^{-1}C^\top)$ | $\operatorname{tr}(\Sigma_k^{-1}\mathcal C_{(k)}\Sigma_{-k}^{-1}\mathcal C_{(k)}^\top)$, any $k$ |
| whitened form | $\lVert\Sigma_1^{-1/2}C\Sigma_2^{-1/2}\rVert_F^2$ | $\lVert\mathcal C\times_1\Sigma_1^{-1/2}\cdots\times_K\Sigma_K^{-1/2}\rVert_F^2$ |
| null distribution | $\chi^2_{p_1p_2}$ | $\chi^2_{\prod p_k}$ |
| flip-flop update | $\frac1{np_2}\sum C_i\Sigma_2^{-1}C_i^\top$ | $\frac1{np_{-k}}\sum\mathcal C_{i(k)}\Sigma_{-k}^{-1}\mathcal C_{i(k)}^\top$ |
| min. (sub)sample size | $p_1/p_2+p_2/p_1+2$ | $\lceil\sum_k p_k^2/p\rceil+2$ (Obukhov) |
| cellwise Shapley | $C\odot\Sigma_1^{-1}C\Sigma_2^{-1}$ | $\mathcal C\odot(\mathcal C\times_1\Sigma_1^{-1}\cdots\times_K\Sigma_K^{-1})$ |
| equivariance | $X\mapsto AXB^\top+D$ | $\mathcal X\mapsto\mathcal X\times_1A_1\cdots\times_KA_K+\mathcal D$ |
| robust estimator | MMCD | TMCD |

---

## 9. What is proven, what is open

| property | MMCD (matrix) | TMCD (tensor), status in Obukhov 2025 |
|---|---|---|
| matrix/tensor affine equivariance | proven (Lemma 3.0.1) | not stated; expected by the same argument, since the flip-flop MLE is equivariant under mode-wise linear maps |
| breakdown point | proven: $\min(n-h+1,\,h-(d+1))/n$, which beats the vectorised MCD | **open**; only checked empirically (contamination up to ~50 %) |
| consistency (elliptical family) | proven (Thm 3.0.2) | **open**; the factor is borrowed from MCD/MMCD |
| C-step monotonicity | proven (eq. 23) | called heuristic in the thesis; the MMCD argument transfers via §5.3 |
| reweighting keeps breakdown | proven (Thm 3.0.3) | open |
| explanations | Shapley values (cell, row, column) | squared whitened cells (not Shapley; see §4.3) |
| implementation | R package `robustmatrix` | R + PyTorch code (GitHub, see thesis ref. [Obu25]) |
| experiments | simulations + real data | simulations, $3\times3\times3$ to $30\times30\times30$; colour surveillance video $633\times128\times160\times3$ (same video as the MMCD grayscale example) |

Each "open" row in the TMCD column is a potential theoretical contribution: lifting the MMCD proofs (Supplement B of Mayrhofer et al.) from 2 to $K$ modes.

---

## 10. Minimal R code to check everything

The following script verifies the vec/Kronecker identity, the three forms of $\mathrm{TMD}^2$, and the tensor Shapley formula. A longer version also checked the mode-2 unfolding identity and ran a 20-sweep flip-flop on $n=50$ simulated $3\times4\times2$ tensors. All checks passed: differences were $\sim10^{-14}$, and after the flip-flop $\sum_i\mathrm{TMD}^2_i = 1200 = np$ exactly.

```r
set.seed(1)
unfold <- function(A, k) { d <- dim(A); perm <- c(k, seq_along(d)[-k]); matrix(aperm(A, perm), d[k]) }
fold   <- function(M, k, d) { perm <- c(k, seq_along(d)[-k]); aperm(array(M, d[perm]), order(perm)) }
mode_prod <- function(A, B, k) { d <- dim(A); d[k] <- nrow(B); fold(B %*% unfold(A, k), k, d) }
msqrt <- function(S, pow = 1/2) { e <- eigen(S, symmetric = TRUE)
  e$vectors %*% diag(e$values^pow, length(e$values)) %*% t(e$vectors) }
rpd <- function(p) { A <- matrix(rnorm(p * p), p); crossprod(A) + diag(p) }

p <- c(3, 4, 2); S <- lapply(p, rpd); A <- array(rnorm(prod(p)), p)

# vec(A x1 S1 x2 S2 x3 S3) == (S3 %x% S2 %x% S1) vec(A)
B <- mode_prod(mode_prod(mode_prod(A, S[[1]], 1), S[[2]], 2), S[[3]], 3)
max(abs(c(B) - kronecker(S[[3]], kronecker(S[[2]], S[[1]])) %*% c(A)))

# TMD^2: Kronecker form, whitened form, mode-1 trace form
Sig <- kronecker(S[[3]], kronecker(S[[2]], S[[1]]))
W <- A; for (k in 1:3) W <- mode_prod(W, msqrt(S[[k]], -1/2), k)
c(drop(t(c(A)) %*% solve(Sig) %*% c(A)),
  sum(W^2),
  sum(diag(solve(S[[1]]) %*% unfold(A, 1) %*% kronecker(solve(S[[3]]), solve(S[[2]])) %*% t(unfold(A, 1)))))

# tensor cellwise Shapley values: sum to TMD^2, differ from W^2
P <- A; for (k in 1:3) P <- mode_prod(P, solve(S[[k]]), k)
phi <- A * P; c(sum(phi), max(abs(phi - W^2)))
```

---

## 11. Suggested reading order

1. Kolda & Bader (2009), *Tensor decompositions and applications*, SIAM Review, §2 only: unfoldings and $n$-mode products.
2. Ohlson, Ahmad & von Rosen (2013) (`ohlson2013.pdf`): the multilinear normal distribution and its properties.
3. Manceur & Dutilleul (2013): tensor flip-flop MLE and sample-size conditions.
4. Obukhov (2025), Ch. 2–3 (`obukhov2025.pdf`): notation, TMCD algorithm and implementation. The thesis is written for readers new to tensors.
5. Mayrhofer, Radojičić & Filzmoser (2025) (`mayrhofer2025a.pdf`), §3 + Supplement B: the proofs to generalise.
6. Llosa-Vite & Maitra (2024) (`llosavite2024.pdf`): elliptically contoured tensor-variate distributions, the natural setting for a TMCD consistency proof.
