# Setting: cellwise outliers in matrix-variate data

Working dossier, 2026-11-04. Scope: **matrix-variate only** (the tensor and functional parts of the proposal are out of scope).

> **Source caveat.** The cloud session that wrote this could not reach arXiv, CRAN, repositum, or Lirias (network policy).
> Paper content below comes from abstracts and search snippets, plus the **actual source code** of
> `robustmatrix` 0.1.5 and `cellWise` 2.5.7 (cloned from GitHub's CRAN mirror). Claims marked **[verified]** were
> checked numerically (`R/check_kronecker_cellwise_identities.R`). Everything else must be re-checked against full texts.

---

## 1. Notation (follows Mayrhofer, Radojičić & Filzmoser 2025)

- Observations $X_1,\dots,X_n \in \mathbb{R}^{p\times q}$, cells $x_{i,jk}$.
- Matrix normal / elliptical: $X \sim \mathcal{MN}(M, \Sigma^{row}, \Sigma^{col})$ $\iff$ $\mathrm{vec}(X)\sim N(\mathrm{vec} M, \Sigma^{col}\otimes\Sigma^{row})$.
  Precisions $\Omega^{row}=(\Sigma^{row})^{-1}$, $\Omega^{col}=(\Sigma^{col})^{-1}$. Identifiability: the scale is fixed via $\Sigma^{col}_{11}=1$ (the code divides by `cov_col(0,0)`).
- Squared matrix Mahalanobis distance:
  $\mathrm{MMD}^2(X)=\mathrm{tr}\big(\Omega^{col}(X-M)'\Omega^{row}(X-M)\big)$.
- Cellwise Shapley values: $\Phi(X)=(X-M)\circ\Omega^{row}(X-M)\Omega^{col}$, with $\sum_{jk}\phi_{jk}=\mathrm{MMD}^2$.
  Row and column Shapley values: $\mathrm{diag}(\Omega^{row}(X-M)\Omega^{col}(X-M)')$ and the column analogue.
- Define the **precision-weighted residual matrix** $R(X) := \Omega^{row}(X-M)\Omega^{col}$, which is $Q\,\mathrm{vec}(X-M)$ with $Q=\Omega^{col}\otimes\Omega^{row}$.

## 2. What MMCD actually does (from the `robustmatrix` source)

- **Objective:** choose $h$ of $n$ matrices that minimize $q\log|\Sigma^{row}|+p\log|\Sigma^{col}|$ of the flip-flop MLE on the subset.
- **C-step:** flip-flop MLE on the h-subset → $\mathrm{MMD}^2$ for all $i$ → keep the $h$ smallest → repeat until the log-det change is below 1e-3.
  Initial subsets start at size $h_0=\lceil p/q+q/p\rceil+2$ and grow until the MLE exists. `nsamp=500` starts, split into subpopulations for large $n$ (FastMCD style).
- **Minimum sample size:** existence needs $n \ge \lceil\max(p/q,q/p)\rceil+1$, and uniqueness needs $\lceil p/q+q/p\rceil+2$. This is why MMCD works with very few observations when $p\approx q$.
- `adapt_alpha=TRUE` maps $\alpha$ to the $h$ with maximal breakdown, $h_{bp}=\lfloor (n+h_{smallest})/2\rfloor$.
- **Consistency:** factor $\alpha / F_{\Gamma(pq/2+1)}(\chi^2_{pq,\alpha}/2)$, or a median rule (`mmd_med`).
  Reweighting is a hard rejection at $\chi^2_{pq,0.975}$, followed by a re-MLE and re-scaling.
- **Regularization hooks:** `lambda` gives shrinkage to the identity inside the flip-flop. `weight_mat_row/col` is a Schur-product tapering or banding mask (a hook for Zhang et al. 2023-style bandable covariances). `diag=` gives diagonal factors.
- **Classification:** `matrix_qda.R` exists but is commented out (quadratic discriminant function with MMD and the log-dets). Matrix QDA is a planned but unfinished item.
- **Bundled data:** `weather` (Hoher Sonnblick, 5 variables × 12 months × years), `darwin` (handwriting, Alzheimer diagnosis). The paper also uses the beach-video example.
- **Robustness is purely casewise.** One bad cell makes the whole matrix an outlier.

## 3. Literature map (Leuven + TU Wien, plus the closest outside work)

### Cellwise core (vector data)
| Work | What it is | Relevance for matrix data |
|---|---|---|
| Rousseeuw & Van den Bossche 2018, *DDC* (Technometrics) | First correlation-aware cell detector; predicts each cell from correlated columns | Template for "matrix DDC" |
| Raymaekers & Rousseeuw 2024, *cellHandler* (Technometrics, arXiv 1912.12446) | Sparse-regression cell flagging given a robust Σ | Plug in Kronecker Σ, see §4 |
| Raymaekers & Rousseeuw 2024, *cellMCD* (JASA 119:2610, arXiv 2207.13493) | Observed likelihood + ℓ0 penalty on flagged cells per column, at least h unflagged per column, eigenvalue floor; C-step decreases the objective | **Direct template for a cellwise MMCD** |
| Rousseeuw 2026, *Analyzing cellwise weighted data* (EcoSta 38:31, arXiv 2209.12697) | Cellwise-weighted likelihood ("unpacking"), EM; `cwLocScat` | Soft-weight alternative to hard flags |
| Raymaekers & Rousseeuw 2024, *Challenges of cellwise outliers* (EcoSta 38:6, arXiv 2302.02156) | Overview, cellwise breakdown results, open questions | Breakdown definitions to adapt |
| Hubert, Raymaekers & Rousseeuw 2026, *Cellwise Outliers* (review, arXiv 2604.14182) | Survey: location/covariance, regression, PCA, **tensor data** | Check its tensor section for any matrix-normal work |
| Centofanti, Hubert & Rousseeuw, *cellPCA* (arXiv 2408.13596) | Single objective with casewise ρ + cellwise ρ, IRLS; residual cellmaps, enhanced outlier maps; casewise and cellwise IF | PCA-based, not likelihood-based |
| Centofanti, Hubert & Rousseeuw 2025, *cellRCov* (arXiv 2505.19925) | Casewise + cellwise + NA covariance in high dimensions via the cellPCA subspace + ridge; also cellRCCA | Vector analogue of what we want; competitor in high dimensions |
| Centofanti, Dağıdır, Hubert & Rousseeuw 2026, *Cellwise Robust DA* (arXiv 2605.30178) | cellLDA and cellQDA; **prediction on test cases with outlying cells** | Vector version of "cellwise matrix DA" |
| Leuven 2026, cellwise/casewise robust multivariate regression with inference (arXiv 2605.08002) | — | Lower priority |

### Matrix and tensor robust work
| Work | What it is |
|---|---|
| Mayrhofer, Radojičić & Filzmoser 2025, *MMCD* (Technometrics 67:516, arXiv 2403.03975) | Casewise robust Kronecker covariance; matrix affine equivariance; breakdown above the vectorized-AE maximum; Shapley explanations |
| Mayrhofer & Filzmoser 2023 (EcoSta, arXiv 2210.10063) | Shapley decomposition of MD²; linear cost in p; relation to cellwise outlyingness |
| Mayrhofer PhD thesis 2024 (TU Wien) | MMCD, Shapley, functional extensions |
| Hirari, Centofanti, Hubert & Van Aelst, *ROMPCA* (arXiv 2503.07327; JCGS 2026, doi 10.1080/10618600.2026.2637632) | Casewise + cellwise robust **multilinear PCA** (Tucker), IRLS, handles NA. **Order 2 = matrix data, so this is the closest existing cellwise method for matrices**, but it is PCA/subspace-based, not a separable covariance model |
| Leuven 2026, *Casewise and cellwise robust tensor-on-tensor regression* (ROTOT, arXiv 2603.25911) | Uses ROMPCA as preprocessing |
| Obukhov 2025, MSc thesis TU Wien | Tensor MCD (out of scope; shows the flip-flop without explicit Kronecker products) |
| Ma, Zhao, Wang et al. 2023 (CSDA, arXiv 2112.06760), and the 2025 HRFPCA via MMCD (arXiv 2509.25957) | Robust factored PCA for matrices (matrix-t); HRFPCA plugs in MMCD. Casewise only |
| He, Kong et al. (arXiv 2112.04186, 2306.03317; R package HDMFA) | Huber matrix factor models, **element-wise Huber loss** (IHR). Econometrics heavy tails, not outlier detection |

### Cellwise mixtures and multi-group (TU Wien and others)
| Work | What it is |
|---|---|
| Puchhammer, Wilms & Filzmoser 2025 (arXiv 2504.02547; R package ssMRCD) | Cellwise robust multi-group GMM; **first cellwise breakdown theory for mixtures** (Hennig-type) |
| Zaccaria, García-Escudero, Greselin & Mayo-Iscar 2025, *cellGMM* (Technometrics 67:643, arXiv 2409.07881) | Cellwise TCLUST-type mixture with cell imputation |
| Aerts & Wilms 2017 (arXiv 1612.07971) | Cellwise robust regularized DA (vector) |

### Matrix normal with missing cells (non-robust, but this is the M-step we need)
- Glanz & Carvalho 2018, JMVA 167:31 (arXiv 1309.6609): EM for the matrix normal with missing entries. The Kronecker structure is lost under arbitrary missingness, so they project the conditional covariance onto the row and column spaces.
- arXiv 2609.00616 (2026): EM for element-wise and structural missingness in matrix-normal **mixtures**. Lachos et al. 2025 give a Cholesky-based M-step.
- Allen & Tibshirani 2010 (AOAS): transposable regularized covariance, penalized EM.

## 4. Gap analysis: cellwise outliers × separable covariance

Searches found no method that combines **cellwise robustness** with a **matrix-normal / Kronecker model**:
- MMCD and HRFPCA: separable, but casewise only.
- cellMCD and cellRCov: cellwise, but unstructured covariance (cellMCD is practical only up to about 20 dimensions).
- ROMPCA: cellwise and matrix-aware, but a low-rank subspace model with no row/column covariance, MMD, or Shapley explanations.
- Matrix-normal EM handles missing cells, but is not robust.

**Why it matters.** Under independent cell contamination at rate ε, a matrix is fully clean with probability $(1-\varepsilon)^{pq}$ **[verified]**:

| pq | ε = 1% | ε = 0.1% |
|---|---|---|
| 25 (5×5) | 0.78 | 0.98 |
| 100 (10×10) | 0.37 | 0.90 |
| 400 (20×20) | 0.02 | 0.67 |
| 10⁴ (100×100 image) | ≈0 | ≈0 |

So MMCD (needs ≥50% clean matrices) breaks down already at 1% cell contamination for 10×20 matrices.

### Key structural identities (all **[verified]**)
Let $a=(j,k)$ be a cell and $\omega_{jk} := \omega^{row}_{jj}\,\omega^{col}_{kk}$.
1. **Single-cell conditional:** $\mathrm{Var}(x_{jk}\mid \text{rest}) = 1/\omega_{jk}$ and $x_{jk}-E(x_{jk}\mid\text{rest}) = R_{jk}/\omega_{jk}$.
   So the **standardized cellwise residual** is $r_{jk}=R_{jk}/\sqrt{\omega_{jk}}$, and all $pq$ residuals cost one $O(pq(p+q))$ product.
2. **Link to Shapley:** $\phi_{jk} = (x_{jk}-m_{jk})\,R_{jk} = (x_{jk}-m_{jk})\, r_{jk}\sqrt{\omega_{jk}}$.
   Shapley values are "deviation × conditional residual". This makes precise the relation to cellwise outlyingness noted in Mayrhofer & Filzmoser (2023).
3. **Cheap E-step for a flagged set $m$:** $E[x_m\mid x_o] = \mu_m - Q_{mm}^{-1}\,[R(\tilde X)]_m$, where $\tilde X$ has the flagged cells set to $M$ and $Q_{mm}[(j,k),(j',k')] = \omega^{row}_{jj'}\omega^{col}_{kk'}$.
   Cost: $O(pq(p+q) + |m|^3)$, with no $(pq)^3$ inversion.
4. **Cheap log-det:** $\log|\Sigma_{oo}| = q\log|\Sigma^{row}|+p\log|\Sigma^{col}| + \log|Q_{mm}|$ (Schur complement).
5. **Observed-data MD:** $\mathrm{MD}^2(x_o)=\mathrm{MMD}^2(\hat X)$, where $\hat X$ is the conditionally imputed matrix.

**Consequence.** A cellMCD-type observed-likelihood objective under a Kronecker model can be evaluated and C-stepped at cost **linear in the number of cells** (times $p+q$), plus small $|m_i|^3$ terms. That makes cellwise robustness feasible for image-sized matrices. This is the computational selling point over cellMCD and cellRCov.

## 5. Research directions (ranked)

**A. cellMMCD: cellwise + casewise robust Kronecker covariance (flagship).**
- **Objective:**
  $\sum_i [\log|(\Sigma^{col}\otimes\Sigma^{row})_{o_i o_i}| + \mathrm{MD}^2(x_{i,o_i})] + \sum_{jk} q_{jk}\,\#\{i: w_{i,jk}=0\}$,
  subject to at least $h$ unflagged entries per cell position $(j,k)$ and eigenvalue floors on both factors.
- **W-step:** flag cell $a$ if $r_a^2 + \log(2\pi/\omega_a^{(i)}) > q_a$, where $\omega_a^{(i)}$ is computed from the current conditional set (identities 1 and 4).
- **Σ-step:** flip-flop on expected sufficient statistics, e.g.
  $\Sigma^{row} \propto \sum_i [(\hat X_i-M)\Omega^{col}(\hat X_i-M)' + \text{corr}_i]$, where $\text{corr}_i$ comes from $Q_{mm}^{-1}$ mapped through $\Omega^{col}$.
  It is an open question whether this flip-flop step monotonically decreases the objective (in cellMCD the C-step provably does).
- **Casewise layer:** keep MMCD's h-subset of matrices, or use a casewise ρ as in cellPCA.
- **Initialization:** Kronecker-aware DDC. Each column $X_{i,\cdot k}\sim N(m_k,\sigma^{col}_{kk}\Sigma^{row})$, so pool the $nq$ column vectors (dimension p) and the $np$ row vectors (dimension q) and run DDC or cellMCD on those low-dimensional pools.
- **Theory:** cellwise breakdown for Kronecker estimators; equivariance (matrix affine equivariance is lost, but row/column permutation and diagonal-scale equivariance should hold); consistency factors.
- **Code:** extend `robustmatrix` (Rcpp/Armadillo, same style).

**B. Structured cellwise contamination (unique to matrix data).**
A matrix observation has four outlier types: cell, within-matrix row (e.g. one broken sensor across time), within-matrix column (one bad time point across variables), and whole matrix (casewise).
Group-ℓ0 penalties on the flag matrix $W_i$ (rows and columns) would flag a whole row when several of its cells are bad.
This includes detection and diagnostics: a cellmap per observation plus row/column outlier maps.

**C. Shapley values under cellwise contamination ("smearing").**
Mean-imputation Shapley values spread a single bad cell to correlated cells via $R=\Omega^{row}(X-M)\Omega^{col}$.
Identity 2 shows exactly how. Study:
- (i) when Shapley flags are faithful;
- (ii) Shapley values computed **after** cellwise imputation (an explain-after-cleaning pipeline);
- (iii) a comparison with cellHandler-style sparse attribution.

This sits on the TU Wien core competence.

**D. Cellwise robust matrix discriminant analysis (cellMLDA/cellMQDA).**
- Builds on A and on the unfinished `matrix_qda`, and mirrors the 2026 Leuven cellLDA/cellQDA but with Kronecker pooling.
- Pooled LDA under a Kronecker model: is a common $\Sigma^{row}$, or a common $\Sigma^{col}$, or both, appropriate?
- Includes prediction for test matrices with bad cells (impute via identity 3, then classify), plus Shapley-based explanations of the class assignment. This matches proposal Hyp. I.5–I.6.

**E. Further items.**
- Cellwise robust matrix PCA via a Kronecker covariance, compared with ROMPCA.
- Cellwise robust test of separability (cell outliers can mimic non-separability).
- High-dimensional p, q > n with tapering (the `weight_mat` hook) or shrinkage.
- Matrix elliptical / matrix-t versions.
- Missing data plus outliers in the weather data.

## 6. Benchmarks and data for experiments
- **Competitors:** MMCD (casewise), cellMCD on vec (only feasible for pq ≲ 20–40), cellRCov/cellPCA on vec, ROMPCA (order 2), matrix-normal EM with oracle NA positions, DDC on vec.
- **Contamination schemes:** ICM cells (ε ∈ {1, 2, 5, 10}%); row-block and column-block structured contamination; casewise (mixed); adversarial cells along the smallest eigenvectors of $\Sigma^{row}$ and $\Sigma^{col}$.
- **Metrics:** KL divergence of $\hat\Sigma^{col}\otimes\hat\Sigma^{row}$ to the truth (as in the MMCD paper); precision and recall of cell flags; imputation MSE; runtime against pq.
- **Real data:** `weather` (glacier station), `darwin`, the beach video (frames as matrices, the walking person as structured cellwise contamination), possibly hyperspectral or EEG.

## 7. Blockers for an autonomous overnight run
1. **Network:** allow `arxiv.org`, `export.arxiv.org`, `repositum.tuwien.at`, `wis.kuleuven.be`, `lirias.kuleuven.be`, `cran.r-project.org`, `cloud.r-project.org` (and optionally `api.openalex.org`, `api.semanticscholar.org`) in the environment's network settings. Paywalled journal versions will stay out of reach. A VPN on the local machine does not affect the cloud container.
2. **R:** not preinstalled. `apt-get install r-base-core r-base-dev` works (R 4.3.3). CRAN is blocked, so packages come from `github.com/cran/<pkg>` via `R CMD INSTALL` unless CRAN is allowed. A setup script should handle this.
3. **Persistence:** the container is ephemeral, so results must be committed and pushed periodically.
