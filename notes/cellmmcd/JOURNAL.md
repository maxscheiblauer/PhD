# cellMMCD journal (append-only, dated)

## 2026-11-04, block 1: candidate formulations and the choice

**Why RJ cannot be the cellwise method.** RJ trims whitened columns z_ij = X_i a_j. A single bad cell (k,l) of X_i enters coordinate k of *every* z_ij (j = 1..q), so RJ must trim all q columns of X_i, i.e. the whole matrix. With ε bad cells, P(matrix clean) = (1-ε)^{pq}: 28 % for ε = .1, pq = 12. That is the 10 % breakdown in `output/FINDINGS.md` §4. Any method whose trimming unit contains many cells has the same problem. The unit must be the cell.

**Candidates considered.**
1. *Product-set deletion* (drop whole rows and columns of X_i so the remaining block stays matrix-normal, N(0, Σᶜ_CC ⊗ Σʳ_RR)): closed-form Σʳ step, but each bad cell costs a full row or column; at ε = .1 half the data goes. Rejected.
2. *Trimmed pseudo-likelihood* (product of full conditionals, each cell as response, trimmed per cell): cheap, but a bad covariate cell corrupts the conditional of every good cell in its row and column (propagation), and the pseudo-likelihood is not the likelihood. Kept as fallback.
3. *Pairwise Kronecker estimator* (robust scale per cell position, robust GK correlation per pair of positions in the same row or the same column, then the Kronecker structure gives Σʳ, Σᶜ directly): O(n(p²q+pq²)), no imputation, no iteration, consistent at the Gaussian with the usual GK factors. Low efficiency, not a likelihood. **Chosen as the initial estimator** (`kpair()`).
4. *Trimmed observed-cell likelihood with Kronecker parameter* = the cellMCD objective (Raymaekers & Rousseeuw 2024) with Σ = Σᶜ ⊗ Σʳ. **Chosen as the method, "cellMMCD".**

**The objective.** Cells flagged in X_i form M_i, the rest O_i; x_i = vec(X_i - M). With λ_kl ≥ 0 and a quota of at least h = ⌈αn⌉ unflagged cells per position (k,l):

F(Σʳ, Σᶜ, W) = Σ_i [ log det Σ_{O_iO_i} + x_{iO_i}' Σ_{O_iO_i}^{-1} x_{iO_i} + |O_i| log 2π ] + Σ_{i,k,l} λ_kl (1 - w_ikl).

Special cases: W ≡ 1 gives the matrix-normal likelihood (flip-flop MLE); flagging whole matrices gives a trimmed-observation likelihood (MMCD-type). No cell value is ever replaced: the estimator is the marginal likelihood of the cells that are kept.

**The identity that makes it cheap (proven, elementary).** For Σ = Σᶜ⊗Σʳ, Ω = Σ^{-1} = Ωᶜ⊗Ωʳ, Y_i = X_i - M with flagged cells set to 0, T_i = Ωʳ Y_i Ωᶜ... more precisely T_i = Ωʳ (X_i - M) Ωᶜ computed on the full matrix, and u_i = (Ωʳ Y_i Ωᶜ)[M_i]:
- det Σ_{OO} = det Σ · det Ω_{MM}, so log det Σ_{O_iO_i} = q log|Σʳ| + p log|Σᶜ| + log det Ω_{M_iM_i};
- x_O' Σ_OO^{-1} x_O = x_O' (Ω_OO - Ω_OM Ω_MM^{-1} Ω_MO) x_O = tr(Ωʳ Y_i Ωᶜ Y_i') - u_i' Ω_{M_iM_i}^{-1} u_i;
- Ω_{MM} = [Ωʳ_{k_a k_b} Ωᶜ_{l_a l_b}]_{ab} is an m_i × m_i matrix read off the factors.
Cost per observation O(p²q + pq² + m_i³) instead of O((pq)³). In the cellwise regime m_i ≈ εpq is small.

**Cell step (exact block coordinate descent).** For a cell c of observation i, with S = M_i ∪ {c}, the conditional law of x_c given the other observed cells has variance σ²_c = (Ω_SS^{-1})_cc and residual r_c = (Ω_SS^{-1} T_S)_c / σ_c (T = Ω x on the full vector; this holds whether c is currently flagged or not). The change of F from unflagging c is Δ_c - λ_kl with Δ_c = log(2πσ²_c) + r_c². F is separable over i given the parameters, so for one position (k,l) the exact minimizer over w_{·kl} flags the cells with the largest Δ - λ > 0, at most n - h of them. Sweeping the positions is therefore monotone. With m_i = 0 the rule reads: flag iff the full-conditional standardized residual (Ωʳ(X_i-M)Ωᶜ)_kl / sqrt(Ωʳ_kk Ωᶜ_ll) is large; that is the Kronecker version of the DDC/cellMCD residual. λ_kl = χ²_{1,quant} + log 2π + log(1/(Ωʳ_kk Ωᶜ_ll)) at the initial estimate, as in cellMCD.

**Parameter step and the EM question (owner guideline 2).** Given W, the first-order condition for Σʳ is Σʳ = (nq)^{-1} Σ_i E[X_i Ωᶜ X_i' | x_{iO_i}], and likewise for Σᶜ, because the observed-data score equals the conditional expectation of the complete-data score (Fisher's identity). So *every* first-order method for this objective evaluates the same conditional-mean and conditional-covariance terms; "EM" here is nothing but the plain fixed-point iteration of the exact score, and with the identity above its cost is O(p²q + pq² + m_i³) per observation, not the O((pq)³) that made EM look prohibitive. I therefore use one flip-flop sweep of this fixed point per outer iteration (a GEM step, so F does not increase). This does not violate guideline 1: the estimator is defined by the observed-cell likelihood and never replaces a cell by a predicted value; the conditional means appear only inside the exact score. The non-EM alternative (gradient step with line search on the same F) uses the same quantities and only adds a step-size search; I see no cheaper exact alternative and record that as the reason.

**Convergence (proven, given boundedness).** F never increases, W takes finitely many values, and for fixed W the parameter step is monotone, so the sequence F_t converges; W stabilizes after finitely many steps. Open: a condition on h, p, q under which F is bounded below (RJ needed h ≥ q+1; here I conjecture h ≥ max(p,q)+1 per position suffices, to be checked numerically).

**Plan for this block.** Implement `kpair`, the objective with a brute-force check, the cell step with a brute-force check of Δ, the parameter step, `cellmmcd()`; tests; one experiment on the FINDINGS cellwise setting (p = 3, q = 4, n = 100, 6-sd cells) against FF, RJ, MMCD, kpair and cellMCD on vec(X) + nearest Kronecker.

## 2026-11-04, block 1 (end): results of steps 1–2
- `exp/cm_step1.R` (n = 100, 20 reps, 6-sd shift cells, random sign; `output/cm_step1.log`): cellMMCD flags the bad cells with recall 1.00 at ε = .05/.10/.20 (0.99 at 5×5, ε = .2) and false-flag rate 0.8–1.9 %. Median KL at 5×5: ε = .10: cellMMCD 0.32, cellMCD-on-vec + nearest Kronecker 0.59, kpair 1.5, RJ 89, MMCD 119, MLE 144; ε = .20: 0.44 vs 1.14 (cellMCD) with everything else > 280. Clean data: KL 0.17 (3×4) / 0.37 (5×5) vs MLE 0.07 / 0.16, MMCD 0.09 / 0.18, RJ 0.17 / 0.32; the loss is mostly scale (tr ratio 0.90; shape KL 0.14 / 0.28).
- Run time (pure R, n = 100): 0.12–0.66 s at 3×4 and 5×5; 10.3 s at 10×10, n = 200, ε = .1 (19 iterations, KL 0.97 vs MLE 640, recall 1.00, fpr 1.6 %), where cellMCD on vec(X) would need n ≥ 500.
- `exp/cm_step2_factor.R`: the clean-data scale ratio is 0.895 ± 0.004 at n = 1000 for both 3×4 and 5×5 (0.885–0.899 at n = 100), with 1.7 % cells flagged at n = 1000 against the nominal 1 %. The naive 1-D tail-trimming factor F_{χ²_3}(c)/F_{χ²_1}(c) = 0.925 is not the right one: the shrunk scale inflates the conditional residuals, which flags more cells, which shrinks more. This is the fixed point of that loop. Shape KL falls from 0.136 to 0.027 (3×4) and 0.257 to 0.059 (5×5) between n = 100 and 1000: the shape looks consistent, the scale needs a factor. Next: apply a consistency factor inside the iteration (scale Σ by 1/c(quant) before the cell step) so that the nominal 1 % is flagged, then derive c.
