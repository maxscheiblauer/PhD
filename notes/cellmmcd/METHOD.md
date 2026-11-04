# cellMMCD: a cellwise matrix MCD for Kronecker covariances

_Draft of 2026-11-04 (block 1). Derivations: `JOURNAL.md`. Code: `R/cellmmcd.R`. Tests: `tests/test_cellmmcd.R`. Experiment: `exp/cm_step1.R`, log `output/cm_step1.log`._

## 1. Setting
$X_i \in \mathbb R^{p\times q}$, $i = 1,\dots,n$, $\operatorname{vec}(X_i) \sim \mathcal N(\operatorname{vec}M,\ \Sigma^c\otimes\Sigma^r)$ at the model; some cells are replaced by arbitrary values. $\Omega = \Sigma^{-1} = \Omega^c\otimes\Omega^r$. Cell $(k,l)$ has vec index $c = k + (l-1)p$. $w_{ikl}\in\{0,1\}$ ($1$ = kept), $O_i$ = kept cells, $M_i$ = flagged cells, $m_i = |M_i|$.

## 2. Definition
$$
\hat\theta = \arg\min_{\Sigma^r,\Sigma^c,W}\ F,\qquad
F = \sum_{i=1}^n\Big[\log\det\Sigma_{O_iO_i} + x_{iO_i}'\Sigma_{O_iO_i}^{-1}x_{iO_i} + |O_i|\log 2\pi\Big] + \sum_{i,k,l}\lambda_{kl}(1-w_{ikl}),
$$
subject to $\sum_i w_{ikl} \ge h = \lceil\alpha n\rceil$ for every position $(k,l)$, where $x_i = \operatorname{vec}(X_i - M)$. This is the cellMCD objective of Raymaekers and Rousseeuw (2024) with the Kronecker constraint on $\Sigma$, exactly as MMCD is the MCD with that constraint. $W\equiv 1$ gives the matrix-normal MLE; flagging whole matrices gives a trimmed-observation likelihood of MMCD type. No cell value is replaced: $F$ is the marginal likelihood of the kept cells.

Defaults (as in cellMCD): $\alpha = 0.75$, $\lambda_{kl} = \chi^2_{1,0.99} + \log 2\pi + \log\big(1/(\Omega^r_{kk}\Omega^c_{ll})\big)$ evaluated at the initial estimate, so a cell is flagged when its conditional standardized residual exceeds the $0.99$ quantile. Cells with $|x_{ikl}|/s_{kl} > 3$ (MAD scale) are flagged permanently. $M$ is the coordinatewise median.

## 3. The Schur identity (proven; tested to $10^{-15}$)
With $Y_i$ the centred matrix with flagged cells set to $0$, $u_i = (\Omega^rY_i\Omega^c)[M_i]$ and $\Omega_{M_iM_i} = [\Omega^r_{k_ak_b}\Omega^c_{l_al_b}]_{ab}$:
$$
\log\det\Sigma_{O_iO_i} = q\log|\Sigma^r| + p\log|\Sigma^c| + \log\det\Omega_{M_iM_i},\qquad
x_{iO_i}'\Sigma_{O_iO_i}^{-1}x_{iO_i} = \operatorname{tr}(\Omega^rY_i\Omega^cY_i') - u_i'\Omega_{M_iM_i}^{-1}u_i .
$$
Both follow from $\det\Sigma_{OO} = \det\Sigma\det\Omega_{MM}$ and $\Sigma_{OO}^{-1} = \Omega_{OO} - \Omega_{OM}\Omega_{MM}^{-1}\Omega_{MO}$. Cost per observation $O(p^2q + pq^2 + m_i^3)$ instead of $O((pq)^3)$.

## 4. Algorithm (block coordinate descent on $F$)
Start: `kpair` (§6). Then repeat until $F$ stops changing:
1. **Cell step.** For every position $(k,l)$ in turn, with $T_i = \Omega^r(X_i - M)\Omega^c$ and $S = M_i\cup\{c\}$:
   $\sigma^2_{ic} = (\Omega_{SS}^{-1})_{cc}$, $r_{ic} = (\Omega_{SS}^{-1}T_{i,S})_c/\sigma_{ic}$, $\Delta_{ic} = \log(2\pi\sigma_{ic}^2) + r_{ic}^2$.
   $\Delta_{ic} - \lambda_{kl}$ is exactly $F(w_{ikl}=1) - F(w_{ikl}=0)$ (tested). Flag the observations with the largest positive values, at most $n-h$ of them. Since $F$ is separable over $i$ for fixed parameters, this is the exact minimiser over the block $w_{\cdot kl}$.
2. **Parameter step.** One flip-flop sweep of the score fixed point
   $\Sigma^r \leftarrow \frac1{nq}\sum_i \mathbb E[X_i\Omega^cX_i'\mid x_{iO_i}]$, $\Sigma^c \leftarrow \frac1{np}\sum_i\mathbb E[X_i'\Omega^rX_i\mid x_{iO_i}]$,
   computed from $\hat X_i$ (conditional mean, $\hat x_{iS} = x_{iS} - \Omega_{SS}^{-1}T_{iS}$ on $S = M_i$) and $C_i = \Omega_{M_iM_i}^{-1}$:
   $\mathbb E[X\Omega^cX'\mid\cdot] = \hat X\Omega^c\hat X' + E_k'(C\circ\Omega^c_{LL})E_k$ with $E_k$ the $m_i\times p$ row-indicator matrix. By Fisher's identity this is the fixed point of the exact observed-data score; it is also one GEM step, so $F$ does not increase (tested). Without flags it is the flip-flop step (tested).

**Convergence (proven given boundedness).** $F$ is non-increasing, $W$ ranges over a finite set, and each parameter step is monotone; hence $F_t$ converges and $W$ is eventually constant. With $W$ fixed the parameter sequence is an EM sequence for a missing-data matrix-normal problem, whose limit points are stationary points of the observed-cell likelihood (Wu 1983) if the likelihood is bounded. **Boundedness (conjectured):** $F$ is bounded below when every position keeps $h \ge \max(p,q)+1$ cells; not yet checked numerically.

## 5. Complexity
Per outer iteration: cell step $O\big(n\,pq\,(m+1)^3 + n(p^2q+pq^2)\big)$, parameter step $O\big(n(p^2q + pq^2 + m^3)\big)$, with $m$ the number of flagged cells per matrix ($\approx \varepsilon pq$ under cellwise contamination). The dependence on $pq$ is linear up to the small $m$-solves, against $O(n(pq)^3)$ for cellMCD on $\operatorname{vec}X$. Timing in `output/cm_step1.log`.

## 6. Initial estimator `kpair` (no iteration, no imputation)
$s_{kl} = \operatorname{MAD}_i(X_{ikl} - M_{kl})$; $\log s_{kl}^2 = a_k + b_l$ by median polish (exact at the model because $\operatorname{var}X_{kl} = \Sigma^r_{kk}\Sigma^c_{ll}$); $\rho^r_{kk'} = \operatorname{med}_l\,\mathrm{GK}(U_{\cdot kl}, U_{\cdot k'l})$ on standardized cells, since $\operatorname{corr}(X_{kl},X_{k'l}) = \rho^r_{kk'}$ for every $l$; same for $\rho^c$; eigenvalue clipping for positive definiteness. Fisher-consistent at the Gaussian (MAD and GK with the standard constants; shape errors 0.02–0.03 at $n = 3000$, tested). Cellwise robust by construction (each statistic uses one or two positions).

## 7. Theory status
| claim | status |
|---|---|
| Schur-form objective and cell criterion | proven, tested |
| cell step is exact block minimisation; parameter step is GEM; $F$ monotone | proven, tested |
| MLE recovered at $W\equiv1$ | proven, tested |
| $F$ bounded below under the quota | conjectured |
| Fisher consistency at the Gaussian with $\lambda$ fixed, size of the tail-trimming bias | open (tr ratio measured in step 1) |
| breakdown under cellwise contamination | open |

## 8. Open problems
1. Boundedness and existence under the quota; the role of $h$ relative to $p$, $q$.
2. Consistency factor for the tail trimming at quant $= 0.99$, and a reweighting or quant $\to 1$ step at the end.
3. Joint estimation of $M$ (the median plug-in is used now).
4. Rcpp for the cell step when $pq$ is large; vectorising the $m = 0$ case.
5. Cells that are bad *and* near their conditional mean cannot be detected (as for every cellwise method).

## 9. Results so far (block 1; details in `JOURNAL.md` and `output/cm_step1.log`)
Median KL, n = 100, 20 reps, cells shifted by 6 sd with random sign:

| ε | cellMMCD | cellMCD on vec + Kronecker | kpair | RJ | MMCD | MLE |
|---|---|---|---|---|---|---|
| 0 (3×4 / 5×5) | 0.17 / 0.37 | 0.28 / 0.65 | 0.45 / 1.04 | 0.17 / 0.32 | 0.09 / 0.18 | 0.07 / 0.16 |
| .05 | 0.13 / 0.38 | 0.26 / 0.63 | 0.40 / 1.20 | 1.7 / 8.5 | 14 / 47 | 27 / 63 |
| .10 | 0.14 / 0.32 | 0.26 / 0.59 | 0.78 / 1.55 | 14 / 89 | 50 / 119 | 61 / 144 |
| .20 | 0.19 / 0.44 | 0.32 / 1.14 | 2.3 / 4.7 | 113 / 289 | 123 / 290 | 134 / 308 |

Recall of bad cells 0.99–1.00, false flags 0.8–1.9 %. Clean-data scale ratio 0.895 (needs a factor, §8 item 2). Run time 0.1–0.7 s at n = 100 (pure R); 10 s at 10×10, n = 200.
