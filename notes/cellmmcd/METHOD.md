# cellMMCD: a cellwise matrix MCD for Kronecker covariances

_Draft of 2026-11-04 (blocks 1–4). Derivations: `JOURNAL.md`. Code: `R/cellmmcd.R`. Tests: `tests/test_cellmmcd.R` (13). Experiments: `exp/cm_step1.R` (cellwise benchmark), `cm_step3_calib.R` (scale calibration), `cm_step3b/3c` (efficiency), `cm_step4_mixed.R` (casewise and mixed), `cm_step5_exist.R` (existence grid), `cm_step6_bias.R` / `cm_step6b_biasformula.R` (bias of the raw estimator), `cm_step7_safeguard.R` (full-column safeguard, negative), `cm_step8_default.R` (choice of the default); logs in `output/`._

## 1. Setting
$X_i \in \mathbb R^{p\times q}$, $i = 1,\dots,n$, $\operatorname{vec}(X_i) \sim \mathcal N(\operatorname{vec}M,\ \Sigma^c\otimes\Sigma^r)$ at the model; some cells are replaced by arbitrary values. $\Omega = \Sigma^{-1} = \Omega^c\otimes\Omega^r$. Cell $(k,l)$ has vec index $c = k + (l-1)p$. $w_{ikl}\in\{0,1\}$ ($1$ = kept), $O_i$ = kept cells, $M_i$ = flagged cells, $m_i = |M_i|$.

## 2. Definition
$$
\hat\theta = \arg\min_{\Sigma^r,\Sigma^c,W}\ F,\qquad
F = \sum_{i=1}^n\Big[\log\det\Sigma_{O_iO_i} + x_{iO_i}'\Sigma_{O_iO_i}^{-1}x_{iO_i} + |O_i|\log 2\pi\Big] + \sum_{i,k,l}\lambda_{kl}(1-w_{ikl}),
$$
subject to $\sum_i w_{ikl} \ge h = \lceil\alpha n\rceil$ for every position $(k,l)$, where $x_i = \operatorname{vec}(X_i - M)$. Write $d_i$ for the $i$-th summand (observed-cell contribution plus its penalties). **Observation weights (block 2):** with $v_i\in\{0,1\}$ and $\sum_i v_i \ge h_{obs} = \lceil\alpha_{obs}n\rceil$,
$$F = \sum_i v_i d_i + \mu\sum_i(1-v_i),\qquad \mu = q\log|\Sigma^r_0| + p\log|\Sigma^c_0| + pq\log 2\pi + \chi^2_{pq,0.99}$$
at the initial estimate ($\mu = -\infty$: hard quota, the MMCD-type trimmed likelihood). A whole matrix is trimmed when its contribution exceeds the $0.99$ quantile of a clean one; flagged cells enter $d_i$ through $\lambda$, i.e. as if they sat at their own $0.99$ conditional quantile. The cell quota then counts the kept observations, $h = \lceil\alpha n_v\rceil$. This is the cellMCD objective of Raymaekers and Rousseeuw (2024) with the Kronecker constraint on $\Sigma$, exactly as MMCD is the MCD with that constraint. $W\equiv 1$ gives the matrix-normal MLE; flagging whole matrices gives a trimmed-observation likelihood of MMCD type. No cell value is replaced: $F$ is the marginal likelihood of the kept cells.

Defaults: $\alpha = \alpha_{obs} = 0.75$, quant $= 0.99$ for the flags, reweighting at $.999$ (§4d), and the **calibrated penalty** $\lambda_{kl} = c/t_1 + \log t_1 + \log 2\pi + \log\big(1/(\Omega^r_{kk}\Omega^c_{ll})\big)$ at the initial estimate, $c = \chi^2_{1,0.99}$, $t_1 = F_{\chi^2_3}(c)/F_{\chi^2_1}(c) = 0.925$; the final $\Sigma$ is divided by $t_1$ (and by $t_{obs} = F_{\chi^2_{pq+2}}(c_{obs})/F_{\chi^2_{pq}}(c_{obs})$ for the observation step). Reason (§4b): without calibration the fixed point of the iteration sits at scale $0.90$ with $1.4\%$ of clean cells flagged; with it, exactly $1 - $ quant are flagged and the scale is Fisher-consistent (proven in the 1-D proxy, checked numerically in the Kronecker model: tr ratio 0.98–0.99). Cells with $|x_{ikl}|/s_{kl} > 3$ (MAD scale) are flagged permanently. $M$ is the coordinatewise median.

## 3. The Schur identity (proven; tested to $10^{-15}$)
With $Y_i$ the centred matrix with flagged cells set to $0$, $u_i = (\Omega^rY_i\Omega^c)[M_i]$ and $\Omega_{M_iM_i} = [\Omega^r_{k_ak_b}\Omega^c_{l_al_b}]_{ab}$:
$$
\log\det\Sigma_{O_iO_i} = q\log|\Sigma^r| + p\log|\Sigma^c| + \log\det\Omega_{M_iM_i},\qquad
x_{iO_i}'\Sigma_{O_iO_i}^{-1}x_{iO_i} = \operatorname{tr}(\Omega^rY_i\Omega^cY_i') - u_i'\Omega_{M_iM_i}^{-1}u_i .
$$
Both follow from $\det\Sigma_{OO} = \det\Sigma\det\Omega_{MM}$ and $\Sigma_{OO}^{-1} = \Omega_{OO} - \Omega_{OM}\Omega_{MM}^{-1}\Omega_{MO}$. Cost per observation $O(p^2q + pq^2 + m_i^3)$ instead of $O((pq)^3)$.

## 4. Algorithm (block coordinate descent on $F$)
Start: `kpair` (§6). Then repeat until $F$ stops changing:
0. **Observation step.** With $(\Sigma, W)$ fixed, trim the observations with the largest $d_i - \mu > 0$, at most $n - h_{obs}$ of them: the exact minimiser over $v$ (tested against brute force). The cell and parameter steps below run over the kept observations only.
1. **Cell step.** For every position $(k,l)$ in turn, with $T_i = \Omega^r(X_i - M)\Omega^c$ and $S = M_i\cup\{c\}$:
   $\sigma^2_{ic} = (\Omega_{SS}^{-1})_{cc}$, $r_{ic} = (\Omega_{SS}^{-1}T_{i,S})_c/\sigma_{ic}$, $\Delta_{ic} = \log(2\pi\sigma_{ic}^2) + r_{ic}^2$.
   $\Delta_{ic} - \lambda_{kl}$ is exactly $F(w_{ikl}=1) - F(w_{ikl}=0)$ (tested). Flag the observations with the largest positive values, at most $n-h$ of them. Since $F$ is separable over $i$ for fixed parameters, this is the exact minimiser over the block $w_{\cdot kl}$.
2. **Parameter step.** One flip-flop sweep of the score fixed point
   $\Sigma^r \leftarrow \frac1{nq}\sum_i \mathbb E[X_i\Omega^cX_i'\mid x_{iO_i}]$, $\Sigma^c \leftarrow \frac1{np}\sum_i\mathbb E[X_i'\Omega^rX_i\mid x_{iO_i}]$,
   computed from $\hat X_i$ (conditional mean, $\hat x_{iS} = x_{iS} - \Omega_{SS}^{-1}T_{iS}$ on $S = M_i$) and $C_i = \Omega_{M_iM_i}^{-1}$:
   $\mathbb E[X\Omega^cX'\mid\cdot] = \hat X\Omega^c\hat X' + E_k'(C\circ\Omega^c_{LL})E_k$ with $E_k$ the $m_i\times p$ row-indicator matrix. By Fisher's identity this is the fixed point of the exact observed-data score; it is also one GEM step, so $F$ does not increase (tested). Without flags it is the flip-flop step (tested).

**Convergence (proven given boundedness).** $F$ is non-increasing, $(v, W)$ ranges over a finite set, and each parameter step is monotone; hence $F_t$ converges and $(v,W)$ is eventually constant. With $(v,W)$ fixed the parameter sequence is an EM sequence for a missing-data matrix-normal problem, whose limit points are stationary points of the observed-cell likelihood (Wu 1983) if the likelihood is bounded.

### 4a. Existence (block 2)
**Proposition (proven for data in general position; proof in `JOURNAL.md`).** Let $p,q\ge2$. If $n - h \ge n/\max(p,q)$, then $\inf F = -\infty$: a singular $\Sigma^r$ with null vector $a$ fits one fully kept column exactly, and the quota allows one flag in every other column, which removes every other null direction supported on the kept cells. Hence the objective has no minimiser at $\alpha = 0.75$ once $\max(p,q)\ge4$, and never at $\alpha = 0.5$. The same counting applies to the unstructured cellMCD objective with $d = pq$. What the algorithm computes is therefore a local minimiser reached from the robust start; in all experiments with $n = 100$ it never approached a degenerate configuration, while at $n\le 12$ it does (`output/cm_step5_exist.log`). **Converse (conjectured):** $n - h < n/\max(p,q)$ suffices for boundedness. **Remedy (conjectured sufficient):** additionally require $p$ fully kept columns $(i,l)$ for every $l$ and $q$ fully kept rows for every $k$; this couples the positions in the cell step and is not implemented yet.

### 4c. Observation threshold and the majority rule (block 2)
The trimming threshold $\mu_i$ is evaluated on the kept cells of $X_i$ only: $\mu_i = \log\det\Sigma_{0,O_iO_i} + |O_i|\log2\pi + \chi^2_{|O_i|,0.99}/t_{obs} + |O_i|\log t_{obs} + \sum_{M_i}\lambda_{kl}$ at the initial estimate (penalties cancel against those in $d_i$). Charging the flagged cells instead (a fixed $\mu$ with $pq$ degrees of freedom) trimmed 12.5 % of the matrices under 5 % cellwise contamination (three flagged cells suffice), which is wrong for a cellwise method. The observed-cell threshold alone, however, misses a shifted matrix most of whose cells are already flagged: its few kept cells look fine by themselves (40 % of shifted matrices trimmed at 20 % contamination, KL 1.33). Hence the constraint that a kept observation has at most `max_flag` $= 1/2$ of its cells flagged; matrices beyond it are trimmed first in the observation step. With this constraint the observation step is the exact minimiser over $v$ of $F$ on the feasible set; a cell step can push a kept matrix over the bound, after which the next observation step trims it, so monotonicity of $F$ is guaranteed only between such transitions (tested: still monotone on the test data). The hard quota ($\mu = -\infty$, exactly $n - h_{obs}$ trimmed) is available as `quant_obs = NULL`; it needs the MCD-type consistency factor $\kappa(\alpha, pq)$, not applied yet (tr ratio 0.87 on clean data).

### 4d. Shape bias of the raw estimator and the reweighting step (block 3)
**Proposition (proven to first order; iterated population map checked numerically).** At the Gaussian, with cutoff $c = \chi^2_{1,\text{quant}}$ and $\kappa = F_1(c) - F_3(c)$, the population fixed point of the cell/parameter iteration at the truth is, to first order in $\kappa$, $\Sigma^\* = \Sigma - \kappa\operatorname{diag}(1/\Omega_{cc})$ projected onto the Kronecker set. Reason: a flagged clean cell is a *tail* cell of its full conditional ($E[r_c^2\mid |r_c|>c] = (1-F_3)/(1-F_1) = 8.45$), but the parameter step treats it as missing (conditional variance $\sigma_c^2$); the cross moments are unchanged by symmetry. The $t_1$ factor removes the trace part; the rest is a shape bias of order $\tfrac14\kappa^2\sum_{i\ne j}\mathrm{pcor}_{ij}^2$ in KL that does not vanish with $n$: observed KL $\to 0.0115$ at 3×4 (quant .99, $n = 8000$), predicted 0.0111 by the iterated map (`cm_step6b_biasformula.R`). At quant .999, $\kappa = 0.0117$ and the bias is 40 times smaller (0.0003). This is the "efficiency gap" of block 2.

**Remedy adopted (default `reweight = TRUE`).** The raw iteration at quant .99 provides $(v, W)$; every cell of a kept matrix is then re-decided at the raw fit with the cutoff $\chi^2_{1,.999}$ on its residual given the other kept cells (forced marginal flags kept, no quota, $v$ kept), and the parameters are refitted with $(v,W)$ fixed by the EM/flip-flop steps (monotone; the observed-cell MLE given the flags), with the scale divided by $t_{rw} = F_3(c_{rw})/F_1(c_{rw})$ and $t_{obs}$. This is the analogue of the reweighted MCD. Clean 3×4: KL ratio to the MLE 1.01 / 1.09 / 1.49 at $n$ = 500 / 2000 / 8000 ($n\cdot$KL bounded); recall at 3-sd cells 0.96 (quant .999 alone: 0.94). Alternative tried: a *censored* parameter step (`tail = TRUE`: flagged cells contribute the tail second moment $\tau\sigma_c^2$), Fisher-consistent in the 1-D proxy (proven) but only removing two thirds of the bias in the Kronecker model and not monotone in $F$; kept as an option.

### 4b. Scale calibration (block 2)
In the 1-D proxy (one cell, working scale $t$, initial scale $t_0$), the cell step flags iff $x^2 > a(t) = t\,(c + \log(t_0/t))$ and the parameter step gives $t' = F_3(a) + t\,(1 - F_1(a))$; the fixed point of the uncalibrated iteration is $t^* = 0.905$ with $1.35\%$ flagged (observed in the Kronecker model: $0.895$, $1.7\%$). With $\lambda$ built from $c/t_1 + \log t_1$ the fixed point is $t_1 = F_3(c)/F_1(c)$ with exactly $1-$quant flagged, and $\Sigma/t_1$ is Fisher-consistent for the scale (proven in the proxy; `cm_proxy_fixpt`). The shape is unaffected to first order: trimming the tails of every full conditional scales $\Omega$ uniformly. Checked numerically: tr ratio $0.98$–$0.99$ at $n = 100$ and $1000$, both shapes.

### 4e. Lock-in between the two trimming steps and the shadow cell step (block 4)
The observation step runs before the cell step, so at the first pass only the marginal flags exist ($|x|/\mathrm{MAD} > 3$). With 10 % of the cells at 3–4 sd in every matrix the still-unflagged cells push more than 25 % of the matrices beyond $\mu_i$, the quota binds, and the trimmed matrices were locked out: the cell step only re-decides cells of kept matrices, so $d_i - \mu_i$ stayed positive at every later observation step (block 3 trimmed 25 % of the matrices in this setting and 4 % of the clean ones). Since the $\lambda$ terms cancel in $d_i - \mu_i$, what the trimmed matrix needs is a current $W_i$. **Shadow cell step:** during the first `shadow_it` $= 10$ iterations the cell step also re-decides the cells of the trimmed matrices by the rule $\Delta_{ic} - \lambda_{kl} > 0$ (coordinate descent on their own $d_i$, no quota); the next observation step then re-admits the matrix if its kept cells pass, while the majority rule still trims a shifted matrix (most of its cells get shadow-flagged). A trimmed matrix enters $F$ only through $\mu_i(W_i)$, so $F$ is not guaranteed monotone during the shadow phase; afterwards the iteration is that of §4 and the convergence argument applies after at most 10 non-monotone iterations (monotone in every run checked, test 18). Results (`cm_step9_interplay.log`, p = q = 5, n = 100, median KL / trimmed): clean 0.182 / 0 (block 3: 0.198 / 0.04); 10 % cells at 3 sd 0.332 / 0 (0.335 / 0.25), at 4 sd 0.259 / 0 (0.297 / 0.25); 20 % shifted matrices 0.344 / 0.20, exactly the shifted ones (0.372 / 0.23); mixed 0.222 / 0.10 (0.205 / 0.14). Alternatives rejected: no observation step (KL 1.96 at 20 % shifted matrices), `quant_obs` $= .999$ (under-trims shifted matrices, KL 0.51).

**Selection and the adaptive count rule (block 4).** The shadow step alone let scaled and structural casewise outliers back in: their worst 20–45 % of cells get flagged, and the remaining (selected) cells pass the kept-cell test (scale 10 / 20 %: KL 0.23 / 0.41, struct 0.36 / 0.41; `cm_step4_mixed_shadow_only.log`). Hence a second constraint on $(v, W)$ next to the majority rule: a kept matrix may have at most $\mathrm{qbinom}(.999, pq, \bar\varepsilon)$ flagged cells, $\bar\varepsilon$ the current overall flag rate among kept matrices (at least $1 - $ quant), the count being Binomial$(pq, \bar\varepsilon)$ under independent cellwise contamination. With it (`cm_step4_mixed.log`, soft column; block 3 in parentheses): none 0.20 (0.22), 5 % cells 0.17 (0.17), shift 10 / 20 % 0.17 / 0.22 (0.18 / 0.21), scale 0.19 / 0.22 (0.19 / 0.22), struct 0.18 / 0.32 (0.19 / 0.29), mixed 0.22 (0.23); clean matrices trimmed 0 % (1–3 %). Structural outliers at 20 % remain the weak case (75 % / 60 % of them trimmed at 10 / 20 %).

## 5. Complexity
Per outer iteration: cell step $O\big(n\,pq\,(m+1)^3 + n(p^2q+pq^2)\big)$, parameter step $O\big(n(p^2q + pq^2 + m^3)\big)$, with $m$ the number of flagged cells per matrix ($\approx \varepsilon pq$ under cellwise contamination). The dependence on $pq$ is linear up to the small $m$-solves, against $O(n(pq)^3)$ for cellMCD on $\operatorname{vec}X$. The cell step and the reweighting pass are in C++ since block 4 (`src/cm_wstep.cpp`, RcppArmadillo, identical results; `R/load.R` falls back to the R loops). Timing (`cm_step10_timing.log`, n = 100, 5 % cells at 6 sd, one core, default settings; the R-loop numbers in `cm_step10_timing_Rengine.log`):

| p = q | cellMMCD (C++) | cellMMCD (R loops) | raw / rw iterations | MMCD | MLE | KL cellMMCD / MMCD / MLE |
|---|---|---|---|---|---|---|
| 5 | 0.3 s | 0.8 s | 8 / 5 | 0.2 s | 0.03 s | 0.11 / 36 / 63 |
| 10 | 0.4 s | 2.7 s | 8 / 5 | 0.8 s | 0.04 s | 0.66 / 227 / 280 |
| 15 | 2.2 s | 12.4 s | 14 / 5 | 2.2 s | 0.10 s | 1.60 / 614 / 697 |
| 20 | 8.6 s | 38.8 s | 15 / 5 | 4.7 s | 0.21 s | 2.84 / 1149 / 1234 |

Recall 1.00 and false flags 0.4–0.5 % at every size. What remains at 20×20 is the R-level parameter step and objective (per-observation solves of $m\times m$ systems in R loops), the next candidates for C++ if larger $pq$ matters.

## 6. Initial estimator `kpair` (no iteration, no imputation)
$s_{kl} = \operatorname{MAD}_i(X_{ikl} - M_{kl})$; $\log s_{kl}^2 = a_k + b_l$ by median polish (exact at the model because $\operatorname{var}X_{kl} = \Sigma^r_{kk}\Sigma^c_{ll}$); $\rho^r_{kk'} = \operatorname{med}_l\,\mathrm{GK}(U_{\cdot kl}, U_{\cdot k'l})$ on standardized cells, since $\operatorname{corr}(X_{kl},X_{k'l}) = \rho^r_{kk'}$ for every $l$; same for $\rho^c$; eigenvalue clipping for positive definiteness. Fisher-consistent at the Gaussian (MAD and GK with the standard constants; shape errors 0.02–0.03 at $n = 3000$, tested). Cellwise robust by construction (each statistic uses one or two positions).

## 7. Theory status
| claim | status |
|---|---|
| Schur-form objective and cell criterion | proven, tested |
| cell, observation and parameter steps: exact, exact, GEM; $F$ monotone | proven, tested |
| MLE recovered at $W\equiv1$ | proven, tested |
| $F$ unbounded below when $n - h \ge n/\max(p,q)$ | proven (general position) |
| $F$ bounded below when $n - h < n/\max(p,q)$; full-column constraint as remedy | conjectured |
| scale fixed point and calibration factor $t_1$ | proven in the 1-D proxy, checked numerically |
| raw estimator inconsistent in shape at the Gaussian: fixed point $\Sigma - \kappa\operatorname{diag}(1/\Omega_{cc})$ to first order | proven to first order; iterated map checked numerically (0.0111 vs 0.0115) |
| reweighted estimator: $n\cdot$KL bounded, ratio to the MLE 1.0–1.5 at $n \le 8000$ | checked numerically |
| full-column constraint prevents degenerate descents | refuted numerically as a post-hoc repair (`cm_step7`) |
| influence functions of the raw and the reweighted estimator in the 1-D proxy (§7a) | proven in the proxy; checked numerically to $10^{-4}$ (test 19) |
| shadow cell step: convergence after at most `shadow_it` non-monotone iterations | proven (reduces to §4); monotone throughout in every run checked |
| breakdown under cellwise contamination | open |

### 7a. Influence functions in the 1-D proxy (block 4; derivation in `JOURNAL.md`)
One cell, residual distribution $F$, start $t_0(F)$ = squared standardized MAD, $N_F(a) = \int x^2 1(x^2\le a)\,dF$, $D_F(a) = \int 1(x^2\le a)\,dF$. The raw estimator is $T(F)/t_1$ with $T$ the fixed point of $t \mapsto N_F(a) + t(1 - D_F(a))$, $a = t(c/t_1 + \log t_1 + \log(t_0/t))$; the reweighted one is $T_{rw}(F) = N_F(a_{rw})/D_F(a_{rw})/t_{rw}$, $a_{rw} = c_{rw}T(F)/t_1$. Both are Fisher-consistent at $\Phi$, and
$$\mathrm{IF}_{raw}(x) = -\frac{1(x^2\le c)(x^2 - t_1) + t_1\varphi(\sqrt c)\frac{c - t_1}{\sqrt c}\mathrm{IF}_{t_0}(x)}{t_1\big[(c-t_1)^2\varphi(\sqrt c)/(t_1\sqrt c) - F_1(c)\big]},\qquad
\mathrm{IF}_{rw}(x) = 1(x^2\le c_{rw})\Big[\frac{x^2}{F_3(c_{rw})} - \frac1{F_1(c_{rw})}\Big] + \sqrt{c_{rw}}\varphi(\sqrt{c_{rw}})\frac{c_{rw}F_1 - F_3}{F_1F_3}\,\mathrm{IF}_{raw}(x).$$
At quant .99 / .999: $\mathrm{IF}_{raw}(x) = [1(x^2\le c)(x^2 - 0.925) + 0.030\,\mathrm{IF}_{t_0}(x)]/0.733$, gross-error sensitivity 7.4, and 0.09 beyond the cutoff (only the MAD start transmits); $\mathrm{IF}_{rw}(x) pprox 1(x^2\le 10.83)(1.013x^2 - 1.001) + 0.058\,\mathrm{IF}_{raw}(x)$, gross-error sensitivity 10.0 just below $\sqrt{c_{rw}} = 3.29$, and 0.006 beyond it. The MLE has $\mathrm{IF} = x^2 - 1$. `cm_proxy_if` / `cm_proxy_T` (exact functionals under $(1-\varepsilon)\Phi + \varepsilon\delta_x$). Not covered by the proxy: the conditional variances and the other cells' flags in the Kronecker model, and the observation step.

## 8. Open problems
1. Degenerate descents at small $n$: the full-column repair failed (§7). Block 4 (`cm_step11_smalln.log`, 4×4 clean): 20 % of the runs degenerate at $n = 8$, 5 % at $n = 10$, none at $n \ge 12$; the operating range is stated as $n \ge 3\max(p,q)$. Still open: the right constraint or stopping rule below that, and the converse of §4a.
2. A *monotone* consistent formulation (censored-cell objective with the cutoff region fixed at the initial estimate) instead of the reweighting device; the residual bias of the reweighted estimator at $n \to \infty$ is $\kappa(.999)^2$-small but not zero.
2b. ~~Interplay of the two trimming steps~~: solved by the shadow cell step (§4e); what remains is a monotone formulation of the shadow phase (e.g. an objective in which a trimmed matrix is charged $\min_{W_i}$ of something sensible).
3. Joint estimation of $M$ (the median plug-in is used now).
4. ~~Rcpp for the cell step~~ done (block 4, §5); the parameter step and the objective are next if $pq > 20$ matters.
5. Cells that are bad *and* near their conditional mean cannot be detected (as for every cellwise method).

## 9. Results so far (details in `JOURNAL.md` and `output/cm_step1.log`, rerun in block 3 with the final default: calibration, observation step, reweighting; the raw-.99 run is `cm_step1_raw99.log`)
Median KL, n = 100, 20 reps, cells shifted by 6 sd with random sign:

| ε | cellMMCD (default) | cellMMCD raw .99 | cellMCD on vec + Kronecker | kpair | RJ | MMCD | MLE |
|---|---|---|---|---|---|---|---|
| 0 (3×4 / 5×5) | 0.087 / 0.198 | 0.125 / 0.25 | 0.28 / 0.65 | 0.45 / 1.04 | 0.17 / 0.32 | 0.09 / 0.18 | 0.07 / 0.16 |
| .05 | 0.083 / 0.184 | 0.10 / 0.25 | 0.26 / 0.63 | 0.40 / 1.20 | 1.7 / 8.5 | 14 / 47 | 27 / 63 |
| .10 | 0.094 / 0.205 | 0.11 / 0.26 | 0.26 / 0.59 | 0.78 / 1.55 | 14 / 89 | 50 / 119 | 61 / 144 |
| .20 | 0.139 / 0.320 | 0.18 / 0.39 | 0.32 / 1.14 | 2.3 / 4.7 | 113 / 289 | 123 / 290 | 134 / 308 |

(cellMMCD column rerun in block 4 with the shadow step, before the count rule, which does not act in these settings; block-3 numbers in `cm_step1_block3.log`: 0.095 / 0.205, 0.089 / 0.187, 0.096 / 0.209, 0.145 / 0.353.)

Recall of bad cells 0.97–1.00, false flags 0.1–0.6 % (default) / 0.4–1.7 % (raw), clean-data scale ratio 0.99–1.00. Run time 0.1–0.7 s at n = 100 (pure R); 10 s at 10×10, n = 200.

**Casewise and mixed contamination (block 2, rerun in block 3 with the reweighted default; `output/cm_step4_mixed.log`, block-2 raw numbers in `cm_step4_mixed_block2.log`; p = q = 5, n = 100, 10 reps, median KL).** Casewise types on 10 % / 20 % of the matrices: shift (+5 in every cell), scale (×3), struct (flipped correlation signs). All cellMMCD columns are reweighted.

| contamination | cells only | cellMMCD (default: observation step, majority rule) | hard quota | MMCD | MLE |
|---|---|---|---|---|---|
| none | 0.20 | 0.20 | 0.38 | 0.18 | 0.16 |
| 5 % cells | 0.17 | 0.17 | 0.37 | 45 | 63 |
| shift 10 % / 20 % | 0.24 / 1.95 | 0.17 / 0.22 | 0.24 / 0.22 | 0.18 / 0.21 | 2.0 / 4.0 |
| scale 10 % / 20 % | 0.30 / 0.95 | 0.19 / 0.22 | 0.30 / 0.26 | 0.20 / 0.27 | 2.9 / 9.0 |
| struct 10 % / 20 % | 0.36 / 1.11 | 0.18 / 0.32 | 0.32 / 0.28 | 0.18 / 0.23 | 4.3 / 10.6 |
| shift 10 % + 5 % cells | 0.32 | 0.22 | 0.29 | 51 | 66 |

(default column rerun in block 4 with the shadow step and the count rule; the other columns are block 3's, `cm_step4_mixed_block3.log`.)

The default trims every shifted or scaled matrix and 85–95 % of the structural ones, at 0–3 % of the clean ones; the hard quota always trims 25 % and is not yet scale-corrected. The clean-data gap to MMCD is now 0.22 vs 0.18 (block 2: 0.31). The clean-data gap to MMCD and the MLE is the efficiency question of §8.
