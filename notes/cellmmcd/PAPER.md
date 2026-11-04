# cellMMCD: a cellwise robust estimator of Kronecker-structured covariances

_Paper draft, blocks 7–8 (2026-11-04). Assembled from `METHOD.md` (§0 and the sections it cites), `JOURNAL.md` and the logs in `output/`. Every number below has a log; every theoretical claim carries its status: **proven**, **checked** (numerically) or **conjectured**. Open points for the owner are marked ▸._

## 1. Setting and aim
Matrix observations $X_i\in\mathbb R^{p\times q}$, $i=1,\dots,n$, with $\operatorname{vec}(X_i)\sim N(\operatorname{vec}M,\ \Sigma^c\otimes\Sigma^r)$ at the model. Some cells are replaced by arbitrary values (cellwise contamination); some whole matrices may be replaced (casewise). The aim is an estimator of $(M,\Sigma^r,\Sigma^c)$ that (i) does not impute flagged cells, (ii) keeps the Kronecker structure in the fit and not only in a post-hoc projection, (iii) is robust to cellwise and, as far as possible, casewise outliers, and (iv) runs in time linear in $pq$ per matrix. $\Sigma^r$ is identified by $\det\Sigma^r = 1$; estimators are compared through $\Sigma^c\otimes\Sigma^r$ (KL divergence to the truth), never through the raw factors.

Related work: cellMCD (Raymaekers & Rousseeuw 2024) on unstructured $\Sigma$ of dimension $pq$, cost $O(n(pq)^3)$ per iteration; MMCD (Mayrhofer, Radojičić & Filzmoser 2024) for casewise outliers with the Kronecker structure; the flip-flop MLE (Dutilleul 1999); the owner's trimmed-likelihood estimators RJ and R1 for whitened columns.

## 2. The estimator
Let $x_i = \operatorname{vec}(X_i - M)$, $w_{ikl}\in\{0,1\}$ cell weights (1 = kept), $O_i$ the kept cells of matrix $i$, $M_i$ the flagged ones, $v_i\in\{0,1\}$ observation weights. cellMMCD minimises the **trimmed observed-cell likelihood with a Kronecker constraint**
$$F(M,\Sigma,W,v)=\sum_i v_i\Big[\log\det\Sigma_{O_iO_i}+x_{iO_i}'\Sigma_{O_iO_i}^{-1}x_{iO_i}+|O_i|\log 2\pi+\sum_{(k,l)\in M_i}\lambda_{kl}\Big]+\sum_i(1-v_i)\,\mu_i,\qquad \Sigma=\Sigma^c\otimes\Sigma^r,$$
subject to (a) a per-position quota $\sum_i v_iw_{ikl}\ge h=\lceil\alpha n_v\rceil$, (b) forced flags for cells beyond 3 MAD from the position median (the *marginal cut*), (c) $\sum_i v_i\ge h_{obs}=\lceil\alpha_{obs}n\rceil$, (d) a kept matrix has at most half of its cells flagged (majority rule) and at most $\mathrm{qbinom}(.999,pq,\bar\varepsilon)$ flagged cells (count rule, $\bar\varepsilon$ a robust estimate of the flag rate). Defaults $\alpha=\alpha_{obs}=0.75$. No cell value is ever replaced: $F$ is the marginal likelihood of the kept cells; flagged cells enter only through the penalty $\lambda_{kl}$, i.e. as if they sat at their own $0.99$ conditional quantile. $W\equiv1$, $v\equiv1$ gives the matrix-normal MLE; trimming whole matrices gives an MMCD-type trimmed likelihood.

**Penalty and calibration.** $\lambda_{kl}=c/t_1+\log t_1+\log2\pi-\log(\Omega^r_{kk}\Omega^c_{ll})$ at the initial estimate, $c=\chi^2_{1,.99}$, $t_1=F_{\chi^2_3}(c)/F_{\chi^2_1}(c)=0.925$, final scale divided by $t_1$. Without the calibration the cell and parameter steps form a shrink-and-flag loop with fixed point at scale 0.905 and 1.4 % flagged; with it the fixed point is Fisher-consistent in scale with exactly 1 % flagged (**proven** in the 1-D proxy, **checked**: trace ratio 0.98–0.99).

**Observation rule.** $\mu_i$ is the $0.99$ quantile of a clean matrix's contribution evaluated on its kept cells (soft threshold, default); the alternative is the hard quota $\mu=-\infty$ with exactly $n-h_{obs}$ trimmed, scaled by the MCD consistency factor $\kappa(\alpha_{obs},pq)$.

**Reweighting (the final estimate).** At the raw fit every cell of a kept matrix is re-decided with cutoff $\chi^2_{1,.999}$ on its residual given the other kept cells (marginal cut kept); for the hard quota the observation weights are re-decided by the kept-cell distance $d_i^2\le\chi^2_{|O_i|,.99}$ plus the count rule; then the observed-cell MLE with $(v,W)$ fixed, scale divided by $t_{rw}\,t_{obs}$. Reason (§5): the raw fit has a shape bias that does not vanish with $n$.

## 3. Algorithm
**The Schur identity (proven, tested to $10^{-15}$).** With $Y_i$ the centred matrix with flagged cells set to 0, $u_i=(\Omega^rY_i\Omega^c)[M_i]$ and $\Omega_{M_iM_i}=[\Omega^r_{k_ak_b}\Omega^c_{l_al_b}]$:
$\log\det\Sigma_{O_iO_i}=q\log|\Sigma^r|+p\log|\Sigma^c|+\log\det\Omega_{M_iM_i}$ and $x_{iO_i}'\Sigma_{O_iO_i}^{-1}x_{iO_i}=\operatorname{tr}(\Omega^rY_i\Omega^cY_i')-u_i'\Omega_{M_iM_i}^{-1}u_i$. Cost $O(p^2q+pq^2+m_i^3)$ per matrix instead of $O((pq)^3)$.

Start: `kpair`, a pairwise Kronecker estimate from position medians, MADs (median polish of $\log s^2_{kl}=a_k+b_l$) and Gnanadesikan–Kettenring correlations, eigenvalue-clipped; Fisher-consistent at the Gaussian (**checked**), cellwise robust by construction. Then block coordinate descent:
0. *Observation step*: with $(\Sigma,W)$ fixed, trim the matrices with the largest $d_i-\mu_i>0$ (exact minimiser over $v$ on the feasible set).
1. *Cell step*: for each position, $\Delta_{ic}-\lambda_{kl}$ with $\Delta_{ic}=\log(2\pi\sigma^2_{ic})+r^2_{ic}$ (conditional variance and standardized residual of the cell given the kept cells of its matrix) equals exactly $F(w_{ikl}=1)-F(w_{ikl}=0)$; flag the largest positive values, at most $n-h$ (exact minimiser over the block). In C++ with $\Omega_{M_iM_i}^{-1}$ cached and updated by bordering / downdating: $O(m^2)$ per candidate. During the first 10 iterations trimmed matrices also have their cells re-decided (*shadow step*), which removes a lock-in between the two trimming steps.
2. *Parameter step*: the EM step for $M$ on the kept cells, then one flip-flop sweep of the observed-data score fixed point, both computed from conditional means and $\Omega_{M_iM_i}^{-1}$; a GEM step (**proven**).

**Convergence (proven).** $F$ is non-increasing after the shadow phase, $(v,W)$ ranges over a finite set, so $F_t$ converges and $(v,W)$ is eventually constant; the parameter sequence is then an EM sequence whose limit points are stationary points of the observed-cell likelihood, provided $F$ is bounded below for that pattern (§6).

**Complexity and timing.** Per iteration $O\big(n(p^2q+pq^2+\bar m^2pq)\big)$, $\bar m$ flags per matrix. One core, $n=100$, 5 % cells: 0.4 / 0.4 / 1.0 / 2.4 / 15 / 111 s at $p=q=5/10/15/20/30/40$ ($n=120$ at 40); MMCD 0.2 / 0.8 / 2.2 / 4.7 / 16 / 46 s; cellMCD on $\operatorname{vec}X$ is only feasible for small $pq$. Recall of 6-sd cells 1.00, false flags ≤ 0.6 % at every size. (`cm_step10_timing.log`, `cm_step22_timing_ext.log`.)

## 4. Theory
| claim | status |
|---|---|
| Schur form, exactness of the cell and observation steps, GEM parameter step, monotone $F$ | proven, tested |
| $F$ unbounded below under the quota alone when $n-h\ge(n-p+1)/\max(p,q)$ (adversarial kept patterns); more generally whenever some row set $T$ has between 1 and $\lvert T\rvert-1$ kept $T$-partial columns (or rows) | proven (Prop. A); the converse is **false** ($2\times2$ counterexample with both factors degenerating together, proven; $3\times3$, $n = 3$ non-transversal instance checked numerically; block 8), so the exact existence criterion for the joint problem is open; factor-wise coercivity under $p$ fully kept columns / $q$ fully kept rows proven (Prop. B); the Prop.-A criterion held in every run checked (`cm_check_pattern`) |
| cross-path descent: $F = \nu\log\delta+O(1)$ along $\Sigma^r = P_{a^\perp}+\delta aa'$, $\Sigma^c = P_{b^\perp}+\delta bb'$ iff the data are orthogonal to the null vectors $N_i(a,b)$ supported on the kept cells ($\nu$ their total dimension plus the $\delta^2$ modes); $\dim N_i$ is a component count of the missing-cell graph on the support stratum; generic data admit no such descent when every stratum has $\sum_i\dim N_i\ge\lvert T\rvert+\lvert U\rvert-1$ | proven (Props. C′, D, Corollary; block 9), slopes checked (0.98, 0.98, 4.00 for $\nu$ = 1, 1, 4); converse conjectured on strata with $\lvert T\rvert,\lvert U\rvert\ge2$; the stratum check never fired on an algorithm-produced pattern ($3\times3$ to $5\times5$, $n$ = 6 to 100) |
| random cellwise contamination at rate $\varepsilon$ puts the kept pattern into Prop. A's unbounded regime once $nq(1-\varepsilon)^p<p$: $20\times20$, $n = 100$: $c_{[p]} = 16$ at $\varepsilon = 0.20$ (KL 3.7), $5$ at $0.25$ (KL 11.0) | proven regime, checked numerically; an operating limit at $p,q\ge20$ |
| scale calibration $t_1$ | proven (1-D proxy), checked |
| raw estimator biased in shape at the Gaussian: $\Sigma-\kappa\operatorname{diag}(1/\Omega_{cc})$, $\kappa=F_1(c)-F_3(c)=0.0745$ | proven to first order; iterated map 0.0111 vs observed 0.0115 |
| reweighted estimator consistent: $n\cdot$KL 8.3 / 8.0 / 7.4 at $n$ = 500 / 2000 / 8000 (MLE 8.9 / 7.4 / 5.9) | checked |
| influence function: closed form in the 1-D proxy (GES 7.4 raw, 10.0 reweighted, ≈ 0 beyond the cutoff); in the Kronecker model $J^{-1}[v\,s(z_O)+K\,\mathrm{IF}_{raw}(z)]$, bounded by the marginal cut, finite limit for gross cells and matrices | proven (proxy: checked to $10^{-4}$; Kronecker: under regularity, $K$ not closed) |
| sensitivity curves bounded and flat beyond the cutoff (one cell 4.95 vs MLE 824 at $x=30$; shifted matrix 5.5 vs 3975) | checked |
| explosion breakdown by gross cells impossible below 1/2 contaminated cells per position (marginal cut); checked: no breakdown at 45 % | proven; checked |
| bias breakdown by moderate cells at the quota $1-\alpha$ per position; casewise at $1-\alpha_{obs}$ for all observation rules; beyond it the quota rules degrade gracefully, the soft rule breaks at 30–35 % | checked |
| location by the EM step: error halved against the cellwise median, which breaks at 20 % shifted matrices | proven (GEM), checked |

**Existence, in words.** The objective can be unbounded below for adversarial kept patterns (one fully kept column and a flag in every other column lets a singular $\Sigma^r$ fit it exactly); the same counting applies to cellMCD. For the patterns the algorithm produces the checkable criterion of Prop. A never fired, and the small-$n$ runs once called "degenerate" are convergent finite minima of an overfitted objective. The estimator is defined as the local minimiser reached from `kpair`; operating range $n\ge3\max(p,q)$ (a variance limit) **and** $nq(1-\varepsilon)^p\gtrsim p$, $np(1-\varepsilon)^q\gtrsim q$ for the cellwise rate $\varepsilon$ (an existence limit, block 9: at $20\times20$, $n = 100$ this is $\varepsilon\lesssim0.2$, at $40\times40$ $\varepsilon\lesssim0.1$).

## 5. Simulations (median KL to the truth, 6-sd cells with random sign, 10–20 reps; logs cited)
**Cellwise, $n=100$** (`cm_step1.log`): 3×4 / 5×5, cellMMCD vs cellMCD-on-vec (+ Kronecker projection) vs MMCD vs MLE:
clean 0.08 / 0.19 vs 0.28 / 0.65 vs 0.09 / 0.18 vs 0.07 / 0.16; 5 % cells 0.09 / 0.20 vs 0.26 / 0.63 vs 14 / 47 vs 27 / 63; 10 % 0.09 / 0.18 vs 0.26 / 0.59 vs 50 / 119 vs 61 / 144; 20 % 0.13 / 0.27 vs 0.32 / 1.14 vs 123 / 290 vs 134 / 308. Recall 0.97–1.00, false flags 0.1–0.6 %.

**Casewise and mixed, 5×5** (`cm_step20_default.log`, $n$ = 100 / 200; soft rule, reweighted hard quota in brackets; MMCD):
shift 20 % 0.20 / 0.11 [0.21 / 0.11] (0.21 / 0.15); scale 20 % 0.21 / 0.10 [0.24 / 0.10] (0.25 / 0.18); struct 20 % 0.21 / 0.14 [0.24 / 0.15] (0.24 / 0.17); shift 10 % + 5 % cells 0.22 / 0.09 [0.21 / 0.09] (52 / 49); struct 20 % + 5 % cells 0.96 / 0.52 [0.37 / 0.23] (50 / 46); clean 0.17 / 0.09 [0.17 / 0.09] (0.16 / 0.08). At 8×8 the same pattern with cellMMCD ahead of MMCD on every casewise case and the weak case gone (0.62 [0.60] at $n=100$).

**Weak case.** Structural casewise outliers (sign-flipped correlations) under 5 % cellwise contamination at small $pq$: after their worst cells are flagged they look like cellwise-contaminated clean matrices; only the count rule can catch them. The reweighted hard quota is the recommended fallback there. ▸ Owner's choice of the default: the soft rule (better on clean data by 0.01–0.03) or the reweighted hard quota (robust in the weak case); both are presented.

**Breakdown** (`cm_step13_breakdown.log`, `cm_step19_gross.log`, `cm_step21_casebreak.log`): cellwise bias breakdown at 25 % of 4-sd cells per position (KL 5 at 25 %, 71 at 30 %), no breakdown up to 45 % gross cells per position; casewise at 25 %.

**Location.** With the EM step the squared error of $M$ halves on every case; the cellwise median breaks at 20 % shifted matrices (3.0 vs 0.3).

Figures: `cm_fig_breakdown.png`, `cm_fig_if.png`, `cm_fig_mixed.png`, `cm_fig_timing.png`, `cm_fig_sc.png`, `cm_fig_casebreak.png` (median KL against the casewise fraction for the three observation rules, shifted and structural matrices; `exp/cm_plot_casebreak.R`), `cm_step6_bias.png`.

## 6. Limits and open problems
1. **Existence.** The "only if" of the configuration criterion (multi-rate degenerations, scale migrating between the factors); a worst-case bound on $n-h$ that guarantees it cannot exist (§4a), so the statement stays local.
2. **Consistency.** The reweighted estimator's residual bias at $n\to\infty$ is $\kappa(.999)^2$-small but not zero; a monotone consistent formulation (censored cells with the cutoff region fixed) is open; the censored parameter step removed two thirds of the bias but is not monotone.
3. **Identifiability of structural casewise outliers under cellwise contamination** at small $pq$ (the weak case).
4. **Breakdown** for moderate cells and the general casewise statement are numerical; the quota $1-\alpha$ is low for a cellwise method (cellMCD has the same).
5. **Cells that are bad and near their conditional mean** are invisible, as for every cellwise method.
6. **Speed beyond 30×30**: the $\bar m^2pq$ term of the cell step; `kpair` and the R-level parameter step are the next candidates.

## 7. Reproducibility
`R/cellmmcd.R` (driver `cellmmcd`, `kpair`, `cm_reweight`, pattern check `cm_check_pattern`), `src/cm_wstep.cpp`; tests `tests/test_cellmmcd.R` (26) and `tests/test_core.R` (31); experiments `exp/cm_step*.R` with seeds in the header; logs and figures in `output/`.
