# Slice-trimmed Kronecker covariance estimation: findings

Status 15 Oct 2026. Every number below comes from the reruns of 15 Oct. Steps 1–6 and `RJ_simple.R` reproduce the earlier results exactly (same seeds).

**Setup unless stated:** p = q = 5, n = 25, Σʳ = AR(0.5), Σᶜ = AR(0.7) (non-diagonal), M = 0 known, β = 0.75, outliers shifted by δ = 5.
- MMCD reference: `robustmatrix::mmcd` 0.1.5, α = 0.75.
- KL = Kullback–Leibler divergence of N(0, Σ̂) from N(0, Σ). Shape-KL is the same divergence minimised over the overall scale.

## Bottom line
- **RJ is the method worth pursuing.**
  - It is block-coordinate descent on a single trimmed likelihood, J\*, which equals the matrix-normal likelihood at β = 1.
  - J\* never increased in any run.
  - It exists whenever every column index keeps **h ≥ q + 1** matrices.
  - **κ(β, p) on the product is the correct consistency factor**. The remaining finite-sample bias is 0.6–4.4 % at n = 1000 and decays like n^(-1/2).
  - It beats MMCD by a factor of 3–6 in KL under column-type contamination. It is as efficient as MMCD on clean data.
  - Weak spots: row-wise contamination, many local optima, and dependence on the start.
- **R1 is dropped.** It has no objective, and an exact limit cycle exists on clean Gaussian data, so no convergence guarantee is possible (§6).

## 1. RJ in one paragraph
Theory: `notes/rj_objective.md`. Code: `tff(rule = "RJ")` in `R/tff.R`. A self-contained version is `fit_rj()` in `exp/RJ_simple.R`.

The unit is the whitened column z_ij = X_i a_j. Each column carries its p/q share of the Jacobian:

J\* = Σ_ij a_ij [log|Σʳ| + (p/q) log|Σᶜ| + d_ij], with Σᶜ = (AA′)⁻¹. The root A is a free parameter. The per-column quota keeps h = βn of the n matrices for every column index j.

Each iteration has three exact block steps:
1. C-step per column index;
2. closed-form Σʳ;
3. column-wise coordinate descent for A.

Without the quota (RJu), J\* → −∞ by trimming a whole column index. That happened in 90 % of clean runs at β = .75 (step 4), so the quota is part of the method.

The two implementations are not step-for-step identical: `tff` solves the A block exactly, while `fit_rj` does one sweep per iteration. The test suite checks that `fit_rj`'s end point is a fixed point of `tff`'s block updates.

## 2. Convergence and existence
- **Monotone:** J\* never increased in any run (steps 1, 4, 9; `tests/test_core.R`).
- **Convergence:**
  - Step 1: all 400 runs converge, in 17–18 iterations.
  - Step 9: all 5400 clean-data fits (n ≥ 25) converge.
- **Existence threshold (step 10, new).** The grid was p ∈ {3, 5}, q ∈ {3, 5, 8}, β ∈ {.5, .75}, and h from min(p,q) − 1 up to max(p,q) + 4, giving 2880 runs.
  - Every run converges once **h ≥ q + 1**. This is the exact threshold in 10 of 12 (p, q, β) settings. The other two are off by one in either direction: for p = 5, q = 3, β = .75 it is q + 2; for p = 3, q = 8, β = .75 it is q. With 15 reps per cell the threshold is only known to ±1.
  - The threshold follows q, not p: for p = 5, q = 3, β = .5, h = 4 is enough although h < p.
  - Below the threshold, runs end singular or degenerate (855) or at maxit (99). The maxit runs creep: the median J\* change over the last 50 iterations is −0.002, so there is no minimiser to converge to.
  - In practice: **n ≥ (q + 1)/β**. This explains the earlier step 5 failures at n = 10, β = .5, p = q = 5, where h = 5 = q.
- **Conjecture for the proof:** for each j, V_j = Σ_{kept i} X_i′(Σʳ)⁻¹X_i must stay positive definite under the rotation freedom of A, which needs more than q kept matrices. This is open.

## 3. Consistency factor (step 9, new): κ(β, p) is right
Clean data, three starts with the lowest J\* kept, 50 reps. The table shows the mean of tr(Σ⁻¹ κ Σ̂)/(pq); 1 means the scale is right, and SE ≤ 0.03.

| | β | n = 25 | 100 | 1000 |
|---|---|---|---|---|
| oracle (true Σᶜ), 5×5 | .5 / .75 / .9 | 1.04 / 1.05 / 1.00 | 1.02 / 1.01 / 1.00 | **1.00 / 1.00 / 1.00** |
| RJ 3×3 | .5 / .75 / .9 | 1.04 / 1.02 / 0.96 | 0.97 / 0.97 / 0.98 | 0.99 / 0.99 / 0.99 |
| RJ 5×5 | .5 / .75 / .9 | 0.92 / 0.96 / 0.93 | 0.95 / 0.96 / 0.97 | 0.98 / 0.99 / 0.99 |
| RJ 3×6 | .5 / .75 / .9 | 0.86 / 0.89 / 0.88 | 0.88 / 0.91 / 0.94 | 0.96 / 0.97 / 0.98 |

- **The factor is right.**
  - With Σᶜ known, κ(β, p) gives exactly 1. The factor the data ask for at n = 1000 equals κ(β, p) to 2 digits in all 9 oracle settings: for example 1.91 vs 1.91 (5×5, β = .5).
  - The MMCD-type factor κ(β, pq) is far off: 1.29 at 5×5, β = .5.
- **The remaining bias comes from fitting the rotation.**
  - RJ's plug-in scale is too small by 0.6–4.4 % at n = 1000. The bias shrinks like n^(-0.3…-0.56), median slope −0.5 for n ≥ 100.
  - It is larger for larger q (the q(q−1)/2 free rotation parameters) and for heavier trimming.
  - Reading: the trimmed set is chosen to fit the rotation that makes the kept columns tight, an optimism of order n^(-1/2). It vanishes, so it does not affect consistency, but it matters at small n.
  - Geometric-mean scale ratios agree (0.95–0.99 at n = 1000).
- A small-sample correction factor, fitted like the MCD's, would remove most of it. It is not yet implemented.

## 4. Robustness: RJ vs MLE, R1, MMCD (step 6, `RJ_simple.R`)
Median KL, n = 25, p = q = 5, 8 reps (IQRs are in `output/step6_mmcd_full.log`).

| | FF (MLE) | R1 | R2 | RJ | MMCD |
|---|---|---|---|---|---|
| clean | **0.63** | 1.59 | 1.05 | 1.31 | 1.29 |
| obs .1 / .2 | 1.7 / 4.4 | 1.3 / 2.0 | **0.90 / 1.07** | 1.3 / 1.7 | 1.0 / 1.1 |
| col .1 / .2 | 26 / 56 | 4.1 / 25 | 9.3 / 44 | **1.9 / 11** | 5.2 / 34 |
| row .1 / .2 | 11 / 24 | **2.2 / 5.6** | 5.2 / 17 | 3.3 / 11 | 3.5 / 13 |
| cell .1 / .2 | 84 / 157 | 48 / 147 | 68 / 155 | **45 / 136** | 59 / 141 |
| wcol .1 / .2 | 9.9 / 23 | 4.6 / 11 | 6.7 / 21 | **1.1 / 2.5** | 5.2 / 15 |

- **Column-type contamination (col, wcol):** RJ beats MMCD by a factor of 3–6. With diagonal Σᶜ (`output/step6_mmcd_diag.log`) the result is the same: col .2 gives 3.0 vs 10.4, wcol .2 gives 1.9 vs 16.5.
- **Clean data:** RJ's efficiency is the same as MMCD's (1.31 vs 1.29).
- **Row contamination is RJ's weak spot.** RJ only trims columns. A bad row touches every z_ij of that X_i, so the whole matrix would need trimming, but the per-column quota cannot do that. With diagonal Σᶜ, RJ gets 10.1 / 33 against MMCD's 5.0 / 30.
- **Observation-wise contamination:** MMCD (R2) is best, as expected.
- **Cellwise** (`output/rj_simple.txt`: n = 100, p = 3, q = 4, 50 reps, 6-sd outliers; first value for scale outliers, second for shift outliers):
  - RJ at β = .75 beats MMCD from 5 % bad cells on: 1.7 vs 2.9 and 2.3 vs 17.6.
  - At 10 % all methods break down; RJ at β = .5 degrades most slowly (4.0 / 4.8 vs 28 / 64).

## 5. Known weaknesses of RJ
- **Many local optima.** 10 starts give 10 distinct kept sets (step 1). On clean data the median sKL between starts (0.90) is of the same order as the KL to the truth (2.3). Multi-start with the lowest J\* is mandatory, as in FastMCD.
- **Start and root dependence.** J\* is exactly flat along A → AQ at the model (`notes/rj_objective.md` §3), so the sample picks the rotation. Different starting roots give different optima: the max sKL across 7 roots is 1.76, against a KL to the truth of 1.5 (step 3). This does not affect consistency, because AA′ is identified.
- **The symmetric-root rule does not apply to RJ**, because there the root is estimated. Only the start uses S^(-1/2).
- **Location:** M is the coordinatewise median and is not updated. Its effect has not been checked.

## 6. R1: documented negative result (archived)
Code: `R/archive/r1.R`, `exp/archive/`. Theory: `notes/r1_convergence.md`.
- R1 alternates two different criteria. Without trimming they coincide; with trimming it is a two-player game, not a descent method.
- **Counterexample** (step 8, rerun): clean Gaussian data, p = q = 2, MLE start.
  - On 15 Oct: n = 7, β = .75. R1 enters an exact period-24 cycle, confirmed by the independent `tff` implementation, in which **every half-step strictly improves its own criterion** (minimum gain 0.16).
  - On 5 Oct the search found n = 17, β = .5, period 4.
  - Either way, **no** Lyapunov function exists, and no convergence theorem is possible without extra conditions.
  - Step 8 searches for a fixed 200 s of wall-clock time, so the number of data sets and the showcase change from run to run.
- **Frequency** (step 8, 21 771 clean 2×2 data sets, MLE start), cycle rate for n in (5,10] / (10,20] / (20,40]:
  - β = .5: 4.9 / 2.7 / 1.2 %;
  - β = .75: 1.4 / 0.5 / 0.1 %;
  - β = .9: 0.4 / 0.05 / 0.1 %.
- The archived step 7 (15 120 runs, 45 min, not rerun) found 6.5 % non-converging runs at β = .5 and 1–2 % at β ≥ .75, and up to 28 % for 2×6 shapes at n = 10.
- **Decision:** non-convergence is a deal-breaker, so R1 is not pursued. Its accuracy at β = .75 was good (best on row contamination), which RJ does not yet match.

## 7. Other spec questions, answered (steps 1–4)
- **Q2 normalisation:** the denominator, κ and the scale identification change only the overall scale; the kept sets and shapes are identical to 3·10⁻¹⁵.
- **Q3 root dependence:** for fixed-root rules (R0, R1, R3), the trimmed set depends on the arbitrary root, and the effect is as large as the estimation error. Only R2 is invariant.
- **R0** does not protect Σᶜ (as bad as the MLE under col). **R3** is the least efficient rule and cycles in 8 % of runs.
- **Reweighting** is not pursued: with slice trimming, a reweighting step would be either flawed or a different (EM-type) algorithm.

## 8. Open questions (RJ)
1. **Row contamination:** can RJ also trim rows within one objective (a joint row/column RJ), or is a two-sided version bound to lose the single objective, as R1 did?
2. **Proof of consistency:** sample-to-population (uniform LLN for trimmed sums modulo the rotation orbit), and the existence condition h ≥ q + 1.
3. **Small-sample factor:** fit c(n, p, q, β) to the step 9 bias; check it under contamination.
4. **Taming multimodality:** deterministic starts (for example from MMCD or the MLE), or a canonical rotation at the end. Measure how many starts are enough.
5. **Location:** joint estimation of M under column trimming, and the effect of the median plug-in.
6. **Larger settings:** larger p and q, and a comparison with cellwise-robust matrix estimators from the literature.

## Reproduce
From the project root, with `R = "C:/Program Files/R/R-4.6.1/bin/Rscript.exe"`. Times are wall-clock on this machine (16 cores) on 15 Oct.

```
Rscript tests/test_core.R                    # 31 tests, ~2 s
Rscript exp/step1_convergence.R              # 65 s  -> output/step1.rds, step1_J_traces.png, step1_status.png
Rscript exp/step2_normalization.R            # 4 s   -> output/step2.rds, step2_shape.png
Rscript exp/step3_root.R                     # 57 s  -> output/step3.rds, step3_root_spread.png
Rscript exp/step4_matching.R                 # 99 s, needs step1 + step3 -> output/step4.rds, step4_accuracy.png
Rscript exp/step5_consistency.R              # 87 s  -> output/step5.rds, step5_ratio.png, step5_qq.png
Rscript exp/step6_mmcd.R full                # 108 s -> output/step6_full.rds, step6_full_boxplots.png
Rscript exp/step6_mmcd.R diag                # 95 s  -> output/step6_diag.rds, step6_diag_boxplots.png
Rscript exp/step9_rj_consistency.R           # 11 min (15 cores) -> output/step9.rds, step9_ratio.png
Rscript exp/step10_rj_existence.R            # 3.5 min (15 cores) -> output/step10.rds
Rscript exp/RJ_simple.R > output/rj_simple.txt                 # 3 min
Rscript exp/archive/step8_r1_cycles.R        # 3.5 min -> output/step8.rds, step8_cycle.png
Rscript exp/archive/step7_r1_practical.R     # ~45 min (15 cores), not rerun on 15 Oct
```
Console output of each run is in `output/<script>.log`.
