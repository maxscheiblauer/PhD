# Morning report: cellMMCD (overnight relay, 2026-11-04, blocks 0–10)

_Written by block 10 (final, session_01Fi3wVf4AFqBEhFFyDjtwqz). Sources: STATE.md, PROGRESS_SUMMARY.md, JOURNAL.md, PAPER.md, LEDGER.csv. Nothing new was computed in this block except a rerun of the two test suites._

## 1. What exists now

**Method: cellMMCD.** The cellMCD objective (trimmed observed-cell Gaussian likelihood, no imputation) under the Kronecker constraint vec(X) ~ N(vec M, Σᶜ ⊗ Σʳ), with a per-position cell quota (1 − α flags per cell position), an observation step for casewise outliers, and reweighting. The estimator is the local minimiser reached from the pairwise Kronecker start `kpair`. Full definition: `notes/cellmmcd/METHOD.md` (§0, §2, §4, §4a–4e); paper form: `notes/cellmmcd/PAPER.md`. The estimator and algorithm have not changed since block 6.

**Algorithm.** Exact cell step (Schur identity det Σ_OO = det Σ · det Ω_MM: O(p²q + pq² + m³) per observation instead of O((pq)³)), observation step (kept-cell χ² test, quota, adaptive binomial count rule with a robust flag rate), GEM parameter step (flip-flop on the kept cells plus EM location step), λ-calibration of the scale, shadow cell step for trimmed matrices in the first 10 iterations, then reweighting (flags at quantile .99, re-decided at .999, refit with flags fixed). F decreases monotonically and converges (proven, tested). C++ cell step (RcppArmadillo) with cached inverses.

**Code.** `R/cellmmcd.R`, `src/cm_wstep.cpp`; diagnostics `cm_check_pattern` (Prop. A criterion) and `cm_check_strata` (stratum count, p, q ≤ 5). Tests: `tests/test_cellmmcd.R` (27) and `tests/test_core.R` (31). Experiments `exp/cm_step1.R` … `exp/cm_step26_strata.R`, all seeded; logs and seven figures in `output/`.

**Theory (status labels as in PAPER §4).**
- Proven: Schur form; exactness of the cell and observation steps; monotone GEM; scale calibration (1-D proxy); first-order shape bias of the raw fit (Σ − κ diag(1/Ω_cc), κ = 0.0745, iterated map 0.0111 vs 0.0115 observed); influence function (closed form in the proxy, chain-rule form in the Kronecker model, bounded by the marginal cut); no explosion breakdown by gross cells below 50 % per position; unboundedness of F under adversarial kept patterns (Prop. A: a row set T with between 1 and |T| − 1 kept T-partial columns); factor-wise coercivity under p fully kept columns (Prop. B); the exact rate F = ν log δ + O(1) along the coupled cross degenerations (Prop. C′), the stratum dimension formula dim N_i = r_i + s_i + comp(Γ_i) − 1 (Prop. D), and the generic count corollary (no (a, b)-descent when every stratum has Σ_i dim N_i ≥ |T| + |U| − 1).
- Checked numerically only: consistency of the reweighted default (n·KL 8.3 / 8.0 / 7.4 at n = 500 / 2000 / 8000; MLE 8.9 / 7.4 / 5.9); sensitivity curves; breakdown points (cellwise bias breakdown at the quota 25 % for 4-sd cells, none up to 45 % gross cells; casewise 25 %); slopes 0.98 / 0.98 / 4.00 for ν = 1 / 1 / 4.
- Conjectured: the full converse "F bounded below iff no stratum violates" (true for |T|, |U| ≥ 2 only as a conjecture; false for |T| = 1 or |U| = 1, where the strata are safe for a different reason); breakdown in general.

**Results (median KL to the truth, n = 100, 6-sd cells).** Cellwise 3×4 / 5×5: clean 0.08 / 0.19; 5 / 10 / 20 % cells 0.09 / 0.20, 0.09 / 0.18, 0.13 / 0.27. Competitors at 10 % cells: cellMCD-on-vec with Kronecker projection 0.26 / 0.59, MMCD 50 / 119, MLE 61 / 144. Recall 0.97–1.00, false flags 0.1–0.6 %. Casewise 5×5, 20 % shifted / scaled / structural matrices: 0.20 / 0.21 / 0.21 (MMCD 0.21 / 0.25 / 0.24); shifted 10 % + 5 % cells 0.22 (MMCD 52). Location error halved by the EM step. Timing (1 core): 10×10 0.4 s, 20×20 2.4 s (MMCD 4.7 s), 30×30 15 s (MMCD 16 s), 40×40 111 s (MMCD 46 s).

## 2. What is solid

- The objective, the Schur reduction, the monotone algorithm, and the convergence proof.
- The cellwise benchmark at small p, q: cellMMCD beats cellMCD-on-vec by 3× and the casewise methods by two orders of magnitude under 5–20 % cells, at the cost of 0.01–0.04 KL on clean data.
- The diagnosis and cure of the raw fit's shape bias (reweighting), confirmed by the large-n consistency run.
- The robustness theory that is proven: bounded influence, no explosion by gross cells, and the existence results listed above. Props. A, C′, D are short proofs and were checked numerically (Prop. D 40 / 40, slopes exact).
- Reproducibility: every experiment has a seed and a log; the test suites pass (block 9: 27 / 27 + 31 / 31; block 10 rerun: see `output/b10_test_*.log`).

## 3. What is shaky

1. **Existence at large p, q is an operating limit, not a theorem.** Under random cellwise contamination at rate ε the kept pattern falls into Prop. A's provably unbounded regime once nq(1 − ε)^p < p. At 20×20, n = 100 this happens at ε = 0.20 (16 fully kept columns, KL 3.7) and 0.25 (5 columns, KL 11.0); 40×40 is affected from ε ≈ 0.1. The algorithm still converges to a finite local minimum, but the estimator is then a local minimiser of an unbounded objective. This is the most consequential open item.
2. **The exact existence criterion is open.** Prop. A's converse is false (block 8: a 2×2 cross-mode counterexample, proven). The stratum count is a proven necessary condition for one family of descents and never fired on an algorithm-produced pattern (3×3 to 5×5), but descents outside the (a, b) family (subspaces, unequal rates, rotating eigenvectors) are not covered.
3. **Consistency is numerical.** The reweighted estimator has a residual bias of order κ(.999)², small but nonzero; a monotone consistent formulation is open.
4. **The weak case.** Structural casewise outliers plus 5 % cellwise contamination at 5×5 are not identifiable by the soft rule (KL 0.96; reweighted hard quota 0.37; MMCD 50). The case disappears at 8×8.
5. **Breakdown** values are numerical; the per-position quota 1 − α = 25 % is low, as for cellMCD.
6. The default observation rule (soft threshold vs reweighted hard quota) is the owner's call: they agree within 0.03 KL except in the weak case (`output/cm_step20_default.log`).

## 4. What failed (dead ends, all documented in JOURNAL.md)

- A fixed observation threshold with pq degrees of freedom (over-trims under cellwise contamination) and the observed-cell threshold alone (misses mostly-flagged shifted matrices), block 2.
- The censored parameter step: removes two thirds of the raw bias but is not monotone; kept as an option, block 3.
- The full-column safeguard against degenerate descents: does not prevent them at n ≤ 12 and breaks monotonicity, block 3.
- A monotone formulation of the shadow phase: impossible without charging flagged cells, block 5.
- A full coercivity proof of the joint (Σʳ, Σᶜ) problem (scale migrating between the factors), block 7; Conjecture C (partial-column counts are exact) refuted, block 8; the unrestricted count criterion C″ refuted on single-row / single-column strata, block 9.
- Certifying boundedness by a generic BFGS over log-Cholesky factors: stalls away from the boundary directions, block 8.
- A C++ parameter step: not worth it (5 % of run time), block 6.

## 5. Recommended next steps (in order)

1. **Quantify and remedy the large-p, q operating limit.** Measure KL of cellMMCD vs the clean MLE and MMCD at 20×20, n = 100 for ε = 0.05 … 0.25 and locate the degradation relative to nq(1 − ε)^p = p. Then try a lower bound on the eigenvalues of the shape factors (a ridge or a penalty τ Σ (log λ_j)²): cheap, coercive by construction; measure its bias at the Gaussian. The alternative is to show that the descent from `kpair` keeps ≥ p fully kept columns.
2. **Sufficiency of the stratum count for |T|, |U| ≥ 2 (C″).** Construct Prop.-A-silent patterns with a violating two-factor stratum (random 4×4 patterns, kept probability 0.7–0.8, n = 3–5), run the least-squares finder of step 26 with `full = TRUE`, and check the slope; a residual > 0 with `full = TRUE` is a counterexample.
3. **Descents outside the (a, b) family.** Prop. C′ for subspaces (A, B) follows the same Schur argument; then unequal rates δ vs δ^γ.
4. **Speed beyond 30×30:** profile the 40×40 run (cell step's m̄²pq term, `kpair`, R-level parameter step).
5. **Paper polish:** PAPER §6 item 1 should cite Props. C′ / D; add figure captions; state the operating range nq(1 − ε)^p ≳ p, np(1 − ε)^q ≳ q next to n ≥ 3 max(p, q).
6. **Owner decisions:** the default observation rule; whether the per-position quota should scale with p (fewer flags per column trades cellwise capacity for existence).

## 6. Spend (from LEDGER.csv)

| block | $ | | block | $ |
|---|---|---|---|---|
| 0 (planning) | 8.00 | | 5 | 9.62 |
| 1 | 8.50 | | 6 | 8.84 |
| 2 | 9.55 | | 7 | 9.19 |
| 3 | 7.60 | | 8 | 7.59 |
| 4 | 8.83 | | 9 | 8.77 |

Closed blocks 0–9: **$86.49** by the in-container meter (safety factor 1.15). Block 10 (this report) is appended to the ledger at close, about $2.44. The supervisor session's own launch overhead is booked only inside block 0's figure, so the account total is somewhat higher than the ledger sum; the launch rows record the supervisor's running estimate (about $83 spent before block 10 was launched).

## Files to read first
`notes/cellmmcd/PAPER.md` (whole), `notes/overnight/STATE.md` (next steps), `notes/cellmmcd/JOURNAL.md` block 9, `output/cm_step26_strata.log` (e), `output/cm_step20_default.log`.


## Addendum (block 11, 2026-11-05, bonus block on the owner's request)

**Question worked on:** STATE next-step 1, the large-$p,q$ existence limit (kept pattern in Prop. A's proven-unbounded regime once $nq(1-\varepsilon)^p<p$).

**What was found (all checked numerically, 20×20, $n = 100$, 6 sd random cells, 3 reps; `output/cm_step27_limit.log`, `cm_step27b_oracle.log`, `cm_step27c_countrule.log`):**
- cellMMCD median KL 2.7 / 2.9 / 3.0 / 3.4 / 3.9 / **9.1** at $\varepsilon$ = 0 / .05 / .10 / .15 / .20 / .25; the clean-data MLE is ≈ 2.1 at every $\varepsilon$; MMCD is 2.4 at $\varepsilon = 0$ and ≥ 1100 from 5 % on. Degradation is gradual (×1.8 the clean MLE at 20 %) and severe only at 25 %.
- The fitted factors never degenerate (condition numbers 9–10 against a true 8.6) even though Prop. A fires in 2 of 3 reps at 20 % and in all at 25 %. **The estimator returns a bounded, well-conditioned local minimiser inside the regime where the objective is provably unbounded.** The threshold $nq(1-\varepsilon)^p = p$ marks where the existence proof stops, not where the estimate breaks.
- At 25 % the oracle observed-cell MLE with the true flags gives KL 3.6–4.0, and switching observation trimming off barely changes cellMMCD (8.2 / 5.3 / 9.0 vs 9.1 / 5.5 / 9.1). The ×2.4 excess is the cell step missing 1.2–1.7 % of the 6 sd cells near the quota $1-\alpha = 0.25$. That is the real large-$p,q$ limit.
- **Remedy (i), a trace-proportional ridge on the shape factors, is rejected:** +28 % KL at the Gaussian for $\rho = 0.05$ and ×9 for $\rho = 0.2$ (the true column factor has condition number 28.5), and no gain at 25 %. It stays in the code as an option (`ridge = 0` default, test 28) and is not part of the method.

**Code:** `cm_ridge`, `ridge` argument (default 0, default estimator unchanged), test 28; three experiment scripts. Tests 28/28 + 31/31.

**Recommended next steps:** (1) the cell-step recall at 25 % (where are the missed cells; `quant_rw` 0.99; a second reweighting pass); (2) an existence statement for the algorithm's output (local minimiser reached from `kpair`), now well supported numerically.

**Spend:** block 11 ≈ $4.56 (meter at closing in `LEDGER.csv`); relay total was ≈ $84.8 before this block.
