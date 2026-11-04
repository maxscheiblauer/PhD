# STATE: baton between overnight blocks

_Last updated: 2026-11-04 13:15 UTC by block 1 (session_01Fi3wVbqSCX1omqYNqN5phw). Each block rewrites this file at its end (PROTOCOL §6)._

## Current method candidate: cellMMCD (status: implemented, tested, clearly working)
- **Definition:** trimmed observed-cell likelihood with Kronecker Σ = Σᶜ⊗Σʳ, penalty λ_kl per flagged cell, quota of h = ⌈αn⌉ kept cells per position (the cellMCD objective with the Kronecker constraint). No imputation: the estimator is the marginal likelihood of the kept cells. Write-up `notes/cellmmcd/METHOD.md`, derivations and results `notes/cellmmcd/JOURNAL.md`.
- **Key device (proven, tested to 1e-15):** Schur identity. log det Σ_OO = q log|Σʳ| + p log|Σᶜ| + log det Ω_MM and x_O'Σ_OO⁻¹x_O = tr(ΩʳYΩᶜY') − u'Ω_MM⁻¹u with Ω_MM the m×m matrix [Ωʳ_{k_a k_b}Ωᶜ_{l_a l_b}]. Cost per observation O(p²q + pq² + m³), m = flagged cells of that matrix.
- **Algorithm:** start `kpair` (pairwise Kronecker: MAD scales + median polish on log s², GK correlations within rows/columns) → repeat {cell step: exact block minimisation per position with criterion Δ_ic = log(2πσ²_ic) + r²_ic vs λ_kl, budget n−h; parameter step: one flip-flop sweep of the score fixed point = GEM step}. F monotone (tested), W finite ⇒ converges. The parameter step uses conditional means only as the exact observed-data score (Fisher's identity); the journal records why no cheaper EM-free step exists. **Owner guideline 2 should be re-read in that light**: it is the plain score fixed point, cost O(pq) per observation, not O((pq)³).
- **Results (exp/cm_step1.R, n = 100):** recall 1.00 of 6-sd cells at ε = .05–.20, false flags 1–2 %. KL at 5×5, ε = .1: cellMMCD 0.32, cellMCD-on-vec 0.59, kpair 1.5, RJ 89, MMCD 119, MLE 144. Clean: 0.37 vs MLE 0.16, mostly scale (tr ratio 0.895).
- **Known defect:** clean-data scale ratio 0.895 (both shapes, n = 100 and 1000), 1.7 % flagged instead of 1 %: a shrink-and-flag-more loop. The 1-D factor 0.925 is not the answer (`exp/cm_step2_factor.R`).

## Code and tests (all passing at commit of this file)
- `R/cellmmcd.R`: `kpair`, `cm_obj`, `cm_delta1`, `cm_wstep`, `cm_pstep`, `cellmmcd`. Loaded by `R/load.R`.
- `tests/test_cellmmcd.R`: 8 tests (brute-force objective, brute-force cell criterion, flip-flop reduction, GEM monotonicity, full monotonicity, flags recall/fpr, beats MLE, kpair consistency). `tests/test_core.R`: 31/31.
- `exp/cm_step1.R` (benchmark; log `output/cm_step1.log`, rds has per-rep values incl. times), `exp/cm_step2_factor.R` (`output/cm_step2_factor.log`).
- `tools/overnight/hb.sh <status>` rewrites heartbeat.json with the meter (block number is hard-coded to 1: edit it).

## Prioritized next steps (concrete first actions)
1. **Consistency factor inside the iteration.** In `cellmmcd`, before each cell step use Σ/c as the working covariance for the residuals (c = scale ratio of the kept cells), and at the end multiply by 1/c. First action: set c = 0.895 empirically, rerun `exp/cm_step2_factor.R`, check that flagged → 1 % and tr ratio → 1. Then derive c(quant) (the fixed point of the shrink-flag loop under the Kronecker conditional law), mark proven/checked.
2. **Boundedness / existence under the quota.** Grid like FINDINGS step 10: p,q ∈ {3,5,8}, α ∈ {.5,.75}, small n; record when F → −∞ or singular Ω_MM. Conjecture: h ≥ max(p,q)+1 per position.
3. **Casewise + cellwise contamination** (brief asks "ideally casewise too"): with α = .75 the quota allows 25 % flagged per position, so a fully bad matrix can be flagged entirely; test `gen_data(type = "obs")` and mixed. If it fails, add an observation-level trim (MMCD-style) as a third block of the same objective.
4. **Location M**: joint update from the observed cells (score fixed point for M is the conditional-mean average); check the median plug-in effect under contamination.
5. **Speed**: Rcpp or vectorised m = 0 case for the cell step (10 s at 10×10, n = 200 in pure R); then larger p, q timings for the complexity section.
6. **Theory write-up**: Fisher consistency at the Gaussian for fixed λ (the shape looked consistent: shape KL 0.136 → 0.027 from n = 100 to 1000), breakdown under cellwise contamination (cellMCD paper's argument should transfer position-wise).
7. Plots for `output/` (none yet), then `notes/cellmmcd/METHOD.md` to paper form.

## Open questions and pitfalls
- λ is fixed from the `kpair` start (cellMCD convention); if the start is poor, λ is off. Consider recomputing λ once after the first parameter step (then F changes definition: document).
- `cm_wstep` loops in R over n·pq cells; fine up to pq ≈ 100, n ≈ 200.
- cellMCD on vec(X) prints verbosely; the benchmark now wraps it in `capture.output`.
- Don't compare raw factors; `scale_ident(..., "det")` is applied at the end of `cellmmcd`.
- `tff_multi(X, rule = "RJ", beta = .75)` already passes trace = FALSE; do not pass it again (silent failure of all starts).

## Files to read first
`notes/overnight/RESEARCH_BRIEF.md`, this file, `notes/cellmmcd/JOURNAL.md`, `notes/cellmmcd/METHOD.md`, `R/cellmmcd.R`, `output/cm_step1.log`.
