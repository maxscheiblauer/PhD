# STATE: baton between overnight blocks

_Last updated: 2026-11-04 13:58 UTC by block 3 (session_01Fi3wVcka367SwTSdePTf7w). Each block rewrites this file at its end (PROTOCOL §6)._

## Current method: cellMMCD, reweighted (status: implemented, tested 17/17 + core 31/31, benchmarked cellwise, casewise and mixed with the final default)
- **Definition** (`notes/cellmmcd/METHOD.md` §2, §4, §4c, §4d): trimmed observed-cell likelihood with Kronecker Σ, penalty λ_kl per flagged cell (calibrated, block 2), quota h = ⌈αn_v⌉ per position, observation weights v_i with the observed-cell threshold μ_i and the majority rule (block 2), **plus (block 3) a reweighting step**: at the raw fit every cell of a kept matrix is re-decided at χ²_{1,.999} on its residual given the other kept cells, then the parameters are refitted with (v, W) fixed (EM, monotone), scale divided by t_rw = F₃(c_rw)/F₁(c_rw) and t_obs. `reweight = FALSE` gives the raw estimator (block 2's default), `$raw` holds it.
- **Why (block 3, proven to first order, checked numerically):** the raw fixed point at quant .99 is biased in shape, Σ − κ diag(1/Ω_cc), κ = F₁(c) − F₃(c) = 0.0745: a flagged clean cell is a tail cell but is treated as missing. KL → 0.0115 (3×4) as n → ∞; iterated population map predicts 0.0111 (`exp/cm_step6b_biasformula.R`). Reweighted: KL ratio to the MLE 1.01 / 1.09 / 1.49 at n = 500 / 2000 / 8000 (`output/cm_step6_bias.log`, plot `cm_step6_bias.png`).
- **Results with the final default (p = q = 5, n = 100, median KL):** clean 0.205 (MMCD 0.18, MLE 0.16); 10 % cells at 6 sd 0.209 (cellMCD on vec 0.59, MMCD 119); 20 % cells 0.353; 20 % shifted / scaled / structural matrices 0.21 / 0.22 / 0.29 (MMCD 0.21 / 0.27 / 0.23); 10 % shifted + 5 % cells 0.23 (MMCD 51). 3×4 clean 0.095 (MLE 0.070). Recall 0.97–1.00 at 6 sd, 0.96 at 3 sd; false flags 0.1–0.6 %.
- **Other block-3 options (off by default):** `tail = TRUE` censored parameter step (`cm_pstep(tau)`; Fisher-consistent in the 1-D proxy, removes 2/3 of the bias in the Kronecker model, not monotone in F); `safeguard = TRUE` full-column repair (`cm_fullcol`; **failed** to prevent degenerate descents, breaks monotonicity, `output/cm_step7_safeguard.log`).

## Code and tests (all passing at commit of this file)
- `R/cellmmcd.R`: sections 1 kpair, 2 objective, 3 cell step, 3a `cm_fullcol`, 3b `cm_vstep`, 4 `cm_pstep(tau)`, 5 driver `cellmmcd(..., tail, reweight, quant_rw, safeguard)`, 5b `cm_r2`, `cm_reweight`, 6 `cm_proxy_fixpt(tau)`.
- `tests/test_cellmmcd.R`: 17 tests (block 2's 13 + censored proxy, reweighting, tail runs, safeguard). `tests/test_core.R`: 31/31.
- Experiments: `exp/cm_step1.R` (cellwise benchmark; log rerun with the final default, raw-.99 rerun in `cm_step1_raw99.log`, block 1 in `cm_step1_block1.log`), `cm_step4_mixed.R` (rerun; block 2 in `cm_step4_mixed_block2.log`), `cm_step6_bias.R`, `cm_step6b_biasformula.R`, `cm_step7_safeguard.R`, `cm_step8_default.R`, `cm_plot_bias.R`. Logs in `output/`.
- `tools/overnight/hb.sh`: block number and session id hard-coded (set to 3: edit for block 4).

## Prioritized next steps (concrete first actions)
1. **Interplay of the two trimming steps (pitfall from `cm_step8`):** with 10 % cells at 3–4 sd in every matrix, the observation step trims 25 % of the matrices (all variants). First action: in `exp/cm_step8_default.R` add a variant `alpha_obs = 1` and one with `quant_obs = 0.999`, compare KL and the trimmed fraction; then consider evaluating μ_i's quantile with a per-matrix allowance for unflagged moderate cells, or trimming only when the matrix also fails the majority rule at a lower cutoff.
2. **Monotone consistent formulation** (METHOD §8.2): censored-cell objective with the cutoff region fixed at the initial estimate, F_cens = Σ_kept log-density + Σ_flagged −2 log P_Σ(cell beyond its cutoff | kept cells). Write it down in JOURNAL first; the 1-D proxy identity F₁(c) − F₃(c) = 2c f₁(c) shows its stationary point is at the truth. Cell step and parameter step would both act on F_cens. Only if the reweighting device is judged insufficient for the paper.
3. **Degenerate descents at n ≤ 12:** the full-column repair failed. Check first whether the "degen" indicator in `cm_step5/7` flags legitimate fits (compare with the MLE's eigenvalue ratio at n = 6–12); if real, try a lower bound on the eigenvalues of Σʳ, Σᶜ relative to the start (a constrained parameter step), or simply document n ≥ 3·max(p,q) as the operating range (no degeneracy observed at n ≥ 40).
4. **Location M** jointly from the kept cells (score fixed point); check the median plug-in under 20 % shifted matrices (`cm_step4`: shift already handled by trimming, so low priority).
5. **Speed:** `cm_wstep` and `cm_reweight` loop in R over n·pq cells; vectorise the m = 0 case (closed form σ² = 1/(Ωʳ_kk Ωᶜ_ll)) or Rcpp; timings for 10×10, 20×20 (step-1 log has one 10×10 timing).
6. **Hard-quota variant:** apply κ(α, pq) (`kappa_mcd` in `R/matnorm.R`) to its scale.
7. Breakdown / influence: the reweighted estimator's influence function in the 1-D proxy is cheap to derive (ratio estimator at cutoff c_rw with flags from the .99 fit); write it in METHOD §7.
8. METHOD.md to paper form (the counting proof and the bias proposition deserve clean statements); more plots (`cm_step8` comparison, benchmark bars).

## Open questions and pitfalls
- The reweighted estimator keeps a κ(.999)²-small bias (predicted shape KL 0.0003 at 3×4): fine for practice, not exactly consistent.
- Casewise shifts of 2 sd in every cell are missed by every variant (recall 0.25, `cm_step8`); 5 sd are caught 100 %.
- t_obs is evaluated at d = pq for all matrices (approximation for matrices with flagged cells).
- The cellMCD paper (Raymaekers & Rousseeuw 2024) was still not consulted (network): check whether they discuss the bias of §4d (they re-estimate at a higher quantile? unknown) and the unbounded objective of §4a.
- `cm_step5/7` "degen" indicator: eigenvalue ratio < 1e-6 or F < F₀ − 50n; may be too crude at n ≤ 12.
- Don't compare raw factors; `scale_ident(..., "det")` is applied at the end. `pkill -f <script>` from a Bash call kills the tool's own shell; use `until grep -q ... ; do sleep; done` to wait for background runs.
- Each block's start-up reading costs ≈ $3 of the cap; keep reads to this file, the brief, METHOD §4d, and the latest journal entry.

## Files to read first
`notes/overnight/RESEARCH_BRIEF.md`, this file, `notes/cellmmcd/METHOD.md` (§2, §4d, §7–9), the block-3 entries of `notes/cellmmcd/JOURNAL.md`, `R/cellmmcd.R` (sections 5 and 5b), `output/cm_step8_default.log`, `output/cm_step6_bias.log`.
