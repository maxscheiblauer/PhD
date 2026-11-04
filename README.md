# PhD: robust covariance estimation for matrix-variate data

**Research question.** Data are matrices X_i (p × q), and the covariance has Kronecker structure, vec(X_i) ~ N(M, Σᶜ ⊗ Σʳ). The MMCD trims **whole observations**. Can we get robustness against contamination of **individual columns** (or rows) of X_i by trimming *slices* of the whitened matrices instead?

**Status (15 Oct 2026).** Full results are in [`output/FINDINGS.md`](output/FINDINGS.md).

| Method | What it trims | Status |
|---|---|---|
| **RJ** (joint trimmed likelihood) | whitened columns z_ij, βn kept per column index | **lead method.** It has one objective, converges monotonically, κ(β,p) is the right consistency factor, it exists for h ≥ q + 1, and it beats MMCD under column-type contamination. |
| R1 (two independent C-steps) | column and row slices | **dropped.** There is no objective, and it can cycle even on clean data (proven counterexample). Archived. |
| R0, R3 | columns (R3: flags whole X_i) | rejected: R0 does not protect Σᶜ; R3 is inefficient and cycles |
| R2 = MMCD | whole observations | reference (`robustmatrix::mmcd`) |
| reweighting | — | not pursued: with slice trimming a reweighting step would need a different algorithm (EM-type) |

**Theory:** [`notes/rj_objective.md`](notes/rj_objective.md) (RJ derivation, link to the MLE, Fisher consistency), [`notes/r1_convergence.md`](notes/r1_convergence.md) (why R1 cannot have a convergence theorem). The original spec and plan are kept as historical notes in `notes/`.

## Layout

```
R/            core functions: matnorm.R (sampler, flip-flop MLE, roots, kappa), tff.R (all trimming rules incl. RJ),
              metrics.R (KL, shape-KL, precision/recall), data.R (contamination types), load.R (sources all)
R/archive/    r1.R (fast R1 with cycle detection and diagnostics; supports the negative result)
tests/        test_core.R (31 unit tests)
exp/          step1..step6 (spec questions Q1-Q7), step9 (RJ consistency), step10 (RJ existence),
              RJ_simple.R (self-contained RJ + cellwise study)
exp/archive/  R1-only experiments (step1b, step7, step8, R1_simple, R1_cellwise, r1_playground)
output/       logs (*.log, *.txt), plots (*.png), FINDINGS.md; *.rds results are git-ignored, rerun to get them
notes/        theory notes; references/ literature list + reading reports (PDFs not tracked)
```

## Reproduce

R 4.6.1 with the CRAN package `robustmatrix` (0.1.5). Run everything from the project root. The exact commands and timings are at the end of `output/FINDINGS.md`. The quick check is:

```
Rscript tests/test_core.R
```

The code in `R/`, `tests/` and `exp/step*` was lost on 7 Oct and restored on 15 Oct by replaying the Claude Code session transcripts of 5 Oct. Restored files were checked byte for byte against file dumps in those transcripts where available. All 30 original tests pass, and steps 1–6 reproduce the original numbers exactly. The R1 notebook `notebooks/trimmed_flip_flop.Rmd` was not restored.

## Setup

Open this folder in VS Code with the [R extension](https://marketplace.visualstudio.com/items?itemName=REditorSupport.r) installed (`PhD.Rproj` serves as the project marker).
