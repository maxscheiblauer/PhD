# Research brief: Kronecker covariance estimation under cellwise outliers

Author of the goal: the PhD student (repo owner). Written 2026-11-04. Scope: **matrix-variate only.**

## Goal
Develop a complete method for estimating a **Kronecker-structured covariance**, vec(X_i) ~ N(vec M, Σᶜ ⊗ Σʳ), from matrix observations X_i (p × q) that contain **cellwise outliers** (and ideally casewise outliers too), at **reasonable computational complexity and run time**.

"Complete" means:
- a well-defined estimator (objective or estimating equations);
- an algorithm with a convergence argument;
- a working implementation with tests;
- supporting theory: existence, consistency or consistency factors, and breakdown or other robustness properties, as far as they can be established;
- an honest simulation comparison;
- complexity and timing results.

## Guidelines from the owner (binding)
1. **Avoid imputation of cellwise outliers.** Statistically it is messy. Do not build the method around replacing flagged cells by predicted values.
2. **EM algorithms** can work here but are extremely costly computationally. Try one **only if no better idea is found**, and document why the alternatives failed.
3. **Develop your own code.** Marcus Mayrhofer's CRAN package `robustmatrix` (0.1.5: MMCD, flip-flop MLE, MMD, Shapley values; source at https://github.com/cran/robustmatrix, C++ in `src/functions.cpp`) is **not a priority** and need not be the code base.
   - Use it where it is genuinely useful, e.g. as the MMCD benchmark or for borrowed pieces.
   - Equally use anything useful found elsewhere: other packages, papers, algorithms.
4. **Start from the owner's preliminary work:** the trimmed-likelihood estimators for Kronecker covariances in this repo, in particular **RJ**, the joint trimmed likelihood on whitened columns with a per-column quota.
5. Matrix-variate only. Tensors and functional data are out of scope.

## Starting material in this repo
- `README.md` and `output/FINDINGS.md`: status of RJ, R1 and the other trimming rules; what works and what fails.
- `notes/rj_objective.md`, `notes/r1_convergence.md`, `notes/trimmed_flip_flop_spec.md`, `notes/from_matrix_to_tensor_mcd.md`: theory notes.
- `R/` (matnorm.R, tff.R, metrics.R, data.R), `tests/test_core.R` (31 tests), `exp/step*`, `exp/RJ_simple.R`.
- `references/`: literature list and reading reports.
- `notes/00_setting_matrix_cellwise.md`: literature map of the cellwise and matrix-variate work from Leuven and TU Wien, MMCD implementation notes, gap analysis, and some verified Kronecker identities.
  Note: identities 3 and 5 there are phrased as conditional-mean computations. Under guideline 1, read them only as tools for evaluating observed-data quantities, if at all.

## Deliverables expected by the morning (as far as the budget allows)
- `notes/cellmmcd/`: method write-up (paper-style LaTeX-in-Markdown is fine). Include the definition, the algorithm, complexity, the theory with what is proven vs conjectured, and the open problems.
- `R/` + `tests/`: implementation and tests (Rcpp only if it pays off).
- `exp/` + `output/`: reproducible experiments with logs and plots, and seeds recorded.
- `notes/overnight/PROGRESS_SUMMARY.md` and the final `notes/overnight/MORNING_REPORT.md`: what was achieved, what failed, and what to do next.
