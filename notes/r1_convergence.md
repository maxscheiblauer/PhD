# Convergence of R1: theory

Setting as in `trimmed_flip_flop_spec.md`, rule R1 (independent C-steps on $d_{ij}$ and $e_{ik}$), $h_r p = h_c q$ ($h_r = \beta n q$, $h_c = \beta n p$). Code: `R/r1.R`; experiments `exp/step7_*.R`, `exp/step8_r1_cycles.R`; results in `out/FINDINGS.md` (Q10).

## 1. Two criteria that are one likelihood without trimming

Give each criterion its share of the Jacobian (as in RJ). With the common log-det part
$$
B(\Sigma^r,\Sigma^c) = h_r\log|\Sigma^r| + h_c\log|\Sigma^c| = \beta n \log|\Sigma^c\otimes\Sigma^r|,
$$
define
$$
L_r^*(\Sigma^r,a;\Sigma^c) = B + \sum_{ij} a_{ij}\, d_{ij},\qquad
L_c^*(\Sigma^c,b;\Sigma^r) = B + \sum_{ik} b_{ik}\, e_{ik},
$$
and the profiled criteria $\tilde L_r = \min_a L_r^*$, $\tilde L_c = \min_b L_c^*$, both minimised over the overall scale $c$ of the product as well ($\tilde L = B + H\log(s/H) + H$, $H = h_r p$, $s$ = trimmed distance sum).

- Both are functions of the product $\Sigma^c\otimes\Sigma^r$ only: they do not depend on how the scale is split between the factors, nor on its overall level. (The spec's $L_r$ alone depends on the scale of $\Sigma^c$ through the $d_{ij}$, so it is not well defined as a function of the product.)
- They differ **only** in which distances are trimmed, so $\tilde L_r - \tilde L_c = H \log(s_d / s_e)$.
- **Without trimming** $\sum_{ij} d_{ij} = \sum_{ik} e_{ik} = \sum_i D_i^2$, so $\tilde L_r = \tilde L_c = $ the profiled $-2\log$-likelihood. Flip-flop is block-coordinate descent on that single function.
- The row step is one C-step of $\tilde L_r$ with $\Sigma^c$ fixed, so it never increases $\tilde L_r$. The column step never increases $\tilde L_c$. (The added $\log|\cdot|$ term of the other factor is constant in that step.)

## 2. R1 is a two-player game, not a descent method

Player $r$ chooses $(\Sigma^r, a)$ to minimise $\tilde L_r$, and player $c$ chooses $(\Sigma^c, b)$ to minimise $\tilde L_c$. R1 is alternating (Gauss–Seidel) **improvement dynamics**. A fixed point is a (C-step) Nash equilibrium: neither player can improve given the other. It is not a minimum of anything.

Alternating improvement is guaranteed to converge in a **potential game**, i.e. when there is a $\Phi$ with $\tilde L_r - \tilde L_c = g_c(\Sigma^c) - g_r(\Sigma^r)$, so that both criteria differ from $\Phi$ only by terms the player does not control. Without trimming, $\Phi = -2\log L$. With trimming, $\tilde L_r - \tilde L_c = H\log(s_d/s_e)$ couples both factors and both kept sets. There is no reason for it to separate, and it does not (see §4).

## 3. Local behaviour near a fixed point

Assume no ties at the trimming thresholds at a fixed point $\Sigma^*$. Then $a$ and $b$ are locally constant, and one R1 iteration is a smooth map $T$ of the shape of $\Sigma^c$ (det 1). $\Sigma^r$ is a function of $\Sigma^c$ once $a$ is frozen.

- **Local convergence holds iff $\rho = $ spectral radius of $DT(\Sigma^*) < 1$.** Then the error shrinks like $\rho^k$. `r1_jacobian()` computes $DT$ by central differences in trace-free log-coordinates.
- **Own gain is second order, cross pull is first order.** Let $\Delta$ be the step size. The row step moves $\Sigma^r$ from the minimiser for the previous $\Sigma^c$ to the minimiser for the current one, so its own gain is $O(\Delta^2)$. The column step changes $\tilde L_r$ by $\langle\nabla_{\Sigma^c}\tilde L_r, \Delta\Sigma^c\rangle + O(\Delta^2)$.
  - At a Nash point $\nabla_{\Sigma^c}\tilde L_c = 0$, but $\nabla_{\Sigma^c}\tilde L_r \neq 0$ in general. That is exactly the gap between a Nash point and a joint optimum. Without trimming the two criteria coincide, so this gradient vanishes.
  - Hence near convergence $|\text{cross pull}| / |\text{own gain}| \sim 1/\Delta \to \infty$.
  - So the answer to "do the own steps outpace the pull?" is **no, asymptotically the pull always dominates.** The individual criteria converge only because the steps shrink geometrically ($\rho < 1$), not because each step makes net progress.
  - The sign of the pull follows the direction of approach. With a negative leading eigenvalue of $DT$ (oscillating approach), $\tilde L_r$, $\tilde L_c$ and $\tilde J$ go up and down until the end. With a positive one, they can converge monotonically *from the wrong side*, increasing into the limit.
- **Empirically** (Step 7, 2520 runs from the MLE start):
  - In the late iterations the ratio $|\text{pull}|/|\text{own gain}|$ has median 300 ($\tilde L_r$) and 6600 ($\tilde L_c$), and is > 1 in 99–100 % of runs.
  - $\tilde J$ increases in 59 % of late half-steps.
  - At the reached fixed points, $\rho$ has median 0.13 and maximum 0.96, against 0.16 / 0.86 for the plain flip-flop at the MLE. This is selection-biased: only attracting fixed points are ever reached.

## 4. No Lyapunov function exists (counterexample)

The search in Step 8 finds clean Gaussian data ($p = q = 2$, $n = 17$, $\beta = .5$) on which R1, started at the MLE, enters an exact limit cycle of period 4. On the cycle, **every** half-step strictly decreases its own criterion (by at least 0.03), yet the state returns to itself.

Hence there is no function $V$ of the state that decreases strictly along R1, whether $\tilde J = \tilde L_r + \tilde L_c$, any weighted sum, or anything else. Along the cycle $V$ would have to satisfy $V(x) < V(x)$. In particular R1 has no ordinal potential, and **no convergence theorem for R1 can exist without extra conditions**. (The finite number of kept sets does not help: the cycle runs through four of them.)

## 5. What to use as a convergence summary

A descent objective is impossible (§4), so the right summary is a **residual**, which is zero exactly at fixed points:
$$
G(\Sigma) = \big[\tilde L_r(\Sigma) - \tilde L_r(C_r\Sigma)\big] + \big[\tilde L_c(\Sigma) - \tilde L_c(C_c\Sigma)\big] \;\ge\; 0,
$$
where $C_r$ / $C_c$ is one C-step of that factor from $\Sigma$, with the other factor held fixed (`r1_gap()`).
- $G = 0$ iff $\Sigma$ is a fixed point (ties aside). It is a lower bound on the Nash gap ("exploitability") of game theory.
- It shrinks like $\rho^{2k}$ under convergence and stays bounded away from 0 on a cycle.
- It is not monotone either, so it measures convergence but does not certify it.

$\tilde J$ remains a sensible score for **choosing among the fixed points** reached from different starts, which is how it is used now. It is not a convergence criterion.

Report per run: the status (converged / cycle / maxit), the final $G$, and the estimated rate $\rho$ (from $DT$, or from the ratio of successive step sizes).
