# Column-trimmed flip-flop estimation for matrix-variate covariances

> **Historical (5 Oct 2026).** This is the original spec that rules R0–R3 come from. The answers to its open questions (§8) are in `output/FINDINGS.md`, and the method now pursued, RJ, is in `notes/rj_objective.md`.

## 1. Model

Observations $X_1,\dots,X_n \in \mathbb{R}^{p\times q}$, i.i.d.

$$
\operatorname{vec}(X_i) \sim \mathcal N_{pq}\big(\operatorname{vec}(M),\; \Sigma^c \otimes \Sigma^r\big),
\qquad \Sigma^r \in \mathbb{R}^{p\times p},\ \Sigma^c \in \mathbb{R}^{q\times q}.
$$

Only $\Sigma^c\otimes\Sigma^r$ is identified, since $(\Sigma^r,\Sigma^c) \sim (c\,\Sigma^r, \Sigma^c/c)$. Below, $M = 0$ unless stated otherwise.

## 2. Whitened vectors

Let $A^c$ satisfy $A^c A^{c\prime} = (\Sigma^c)^{-1}$, and $A^r$ satisfy $A^{r\prime} A^r = (\Sigma^r)^{-1}$. The default is symmetric roots.

$$
Z_i = X_i A^c = [\,z_{i1},\dots,z_{iq}\,],\quad z_{ij}\in\mathbb{R}^p,
\qquad
U_i = A^r X_i = [\,u_{i1},\dots,u_{ip}\,]',\quad u_{ik}\in\mathbb{R}^q .
$$

At the true parameters, $z_{ij} \overset{iid}{\sim} \mathcal N_p(0,\Sigma^r)$ and $u_{ik} \overset{iid}{\sim} \mathcal N_q(0,\Sigma^c)$.

## 3. Distances

$$
d_{ij} = z_{ij}'(\Sigma^r)^{-1} z_{ij} \;\sim\; \chi^2_p,
\qquad
e_{ik} = u_{ik}'(\Sigma^c)^{-1} u_{ik} \;\sim\; \chi^2_q \qquad (\text{at true parameters}).
$$

Identity, which holds for any choice of roots:

$$
\sum_{j=1}^q d_{ij} \;=\; \sum_{k=1}^p e_{ik} \;=\; D_i^2 := \operatorname{tr}\!\big((\Sigma^c)^{-1}X_i'(\Sigma^r)^{-1}X_i\big) \sim \chi^2_{pq}.
$$

## 4. Weighted (trimmed) updates

Weights $a_{ij}\in\{0,1\}$ apply to the $\Sigma^r$ step and $b_{ik}\in\{0,1\}$ to the $\Sigma^c$ step. The retained counts are $h_r=\sum a_{ij}$ and $h_c=\sum b_{ik}$.

$$
\hat\Sigma^r \leftarrow \frac{\kappa_r}{N_r}\sum_{i=1}^n\sum_{j=1}^q a_{ij}\, z_{ij} z_{ij}',
\qquad
\hat\Sigma^c \leftarrow \frac{\kappa_c}{N_c}\sum_{i=1}^n\sum_{k=1}^p b_{ik}\, u_{ik} u_{ik}'.
$$

- **Denominator** $N_r \in \{h_r,\ nq\}$ and $N_c \in \{h_c,\ np\}$.
- **Consistency factor.** With retained fraction $\beta = h/N$ and dimension $m$ ($p$ for $\Sigma^r$, $q$ for $\Sigma^c$), the MCD-type factor is
  $$
  \kappa(\beta, m) = \frac{\beta}{F_{\chi^2_{m+2}}\!\big(\chi^2_{m,\beta}\big)}, \qquad \kappa \equiv 1 \text{ if no correction.}
  $$
- **Scale identification**, applied after each iteration: one of none, $\det\hat\Sigma^c = 1$, $\operatorname{tr}\hat\Sigma^c = q$, or $\hat\Sigma^c_{11}=1$.

With $a\equiv b\equiv 1$, $N_r = nq$, $N_c=np$ and $\kappa = 1$, this is the ordinary flip-flop MLE.

## 5. Selection rules

- **Concentration (C-step):** $a_{ij}=1$ for the $h_r$ smallest $d_{ij}$, and $b_{ik}=1$ for the $h_c$ smallest $e_{ik}$, recomputed every iteration.

The matching rules couple $a$ and $b$:

| Rule | $a$ ($\Sigma^r$ step) | $b$ ($\Sigma^c$ step) |
|---|---|---|
| R0 | C-step on $d_{ij}$ | $b\equiv 1$ (no trimming) |
| R1 | C-step on $d_{ij}$ | C-step on $e_{ik}$ (independent) |
| R2 = MMCD | $a_{ij}=s_i$ | $b_{ik}=s_i$, where $s_i=1$ for the $h$ smallest $D_i^2$ |
| R3 | C-step on $d_{ij}$ | $b_{ik}=1$ iff $X_i$ is not flagged by $a$ (e.g. $\sum_j a_{ij}\ge q-t$) |

## 6. Candidate objective (R1)

$$
J(\Sigma^r,\Sigma^c;a,b) = h_r\log|\Sigma^r| + \sum_{ij} a_{ij} d_{ij} \;+\; h_c\log|\Sigma^c| + \sum_{ik} b_{ik} e_{ik}.
$$

Each half is minimized by its own update, given the other factor. The two halves are not one likelihood, so there is no guarantee that $J$ decreases monotonically under alternation.

## 7. Known structural issues

- **Rotation non-invariance.** $Z_iQ$ is equally white for any orthogonal $Q$. With $a\equiv1$, $\hat\Sigma^r$ is invariant to $Q$; with trimming it is not. The set of "columns" $\{z_{ij}\}$ therefore depends on the arbitrary choice of root $A^c$.
- **Mixing.** $z_{ij}$ is a linear combination of all columns of $X_i$. Trimming $z_{ij}$ does not remove an original variable unless $\Sigma^c$ is diagonal.
- **Plug-in.** With $\hat\Sigma^c$ in place of $\Sigma^c$, the $z_{ij}$ are neither independent nor exactly $\chi^2_p$-distributed, so $\kappa$ is only approximate.
- **Scale invariance of the ranks (conjecture).** Multiplying $\hat\Sigma^r$ by $c_1$ and $\hat\Sigma^c$ by $c_2$ scales every $d_{ij}$ and $e_{ik}$ by $1/(c_1c_2)$. Rank-based C-steps are then unchanged, so $N$, $\kappa$ and the scale identification should affect only the scale of $\hat\Sigma^c\otimes\hat\Sigma^r$, not its shape or the selected sets.

## 8. Open questions

1. Does the alternation of updates and C-steps converge, or does it cycle? Under which rule, if any, is $J$ monotone?
2. Is the conjecture in §7 correct, i.e. do the normalization choices affect only scale?
3. How strongly does $\hat\Sigma^r$ depend on the root ($A^c$ symmetric vs. Cholesky vs. a random rotation)? Is there a canonical, data-driven root?
4. Is a matching rule necessary? Does R0 or R1 break in ways that R2 or R3 avoid?
5. Is there a single objective whose block-coordinate descent yields both updates under column-level trimming?
6. Is $\kappa(\beta,p)$ approximately correct for clean data at small $n$? How large is the plug-in bias?
7. Against which contamination type (observation-, column-, row-, cellwise) does column-level trimming beat MMCD, if any?
8. What is the minimal $h_r$ or $h_c$ for existence ($h_r \ge p$? plus conditions on $\Sigma^c$)?
9. How should $M$ be estimated under column-level trimming?