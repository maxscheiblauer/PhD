# RJ: the unified objective, its relation to the MLE, and consistency

Setting and notation as in `trimmed_flip_flop_spec.md`; code `exp/RJ_simple.R`. Throughout $M = 0$ (in the code $M$ is the coordinatewise median, which is consistent at the model by symmetry).

## 1. Derivation

**Step 1: the matrix-normal likelihood, written in whitened columns.** Let $A \in \mathbb{R}^{q\times q}$ be invertible with $AA' = (\Sigma^c)^{-1}$, and $Z_i = X_iA = [z_{i1},\dots,z_{iq}]$. Then
$$
\operatorname{vec} Z_i = (A'\otimes I_p)\operatorname{vec}X_i \sim \mathcal N_{pq}\big(0,\ A'\Sigma^cA \otimes \Sigma^r\big) = \mathcal N_{pq}(0,\ I_q\otimes\Sigma^r),
$$
so the columns $z_{i1},\dots,z_{iq}$ are i.i.d. $\mathcal N_p(0,\Sigma^r)$. This holds for **every** root $A$, not just the symmetric one. The Jacobian of $X_i\mapsto Z_i$ is $|\det A|^p$, so
$$
f(X_i) = |\det A|^{p}\prod_{j=1}^q \varphi_p(z_{ij};\Sigma^r).
$$

**Step 2: give each column its share of the Jacobian.** Spread the factor $|\det A|^p$ evenly over the $q$ columns of $X_i$. Then the matrix density is a product of $nq$ "column likelihoods":
$$
-2\log L = \sum_{i=1}^n\sum_{j=1}^q \ell_{ij},\qquad
\ell_{ij} = \log|\Sigma^r| - \tfrac{2p}{q}\log|\det A| + d_{ij},\qquad d_{ij} = z_{ij}'(\Sigma^r)^{-1}z_{ij}
$$
(up to constants). Because $-2\log|\det A| = \log|\Sigma^c|$, each column carries $\log|\Sigma^r| + \tfrac pq\log|\Sigma^c|$. Summed over all units this is $nq\log|\Sigma^r| + np\log|\Sigma^c|$, the usual MLE log-determinant part. The quadratic part sums to $\sum_j d_{ij} = D_i^2$.

**Step 3: trim units, not observations.** Let $a_{ij}\in\{0,1\}$ mark the kept columns. The **RJ objective** is the trimmed version of the same sum:
$$
\boxed{\;J^*(\Sigma^r, A;\, a) = \sum_{i,j} a_{ij}\Big[\log|\Sigma^r| - \tfrac{2p}{q}\log|\det A| + d_{ij}\Big],\qquad \sum_{i} a_{ij} = h = \lfloor\beta n\rfloor \ \ \forall j\;}
$$
With the quota, $H = hq$ and $J^* = H\log|\Sigma^r| - 2ph\log|\det A| + \sum_{a_{ij}=1} d_{ij}$, which is `rj_crit()`. RJ minimizes $J^*$ jointly over $(a, \Sigma^r, A)$.

**Step 4: the three block minimizers** (each exact, so $J^*$ never increases):
- *C-step.* Given $(\Sigma^r, A)$, $J^*$ is linear in $a$ with coefficients $d_{ij}$ plus terms that are constant under the quota. In each column $j$, keep the $h$ smallest $d_{ij}$.
- *$\Sigma^r$ step.* $J^*$ is the Gaussian $-2\log L$ of the $H$ kept $p$-vectors, so $\Sigma^r = H^{-1}\sum_{a_{ij}=1} z_{ij}z_{ij}'$.
- *A step, one column at a time.* With $V_j = \sum_{i: a_{ij}=1} X_i'(\Sigma^r)^{-1}X_i$,
  $$
  J^* = -2ph\log|\det A| + \sum_j a_j'V_ja_j + \text{const}.
  $$
  With the other columns fixed, $\det A = c\cdot w_j'a_j$, where $w_j$ is column $j$ of $A^{-T}$. Setting the gradient $-2ph\,w_j/(w_j'a_j) + 2V_ja_j$ to zero gives $a_j = \lambda V_j^{-1}w_j$. Minimizing $-2ph\log\lambda + \lambda^2 w_j'V_j^{-1}w_j$ then gives $\lambda = \sqrt{ph / w_j'V_j^{-1}w_j}$.

## 2. How RJ differs from the MLE

| | MLE (flip-flop) | RJ |
|---|---|---|
| Unit of the likelihood | observation $X_i$ ($n$ units) | whitened column $z_{ij}$ ($nq$ units, Jacobian shared) |
| What is trimmed | nothing (MMCD: whole $X_i$) | $(1-\beta)n$ columns **per column index** $j$ |
| Parameter | $(\Sigma^r, \Sigma^c)$ | $(\Sigma^r, A)$, with $q^2$ entries instead of $q(q+1)/2$ |
| Invariance $A \to AQ$, $Q$ orthogonal | yes | only if $a\equiv 1$ |
| $\Sigma^c$ update | $\Sigma^c = \frac{1}{np}\sum_i X_i'(\Sigma^r)^{-1}X_i$ | column-wise $a_j \propto V_j^{-1}w_j$, with a different $V_j$ per column |
| Objective | likelihood, concave in $(\Sigma^r)^{-1}$ and $(\Sigma^c)^{-1}$ separately | combinatorial (C-step), many local minima |
| Scale of $\hat\Sigma$ | unbiased up to $O(1/n)$ | shrunk by $1/\kappa(\beta,p)$, needs correction |

Points to stress:
- **$\beta = 1$ recovers the MLE exactly.** Then all $V_j = V = \sum_iX_i'(\Sigma^r)^{-1}X_i$ and $J^* = -np\log|AA'| + \operatorname{tr}(V AA') + \dots$ depends on $A$ only through $AA'$. Its minimizer $AA' = npV^{-1}$ is the flip-flop column update, and the C-step is void. RJ is therefore block-coordinate descent on a **trimmed likelihood** that reduces to the matrix-normal likelihood at $\beta = 1$. It is not an ad hoc pair of updates.
- **Unified, unlike R1.** R1 has two criteria (one in $d_{ij}$, one in $e_{ik}$) that differ by $H\log(s_d/s_e)$, and it can cycle (`r1_convergence.md` §4). RJ has a single criterion. The information on $\Sigma^c$ enters through the shared Jacobian term $-\tfrac{2p}{q}\log|\det A|$ and through $d_{ij}$, not through a second set of distances. Hence the iteration is monotone and the kept set is finite, so it terminates.
- **The root becomes a parameter.** Once $a \not\equiv 1$, which vectors count as "columns" depends on $A$ itself, not only on $AA'$. The trimmed objective can prefer a rotation whose columns separate the outlying directions. This is the source of RJ's gain under column-like contamination, and also of its start and root dependence (FINDINGS Q3, Q5).
- **The quota is part of the objective, not a tuning detail.** If column $j$ kept $h_j = 0$ units, nothing would bound $\|a_j\|$ and $-\log|\det A| \to -\infty$. With $h_j = h$, scaling $a_j \to ta_j$ changes $J^*$ by $-2ph\log t + (t^2-1)\sum_i a_{ij}d_{ij}$, which is bounded below. The quota is also what makes the population argument below go through.
- **Efficiency.** Clean-data KL is 0.95 for RJ at $\beta = .5$, against 0.15 for the MLE (`output/rj_simple.txt`). This is the usual price of trimming half the units, made somewhat larger by the extra rotation parameters.

## 3. Why it should still be consistent

The argument is Fisher consistency: the population version of $J^*$ is minimized at the truth, up to the rotation orbit and a known scale. Let $X$ follow the model with $(\Sigma^r_0, \Sigma^c_0)$ and $\Sigma_0 = \Sigma^c_0\otimes\Sigma^r_0$. The population objective (per observation) is
$$
J(\Sigma^r, A) = \sum_{j=1}^q \Big\{\beta\log|\Sigma^r| - \tfrac{2p}{q}\beta\log|\det A| + \mathbb E\big[d_j\,\mathbb 1\{d_j \le t_j\}\big]\Big\},
$$
with $d_j = (Xa_j)'(\Sigma^r)^{-1}(Xa_j)$ and $t_j$ the $\beta$-quantile of $d_j$. This is the C-step in the population, applied per column under the quota.

**(i) Every column has the right shape, whatever $A$ is.** $z_j = Xa_j \sim \mathcal N_p(0,\ s_j\Sigma^r_0)$ with $s_j = a_j'\Sigma^c_0a_j$. Write $B = (\Sigma^r_0)^{-1/2}\Sigma^r(\Sigma^r_0)^{-1/2}$ and $w\sim\mathcal N_p(0,I)$. Then $d_j \overset{d}{=} s_j\,w'B^{-1}w$. Scaling by $s_j$ does not change which $\beta$-fraction is smallest, so
$$
J = q\beta\log|B| + \Big(\sum_j s_j\Big)\,T(B) - 2p\beta\log|\det A| + \text{const},\qquad T(B) = \mathbb E\big[w'B^{-1}w\,\mathbb 1\{w'B^{-1}w \le \text{its }\beta\text{-quantile}\}\big].
$$
Only the marginals of the $z_j$ enter, so their correlation for a "wrong" $A$ plays no role.

**(ii) Profile out $\Sigma^r$: this is the MCD at $\mathcal N(0,I)$.** For fixed $\bar s = \sum_j s_j$, minimizing $q\beta\log|B| + \bar s\,T(B)$ over $B$ is the trimmed-likelihood form of the MCD, applied to a spherical normal. Profiling $B$ for a fixed set gives $\log|\operatorname{cov}(\text{set})|$, which is minimized by the centred ball (Butler, Davies & Jhun 1993). Hence $B = bI$, i.e. $\Sigma^r \propto \Sigma^r_0$. With $T(bI) = p\,F_{\chi^2_{p+2}}(\chi^2_{p,\beta})/b =: \tau/b$, the optimal $b$ is $\bar s\tau/(pq\beta)$, and
$$
\min_{\Sigma^r} J = p\beta\Big[q\log\textstyle\sum_j a_j'\Sigma^c_0a_j - 2\log|\det A|\Big] + \text{const}.
$$

**(iii) Minimize over $A$: AM–GM.** Put $G = (\Sigma^c_0)^{1/2}A$. The bracket becomes $q\log\operatorname{tr}(G'G) - \log\det(G'G) + \text{const}$. By AM–GM on the eigenvalues of $G'G$, $\det(G'G)^{1/q} \le \operatorname{tr}(G'G)/q$, with equality iff $G'G = cI$. So the minimizers are exactly
$$
A = c\,(\Sigma^c_0)^{-1/2}Q,\qquad Q \text{ orthogonal}\quad\Longrightarrow\quad \Sigma^c = (AA')^{-1} \propto \Sigma^c_0 .
$$

**(iv) The scale.** At the minimizer, $d_j \sim \chi^2_p/(c_rc_c)$, where $\Sigma^r = c_r\Sigma^r_0$ and $\Sigma^c = c_c\Sigma^c_0$. Both first-order conditions ($\Sigma^r$ step: $\Sigma^r = \beta^{-1}\mathbb E[z_jz_j'\mathbb 1]$; A step: $\mathbb E[d_j\mathbb 1] = p\beta$) give the same product,
$$
c_rc_c = F_{\chi^2_{p+2}}(\chi^2_{p,\beta})/\beta = 1/\kappa(\beta,p).
$$
Hence $\kappa(\beta,p)\,\hat\Sigma^c\otimes\hat\Sigma^r \to \Sigma_0$ at the population level. This is why the code applies $\kappa$ with dimension $p$, the dimension of the trimmed unit, and applies it once to the product.

**Reading of the result.**
- The **product** $\Sigma^c\otimes\Sigma^r$ is identified: the unique population minimizer, after $\kappa$.
- The **root** is not: $J$ is exactly flat along $A_0Q$ at the model. In finite samples, noise picks the rotation, which explains why 10 starts give 10 different kept sets with nearly the same estimate. This does not harm consistency of $\Sigma^c$, because $AA'$ is constant on the orbit.
- Same estimand as the MLE: at the model both target $\Sigma_0$. They differ in efficiency and in robustness, not in what they estimate. The difference from MMCD is $\kappa(\beta,p)$ instead of $\kappa(\beta,pq)$, because RJ's trimmed unit is a $p$-vector, not a $pq$-vector.

## 4. What is still missing for a full proof

1. **Sample to population.** We need $\hat\Sigma \to \Sigma_0$ a.s. This requires a uniform LLN for trimmed sums of quadratic forms over $(\Sigma^r, A)$ in a compact set (Butler et al. 1993; Cuesta-Albertos, Gordaliza & Matrán 2008 for trimmed likelihoods), plus well-separatedness of the minimum **modulo $Q$**. Compactness needs the quota. Empirically (step 10, FINDINGS §2) RJ exists and converges iff every column index keeps $h \ge q+1$ matrices, independent of $p$, so the condition is on $q$, not on $h \ge p$ as first guessed. Proving that $h \ge q+1$ keeps every $V_j$ positive definite along the minimising sequence is open.
2. **Plug-in location.** Above $M$ is known. With the coordinatewise median, $\hat M - M = O_p(n^{-1/2})$ should not affect the limit, but this is unchecked.
3. **Local vs global.** Consistency is a statement about the **global** minimizer of $J^*$. The algorithm finds local minima (FINDINGS Q5), so in practice multiple starts are needed, picking the lowest $J^*$, as for FastMCD.
4. **Finite-sample scale.** Step 9 (FINDINGS §3, n up to 1000) supports the Fisher-consistency argument. With $\Sigma^c$ known, $\kappa(\beta,p)$ gives the exact scale. The factor the data ask for converges to $\kappa(\beta,p)$, not to $\kappa(\beta,pq)$. The plug-in RJ scale is too small by 0.6–4.4 % at $n = 1000$, and the gap shrinks like $n^{-1/2}$. It grows with $q$ and with the trimmed share, which is consistent with the rotation of $A$ being fitted to the trimmed set (optimism of order $n^{-1/2}$). A small-sample correction $c(n,p,q,\beta)$, as for the MCD, is the remaining practical step.
5. **Contaminated model.** Fisher consistency at $F_0$ says nothing about bias under contamination. The robustness claim (columnwise breakdown under the quota) is separate work.
