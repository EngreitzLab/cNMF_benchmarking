# ProgramDE

## Overview
A pipeline for testing differential effect (DE) of gene targets on gene programs. This pipeline implements a Conditional Randomization Test (CRT) to assess the statistical significance of perturbations.

The core analysis pipeline lives in the `src.sceptre` module. Github for this test: https://github.com/edtnguyen/programDE/tree/main (Written by Tri Nguyen) 

---

## Method overview

This pipeline tests **target gene → program usage** effects using a Conditional Randomization Test (CRT) with optional skew-normal calibration.

### Notation

* Cells: $i = 1,\dots,N$
* Programs: $k = 1,\dots,K$ (here $K \approx 70$)
* Covariates: $C \in \mathbb{R}^{N \times p}$ (includes an intercept column)
* cNMF usage (composition): $u_i \in \Delta^{K-1}$ with $\sum_{k=1}^K u_{ik}=1$ and $u_{ik}\ge 0$

## 1) Linear model for gene effect on program usage

### 1.1 Gene-level perturbation indicator (union)

For each target gene $g$, define a **gene-level union indicator**:

```math
x_i \in \{0,1\},\qquad
x_i = \mathbf{1}\{ \exists\ \mathrm{guide\ targeting}\ g\ \mathrm{in\ cell}\ i \}.
```

### 1.2 CLR transform of program usage

Because program usages are compositional (each row sums to 1), we apply a **centered log-ratio (CLR)** transform. After flooring and renormalization (to avoid $\log 0$), define:

```math
u'_{ik}=\frac{\max(u_{ik},\varepsilon)}{\sum_{j=1}^K \max(u_{ij},\varepsilon)}.
```

Then the CLR-transformed outcome is:

```math
Y_{ik}=\mathrm{CLR}(u'_i)_k
= \log u'_{ik} - \frac{1}{K}\sum_{j=1}^K \log u'_{ij}.
```

Collect outcomes in $Y \in \mathbb{R}^{N \times K}$.

### 1.3 Per-program linear regression

For each program $k$, fit:

```math
Y_{ik} = \beta_k\, x_i + C_i^\top \gamma_k + \varepsilon_{ik}.
```

where $\beta_k$ is the gene effect on program $k$ (on the CLR scale), and $\gamma_k$ are covariate effects.

**Test statistic:** the OLS coefficient $\hat\beta_k$.

## 2) Conditional Randomization Test (CRT)

The CRT tests $H_0: Y \perp x \mid C$ by resampling $x$ from its conditional distribution given covariates.

### 2.1 Propensity model

Fit a regularized logistic regression for the union indicator:

```math
p_i = \mathbb{P}(x_i = 1 \mid C_i)
= \sigma\!\left(C_i^\top \theta\right),
\qquad
\sigma(t)=\frac{1}{1+e^{-t}}.
```

### 2.2 Null resampling of $x$

Generate $B$ synthetic perturbation vectors:

```math
\tilde x_i^{(b)} \sim \mathrm{Bernoulli}(p_i),
\qquad b=1,\dots,B,\ i=1,\dots,N.
```

**Efficient Bernoulli resampling (index sampler).** Instead of drawing $B$ Bernoullis for every cell, we use an equivalent two-stage procedure per cell $i$:

1. Draw how many resamples include cell $i$ as treated:

```math
M_i \sim \mathrm{Binomial}(B, p_i).
```

2. Sample $M_i$ distinct resample indices uniformly without replacement:

```math
S_i \subset \{1,\dots,B\},\quad |S_i|=M_i,\quad S_i\ \mathrm{uniform}.
```

3. Set $\tilde x_i^{(b)}=1$ for $b\in S_i$ and $\tilde x_i^{(b)}=0$ otherwise.

This yields exactly the same distribution as i.i.d. Bernoulli draws across $b$, but is faster when $p_i$ is small (sparse perturbations).

### 2.3 Recompute the test statistic under the null

For each resample $b$, compute OLS coefficients $\hat\beta_k^{(b)}$ for all programs using $\tilde x^{(b)}$:

```math
Y_{ik} = \beta_k^{(b)}\, \tilde x_i^{(b)} + C_i^\top \gamma_k^{(b)} + \varepsilon_{ik}^{(b)}.
```

In the implementation, $\hat\beta_k^{(b)}$ is computed using precomputed OLS summary quantities (no repeated least-squares solves).

### 2.4 Empirical CRT p-values

Let $\hat\beta_k^{(\mathrm{obs})}$ be from observed $x$, and $\hat\beta_k^{(b)}$ from resample $b$. The two-sided CRT p-value is:

```math
p_k
=
\frac{
1 + \sum_{b=1}^B \mathbf{1}\!\left(\left|\hat\beta_k^{(b)}\right| \ge \left|\hat\beta_k^{(\mathrm{obs})}\right|\right)
}{
B+1
}.
```

### 2.5 Null Model Testing
Using the beta_hat matrix (program p × gene g), perform leave-one-out validation:

Compare beta_hat[i] vs beta_hat[b-i] for each gene
This generates g × (g-1) p-values for null calibration