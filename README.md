# simdatacopulas

`simdatacopulas` is an extension of `simdata`. It wraps different copula types (dependence structures) into `simdata::simdesign()` templates, so you can generate multivariate data in a unified way with "arbitrary marginals + controllable dependence".

The core idea of this package is a two-step workflow:

1. Use a copula to generate dependent \(U\in(0,1)^d\) (joint dependence structure)
2. Apply a marginal quantile function `dist[[j]](u)` to each column \(U_j\), mapping \(U\) to the target variable scales (marginal distributions)

In short, the copula controls correlation/tail dependence/high-dimensional structure, while the marginal quantile functions control each univariate distribution shape (continuous/discrete/truncated, etc.).

---

## Installation and dependencies

- Runtime dependencies: `simdata`, `copula`
- Optional dependency: `VineCopula` (only required for the vine-copula module)

---

## Quick start (generic workflow)

```r
library(simdata)
library(simdatacopulas)

dist <- list(
  function(u) qnorm(u, mean = 0, sd = 1),
  function(u) qexp(u, rate = 0.2),
  function(u) qbeta(u, 2, 8)
)

design <- simdesign_gaussian_copula(
  dist = dist,
  dim = 3,
  rho = 0.4,
  names_final = c("x1", "t", "b")
)

x <- simulate_data(design, n_obs = 1000, seed = 1)
head(x)
```

---

## Model overview (what each copula family is for)

This package currently provides six joint-generation modules:

- Elliptical copulas: `simdesign_elliptical_copula()` / `simdesign_gaussian_copula()` / `simdesign_elliptical_copula_from_data()`
- Archimedean copulas: `simdesign_archimedean_copula()` plus convenience constructors for Clayton/Gumbel/Frank/Joe and `*_from_data()`
- Vine copulas (pair-copula construction): `simdesign_vine_copula()` / `simdesign_vine_copula_from_data()`
- Empirical copulas (resample the joint rank structure): `simdesign_empirical_copula()` / `simdesign_empirical_copula_from_data()`
- Factor copulas (low-rank dependence, suitable for high dimensions): `simdesign_factor_copula()` / `simdesign_factor_copula_from_data()`
- Latent threshold Gaussian model (joint generation for discrete/ordinal endpoints): `simdesign_latent_threshold_gaussian()` / `simdesign_latent_threshold_gaussian_from_data()`

Below we describe each model's purpose, when to use it, and example usage.

---

## 1) Elliptical copulas (Gaussian / t)

### Purpose

- Use a correlation matrix (or an exchangeable correlation coefficient) to define the dependence structure; suitable for simulating multivariate endpoints with **linear correlation / symmetric dependence**.
- A `Gaussian copula` is appropriate when dependence is largely captured by correlations.
- A `t copula` adds degrees of freedom (`df`) on top of the Gaussian case and can represent stronger **tail dependence** (extreme events co-occur more often).

### When to use

- You have a target correlation matrix (or want an exchangeable structure) and want to couple multiple endpoints (e.g., continuous PFS/OS scores, biomarker measures, etc.).
- Dimension is moderate and you want a robust baseline dependence structure.

### Key interfaces

- `simdesign_elliptical_copula(copula, dist, ...)`: pass a `copula::normalCopula()` or `copula::tCopula()` object, etc.
- `simdesign_gaussian_copula(dist, structure = "ex"/"un", ...)`: convenience entry point for the Gaussian copula
- `simdesign_elliptical_copula_from_data(data, vars, family = "gaussian"/"t", ...)`: fit dependence from data; marginals default to empirical quantiles

### Example: fit a Gaussian copula from data

```r
library(simdata)
library(simdatacopulas)

set.seed(1)
pilot <- data.frame(
  pfs = rexp(500, 0.15),
  os  = rexp(500, 0.10),
  biomarker = rnorm(500, 0, 1)
)

design_fit <- simdesign_elliptical_copula_from_data(
  data = pilot,
  vars = c("pfs", "os", "biomarker"),
  family = "gaussian",
  structure = "un"
)

sim <- simulate_data(design_fit, n_obs = 1000, seed = 2)
head(sim)
```

---

## 2) Archimedean copulas (Clayton / Gumbel / Frank / Joe)

### Purpose

Archimedean copulas parameterize dependence strength using a single parameter (often denoted `theta`) and are commonly used to model **asymmetric tail dependence**:

- Clayton: stronger **lower-tail dependence** (low values / worse outcomes co-occur more often)
- Gumbel: stronger **upper-tail dependence** (high values / better outcomes or extreme highs co-occur more often)
- Frank: symmetric dependence, typically **no tail dependence** (dependence is more pronounced in the middle)
- Joe: upper-tail dependence (often "sharper" than Gumbel)

### When to use

- You care about co-movement in extreme/tail events, e.g., synchrony between severe AEs and extreme efficacy responses.
- You want a lightweight dependence model with few parameters while retaining different tail behaviors.

### Key interfaces

- `simdesign_archimedean_copula(copula, dist, ...)`
- Convenience constructors: `simdesign_clayton_copula()`, `simdesign_gumbel_copula()`, `simdesign_frank_copula()`, `simdesign_joe_copula()`
- `simdesign_archimedean_copula_from_data(data, vars, family = ..., ...)`: fit `theta` from data (with empirical margins)

### Example: simulate "upper-tail dependence" with a Gumbel copula

```r
library(simdata)
library(simdatacopulas)

dist <- list(
  function(u) qlnorm(u, meanlog = 0, sdlog = 0.6),
  function(u) qlnorm(u, meanlog = 0.2, sdlog = 0.6)
)

design <- simdesign_gumbel_copula(
  dist = dist,
  theta = 2,
  dim = 2,
  names_final = c("endpoint1", "endpoint2")
)

x <- simulate_data(design, n_obs = 2000, seed = 1)
head(x)
```

---

## 3) Vine copula（C-vine / D-vine / R-vine）

### Purpose

Vine copulas build high-dimensional dependence by stitching together a tree of pair-copulas (pair-copula construction):

- More flexible than a single "global" copula in high dimensions: different variable pairs can use different families and tail behaviors
- Can express complex conditional dependence: e.g., dependence between \(X_1\) and \(X_3\) may be largely mediated through \(X_2\)

### When to use

- Dimension is high (e.g., 6+ endpoints/covariates) and you do not want to be constrained by a single correlation matrix / single parameter.
- You want data-driven selection of pair-copula families (AIC/BIC, etc.).

### Dependency note

The vine module depends on the `VineCopula` package; if it is not installed, an error will prompt you to install it.

### Key interfaces

- `simdesign_vine_copula(vine, dist, ...)`: pass a `VineCopula` vine object (e.g., `RVineMatrix`)
- `simdesign_vine_copula_from_data(data, vars, vine_type = "rvine"/"cvine"/"dvine", ...)`: fit a vine structure and pair-copulas from data

### Example: fit an R-vine from data (automatic structure selection)

```r
library(simdata)
library(simdatacopulas)

set.seed(1)
pilot <- data.frame(
  x1 = rnorm(600),
  x2 = rnorm(600),
  x3 = rnorm(600),
  x4 = rnorm(600)
)

design_fit <- simdesign_vine_copula_from_data(
  data = pilot,
  vars = c("x1", "x2", "x3", "x4"),
  vine_type = "rvine",
  selectioncrit = "AIC",
  method = "itau"
)

sim <- simulate_data(design_fit, n_obs = 1000, seed = 2)
head(sim)
```

---

## 4) Empirical copula (rank resampling)

### Purpose

An empirical copula does not assume a parametric copula family explicitly. Instead, it:

- **Resamples** the joint rank structure from existing `U` (pseudo-observations)
- Combines it with your chosen marginal quantile functions to generate new samples

Its core value is that, when sample size is sufficient and you believe the observed joint rank structure is representative, you can very directly "replicate" the dependence structure.

### When to use

- You have reliable historical/real data and want to preserve its joint dependence structure (especially complex dependence that is hard to parameterize).
- You want to avoid modeling bias introduced by choosing a copula family.

### Notes

- Extrapolation is limited: it is good at reproducing observed dependence, but not at inferring unobserved tail structure.
- Optional `jitter` can help mitigate repeated points due to many ties or discretization.

### Key interfaces

- `simdesign_empirical_copula(u_data, dist, replace = TRUE, jitter = 0, ...)`
- `simdesign_empirical_copula_from_data(data, vars, ...)`: automatically compute `u_data` via `pobs()`

### Example: build an empirical-copula design directly from data

```r
library(simdata)
library(simdatacopulas)

set.seed(1)
pilot <- data.frame(
  a = rexp(800, 0.3),
  b = rlnorm(800, 0, 0.5),
  c = rnorm(800)
)

design_fit <- simdesign_empirical_copula_from_data(
  data = pilot,
  vars = c("a", "b", "c"),
  jitter = 0
)

sim <- simulate_data(design_fit, n_obs = 1000, seed = 2)
head(sim)
```

---

## 5) Factor copula (Factor Gaussian copula)

### Purpose

Factor copulas generate high-dimensional dependence using a structure of "a few common factors + idiosyncratic noise":

\[
Z = F \Lambda^\top + E,\quad U = \Phi(\tilde Z)
\]

- \(\Lambda\) (`loadings`) controls each variable's sensitivity to the common factors
- `uniq_var` controls each variable's idiosyncratic variance

This effectively constrains the covariance/correlation structure to "low-rank + diagonal", making it more stable and parsimonious in high dimensions.

### When to use

- Dimension is high (e.g., multiple endpoints + covariates) and you want a small-parameter representation of shared drivers of correlation (e.g., disease severity, exposure, treatment sensitivity).
- You need a correlation structure that is more stable and interpretable than a full correlation matrix.

### Key interfaces

- `simdesign_factor_copula(loadings, uniq_var, dist, ...)`: provide the factor structure and marginals directly
- `simdesign_factor_copula_from_data(data, vars, n_factors = 1, ...)`: convert data to pseudo-observations and map to Gaussian space, then fit a factor model

### Example: fit a 1-factor copula from data

```r
library(simdata)
library(simdatacopulas)

set.seed(1)
pilot <- data.frame(
  x1 = rexp(2000, 0.2),
  x2 = rlnorm(2000, 0, 0.6),
  x3 = rnorm(2000)
)

design_fit <- simdesign_factor_copula_from_data(
  data = pilot,
  vars = c("x1", "x2", "x3"),
  n_factors = 1
)

sim <- simulate_data(design_fit, n_obs = 1000, seed = 2)
head(sim)
```

---

## 6) Latent threshold Gaussian model

### Purpose

This module jointly generates **discrete/ordinal endpoints** (e.g., AE grades, binary ORR, toxicity grades):

- First generate continuous latent variables \(Z\) using a factor structure
- Each discrete variable is formed by thresholding \(Z_j\) using a threshold vector `thresholds[[j]]` into levels 1/2/… etc.
- Optionally provide `levels` to set ordered factor labels for each column

The value of this model is that it expresses correlation among discrete endpoints via a Gaussian factor (or multivariate Gaussian) structure while preserving each discrete variable's marginal distribution (determined by the thresholds).

### When to use

- Endpoints include ordinal/binary variables, and you want them to share latent correlation sources (common factors) with other endpoints.
- You want explicit control of category proportions (thresholds) and the correlation structure across levels.

### Key interfaces

- `simdesign_latent_threshold_gaussian(loadings, uniq_var, thresholds, levels = NULL, ...)`
- `simdesign_latent_threshold_gaussian_from_data(data, vars, n_factors = 1, ...)`:
  - Automatically estimate thresholds from each column (map cumulative category probabilities to `qnorm` cut points)
  - Convert data to pseudo-observations, map to Gaussian space, then fit a factor structure

### Example: manually specify thresholds to generate a binary + three-level joint distribution

```r
library(simdata)
library(simdatacopulas)

loadings <- matrix(c(0.9, 0.6, 0.7), nrow = 3, ncol = 1)
uniq <- c(0.4, 0.7, 0.6)

thresholds <- list(
  qnorm(0.30),
  qnorm(0.55),
  qnorm(c(0.20, 0.70))
)

design <- simdesign_latent_threshold_gaussian(
  loadings = loadings,
  uniq_var = uniq,
  thresholds = thresholds,
  names_final = c("y1", "y2", "y3")
)

sim <- simulate_data(design, n_obs = 2000, seed = 1)
head(sim)
```

---

## Integration with TrialSimulator (recommended)

`TrialSimulator` expects a data generator of the form "a `generator(n)` that returns a data.frame of endpoint columns". This package returns `simdata` `simdesign` objects, so a typical bridge looks like:

1. Build a `design` using `simdatacopulas::*`
2. Generate a batch of jointly dependent base variables using `simdata::simulate_data(design, n_obs = n, seed = ...)`
3. Inside TrialSimulator's endpoint generator, map those base variables into the endpoint format TrialSimulator needs (e.g., `*_event`, censoring rules, etc.)

The benefit is that dependence/marginal modeling is decoupled from the trial process (enrollment/censoring/analysis), making it easier to reuse and iterate.
