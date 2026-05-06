# fastsemR

<!-- badges: start -->
<!--
  Drop a hex sticker at man/figures/logo.png (square PNG, ~600×600) and
  uncomment the line below to get the BGmisc-style top-right logo.
-->
<!-- <a href="https://lf-araujo.github.io/fastsemR/"><img src="man/figures/logo.png" align="right" height="139" alt="fastsemR logo" /></a> -->

[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![R-CMD-check](https://github.com/lf-araujo/fastsemR/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/lf-araujo/fastsemR/actions/workflows/R-CMD-check.yml)
[![pkgdown](https://github.com/lf-araujo/fastsemR/actions/workflows/pkgdown.yml/badge.svg)](https://github.com/lf-araujo/fastsemR/actions/workflows/pkgdown.yml)
[![Codecov test coverage](https://codecov.io/gh/lf-araujo/fastsemR/branch/main/graph/badge.svg)](https://app.codecov.io/gh/lf-araujo/fastsemR)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

**fastsemR** is an R interface to [fastsem](https://github.com/lf-araujo/fastsem), a compiled Nim SEM engine that supports:

- **ML** (full-data maximum likelihood) with analytical gradient and optional GPU acceleration
- **FIML** (full-information ML for missing data and definition variables)
- **DWLS / WLSMV** (diagonally weighted least squares for ordinal indicators)
- **Multi-group models** with configural or constrained parameters
- **Equality constraints**, parameter bounds, and derived parameters (delta-method SEs)
- **Cluster-robust standard errors**
- **Standardised estimates** (StdAll) with SEs via the delta method

The package translates `umx` / OpenMx RAM model objects into fastsem lavaan syntax,
fits the model, and injects estimates back into the model object so that
`summary()`, `omxGetParameters()`, and `umxCompare()` continue to work.

---

## Installation

```r
# Install from GitHub (requires devtools)
devtools::install_github("lf-araujo/fastsemR")
```

The compiled shared library is downloaded automatically from GitHub Releases
the first time the package is loaded.  To force a re-download (e.g. after a
fastsem update):

```r
library(fastsemR)
fastsem_update()
```

If you have built the library yourself from the fastsem source tree
(`nimble buildRnim`), point the package to your local build:

```r
fastsem_load("~/fastsem/libfastsem_r.so")
```

---

## Example 1 — General statistics: mediation analysis

Decompose the effect of engine displacement (`disp`) on fuel economy
(`mpg`) through vehicle weight (`wt`), labelling paths `a`, `b`, and
`c_prime` so that the indirect effect `ab := a * b` is computed via the
delta method.

```r
library(fastsemR)
library(umx)

df <- as.data.frame(scale(mtcars[, c("mpg", "disp", "wt")]))

m <- umxRAM("Mediation",
  umxPath("disp", to = "wt",  labels = "a"),
  umxPath("wt",   to = "mpg", labels = "b"),
  umxPath("disp", to = "mpg", labels = "c_prime"),
  umxPath(var   = c("disp", "wt", "mpg")),
  umxPath(means = c("disp", "wt", "mpg")),
  data = df, autoRun = FALSE
)

m_fit <- run_fastsem(m)
summary(m_fit)
```

```
fastsem  --  FIML  [SE: OIM]
  N=32      free params=9    df=0

  #    Parameter                  Estimate          SE         z
  ─────────────────────────────────────────────────────────────
  1    a (disp->wt)               0.888       0.038     23.36
  2    b (wt->mpg)               -0.541       0.112     -4.83
  3    c_prime (disp->mpg)       -0.408       0.112     -3.64
  ...
```

The indirect effect and its SE are computed automatically when a
`derived_params` line is added to the syntax, e.g.:

```r
# Append to the umxPath-generated syntax before fitting
syntax <- paste(umx_to_lavaan(m), "ab := a * b", sep = "\n")
res    <- fastsem_fit(syntax, df)
cat("Indirect effect:", res$estimates[res$paramNames == "ab"],
    "  SE:", res$se[res$paramNames == "ab"], "\n")
```

---

## Example 2 — Behaviour genetics: univariate ACE twin model

The **ACE model** partitions phenotypic variance into additive genetic
(*A*), shared-environment (*C*), and unique-environment (*E*) components.
It is estimated as a two-group (MZ / DZ) SEM with **separate** A factors
per twin (`A1`, `A2`) and a per-group cross-twin genetic covariance
(`c(1, 0.5)` for MZ / DZ).

```r
library(fastsemR)

# ── Simulate twin data ─────────────────────────────────────────────────────────
set.seed(42)
sim_twins <- function(n, h2 = 0.5, c2 = 0.2) {
  e2 <- 1 - h2 - c2
  a  <- sqrt(h2);  cc <- sqrt(c2);  e <- sqrt(e2)
  A_mz <- rnorm(n);  C_mz <- rnorm(n)
  A_dz1 <- rnorm(n)
  A_dz2 <- 0.5 * A_dz1 + sqrt(1 - 0.5^2) * rnorm(n)
  C_dz  <- rnorm(n)
  mz <- data.frame(
    T1 = a*A_mz  + cc*C_mz  + e*rnorm(n),
    T2 = a*A_mz  + cc*C_mz  + e*rnorm(n), grp = 1L
  )
  dz <- data.frame(
    T1 = a*A_dz1 + cc*C_dz + e*rnorm(n),
    T2 = a*A_dz2 + cc*C_dz + e*rnorm(n), grp = 2L
  )
  rbind(mz, dz)
}

all_dat <- sim_twins(500)   # 500 pairs per zygosity

# ── ACE syntax ────────────────────────────────────────────────────────────────
# Separate A factors per twin with cross-twin genetic covariance c(1, 0.5).
# Setting Var(A1)=Var(A2)=1 means a²+c²+e² = 1 in both groups.

ace_syntax <- "
group: grp

A1 =~ a*T1
A2 =~ a*T2
C  =~ c*T1 + c*T2
E1 =~ e*T1
E2 =~ e*T2

A1 ~~ 1*A1;  A2 ~~ 1*A2
C  ~~ 1*C
E1 ~~ 1*E1;  E2 ~~ 1*E2

A1 ~~ c(1, 0.5)*A2   # genetic covariance: 1 for MZ, 0.5 for DZ

A1 ~~ 0*C;   A2 ~~ 0*C
A1 ~~ 0*E1;  A1 ~~ 0*E2
A2 ~~ 0*E1;  A2 ~~ 0*E2
C  ~~ 0*E1;  C  ~~ 0*E2
E1 ~~ 0*E2

T1 ~~ 0*T1;  T2 ~~ 0*T2
T1 ~ 1;      T2 ~ 1
"

res <- fastsem_fit(ace_syntax, all_dat)
print_fastsem(res)

# Labeled parameters appear in paramNames under their label
a_hat <- res$estimates[res$paramNames == "a"]
c_hat <- res$estimates[res$paramNames == "c"]
e_hat <- res$estimates[res$paramNames == "e"]
cat(sprintf("\nh²=%.3f  c²=%.3f  e²=%.3f  (total=%.3f)\n",
            a_hat^2, c_hat^2, e_hat^2, a_hat^2 + c_hat^2 + e_hat^2))
```

```
h²=0.349  c²=0.293  e²=0.339  (total=0.981)
```

Note: ML estimates can differ from method-of-moments (2*(MZ-DZ covariance) for
h²) because FIML accounts for the full variance structure — including that
Var(T1) ≠ Var(T2) in finite samples — while MOM uses only covariance
differences.  The chi-square test (df=4, p=0.42) confirms good model fit.

---

## Supported features at a glance

| Feature | fastsemR |
|---|:---:|
| ML (complete data) | ✓ |
| FIML (missing data) | ✓ |
| FIML with definition variables | ✓ |
| DWLS / WLSMV (ordinal data) | ✓ |
| Multi-group models | ✓ |
| Equality constraints | ✓ |
| Parameter bounds | ✓ |
| Derived parameters + delta-method SEs | ✓ |
| Cluster-robust SEs | ✓ |
| Standardised estimates (StdAll) | ✓ |
| GPU acceleration (OpenCL) | ✓ (engine only) |
| umx / OpenMx bridge (`run_fastsem`) | ✓ |

---

## Vignettes

- [Getting started with fastsemR](vignettes/getting-started.Rmd) —
  path models, CFA, FIML, labels, bounds, derived parameters.
- [Behaviour genetics](vignettes/behavior-genetics.Rmd) —
  univariate ACE twin model, bivariate Cholesky, RI-CLPM.
- [State-space models](vignettes/state-space-models.Rmd) —
  Kalman filtering, latent trajectories, longitudinal data.

See also the website articles for [comparison with other SEM
software](https://lf-araujo.github.io/fastsemR/articles/comparison.html),
[benchmarks](https://lf-araujo.github.io/fastsemR/articles/benchmarks.html),
and the [roadmap](https://lf-araujo.github.io/fastsemR/articles/roadmap.html).

---

## Citation

If you use **fastsemR** in your research, please cite it.  A BibTeX
entry is generated by:

``` r
citation("fastsemR")
```

A peer-reviewed methods paper is in preparation; this section will be
updated when the preprint is available.

---

## Contributing

Contributions are welcome.  Please open an issue at
<https://github.com/lf-araujo/fastsemR/issues> for bug reports,
feature requests, or questions about the R interface.  For issues
specific to the underlying Nim engine (numerical results, GPU
kernels, performance), file at
<https://github.com/lf-araujo/fastsem/issues>.

Pull requests targeting `main` should pass `R CMD check` with no
errors or warnings; the GitHub Actions workflow runs the same checks
on Linux, macOS, and Windows.

---

## License

fastsemR is released under the MIT License.  See
[LICENSE](https://github.com/lf-araujo/fastsemR/blob/main/LICENSE)
for details.  The compiled **fastsem** engine is distributed under
its own license; see the
[fastsem repository](https://github.com/lf-araujo/fastsem) for
terms.
