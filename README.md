# fastsemR

<p align="center">
  <strong>📚 Escolha o idioma · 选择语言 · Choose language</strong><br>
  <a href="https://lf-araujo.github.io/fastsemR/articles/getting-started_pt.html">🇧🇷 Português</a> ·
  <a href="https://lf-araujo.github.io/fastsemR/articles/getting-started_zh.html">🇨🇳 中文</a> ·
  <a href="https://lf-araujo.github.io/fastsemR/articles/getting-started.html">🇬🇧 English</a>
</p>

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

**fastsemR** is an R interface to fastsem, a compiled Nim SEM engine that supports:

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

## Example: mediation

Decompose the effect of engine displacement (`disp`) on fuel economy
(`mpg`) through vehicle weight (`wt`).  Specify the model in umx (or
OpenMx) as usual, then let fastsem fit it: `run_fastsem()` returns the
same `MxModel`, so every umx/OpenMx tool keeps working — including
algebras such as the indirect effect `ab`, with delta-method SEs.

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
  mxAlgebra(a * b, name = "ab"),
  data = df, autoRun = FALSE
)

m_fit <- run_fastsem(m)
umxSummary(m_fit)        # Algebra 'ab' = -0.483 CI95[-0.808, -0.158]
mxSE(ab, m_fit)          # 0.166, identical to mxRun()
```

The same works from lavaan: build the model with `do.fit = FALSE` and
`run_fastsem()` returns a fitted `lavaan` object.

```r
library(lavaan)
fit <- run_fastsem(sem("wt ~ a*disp\nmpg ~ b*wt + c*disp\nab := a*b",
                       data = df, do.fit = FALSE))
summary(fit)
```

`run_fastsem()` checks its own work: OpenMx / lavaan re-evaluate the
model at the fastsem estimates, and a −2lnL that does not match the
engine's is reported.



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
for details.  