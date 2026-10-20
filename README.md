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

## Roadmap / known limitations

- [x] **`fastsem_fit()` drops `:=` derived parameters.** Found by the
  spec-equivalence suite (`2026-fastsem/2026-spec-equivalence.org`, M3 —
  mediation `ab := a * b`): the engine computes them (the CLI prints
  `Derived Parameters (:=)` with delta-method SEs — `ab = −0.4147, SE
  0.1073`, matching lavaan exactly), but the R result list had no
  `derived*` field, so `ab` silently vanished from the result.
  *Fixed 2026-10-19:* the result list now carries `derivedNames`,
  `derivedEstimates`, `derivedSEs` (kept separate from
  `paramNames`/`estimates` so `nFreeParams` stays honest).
- [x] **Labeled parameters are reported by label only.** Same suite: a
  path written `wt ~ a*hp` came back as `paramNames == "a"` with no
  record of the `hp->wt` mapping. *Fixed 2026-10-19:* the result list
  now carries `paramPaths`, position-aligned with `paramNames`, holding
  the label-independent structural name.

---

- [ ] **Multi-group *ordinal twin* models don't fit through either R route.**
  Found while trying to reproduce an OpenMx three-stage CCC (ordinal
  puff→use→rep, MZ/DZ) to get bootstrap CIs. Two failures:
  - `run_fastsem()` on an ordinal `umx` model (e.g. `umxACE(..., ordinal)`)
    **errors** with `no applicable method for '@' applied to an object of
    class "NULL"`. The OpenMx→fastsem bridge (`umx_to_lavaan` / RAM path)
    doesn't carry ordinal **thresholds**, so the ordinal expectation
    resolves to `NULL`. It should either translate thresholds or fail with a
    clear, actionable message (currently the RAM helper only notes "no
    ordinal thresholds" in a comment).
  - Hand-written multi-group `fastsem_fit()` syntax mis-handles **cross-twin
    loading equality**: giving both twins' loadings the same label
    (`A1 =~ a*puff_T1`, `A2 =~ a*puff_T2`) does **not** equate them — a
    separate `A2->puff_T2` parameter is created — so the ACE structure
    collapses (std → 0/0/1, `-2lnL` blows up, optimiser stops in 1 step).
    Also `ordered:` here **frees** the thresholds ("N free thresholds added
    to optimisation") rather than fixing them at empirical probit as
    `benchmark/compare_accuracy_gl3.R` assumes.
  - *Needs:* a documented multi-group convention (how to equate a loading
    across twins; how to fix a cross-group covariance to `c(1, 0.5)`; whether
    `ordered:` fixes or frees thresholds and how to force fixed) **plus a
    worked twin ACE / CCC vignette or test**. The `behaviour-genetics`
    vignette currently covers only *continuous* ACE; the only ordinal example
    (`run_ordinal_riclpm`) is single-group. NEWS.md claims twin ACE/CCC
    ordinal "matches OpenMx" — that path needs an R-level example.
  - *Repro:* `fastsem_load()` OK; GPU OpenCL FIML runs; OpenMx `umxACE` on the
    same binary `puff_ecig` data fits fine (−2lnL 1594.87, std A/C/E
    .37/.44/.19) — so the blocker is purely the fastsem multi-group ordinal
    twin interface, not the data or the engine's single-group ordinal FIML.

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
