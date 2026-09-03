# fastsemR (development version)

## Engine update (2026-09)

* Refreshed the downloaded engine binary (release `0.1`) with the fixed fastsem
  engine: **multi-group ordinal FIML** now works (twin ACE / CCC models were
  previously fit as continuous-on-codes, collapsing the genetic contrast into
  C); **Newton-TR** is the primary optimiser for ordinal and SSM/Kalman fits;
  and the Linux R library is now **GPU-accelerated by default** (OpenCL, with
  graceful CPU fallback). Estimates and standard errors match OpenMx.
  Run `fastsem_update()` to pull the new binary.

## Initial development

* R bindings to the fastsem Nim SEM engine: `fastsem_fit()`,
  `fastsem_sem_vars()`, `print_fastsem()`.
* umx / OpenMx bridge: `run_fastsem()`, `run_fastsem_ram()`,
  `umx_to_lavaan()`, `umx_to_fastsem()`.
* Library management: `fastsem_install()`, `fastsem_update()`,
  `fastsem_load()` — the compiled shared library is downloaded
  automatically from GitHub Releases on first use.
* Vignettes: getting started, behaviour genetics (univariate ACE,
  bivariate Cholesky, RI-CLPM), state-space models.
* pkgdown site at <https://lf-araujo.github.io/fastsemR/>.
