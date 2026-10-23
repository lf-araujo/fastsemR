# fastsemR 0.2.1

## Engine update (2026-10, release `0.2.1`)

* Multi-group configural fits fixed: the per-group DWLS warm start omitted
  the variance parameterisation, so under the now-default `vartype: direct`
  it returned log-scale variances that the FIML stage read as direct values —
  group 2 of a configural fit diverged (negative variances, `-2lnL` ~1e10,
  start values reported as estimates).  Both groups now match OpenMx exactly.
* The R process is no longer pinned to a single CPU core after a fit.
  Weave's thread pool pins the calling (R main) thread to CPU 0 on init and
  never restores the mask; threads created afterwards (OpenBLAS's lazily-
  initialised pool) inherited the one-core mask and spin-waited, making
  lavaan's `chol2inv`-heavy fits up to ~1000× slower in-process.  The engine
  now saves and restores the affinity mask around the pool.
* The multi-group ordinal twin conventions — shared-label equality across
  twins, `c()` per-group fixed covariances, labelled (free, pooled)
  thresholds — are pinned by engine tests verified against OpenMx joint
  multigroup FIML.
* The free-tier license gate is disabled while the payment placeholders are
  unconfigured (models of any size fit without a key).

## R package fixes

* `run_fastsem()` / `umx_to_lavaan()` on umx twin supermodels (`umxACE`,
  `umxCP`, ...) now fail with a clear, actionable message that spells out
  the hand-written multi-group `fastsem_fit()` route (previously an opaque
  `no applicable method for '@'` error).
* Estimates are now injected back into the OpenMx matrices for *labelled*
  parameters — previously only unlabelled ones were written, so umx start
  values could survive as "estimates" in `omxGetParameters()`.
* Multi-group `engine = "auto"` routes through the lavaan engine (the native
  multi-group RAM API is not implemented engine-side; the dispatch picks it
  up automatically if a future engine ships it).
* New worked example: multi-group ordinal twin ACE verified against OpenMx
  (`tests/testthat/test-ordinal-twin-ace.R`).

*Note:* the 0.2.0 package pointed at an engine release that was never
published; 0.2.1 supersedes it.  Run `fastsem_update()` after upgrading to
pull the new binary.

# fastsemR 0.1.4

## Engine update (2026-10)

* Refreshed the downloaded engine binary (new release `0.1.4`). Adds lavaan's
  `y ~*~ y` (delta / unit-liability) operator for ordinal models: the
  model-implied covariance is diagonal-standardised on the marked variables so
  their thresholds become free — the smooth equivalent of OpenMx's
  `mxConstraint(diag2vec(Vb) == 1)`, compatible with the trust-region optimiser
  and both variance parameterisations. This makes the strictly-contingent
  multi-group ordinal cascade (e.g. a three-stage CCC twin model) identifiable
  where fixed thresholds left the liability scale free. Standardisation runs on
  the F64 GPU kernel (−2logL matches CPU exactly); the analytic gradient carries
  the exact VJP (verified against finite differences).
* Multi-group ordinal fits now report standardised-estimate SEs (`stdSE`,
  previously returned as zeros) via the delta method over the identified
  subspace, so scale-invariant standardised estimates get finite, meaningful SEs
  even when raw-scale parameters sit on an unidentified scale ray; and the
  scale-free Newton-decrement convergence value is now reported on this path
  (previously `NaN`).
* Standard errors under the direct-symmetric variance parameterisation are no
  longer `NaN` on ill-conditioned-but-positive-definite information matrices
  (nearest-SPD Cholesky inverse). Estimates and SEs continue to match OpenMx.
  Run `fastsem_update()` after upgrading to pull the new binary.

# fastsemR 0.1.3

## Engine update (2026-09)

* Refreshed the downloaded engine binary (new release `0.1.3`). Fixes an
  intermittent crash in the definition-variable algebra-FIML path (a
  thread-safety bug in the parallel gradient) and improves the optimiser's
  convergence test: it now stops on a scale-free Newton decrement (√λ² ≈
  remaining move in SE units, invariant to how variables are scaled) rather than
  a unit-dependent gradient norm, so well-conditioned models converge in fewer
  iterations while weakly-identified directions are reported honestly. Fit
  objects now carry `newtonDecrement`, a convergence-quality value. Estimates and
  standard errors continue to match OpenMx. Run `fastsem_update()` after
  upgrading to pull the new binary.

# fastsemR 0.1.2

## Engine update (2026-09)

* Refreshed the downloaded engine binary (new release `0.1.2`) with the latest
  fastsem engine. Ordinal FIML now runs through a single estimation path shared
  by single- and multi-group models, so features compose: `vartype: direct`,
  definition variables, and grouping now all run on the GPU for ordinal models
  (multi-group ordinal previously forced definition variables onto the CPU and
  rejected `vartype: direct`). Estimates and standard errors continue to match
  OpenMx. Run `fastsem_update()` after upgrading the package to pull the new
  binary.

# fastsemR 0.1.1

## Engine update (2026-09)

* Refreshed the downloaded engine binary (new release `0.1.1`) with the latest
  fastsem engine. High-dimensional FIML now parallelises to models with up to
  64 observed variables (previously 32) for both the fit and the standard-error
  Hessian; the per-pattern gradient's memory is bounded, so large models that
  previously exhausted memory now complete; and a latent crash in the
  standard-error step for models above 32 observed variables is fixed. The
  optimiser dispatches a second-order trust-region step scaled to the parameter
  count. Estimates and standard errors continue to match OpenMx.
  Run `fastsem_update()` after upgrading the package to pull the new binary.
* The downloaded engine is now cached under a release-versioned filename, so a
  future package upgrade fetches the matching binary automatically on first
  use — no manual `fastsem_update()` needed. Any legacy version-less cache file
  is superseded.

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
