# fastsemR 0.3.0

## One entry point for umx / OpenMx and lavaan

* `run_fastsem()` now takes an `MxModel` (umxRAM, `mxModel(type = "RAM")`,
  or a multi-group container of RAM models) **or** an unfitted lavaan object
  (`do.fit = FALSE`), and returns the same kind of object.  The model is read
  into an internal parameter table (one row per RAM cell and group, like
  lavaan's `parTable()`) and sent to the engine as typed columns; the engine
  lowers it with internal names and round-trip numbers and returns estimates
  and SEs per table row.  The column layout is checked against the engine's
  at fit time.
* OpenMx / lavaan then evaluate their own model at the fastsem estimates
  (`mxRun(useOptimizer = FALSE)`; lavaan with `optim.method = "none"`), so fit
  statistics, algebras (`mxSE()`, `umxSummary()` CIs), standardised solutions
  and lavaan's robust SEs/tests are computed natively — no more hand-filled
  `@output`.  For OpenMx, fastsem's SEs and parameter covariance are attached
  (`vcov()` works).  The fastsem result list is in `model@output$fastsem` /
  `fit@external$fastsem`.
* Self-check: the −2lnL that OpenMx/lavaan compute at the fastsem estimates
  must equal the engine's; `check = "warn" | "error" | "off"`
  (`options(fastsemR.check = )`).
* Bounds (OpenMx `lbound`/`ubound`, lavaan `lower`/`upper` and `a > 0`
  lines) are enforced by the engine as hard bounds: a Heywood case stops at
  `lbound = 0` exactly as in OpenMx.  Inequalities between two parameters
  cannot be sent and are only checked after the fit.
* Newly supported through the bridge: labelled and shared means; multi-group
  models with shared (same label) and group-specific (unlabelled) parameters,
  as OpenMx and lavaan define them; `group.equal`; lavaan `:=`, MLR and other
  ML-based estimators, `fixed.x` covariates; `mxConstraint` equalities between
  labels; plain threshold matrices and umx's deviation thresholds
  (`umxThresholdMatrix`, incl. Mehta's default for 3+ categories); ordinal
  models mixed with continuous indicators; means fixed at non-zero values.
  Ordinal variables must be identified (theta, Mehta, binary conventions);
  the delta parameterisation (a total-variance mxConstraint) is not
  translated.
  Unsupported constructs fail with a clear error instead of being dropped:
  umx twin supermodels built from algebras (`umxACE`, `umxCP`, ...), cells
  bound to algebras, definition variables on variances/means, lavaan ordered
  indicators, and umxRAM's deviation-based threshold algebras.
* `umx_to_lavaan()`, `umx_to_fastsem()` and `run_fastsem_ram()` are
  deprecated; `run_fastsem(engine = )` is ignored.

## Engine (release `0.3.0`)

The bridge requires these engine fixes; `run_fastsem()` refuses to run on an
older engine (`fastsem_update()`).

* Labelled, shared-label and start-valued intercepts (`x ~ m*1`) were
  silently dropped — the model became covariance-only with a meaningless
  −2lnL.  Fixed non-zero intercepts (`x ~ 5*1`, `c(NA, 5)*1`) were dropped
  too.  Both are now fitted (a fixed non-zero mean needs at least one free
  mean in the model).
* Fixed values are kept in double precision (they were stored as float32).
* Multi-group fits now follow lavaan/OpenMx semantics: a labelled parameter is
  shared across groups, an unlabelled one is group-specific, `c(NA, ...)` is a
  separate parameter per group, `c()` may carry different labels per group,
  and a parameter fixed by `c()` in some groups no longer receives gradient
  from them.  Previously `c(NA, NA)` equated across groups whenever any `c()`
  was present, and labels were ignored when none was.  The joint fit's SEs
  now come from the joint Hessian.  Verified against OpenMx (−2lnL, estimates
  and SEs).
* Parameter-table input (`fitSemTableR`, `tableSchemaR`, `tableSyntaxR`;
  `partable.nim`): the row type is the single definition of the columns
  (schema and decoder are derived from its fields).
* `:=` derived parameters are returned for multi-group models (they were
  dropped), with delta-method SEs from the joint covariance.
* Engine programs no longer read stdin at start-up: the CLI's pre-pass
  (skip the OpenCL compile for `expected:` models) moved into `main`; the
  test binaries hung when stdin was an open pipe.
* **Ordinal FIML**, the engine's focus:
  - Models with continuous *and* ordinal observed variables are fitted by the
    joint likelihood, as in OpenMx: per row the continuous values by their
    normal density and the ordinal ones by the rectangle probability of their
    latent responses conditional on the continuous values.  Previously the
    continuous variables were silently integrated out, so their data did not
    enter the fit.  Runs on the F64 GPU kernel (CPU otherwise); integrating
    only the ordinal dimensions also makes these rows cheaper.
  - Mehta's parameterisation (two thresholds fixed, latent-response means and
    residual variances free) is fitted by FIML; it was routed to DWLS.  An
    ordinal variable is FIML-identified when two of {mean, residual variance,
    thresholds} are fixed; an under-identified one is now an **error** instead
    of a silent switch to DWLS, which is available on request with
    `estimator: dwls`.
  - A free threshold above fixed ones is kept above them (and started one
    fixed-threshold spacing higher), so a category cannot collapse.
  - The covariance of (structural parameters, thresholds) is returned
    (`vcovAug`), e.g. for SEs of umx threshold deviations.
* `likelihood: normal` directive: complete-data ML on the N-denominator
  likelihood (OpenMx and lavaan's default); the default remains Wishart.
* Bounds (`a > lo`, `a < hi`) are hard bounds, enforced by an active-set loop
  around every optimiser; an estimate sits exactly on a bound when it binds.
  They were a log-barrier penalty (weight 1000) that also moved solutions
  whose bound did not bind, and under `vartype: direct` (the default) a
  variance bound was applied to exp(σ²), so `v > 0` never bound and pushed
  the variance up.
* **`expected:` (algebra / continuous-time) models fit on the GPU** by
  default: the engine compiles each model's algebra into its own OpenCL
  kernels (objective and analytical gradient), with a note that this can be
  slower on some hardware.  A bivariate CT-CLPM with per-person intervals:
  N = 10⁴ 7.6 s → 1.3 s, N = 10⁵ 87 s → 11 s, identical results.
  `algebragpu: off` or `Sys.setenv(FASTSEM_ALGEBRA_GPU = "off")` keeps the
  CPU.  The CPU algebra evaluator also gained `/`, `%x%`, `solve()`, `tr()`
  and `diag()` for models with definition variables or missing data.
* **Ordinal models above 8 observed variables follow the GPU type**: on a GPU
  that shares memory with the host (APUs, integrated graphics, Apple
  Silicon) they run on the CPU, because long GPU jobs there can freeze the
  display; on a dedicated GPU they run on the GPU, with a note.
  `options(fastsemR.ordinal_gpu = TRUE / FALSE)` overrides either way; the
  session question is gone.
* Complete-data ML without means follows lavaan's `meanstructure = FALSE`
  throughout: −2logL, AIC and BIC are evaluated at the sample means (they
  were evaluated at μ = 0, adding N·ȳᵀΣ⁻¹ȳ on uncentred data), and
  `expected:` models no longer fit the uncentred second moment.  Models from
  `run_fastsem()` always carry means and are unaffected.
* The CPU and GPU paths share one source for their numerics (Cholesky and
  solves, the Woodbury update for definition variables, the matrix
  exponential, the Genz rectangle probability); the GPU Woodbury kernels no
  longer add a 1e-8 diagonal jitter the CPU did not.

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
