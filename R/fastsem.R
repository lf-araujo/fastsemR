## fastsem.R — core SEM fitting, result printing, deprecated bridge wrappers

# ── Run-control helpers ──────────────────────────────────────────────────────
#
# `control` is a named list with optional fields (absent or 0 = engine default):
#   maxIter    integer  optimiser iteration cap
#   tol        numeric  gradient tolerance
#   ftol       numeric  objective-change tolerance
#   tryHard    integer  extra rounds while not converged; each restarts from the
#                       last solution with 3× maxIter and tol/10
#   quadrature integer  ordinal FIML Gauss-Legendre points per dimension (1–15)
# Each is applied as the matching model-syntax directive (`maxiter:`, `tol:`,
# `ftol:`, `tryhard:`, `quadrature:`), appended after the user's syntax so the
# control list takes precedence over a directive written in the syntax.
.fastsem_control_keys <- c(maxIter = "maxiter", tol = "tol", ftol = "ftol",
                           tryHard = "tryhard", quadrature = "quadrature")

.fastsem_control_directives <- function(control) {
  if (is.null(control) || length(control) == 0L) return(character(0))
  if (!is.list(control) || is.null(names(control)))
    stop("fastsem control must be a named list")
  unknown <- setdiff(names(control), names(.fastsem_control_keys))
  if (length(unknown))
    stop("unknown fastsem control field(s): ", paste(unknown, collapse = ", "),
         "; recognised: ", paste(names(.fastsem_control_keys), collapse = ", "))
  out <- character(0)
  for (nm in names(control)) {
    v <- control[[nm]]
    if (is.null(v) || isTRUE(as.numeric(v[1]) == 0)) next
    val <- if (nm %in% c("tol", "ftol")) format(as.numeric(v[1]), digits = 15) else as.integer(v[1])
    out <- c(out, paste0(.fastsem_control_keys[[nm]], ": ", val))
  }
  out
}

# Legacy numeric control vector still passed to the native entry points (the
# current engine ignores it; run control travels as syntax directives).
.fastsem_control_vec <- function(control) {
  if (is.null(control) || length(control) == 0L) return(numeric(0))
  if (!is.list(control))
    stop("fastsem control must be a named list")
  pick <- function(nm, default) {
    v <- control[[nm]]
    if (is.null(v)) return(default)
    as.numeric(v[1])
  }
  c(pick("maxIter", 0),
    pick("tol",     0),
    pick("ftol",    0),
    pick("tryHard", 0))
}

# ── Low-level fit wrapper ──────────────────────────────────────────────────────

#' Fit a structural equation model from lavaan syntax.
#'
#' The workhorse function of fastsemR.  Passes lavaan model syntax and a
#' data matrix to the fastsem compiled engine and returns a named list of fit
#' statistics and parameter estimates.
#'
#' @section Estimators:
#' The estimator is chosen automatically based on the data and model:
#' * **ML** — complete data, no definition variables, no ordinal indicators.
#' * **FIML** — missing data (`NA`) present, or definition variables used.
#' * **DWLS / WLSMV** — ordinal indicators declared with `ordered:` or
#'   threshold syntax (`y | t1 + t2`).
#'
#' @section Extended lavaan syntax:
#' Beyond standard lavaan, fastsem supports:
#' * **Parameter labels** — `b1*var`: assign the label `b1` to the path.
#' * **Equality constraints** — two params with the same label are equated,
#'   or `b1 == b2` on a standalone line.
#' * **Bounds** — `b1 > 0.5` / `b1 < 2.0`.
#' * **Starting values** — `start(0.8)*var`.
#' * **Derived parameters** — `indirect := a * b` (delta-method SE).
#' * **Definition variables** — `data.colname * var` on the RHS of `=~`.
#' * **Multi-group** — `group: colname` directive plus `c(v1, v2)*var`
#'   per-group annotations.  Each element is a number (fixed in that group),
#'   `NA` (free), `NaN`/`_` (variable blanked in that group) or a label: in
#'   `A1 ~~ c(vA, 0)*A2` the covariance equals the parameter labelled `vA` in
#'   group 1 and is fixed at 0 in group 2.  A parameter may carry one label
#'   across groups (mixed with fixed values, not with `NA`).
#' * **Weighted observations** — `weight: colname`.
#' * **Cluster-robust SEs** — `cluster: colname`.
#' * **Ordinal thresholds** — `y | t1 + t2`; each term may be `lab*t1`
#'   (equality label), `0.5*t1` (fixed) or `start(0.5)*t1`, and modifiers
#'   combine: `start(0.5)*lab*t1`.
#'
#' @section Model directives:
#' Lines of the form `key: value` set model-wide options:
#' * `group: col`, `weight: col`, `cluster: col` — see above.
#' * `vartype: log | direct` — variance parameterisation.  `log` (default)
#'   keeps variances positive; `direct` estimates them unconstrained, so they
#'   can go negative (as OpenMx's direct-symmetric parameterisation).
#' * `warmstart: on | off` — the polychoric/DWLS warm start (default on); turn
#'   off for sparse or contingent ordinal cells.
#' * `y ~*~ y` — unit-liability (delta) scaling of observed `y`: the implied
#'   variance of `y` is standardised to 1, so its thresholds are free (OpenMx's
#'   `diag(Vb) = 1` constraint).
#' * `quadrature: N` — Gauss–Legendre points per integrated dimension for
#'   ordinal FIML, 1–15 (default 3, or the `FASTSEM_ORD_NQ` environment
#'   variable).  Cost per row grows as \eqn{N^{k-1}}; 3 is fast but biased
#'   (about 6 −2lnL units on a 6-variable twin model), 7–9 matches an
#'   accurately integrated OpenMx fit.
#' * `maxiter: N` — optimiser iteration cap.
#' * `tol: x` — gradient tolerance.
#' * `ftol: x` — objective-change tolerance (L-BFGS: absolute
#'   \eqn{|\Delta f|}; trust-region optimisers: relative
#'   \eqn{|\Delta f|/(1+|f|)}).
#' * `tryhard: N` — up to `N` extra rounds while the gradient is above
#'   tolerance, each restarting from the last solution with 3× `maxiter` and
#'   `tol`/10 (cf. OpenMx's `mxTryHard`).
#'
#' `maxiter`, `tol`, `ftol` and `tryhard` apply to every optimiser used by the
#' fit; when omitted, each fitting regime keeps its own default.  The same
#' settings can be passed as `control`.
#'
#' @param syntax Character string of lavaan model syntax.
#' @param data Data frame or numeric matrix.  Non-numeric columns are coerced
#'   to `NA`.  Row order is preserved.
#' @param col_names Character vector of column names.  Defaults to
#'   `colnames(data)`.  Must include every observed variable named in
#'   `syntax`.
#' @param control Optional named list of run-control settings: `maxIter`,
#'   `tol`, `ftol`, `tryHard`, `quadrature`.  Each is applied as the matching
#'   model directive (`maxiter:`, `tol:`, `ftol:`, `tryhard:`, `quadrature:`;
#'   see *Model directives*) and takes precedence over a directive in
#'   `syntax`.  Absent or `0` keeps the default; unknown names are an error.
#'
#' @return A named list with the following elements (errors raised inside
#'   the engine, e.g. a syntax error, are reported as R errors):
#'   \describe{
#'     \item{`estimator`}{Character. `"ML"`, `"FIML"`, or `"DWLS"`.}
#'     \item{`seType`}{Character. SE method, e.g. `"OIM"` or `"cluster-robust"`.}
#'     \item{`nobs`}{Integer. Number of observations (or effective N for weighted
#'       models).}
#'     \item{`nFreeParams`}{Integer. Number of free parameters.}
#'     \item{`df`}{Integer. Model degrees of freedom.}
#'     \item{`chi2`}{Numeric. Chi-square test statistic.}
#'     \item{`pvalue`}{Numeric. Chi-square p-value.}
#'     \item{`neg2logL`}{Numeric. \eqn{-2 \ln L} value.}
#'     \item{`aic`}{Numeric. Akaike information criterion.}
#'     \item{`bic`}{Numeric. Bayesian information criterion.}
#'     \item{`paramNames`}{Character vector. Structural names of all free
#'       parameters (e.g. `"f->y1"`, `"Var(y1)"`, `"y1~1"`).}
#'     \item{`estimates`}{Numeric vector. Parameter estimates in the same
#'       order as `paramNames`.  Variance parameters are on the original
#'       (not log) scale.}
#'     \item{`se`}{Numeric vector. Standard errors.}
#'     \item{`z`}{Numeric vector. z-statistics (`estimates / se`).}
#'     \item{`stdEst`}{Numeric vector. StdAll standardised estimates.}
#'     \item{`stdSE`}{Numeric vector. Delta-method SEs for `stdEst`.}
#'     \item{`groupEstimates`}{List of numeric vectors (multi-group models).
#'       `groupEstimates[[g]]` contains the per-group estimates for group
#'       `g` in the shared parameter space.}
#'     \item{`groupSEs`}{List of numeric vectors (multi-group models).}
#'     \item{`newtonDecrement`}{Numeric. Scale-free convergence measure
#'       \eqn{g^T H^+ g} (≈ remaining move in SE units).}
#'     \item{`thresholdNames`, `thresholds`, `thresholdSE`}{Ordinal thresholds
#'       (`"y|t1"`, …), their estimates and SEs (`NaN` for fixed thresholds).
#'       With free thresholds, the structural SEs include the threshold
#'       uncertainty.}
#'     \item{`freeNames`, `vcov`}{Names of the free structural parameters and
#'       their covariance matrix (raw scale; pseudo-inverse over the
#'       identified subspace), for delta-method CIs of functions of them,
#'       e.g. standardised variance components.}
#'   }
#'
#' @seealso [print_fastsem()] to display the result;
#'   [run_fastsem()] to fit OpenMx/umx or lavaan models.
#'
#' @export
#' @examples
#' \dontrun{
#' # Single-factor CFA on simulated data
#' set.seed(1)
#' n  <- 200
#' f  <- rnorm(n)
#' df <- data.frame(y1 = 0.8*f + rnorm(n, sd=0.6),
#'                  y2 = 0.9*f + rnorm(n, sd=0.5),
#'                  y3 = 0.7*f + rnorm(n, sd=0.7))
#'
#' syntax <- "
#'   f =~ 1*y1 + y2 + y3
#'   f ~~ f
#'   y1 ~~ y1;  y2 ~~ y2;  y3 ~~ y3
#' "
#' res <- fastsem_fit(syntax, df)
#' print_fastsem(res)
#'
#' # Mediation with derived indirect effect
#' df2 <- as.data.frame(scale(mtcars[, c("mpg","hp","wt")]))
#' syntax2 <- "
#'   wt  ~ a*hp
#'   mpg ~ b*wt + c_prime*hp
#'   hp ~~ hp;  wt ~~ wt;  mpg ~~ mpg
#'   hp ~ 1;  wt ~ 1;  mpg ~ 1
#'   ab := a * b
#' "
#' res2 <- fastsem_fit(syntax2, df2)
#' cat("Indirect hp->wt->mpg:", res2$estimates[res2$paramNames == "ab"],
#'     "+/-", res2$se[res2$paramNames == "ab"], "\n")
#' }
fastsem_fit <- function(syntax, data, col_names = NULL, control = list()) {
  if (!isTRUE(.fastsem_env$loaded))
    stop("fastsem native library is not loaded. ",
         "Run fastsem_install() or fastsem_load().")

  mat <- .df_to_numeric_matrix(data)
  if (is.null(col_names)) col_names <- colnames(mat)
  if (is.null(col_names) || length(col_names) != ncol(mat))
    stop("fastsem_fit: col_names length must equal the number of data columns")

  # Run control travels as syntax directives (see .fastsem_control_directives).
  dirs <- .fastsem_control_directives(control)
  if (length(dirs))
    syntax <- paste(c(as.character(syntax)[1], dirs), collapse = "\n")
  ctl <- .fastsem_control_vec(control)
  sym <- getNativeSymbolInfo("fitSemR", .fastsem_env$dll_handle)
  res <- .Call(sym, as.character(syntax)[1], mat, as.character(col_names), ctl)
  if (is.null(res))
    stop("fastsem_fit: the engine returned no result (an error inside an older ",
         "engine build was not reported); rebuild/update the native library")
  if (!is.null(res$error)) stop(res$error, call. = FALSE)
  # Engines that return the free-parameter covariance send it flattened
  # (row-major); reshape to a named matrix for delta-method work.
  nf <- length(res$freeNames)
  if (nf > 0 && length(res$vcov) == nf * nf)
    res$vcov <- matrix(res$vcov, nf, nf, byrow = TRUE, dimnames = list(res$freeNames, res$freeNames))
  res
}

#' Extract observed variable names referenced by a lavaan syntax string.
#'
#' Parses the model syntax and returns only those variable names that are
#' also present in `col_names`.  Useful for subsetting a data frame before
#' passing it to [fastsem_fit()].
#'
#' @param syntax Character. lavaan model syntax string.
#' @param col_names Character vector of all column names in the data set.
#'
#' @return Character vector of variable names that appear in both `syntax`
#'   and `col_names`.
#'
#' @export
#' @examples
#' \dontrun{
#' syntax   <- "f =~ y1 + y2 + y3\nf ~~ f\ny1 ~~ y1\ny2 ~~ y2\ny3 ~~ y3"
#' all_cols <- c("id", "y1", "y2", "y3", "y4", "group")
#' fastsem_sem_vars(syntax, all_cols)
#' # [1] "y1" "y2" "y3"
#' }
fastsem_sem_vars <- function(syntax, col_names) {
  if (!isTRUE(.fastsem_env$loaded))
    stop("fastsem native library is not loaded. ",
         "Run fastsem_install() or fastsem_load().")
  sym <- getNativeSymbolInfo("extractSemVarsR", .fastsem_env$dll_handle)
  .Call(sym, as.character(syntax)[1], as.character(col_names))
}

# ── Deprecated umx bridge entry points ────────────────────────────────────────

#' Deprecated: translate an OpenMx RAM model to fastsem syntax.
#'
#' Use [run_fastsem()], which fits OpenMx/umx and lavaan models directly and
#' returns the same kind of object.  This now returns the syntax the engine
#' lowers the model to, which uses internal names (variables `v1, v2, ...`,
#' parameters `p1, p2, ...`; see `var_map` / `param_map`).
#'
#' @param model An `MxModel` (RAM, or a multi-group container of RAM models).
#' @param include_map Logical.  When `TRUE`, returns a list with `syntax`,
#'   `var_map` (model name → internal name) and `param_map` (internal label →
#'   OpenMx parameter name).
#' @return Character string, or a list when `include_map = TRUE`.
#' @export
umx_to_lavaan <- function(model, include_map = FALSE) {
  .Deprecated("run_fastsem", package = "fastsemR")
  pt  <- .ptab_from_mx(model)
  syn <- .ptab_syntax(pt)
  if (!include_map) return(syn$syntax)
  fr <- pt[pt$free & !duplicated(pt$plabel), , drop = FALSE]
  list(syntax = syn$syntax, var_map = syn$var_map,
       param_map = setNames(fr$omx, syn$plabel[fr$id]))
}

#' Deprecated: fit an OpenMx RAM model and return the raw fastsem result list.
#'
#' Use [run_fastsem()]; its fastsem result list is in `model@output$fastsem`.
#'
#' @param model An `MxModel`.
#' @param data Ignored (the model's own data are used).
#' @param ... Ignored.
#' @return The fastsem result list (see [fastsem_fit()]), with `lavaan_syntax`
#'   holding the syntax sent to the engine.
#' @export
umx_to_fastsem <- function(model, data = NULL, ...) {
  .Deprecated("run_fastsem", package = "fastsemR")
  if (!is.null(data))
    warning("umx_to_fastsem(data = ) is ignored; the model's own data are used.",
            call. = FALSE)
  fs <- .ptab_fit(.ptab_from_mx(model), list())
  fs$lavaan_syntax <- fs$syntax
  fs
}

#' Deprecated alias for [run_fastsem()].
#'
#' @inheritParams run_fastsem
#' @return The fitted `MxModel`.
#' @export
run_fastsem_ram <- function(model, control = list()) {
  .Deprecated("run_fastsem", package = "fastsemR")
  run_fastsem(model, control = control)
}

#' Print a fastsem result list in a human-readable table.
#'
#' Displays the estimator, fit statistics, parameter table (with estimates,
#' SEs, and z-scores), and — when present — the StdAll standardised
#' solution.
#'
#' @param res Named list as returned by [fastsem_fit()] or
#'   [umx_to_fastsem()].
#'
#' @return Invisibly, `res` (the original list, unchanged).
#'
#' @export
#' @examples
#' \dontrun{
#' res <- fastsem_fit(syntax, df)
#' print_fastsem(res)
#' }
print_fastsem <- function(res) {
  cat(sprintf("fastsem  --  %s  [SE: %s]\n", res$estimator, res$seType))
  cat(sprintf("  N=%-6d  free params=%-4d  df=%-4d\n",
              res$nobs, res$nFreeParams, res$df))
  cat(sprintf("  chi2=%8.4f  p=%6.4f  AIC=%10.4f  BIC=%10.4f\n\n",
              res$chi2, res$pvalue, res$aic, res$bic))

  # Equality-linked duplicates carry the same value+SE as their canonical
  # sibling (see Nim's `paramNamesResolved` and rnim_api's `isCanonical`
  # flag).  Hide them in the printed table — only the canonical row.
  show_idx <- if (!is.null(res$isCanonical) &&
                  length(res$isCanonical) == length(res$paramNames))
                which(res$isCanonical) else seq_along(res$paramNames)

  cat(sprintf("  %-3s  %-26s  %12s  %10s  %8s\n",
              "#", "Parameter", "Estimate", "SE", "z"))
  cat("  ", strrep("-", 65), "\n", sep = "")
  display_n <- 0L
  for (i in show_idx) {
    display_n <- display_n + 1L
    cat(sprintf("  %-3d  %-26s  %12.6f  %10.6f  %8.3f\n",
                display_n, res$paramNames[i],
                res$estimates[i], res$se[i], res$z[i]))
  }

  if (length(res$stdEst) > 0) {
    cat("\n  Standardized Estimates (StdAll)\n")
    cat("  ", strrep("-", 50), "\n", sep = "")
    display_n <- 0L
    for (i in show_idx) {
      display_n <- display_n + 1L
      cat(sprintf("  %-3d  %-26s  %10.6f\n",
                  display_n, res$paramNames[i], res$stdEst[i]))
    }
  }
  invisible(res)
}

# ── Internal helpers ───────────────────────────────────────────────────────────

`%||%` <- function(a, b) if (!is.null(a)) a else b

.df_to_numeric_matrix <- function(x) {
  if (is.matrix(x) && is.numeric(x)) return(x)
  if (is.data.frame(x)) {
    nm  <- colnames(x)
    mat <- matrix(NA_real_, nrow = nrow(x), ncol = ncol(x),
                  dimnames = list(NULL, nm))
    for (j in seq_len(ncol(x))) {
      col <- x[[j]]
      if (is.numeric(col))       mat[, j] <- col
      else if (is.logical(col))  mat[, j] <- as.numeric(col)
      else warning("fastsem: column '", nm[j], "' is non-numeric; treated as NA")
    }
    return(mat)
  }
  tryCatch(
    matrix(as.numeric(x), nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x)),
    error = function(e)
      stop("fastsem: cannot coerce data to numeric matrix: ", e$message)
  )
}
