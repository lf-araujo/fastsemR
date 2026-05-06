## fastsem.R — core SEM fitting and umx/OpenMx bridge

# ── FIML LBFGS control helper ────────────────────────────────────────────────
#
# `control` is a named list with optional fields:
#   maxIter   integer   max LBFGS iterations (default 2000 in core)
#   tol       numeric   gradient-norm tolerance (default 1e-6)
#   ftol      numeric   objective-change tolerance (default 1e-5)
#   tryHard   integer   extra retry rounds on non-convergence; each round
#                       multiplies maxIter by 3 and divides tol by 10.
#                       0 = no retry.
#
# Use `tryHard >= 1` for behavior analogous to OpenMx's `mxTryHard`.
# Returned: numeric(4) — the wire format consumed by the Nim ingester
# (`applyFimlControl` in rnim_api.nim).  Length-0 means "all defaults".
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
#'   per-group annotations.
#' * **Weighted observations** — `weight: colname`.
#' * **Cluster-robust SEs** — `cluster: colname`.
#'
#' @param syntax Character string of lavaan model syntax.
#' @param data Data frame or numeric matrix.  Non-numeric columns are coerced
#'   to `NA`.  Row order is preserved.
#' @param col_names Character vector of column names.  Defaults to
#'   `colnames(data)`.  Must include every observed variable named in
#'   `syntax`.
#'
#' @return A named list with the following elements:
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
#'   }
#'
#' @seealso [print_fastsem()] to display the result;
#'   [run_fastsem()] for the umx/OpenMx bridge.
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

  ctl <- .fastsem_control_vec(control)
  sym <- getNativeSymbolInfo("fitSemR", .fastsem_env$dll_handle)
  .Call(sym, as.character(syntax)[1], mat, as.character(col_names), ctl)
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

# ── umx / OpenMx bridge ───────────────────────────────────────────────────────

#' Convert a umx / OpenMx RAM model to lavaan syntax.
#'
#' Inspects the RAM matrices (A, S, M) of an `MxModel` and generates the
#' extended lavaan syntax string that [fastsem_fit()] expects.  This is
#' called internally by [run_fastsem()] but is also useful for inspection
#' or for further hand-editing before fitting.
#'
#' @section Translations performed:
#' \describe{
#'   \item{A matrix (directional)}{
#'     Latent → manifest paths become `factor =~ indicator` lines.
#'     Latent → latent become `~` regression lines.
#'     Manifest → manifest become `~` regression lines.
#'     A label starting with `"data."` is treated as a definition variable
#'     and emitted as `data.colname * indicator`.
#'   }
#'   \item{S matrix (symmetric)}{
#'     Free diagonal entries → `var ~~ var` (with `start()` or label if set).
#'     Fixed non-zero diagonal entries → `var ~~ value*var`.
#'     Free off-diagonal → `a ~~ b` covariance.
#'   }
#'   \item{M matrix (means)}{
#'     Free entries → `var ~ 1` (with `start()` if the starting value is
#'     non-zero).
#'     Fixed non-zero entries → `var ~ value*1`.
#'   }
#'   \item{Bounds}{
#'     `S@lbound` / `S@ubound` / `A@lbound` / `A@ubound` → `label > lo` /
#'     `label < hi` lines.  Variance lower bounds of 0 are suppressed
#'     because fastsem's log-variance parameterisation already prevents
#'     negative variances.
#'   }
#'   \item{Equality constraints}{
#'     Labels appearing on more than one free parameter are emitted as
#'     `label*var`, which fastsem treats as an equality constraint.
#'   }
#' }
#'
#' @param model An `MxModel` as returned by `umxRAM()` or `mxModel()`.
#' @param include_map Logical (default `FALSE`).  When `TRUE` returns a list
#'   with two elements:
#'   \describe{
#'     \item{`syntax`}{The lavaan syntax string.}
#'     \item{`param_map`}{Named character vector mapping fastsem structural
#'       names (e.g. `"f->y1"`) to OpenMx free-parameter labels
#'       (e.g. `"g_to_y1"`).  Used by [run_fastsem()] for inject-back.}
#'   }
#'
#' @return Character string (default) or named list when `include_map = TRUE`.
#'
#' @seealso [run_fastsem()] to fit and inject back in one step.
#'
#' @export
#' @examples
#' \dontrun{
#' library(umx)
#' df <- as.data.frame(scale(iris[, c(1, 3, 4)]))
#' names(df) <- c("sl", "pl", "pw")
#'
#' m <- umxRAM("CFA",
#'   umxPath(from = "g", to = c("sl", "pl", "pw")),
#'   umxPath(v1m0 = "g"),
#'   umxPath(var  = c("sl", "pl", "pw")),
#'   data = df, autoRun = FALSE
#' )
#'
#' # Inspect the generated syntax
#' cat(umx_to_lavaan(m))
#'
#' # Include the parameter map for custom inject-back
#' out <- umx_to_lavaan(m, include_map = TRUE)
#' str(out$param_map)
#' }
umx_to_lavaan <- function(model, include_map = FALSE) {
  if (!requireNamespace("OpenMx", quietly = TRUE))
    stop("The 'OpenMx' package is required for umx_to_lavaan().")

  man_vars <- model@manifestVars
  lat_vars <- model@latentVars
  A_mat    <- model@matrices$A
  S_mat    <- model@matrices$S

  all_vars <- rownames(A_mat@values)
  if (is.null(all_vars)) all_vars <- c(man_vars, lat_vars)

  A_vals <- A_mat@values
  A_free <- A_mat@free
  A_labs <- if (!is.null(A_mat@labels)) A_mat@labels else
              matrix(NA_character_, nrow = nrow(A_vals), ncol = ncol(A_vals),
                     dimnames = dimnames(A_vals))
  S_vals <- S_mat@values
  S_free <- S_mat@free
  S_labs <- if (!is.null(S_mat@labels)) S_mat@labels else
              matrix(NA_character_, nrow = nrow(S_vals), ncol = ncol(S_vals),
                     dimnames = dimnames(S_vals))

  lines     <- character(0)
  param_map <- character(0)

  # Pre-scan: labels appearing on >1 free parameter are equality constraints
  all_free_labs <- c(A_labs[A_free], S_labs[S_free])
  all_free_labs <- all_free_labs[!is.na(all_free_labs) & nzchar(all_free_labs)]
  equality_labels <- names(which(table(all_free_labs) > 1))

  emit_path_term <- function(from_var, to_var, term_var, val, free_flag, lab) {
    struct_name <- paste0(from_var, "->", to_var)
    if (isTRUE(free_flag)) {
      if (!is.na(lab) && nzchar(lab) && !startsWith(lab, "data.")) {
        is_eq <- lab %in% equality_labels
        if (is_eq) {
          param_map[struct_name] <<- lab
          return(paste0(lab, "*", term_var))
        }
        param_map[struct_name] <<- lab
        has_start <- !is.na(val) && is.finite(val) && val != 0 && val != 1
        if (has_start)
          return(sprintf("start(%g)*%s*%s", val, lab, term_var))
        return(paste0(lab, "*", term_var))
      }
      if (!is.na(val) && val != 0 && val != 1)
        return(sprintf("start(%g)*%s", val, term_var))
      return(term_var)
    } else {
      if (!is.na(lab) && startsWith(lab, "data."))
        return(paste0(lab, "*", term_var))
      if (!is.na(val) && val != 0)
        return(sprintf("%g*%s", val, term_var))
      return(NULL)
    }
  }

  # A matrix: latent factor loadings
  for (lv in lat_vars) {
    if (!(lv %in% colnames(A_vals))) next
    terms <- character(0)
    for (mv in man_vars) {
      if (!(mv %in% rownames(A_vals))) next
      tok <- emit_path_term(lv, mv, mv,
                            A_vals[mv, lv], A_free[mv, lv], A_labs[mv, lv])
      if (!is.null(tok)) terms <- c(terms, tok)
    }
    if (length(terms) > 0)
      lines <- c(lines, sprintf("%s =~ %s", lv, paste(terms, collapse = " + ")))
  }

  # A matrix: structural paths between latents
  for (lv_to in lat_vars) {
    if (!(lv_to %in% rownames(A_vals))) next
    preds <- character(0)
    for (lv_from in lat_vars) {
      if (lv_from == lv_to || !(lv_from %in% colnames(A_vals))) next
      tok <- emit_path_term(lv_from, lv_to, lv_from,
                            A_vals[lv_to, lv_from], A_free[lv_to, lv_from],
                            A_labs[lv_to, lv_from])
      if (!is.null(tok)) preds <- c(preds, tok)
    }
    if (length(preds) > 0)
      lines <- c(lines, sprintf("%s ~ %s", lv_to, paste(preds, collapse = " + ")))
  }

  # A matrix: regressions (manifest ~ manifest predictors)
  for (mv_to in man_vars) {
    if (!(mv_to %in% rownames(A_vals))) next
    preds <- character(0)
    for (iv in all_vars) {
      if (iv == mv_to || iv %in% lat_vars || !(iv %in% colnames(A_vals))) next
      tok <- emit_path_term(iv, mv_to, iv,
                            A_vals[mv_to, iv], A_free[mv_to, iv],
                            A_labs[mv_to, iv])
      if (!is.null(tok)) preds <- c(preds, tok)
    }
    if (length(preds) > 0)
      lines <- c(lines, sprintf("%s ~ %s", mv_to, paste(preds, collapse = " + ")))
  }

  # S matrix: variances and covariances
  n_vars <- length(all_vars)
  for (i in seq_len(n_vars)) {
    vi <- all_vars[i]
    if (!(vi %in% rownames(S_vals))) next

    if (isTRUE(S_free[vi, vi])) {
      struct_name <- paste0("Var(", vi, ")")
      lab <- S_labs[vi, vi];  sv <- S_vals[vi, vi]
      if (!is.na(lab) && nzchar(lab)) {
        param_map[struct_name] <- lab
        lines <- c(lines, sprintf("%s ~~ %s*%s", vi, lab, vi))
      } else if (!is.na(sv) && is.finite(sv) && sv > 0) {
        lines <- c(lines, sprintf("%s ~~ start(%g)*%s", vi, sv, vi))
      } else {
        lines <- c(lines, sprintf("%s ~~ %s", vi, vi))
      }
    } else {
      v <- S_vals[vi, vi]
      if (!is.na(v) && v != 0)
        lines <- c(lines, sprintf("%s ~~ %g*%s", vi, v, vi))
    }

    if (i < n_vars) {
      for (j in seq(i + 1L, n_vars)) {
        vj <- all_vars[j]
        if (!(vj %in% rownames(S_vals))) next
        if (isTRUE(S_free[vi, vj]) || isTRUE(S_free[vj, vi])) {
          struct_name <- paste0(vi, "~~", vj)
          lab <- S_labs[vi, vj]
          if (is.na(lab) || !nzchar(lab)) lab <- S_labs[vj, vi]
          sv  <- S_vals[vi, vj]
          if (!is.na(lab) && nzchar(lab)) {
            param_map[struct_name] <- lab
            lines <- c(lines, sprintf("%s ~~ %s*%s", vi, lab, vj))
          } else if (!is.na(sv) && is.finite(sv) && sv != 0) {
            lines <- c(lines, sprintf("%s ~~ start(%g)*%s", vi, sv, vj))
          } else {
            lines <- c(lines, sprintf("%s ~~ %s", vi, vj))
          }
        } else {
          v <- S_vals[vi, vj]
          if (!is.na(v) && v != 0)
            lines <- c(lines, sprintf("%s ~~ %g*%s", vi, v, vj))
        }
      }
    }
  }

  # M matrix: intercepts
  has_M <- "M" %in% names(model@matrices)
  if (has_M) {
    M_mat  <- model@matrices$M
    M_vals <- M_mat@values
    M_free <- M_mat@free
    M_labs <- if (!is.null(M_mat@labels)) M_mat@labels else
                matrix(NA_character_, nrow = 1L, ncol = ncol(M_vals),
                       dimnames = dimnames(M_vals))
    for (vn in colnames(M_vals)) {
      if (isTRUE(M_free[1, vn])) {
        struct_name <- paste0(vn, "~1")
        lab <- M_labs[1, vn];  sv <- M_vals[1, vn]
        if (!is.na(lab) && nzchar(lab)) param_map[struct_name] <- lab
        if (!is.na(sv) && is.finite(sv) && sv != 0)
          lines <- c(lines, sprintf("%s ~ start(%g)*1", vn, sv))
        else
          lines <- c(lines, sprintf("%s ~ 1", vn))
      } else {
        mv <- M_vals[1, vn]
        if (!is.na(mv) && is.finite(mv) && mv != 0)
          lines <- c(lines, sprintf("%s ~ %g*1", vn, mv))
      }
    }
  }

  # Bounds from S@lbound/ubound and A@lbound/ubound
  .emit_matrix_bounds <- function(mat_obj, acc) {
    if (is.null(mat_obj)) return(acc)
    lb <- mat_obj@lbound;  ub <- mat_obj@ubound
    labs <- mat_obj@labels;  fr <- mat_obj@free
    if (is.null(lb) && is.null(ub)) return(acc)
    for (i in seq_len(nrow(labs))) {
      for (j in seq_len(ncol(labs))) {
        if (!isTRUE(fr[i, j])) next
        lab <- labs[i, j]
        if (is.na(lab) || !nzchar(lab)) next
        lo <- if (!is.null(lb)) lb[i, j] else NA_real_
        hi <- if (!is.null(ub)) ub[i, j] else NA_real_
        is_diag <- (i == j)
        if (!is_diag && !is.na(lo) && is.finite(lo) && lo > 0)
          acc <- c(acc, sprintf("%s > %g", lab, lo))
        if (!is.na(hi) && is.finite(hi))
          acc <- c(acc, sprintf("%s < %g", lab, hi))
      }
    }
    acc
  }
  lines <- .emit_matrix_bounds(model@matrices$S, lines)
  lines <- .emit_matrix_bounds(model@matrices$A, lines)

  syntax <- paste(lines, collapse = "\n")
  if (include_map) list(syntax = syntax, param_map = param_map)
  else syntax
}

#' Fit a umx / OpenMx RAM model with fastsem and return the updated model.
#'
#' The primary high-level function for users who work with `umxRAM()` models.
#' It:
#' 1. Translates the model to lavaan syntax via [umx_to_lavaan()].
#' 2. Extracts raw data from `model@data@observed`.
#' 3. Fits the model with [fastsem_fit()].
#' 4. Writes estimates back into the model matrices and populates
#'    `model@output` so that `summary()`, `omxGetParameters()`, and
#'    `umxCompare()` work as if `mxRun()` had been called.
#'
#' @section Multi-group models:
#' When `model@submodels` is non-empty, [run_fastsem()] fits the joint
#' multi-group model and injects per-group estimates back into each submodel.
#' Two engines are available (selected via `engine`):
#' * **`"ram"`** (default for multi-group) sends the OpenMx RAM matrices
#'   (A, S, optionally M) directly to fastsem's Nim core via
#'   `fitSemFromMatricesMultigroupR`.  This bypasses the lavaan-syntax
#'   round-trip and honors cross-group fixed covariances exactly (e.g.
#'   cross-twin A=1.0 in MZ vs 0.5 in DZ).
#' * **`"lavaan"`** stacks data with a `fastsem_grp` indicator, builds a
#'   combined lavaan syntax with `c()` per-group annotations, and fits via
#'   [fastsem_fit()].  Use this when the RAM engine errors on a model that
#'   uses features it does not yet support (definition variables, ordinal
#'   thresholds, configural cross-group different-label free parameters).
#' Parameters that carry the same OpenMx label in all submodels are equated;
#' others are group-specific.
#'
#' @param model An `MxModel` as returned by `umxRAM()` (or a bare `mxModel`
#'   container holding submodels for multi-group use).  The model must have
#'   raw (row-level) data in `model@data@observed`.
#' @param control Optional named list passed to the optimizer.  Recognized
#'   fields: `maxIter` (int), `tol` (gradient-norm tolerance),
#'   `ftol` (objective-change tolerance), `tryHard` (extra retry rounds on
#'   non-convergence; each round multiplies `maxIter` by 3 and divides `tol`
#'   by 10).  Any field absent or `0`/`0.0` uses the core default.
#' @param engine One of `"auto"`, `"ram"`, `"lavaan"`.  For multi-group
#'   models, `"auto"` selects the RAM engine.  For single-group models, the
#'   lavaan engine is always used (RAM single-group is not yet wired in the
#'   R wrapper).
#'
#' @return The same `MxModel` object with:
#' \describe{
#'   \item{Matrix values updated}{A, S, and M matrices hold the fastsem
#'     point estimates.}
#'   \item{`@output` populated}{`$estimate`, `$standardErrors`, `$fit`
#'     (`-2lnL`), `$Chi`, `$ChiDoF`, plus `$.fastsem_*` slots for the
#'     raw fastsem result list.}
#'   \item{`@.wasRun` set to `TRUE`}{Makes `summary()` and `umxCompare()`
#'     work.}
#' }
#'
#' @seealso [umx_to_lavaan()] for syntax inspection;
#'   [umx_to_fastsem()] for the raw result list without inject-back;
#'   [fastsem_fit()] for fitting from raw lavaan syntax.
#'
#' @export
#' @examples
#' \dontrun{
#' library(umx)
#'
#' # ── CFA ──────────────────────────────────────────────────────────────────
#' df <- as.data.frame(scale(iris[, c(1, 3, 4)]))
#' names(df) <- c("sl", "pl", "pw")
#'
#' m <- umxRAM("CFA",
#'   umxPath(from = "g", to = c("sl", "pl", "pw")),
#'   umxPath(v1m0 = "g"),
#'   umxPath(var  = c("sl", "pl", "pw")),
#'   data = df, autoRun = FALSE
#' )
#' m_fit <- run_fastsem(m)
#' summary(m_fit)
#' omxGetParameters(m_fit)
#'
#' # ── Multi-group CFA ───────────────────────────────────────────────────────
#' data(HolzingerSwineford1939, package = "lavaan")
#' hs <- HolzingerSwineford1939
#'
#' mg_p <- umxRAM("mg_p",
#'   umxPath(from = "g", to = c("x1","x2","x3")),
#'   umxPath(v1m0 = "g"), umxPath(var = c("x1","x2","x3")),
#'   umxPath(means = c("x1","x2","x3")),
#'   data = hs[hs$school == "Pasteur", c("x1","x2","x3")]
#' )
#' mg_gw <- umxRAM("mg_gw",
#'   umxPath(from = "g", to = c("x1","x2","x3")),
#'   umxPath(v1m0 = "g"), umxPath(var = c("x1","x2","x3")),
#'   umxPath(means = c("x1","x2","x3")),
#'   data = hs[hs$school == "Grant-White", c("x1","x2","x3")]
#' )
#' r <- run_fastsem(mxModel("HS_mg", mg_p, mg_gw))
#' omxGetParameters(r@submodels$mg_p)
#' }
run_fastsem <- function(model, control = list(),
                        engine = c("auto", "ram", "lavaan")) {
  if (!requireNamespace("OpenMx", quietly = TRUE))
    stop("The 'OpenMx' package is required for run_fastsem().")
  if (!isTRUE(.fastsem_env$loaded))
    stop("fastsem native library is not loaded. ",
         "Run fastsem_install() or fastsem_load().")
  engine <- match.arg(engine)

  is_multigroup <- length(model@submodels) > 0

  if (!is_multigroup) {
    if (engine == "ram")
      stop("run_fastsem(engine = 'ram'): single-group RAM dispatch is not ",
           "yet wired in the R wrapper. Use engine = 'auto' or 'lavaan'.")
    return(.run_fastsem_lavaan_sg(model, control))
  }

  # Multi-group: default to RAM unless caller forced lavaan.
  if (engine == "lavaan")
    return(.run_fastsem_lavaan_mg(model, control))
  .run_fastsem_ram_mg(model, control)
}

# ── Multi-group lavaan engine ────────────────────────────────────────────────
.run_fastsem_lavaan_mg <- function(model, control) {
  sub_names    <- names(model@submodels)
  nGroups      <- length(sub_names)
  grp_col_name <- "fastsem_grp"

  grp_data_list <- vector("list", nGroups)
  for (gi in seq_along(sub_names)) {
    gdat <- .extract_model_data(model@submodels[[sub_names[gi]]])
    gdat[[grp_col_name]] <- gi
    grp_data_list[[gi]] <- gdat
  }
  stacked_data <- do.call(rbind, grp_data_list)

  lav_syntax <- .umx_multigroup_lavaan(model@submodels, sub_names, grp_col_name)
  col_names  <- colnames(stacked_data)
  fs_res     <- fastsem_fit(lav_syntax, stacked_data, col_names, control = control)

  for (gi in seq_along(sub_names)) {
    sub      <- model@submodels[[sub_names[gi]]]
    lav_info <- umx_to_lavaan(sub, include_map = TRUE)
    g_res    <- fs_res
    if (length(fs_res$groupEstimates) >= gi)
      g_res$estimates <- fs_res$groupEstimates[[gi]]
    if (length(fs_res$groupSEs) >= gi)
      g_res$se <- fs_res$groupSEs[[gi]]
    model@submodels[[sub_names[gi]]] <-
      .inject_fastsem_into_model(sub, g_res, lav_info$param_map)
  }

  model <- .populate_model_output(model,
    neg2logL = fs_res$neg2logL %||% NA_real_,
    chi2     = fs_res$chi2     %||% NA_real_,
    df       = fs_res$df       %||% NA_integer_,
    fs_res   = fs_res)
  model@.wasRun           <- TRUE
  model@.modifiedSinceRun <- FALSE
  model
}

# ── Single-group lavaan engine ───────────────────────────────────────────────
.run_fastsem_lavaan_sg <- function(model, control) {
  data      <- .extract_model_data(model)
  lav_info  <- umx_to_lavaan(model, include_map = TRUE)
  col_names <- colnames(data)
  if (is.null(col_names))
    stop("run_fastsem: data must have column names")
  fs_res <- fastsem_fit(lav_info$syntax, data, col_names, control = control)
  .inject_fastsem_into_model(model, fs_res, lav_info$param_map)
}

# ── Multi-group RAM engine ───────────────────────────────────────────────────
# Bypasses the lavaan syntax translation: extracts each submodel's RAM
# matrices (A, S, optionally M) with @free flags and @labels, packs them into
# flat group-major vectors, and calls fitSemFromMatricesMultigroupR.
# Cross-group equality is conveyed by shared OpenMx labels (auto-detected);
# per-group fixed values that differ across groups (e.g. the cross-twin A
# correlation fixed at 1.0 in MZ vs 0.5 in DZ for ACE/Cp twin models) flow
# through directly without `c()` annotations.
#
# v1 restrictions (mirroring the Nim ingester):
# • Same RAM structure across groups: a position is free in every group or
#   fixed in every group.
# • Free parameters equated across groups must share the OpenMx label.
# • No definition variables, no ordinal thresholds.
.run_fastsem_ram_mg <- function(model, control) {
  sub_names <- names(model@submodels)
  if (length(sub_names) < 2)
    stop("RAM engine requires ≥2 submodels.")
  nG <- length(sub_names)

  ref_sub  <- model@submodels[[sub_names[1]]]
  lat_vars <- ref_sub@latentVars
  man_vars <- ref_sub@manifestVars
  nLat     <- length(lat_vars)
  nObs     <- length(man_vars)
  p        <- nLat + nObs
  canonical_vars <- c(lat_vars, man_vars)

  for (gi in seq_along(sub_names)) {
    sub <- model@submodels[[sub_names[gi]]]
    if (!setequal(sub@latentVars,   lat_vars) ||
        !setequal(sub@manifestVars, man_vars))
      stop("run_fastsem_ram(): submodel '", sub_names[gi],
           "' has different latent/manifest vars from reference submodel '",
           sub_names[1], "'")
  }

  has_M <- "M" %in% names(ref_sub@matrices)

  aVals <- numeric(p * p * nG); aFree <- logical(p * p * nG); aLabs <- character(p * p * nG)
  sVals <- numeric(p * p * nG); sFree <- logical(p * p * nG); sLabs <- character(p * p * nG)
  if (has_M) {
    mVals <- numeric(p * nG); mFree <- logical(p * nG); mLabs <- character(p * nG)
  } else {
    mVals <- numeric(0); mFree <- logical(0); mLabs <- character(0)
  }

  # Bounds — group-0 only (umx applies bounds identically across groups).
  # NaN in any cell means "unbounded on that side"; both sides NaN = no bound.
  A_ref <- ref_sub@matrices$A
  S_ref <- ref_sub@matrices$S
  aLbound <- as.numeric(A_ref@lbound[canonical_vars, canonical_vars, drop = FALSE])
  aUbound <- as.numeric(A_ref@ubound[canonical_vars, canonical_vars, drop = FALSE])
  sLbound <- as.numeric(S_ref@lbound[canonical_vars, canonical_vars, drop = FALSE])
  sUbound <- as.numeric(S_ref@ubound[canonical_vars, canonical_vars, drop = FALSE])
  if (has_M) {
    M_ref   <- ref_sub@matrices$M
    mLbound <- as.numeric(M_ref@lbound[1, canonical_vars, drop = TRUE])
    mUbound <- as.numeric(M_ref@ubound[1, canonical_vars, drop = TRUE])
  } else {
    mLbound <- numeric(0); mUbound <- numeric(0)
  }

  for (gi in seq_along(sub_names)) {
    sub <- model@submodels[[sub_names[gi]]]
    A   <- sub@matrices$A
    S   <- sub@matrices$S

    A_vals <- A@values [canonical_vars, canonical_vars, drop = FALSE]
    A_free <- A@free   [canonical_vars, canonical_vars, drop = FALSE]
    A_labs <- A@labels [canonical_vars, canonical_vars, drop = FALSE]
    A_labs[is.na(A_labs)] <- ""

    S_vals <- S@values [canonical_vars, canonical_vars, drop = FALSE]
    S_free <- S@free   [canonical_vars, canonical_vars, drop = FALSE]
    S_labs <- S@labels [canonical_vars, canonical_vars, drop = FALSE]
    S_labs[is.na(S_labs)] <- ""

    g_off <- (gi - 1L) * p * p
    aVals[g_off + seq_len(p * p)] <- as.numeric(A_vals)
    aFree[g_off + seq_len(p * p)] <- as.logical(A_free)
    aLabs[g_off + seq_len(p * p)] <- as.character(A_labs)
    sVals[g_off + seq_len(p * p)] <- as.numeric(S_vals)
    sFree[g_off + seq_len(p * p)] <- as.logical(S_free)
    sLabs[g_off + seq_len(p * p)] <- as.character(S_labs)

    if (has_M) {
      M <- sub@matrices$M
      M_vals <- M@values[1, canonical_vars, drop = TRUE]
      M_free <- M@free  [1, canonical_vars, drop = TRUE]
      M_labs <- M@labels[1, canonical_vars, drop = TRUE]
      M_labs[is.na(M_labs)] <- ""
      m_off <- (gi - 1L) * p
      mVals[m_off + seq_len(p)] <- as.numeric(M_vals)
      mFree[m_off + seq_len(p)] <- as.logical(M_free)
      mLabs[m_off + seq_len(p)] <- as.character(M_labs)
    }
  }

  grp_data_list <- vector("list", nG)
  for (gi in seq_along(sub_names)) {
    gdat <- .extract_model_data(model@submodels[[sub_names[gi]]])
    gdat[["fastsem_grp"]] <- gi
    grp_data_list[[gi]] <- gdat
  }
  stacked <- do.call(rbind, grp_data_list)
  stacked_mat <- data.matrix(stacked)
  storage.mode(stacked_mat) <- "double"
  col_names <- colnames(stacked_mat)
  if (is.null(col_names))
    stop("run_fastsem_ram(): stacked data is missing column names")

  ctl <- .fastsem_control_vec(control)
  sym <- getNativeSymbolInfo("fitSemFromMatricesMultigroupR",
                              .fastsem_env$dll_handle)
  fs_res <- .Call(sym,
                  aVals, aFree, aLabs,
                  sVals, sFree, sLabs,
                  mVals, mFree, mLabs,
                  aLbound, aUbound,
                  sLbound, sUbound,
                  mLbound, mUbound,
                  lat_vars, man_vars,
                  stacked_mat, col_names, ctl)

  for (gi in seq_along(sub_names)) {
    sub   <- model@submodels[[sub_names[gi]]]
    g_res <- fs_res
    if (length(fs_res$groupEstimates) >= gi)
      g_res$estimates <- fs_res$groupEstimates[[gi]]
    if (length(fs_res$groupSEs) >= gi)
      g_res$se <- fs_res$groupSEs[[gi]]
    # Build a structural-name → OpenMx-label map so SEs land in the right
    # @output$standardErrors slot.  Keys: "from->to", "Var(x)", "x~~y", "x~1".
    param_map <- .ram_build_param_map(sub)
    model@submodels[[sub_names[gi]]] <-
      .inject_fastsem_into_model(sub, g_res, param_map)
  }

  model <- .populate_model_output(model,
    neg2logL = fs_res$neg2logL %||% NA_real_,
    chi2     = fs_res$chi2     %||% NA_real_,
    df       = fs_res$df       %||% NA_integer_,
    fs_res   = fs_res)

  # Aggregate parameters across submodels into the parent model's @output so
  # summary() can render the parameter table.  omxGetParameters() crawls the
  # submodels and merges by label; submodel @output$standardErrors are pulled
  # in by matching label.
  parent_estimate <- tryCatch(OpenMx::omxGetParameters(model),
                              error = function(e) setNames(numeric(0), character(0)))
  fp_names <- names(parent_estimate)
  parent_se <- setNames(rep(NA_real_, length(fp_names)), fp_names)
  for (sn in sub_names) {
    sub_se_mat <- model@submodels[[sn]]@output$standardErrors
    if (is.null(sub_se_mat) || !is.matrix(sub_se_mat)) next
    for (lbl in rownames(sub_se_mat)) {
      if (lbl %in% fp_names && is.na(parent_se[[lbl]]))
        parent_se[[lbl]] <- sub_se_mat[lbl, 1]
    }
  }
  model@output$estimate       <- parent_estimate
  model@output$standardErrors <- matrix(parent_se, ncol = 1L,
                                         dimnames = list(names(parent_se), "SE"))
  model@.wasRun           <- TRUE
  model@.modifiedSinceRun <- FALSE

  model
}

#' Deprecated alias for `run_fastsem(..., engine = "ram")`.
#'
#' Retained for backwards compatibility.  New code should call
#' [run_fastsem()] with `engine = "ram"` (the default for multi-group
#' models) or `engine = "auto"`.
#'
#' @inheritParams run_fastsem
#' @return The fitted `MxModel`.
#' @export
run_fastsem_ram <- function(model, control = list()) {
  .Deprecated("run_fastsem(model, engine = \"ram\")", package = "fastsemR")
  run_fastsem(model, control = control, engine = "ram")
}

#' Fit a umx model and return the raw fastsem result list.
#'
#' A lower-level alternative to [run_fastsem()] that returns the raw fastsem
#' result list without injecting estimates back into the model object.  Useful
#' for inspecting fit statistics or parameter tables before committing to the
#' OpenMx model representation, or when you do not need the full inject-back
#' machinery.
#'
#' @param model An `MxModel` object (from `umxRAM()` etc.).
#' @param data Data frame of observed data.  Defaults to
#'   `model@data@observed` when `NULL`.
#' @param ... Reserved for future use.
#'
#' @return Named list as returned by [fastsem_fit()], with an additional
#'   `$lavaan_syntax` element containing the translated lavaan syntax string
#'   for inspection.
#'
#' @seealso [run_fastsem()] for the full inject-back workflow;
#'   [fastsem_fit()] for fitting from raw lavaan syntax.
#'
#' @export
#' @examples
#' \dontrun{
#' library(umx)
#' df <- as.data.frame(scale(mtcars[, c("mpg","hp","wt")]))
#'
#' m <- umxRAM("path",
#'   umxPath(from = c("hp","wt"), to = "mpg"),
#'   umxPath(cov   = c("hp","wt")),
#'   umxPath(var   = c("hp","wt","mpg")),
#'   umxPath(means = c("hp","wt","mpg")),
#'   data = df, autoRun = FALSE
#' )
#' res <- umx_to_fastsem(m)
#' cat("lavaan syntax used:\n", res$lavaan_syntax, "\n")
#' cat("chi2 =", res$chi2, "  df =", res$df, "\n")
#' }
umx_to_fastsem <- function(model, data = NULL, ...) {
  if (!requireNamespace("OpenMx", quietly = TRUE))
    stop("The 'OpenMx' package is required for umx_to_fastsem().")
  if (!isTRUE(.fastsem_env$loaded))
    stop("fastsem native library is not loaded. ",
         "Run fastsem_install() or fastsem_load().")

  if (is.null(data)) data <- .extract_model_data(model)
  mat       <- .df_to_numeric_matrix(data)
  col_names <- colnames(mat)
  if (is.null(col_names))
    stop("umx_to_fastsem: data must have column names")

  lav_info <- umx_to_lavaan(model, include_map = TRUE)
  res      <- fastsem_fit(lav_info$syntax, mat, col_names, ...)
  res$lavaan_syntax <- lav_info$syntax
  res
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

.extract_model_data <- function(model) {
  dat <- model@data
  if (is.null(dat))
    stop("run_fastsem: model has no @data slot.")
  if (!is.null(dat@type) && dat@type != "raw")
    stop("run_fastsem: only raw (row-level) data is supported; got type = '",
         dat@type, "'.")
  obs <- dat@observed
  if (is.null(obs)) stop("run_fastsem: model@data@observed is NULL.")
  obs
}

.inject_fastsem_into_model <- function(model, fs_res, param_map = NULL) {
  n_fs           <- length(fs_res$estimates)
  struct_to_fidx <- setNames(seq_len(n_fs), fs_res$paramNames)

  label_to_fidx <- list()
  if (!is.null(param_map) && length(param_map) > 0) {
    for (sname in names(param_map)) {
      omx_lab <- param_map[[sname]]
      fidx_v  <- struct_to_fidx[sname]
      if (!is.na(fidx_v) && !is.na(omx_lab) && nzchar(omx_lab))
        label_to_fidx[[omx_lab]] <- fidx_v[[1]]
    }
  }

  lat_vars   <- model@latentVars
  man_vars   <- model@manifestVars
  model_vars <- c(lat_vars, man_vars)
  has_M      <- "M" %in% names(model@matrices)

  for (sname in names(struct_to_fidx)) {
    fidx <- struct_to_fidx[[sname]]
    est  <- fs_res$estimates[[fidx]]
    if (is.na(est)) next

    if (grepl("->", sname, fixed = TRUE)) {
      parts <- strsplit(sname, "->", fixed = TRUE)[[1]]
      if (length(parts) == 2 &&
          parts[1] %in% model_vars && parts[2] %in% model_vars)
        model@matrices$A@values[parts[2], parts[1]] <- est

    } else if (startsWith(sname, "Var(") && endsWith(sname, ")")) {
      vn <- substr(sname, 5L, nchar(sname) - 1L)
      if (vn %in% model_vars) model@matrices$S@values[vn, vn] <- est

    } else if (grepl("~~", sname, fixed = TRUE)) {
      parts <- strsplit(sname, "~~", fixed = TRUE)[[1]]
      if (length(parts) == 2 &&
          parts[1] %in% model_vars && parts[2] %in% model_vars) {
        model@matrices$S@values[parts[1], parts[2]] <- est
        model@matrices$S@values[parts[2], parts[1]] <- est
      }

    } else if (endsWith(sname, "~1")) {
      vn <- substr(sname, 1L, nchar(sname) - 2L)
      if (has_M && vn %in% colnames(model@matrices$M@values))
        model@matrices$M@values[1L, vn] <- est
    }
  }

  out_estimates <- tryCatch(OpenMx::omxGetParameters(model),
                            error = function(e) setNames(numeric(0), character(0)))
  fp_names <- names(out_estimates)

  out_se <- setNames(rep(NA_real_, length(fp_names)), fp_names)
  for (fp in fp_names) {
    fidx <- label_to_fidx[[fp]]
    if (is.null(fidx)) { v <- struct_to_fidx[fp]; if (!is.na(v)) fidx <- v[[1]] }
    if (!is.null(fidx) && !is.na(fidx) && fidx <= length(fs_res$se))
      out_se[[fp]] <- fs_res$se[[fidx]]
  }

  se_mat <- matrix(out_se, ncol = 1L, dimnames = list(names(out_se), "SE"))

  model <- .populate_model_output(model,
    neg2logL = fs_res$neg2logL, chi2 = fs_res$chi2, df = fs_res$df,
    fs_res = fs_res)
  model@output$estimate       <- out_estimates
  model@output$standardErrors <- se_mat

  # Populate output$expected.means so that summary() displays the fitted means.
  # fs_res$muHat is a named numeric vector (observed-variable order) returned
  # by fitSemR; it holds (I-A)^{-1} alpha filtered to observed variables.
  if (!is.null(fs_res$muHat) && length(fs_res$muHat) > 0) {
    exp_means <- matrix(fs_res$muHat, nrow = 1L,
                        dimnames = list(NULL, names(fs_res$muHat)))
    model@output$expected.means <- exp_means
  }

  model@.wasRun              <- TRUE
  model@.modifiedSinceRun    <- FALSE
  model
}

.populate_model_output <- function(model, neg2logL, chi2, df, fs_res = NULL) {
  model@output$fit                 <- neg2logL
  model@output$Minus2LogLikelihood <- neg2logL
  model@output$status       <- list(code = 0L, status = "OK", message = "",
                                    iterations = 1L, evaluations = 1L,
                                    gradient.norm = 0)
  model@output$infoDefinite <- TRUE
  model@output$iterations   <- 1L
  model@output$timestamp    <- format(Sys.time(), "%a %b %d %H:%M:%S %Z %Y")
  model@output$wallTime     <- 0
  model@output$mxVersion    <- tryCatch(
    as.character(utils::packageVersion("OpenMx")), error = function(e) "unknown")
  model@output$fitUnits     <- "-2lnL"

  if (!is.null(chi2) && !is.na(chi2)) {
    model@output$saturatedLikelihood <- neg2logL - chi2
    model@output$Chi                 <- chi2
    model@output$ChiDoF              <- df
  }

  if (!is.null(fs_res)) {
    model@output$.fastsem_chi2   <- fs_res$chi2
    model@output$.fastsem_df     <- fs_res$df
    model@output$.fastsem_pvalue <- fs_res$pvalue
    model@output$.fastsem_aic    <- fs_res$aic
    model@output$.fastsem_bic    <- fs_res$bic
    model@output$.fastsem_result <- fs_res
  }
  model
}

.ram_build_param_map <- function(sub) {
  # Build a {fastsem_struct_name → OpenMx_label} map for a single submodel,
  # mirroring the structural-name conventions emitted by the multi-group RAM
  # ingester in rnim_api.nim:
  #   A[row, col] free  →  "<col_name>-><row_name>"
  #   S[i,i]    free    →  "Var(<i_name>)"
  #   S[i,j] i<j free   →  "<i_name>~~<j_name>"
  #   M[1, i]   free    →  "<i_name>~1"
  pm <- character(0)
  A  <- sub@matrices$A
  S  <- sub@matrices$S
  has_M <- "M" %in% names(sub@matrices)

  vn <- rownames(A@values)
  for (col in seq_along(vn)) for (row in seq_along(vn)) {
    if (row == col) next
    if (isTRUE(A@free[row, col])) {
      lab <- A@labels[row, col]
      if (!is.na(lab) && nzchar(lab)) {
        struct <- paste0(vn[col], "->", vn[row])
        pm[struct] <- lab
      }
    }
  }
  for (i in seq_along(vn)) {
    if (isTRUE(S@free[i, i])) {
      lab <- S@labels[i, i]
      if (!is.na(lab) && nzchar(lab))
        pm[paste0("Var(", vn[i], ")")] <- lab
    }
    for (j in seq(i + 1L, length(vn))) {
      if (j > length(vn)) break
      if (isTRUE(S@free[i, j]) || isTRUE(S@free[j, i])) {
        lab <- S@labels[i, j]
        if (is.na(lab) || !nzchar(lab)) lab <- S@labels[j, i]
        if (!is.na(lab) && nzchar(lab))
          pm[paste0(vn[i], "~~", vn[j])] <- lab
      }
    }
  }
  if (has_M) {
    M  <- sub@matrices$M
    mvn <- colnames(M@values)
    for (i in seq_along(mvn)) {
      if (isTRUE(M@free[1, i])) {
        lab <- M@labels[1, i]
        if (!is.na(lab) && nzchar(lab))
          pm[paste0(mvn[i], "~1")] <- lab
      }
    }
  }
  pm
}

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

# ── Multi-group lavaan builder ─────────────────────────────────────────────────

.umx_multigroup_lavaan <- function(submodels, sub_names, grp_col_name) {
  nG    <- length(sub_names)
  ref   <- submodels[[sub_names[1]]]
  lines <- c(paste0("group: ", grp_col_name))

  A_ref    <- ref@matrices$A
  S_ref    <- ref@matrices$S
  has_M    <- "M" %in% names(ref@matrices)
  all_vars <- rownames(A_ref@values)
  lat_vars <- ref@latentVars
  man_vars <- ref@manifestVars

  a_lab <- function(g, row, col) {
    mat <- submodels[[sub_names[g]]]@matrices$A
    if (!isTRUE(mat@free[row, col])) return(NA_character_)
    lab <- mat@labels[row, col]
    if (is.null(lab)) NA_character_ else lab
  }
  s_lab <- function(g, row, col) {
    mat <- submodels[[sub_names[g]]]@matrices$S
    if (!isTRUE(mat@free[row, col]) && !isTRUE(mat@free[col, row]))
      return(NA_character_)
    lab <- mat@labels[row, col]
    if (is.null(lab) || is.na(lab)) lab <- mat@labels[col, row]
    if (is.null(lab)) NA_character_ else lab
  }

  c_term <- function(dest_var, labels_g, vals_g = NULL) {
    non_na_labs <- labels_g[!is.na(labels_g)]
    if (length(non_na_labs) == nG && length(unique(non_na_labs)) == 1)
      return(paste0(non_na_labs[1], "*", dest_var))
    lab_counts <- if (length(non_na_labs) > 0) table(non_na_labs) else integer(0)
    c_vals <- vapply(seq_len(nG), function(g) {
      lab <- labels_g[g]
      if (!is.null(vals_g) && !is.na(vals_g[g])) return(as.character(vals_g[g]))
      if (is.na(lab)) return("NA")
      if (lab_counts[lab] > 1L) return(lab)
      "NA"
    }, character(1))
    paste0("c(", paste(c_vals, collapse = ", "), ")*", dest_var)
  }

  for (lv in lat_vars) {
    if (!(lv %in% colnames(A_ref@values))) next
    terms <- character(0)
    for (mv in man_vars) {
      if (!(mv %in% rownames(A_ref@values))) next
      free_any <- any(vapply(seq_len(nG),
        function(g) isTRUE(submodels[[sub_names[g]]]@matrices$A@free[mv, lv]),
        logical(1)))
      val_ref <- A_ref@values[mv, lv];  lab_ref <- A_ref@labels[mv, lv]
      if (!free_any) {
        if (!is.na(lab_ref) && startsWith(lab_ref, "data."))
          terms <- c(terms, paste0(lab_ref, "*", mv))
        else if (!is.na(val_ref) && val_ref != 0)
          terms <- c(terms, sprintf("%g*%s", val_ref, mv))
      } else {
        labs_g <- vapply(seq_len(nG), function(g) a_lab(g, mv, lv), character(1))
        terms  <- c(terms, c_term(mv, labs_g, NULL))
      }
    }
    if (length(terms) > 0)
      lines <- c(lines, sprintf("%s =~ %s", lv, paste(terms, collapse = " + ")))
  }

  for (lv_to in lat_vars) {
    if (!(lv_to %in% rownames(A_ref@values))) next
    preds <- character(0)
    for (lv_from in lat_vars) {
      if (lv_from == lv_to || !(lv_from %in% colnames(A_ref@values))) next
      free_any <- any(vapply(seq_len(nG),
        function(g) isTRUE(submodels[[sub_names[g]]]@matrices$A@free[lv_to, lv_from]),
        logical(1)))
      if (!free_any) {
        v <- A_ref@values[lv_to, lv_from]
        if (!is.na(v) && v != 0) preds <- c(preds, sprintf("%g*%s", v, lv_from))
      } else {
        labs_g <- vapply(seq_len(nG), function(g) a_lab(g, lv_to, lv_from), character(1))
        preds  <- c(preds, c_term(lv_from, labs_g, NULL))
      }
    }
    if (length(preds) > 0)
      lines <- c(lines, sprintf("%s ~ %s", lv_to, paste(preds, collapse = " + ")))
  }

  for (mv_to in man_vars) {
    if (!(mv_to %in% rownames(A_ref@values))) next
    preds <- character(0)
    for (iv in all_vars) {
      if (iv == mv_to || iv %in% lat_vars || !(iv %in% colnames(A_ref@values))) next
      free_any <- any(vapply(seq_len(nG),
        function(g) isTRUE(submodels[[sub_names[g]]]@matrices$A@free[mv_to, iv]),
        logical(1)))
      if (!free_any) {
        v <- A_ref@values[mv_to, iv]
        if (!is.na(v) && v != 0) preds <- c(preds, sprintf("%g*%s", v, iv))
      } else {
        labs_g <- vapply(seq_len(nG), function(g) a_lab(g, mv_to, iv), character(1))
        preds  <- c(preds, c_term(iv, labs_g, NULL))
      }
    }
    if (length(preds) > 0)
      lines <- c(lines, sprintf("%s ~ %s", mv_to, paste(preds, collapse = " + ")))
  }

  nv <- length(all_vars)
  for (i in seq_len(nv)) {
    vi <- all_vars[i]
    if (!(vi %in% rownames(S_ref@values))) next
    free_any <- any(vapply(seq_len(nG),
      function(g) isTRUE(submodels[[sub_names[g]]]@matrices$S@free[vi, vi]),
      logical(1)))
    if (free_any) {
      labs_g <- vapply(seq_len(nG), function(g) s_lab(g, vi, vi), character(1))
      lines  <- c(lines, sprintf("%s ~~ %s", vi, c_term(vi, labs_g, NULL)))
    } else {
      fv <- S_ref@values[vi, vi]
      if (!is.na(fv) && fv != 0)
        lines <- c(lines, sprintf("%s ~~ %g*%s", vi, fv, vi))
    }
    if (i < nv) {
      for (j in seq(i + 1L, nv)) {
        vj <- all_vars[j]
        if (!(vj %in% rownames(S_ref@values))) next
        free_any <- any(vapply(seq_len(nG), function(g) {
          m <- submodels[[sub_names[g]]]@matrices$S
          isTRUE(m@free[vi, vj]) || isTRUE(m@free[vj, vi])
        }, logical(1)))
        if (free_any) {
          labs_g <- vapply(seq_len(nG), function(g) s_lab(g, vi, vj), character(1))
          lines  <- c(lines, sprintf("%s ~~ %s", vi, c_term(vj, labs_g, NULL)))
        }
      }
    }
  }

  if (has_M) {
    M_ref <- ref@matrices$M
    for (vn in colnames(M_ref@values)) {
      free_any <- any(vapply(seq_len(nG),
        function(g) isTRUE(submodels[[sub_names[g]]]@matrices$M@free[1, vn]),
        logical(1)))
      if (free_any) {
        sv <- M_ref@values[1, vn]
        if (!is.na(sv) && is.finite(sv) && sv != 0)
          lines <- c(lines, sprintf("%s ~ start(%g)*1", vn, sv))
        else
          lines <- c(lines, sprintf("%s ~ 1", vn))
      }
    }
  }

  paste(lines, collapse = "\n")
}
