## fastsem.R — core SEM fitting and umx/OpenMx bridge

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
fastsem_fit <- function(syntax, data, col_names = NULL) {
  if (!isTRUE(.fastsem_env$loaded))
    stop("fastsem native library is not loaded. ",
         "Run fastsem_install() or fastsem_load().")

  mat <- .df_to_numeric_matrix(data)
  if (is.null(col_names)) col_names <- colnames(mat)
  if (is.null(col_names) || length(col_names) != ncol(mat))
    stop("fastsem_fit: col_names length must equal the number of data columns")

  sym <- getNativeSymbolInfo("fitSemR", .fastsem_env$dll_handle)
  .Call(sym, as.character(syntax)[1], mat, as.character(col_names))
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
#' When `model@submodels` is non-empty, [run_fastsem()] stacks the data from
#' all submodels with a `fastsem_grp` group indicator, builds a combined
#' lavaan syntax with `c()` per-group annotations (see [fastsem_fit()]),
#' fits the joint model, and injects per-group estimates back into each
#' submodel.  Parameters that carry the same OpenMx label in all submodels
#' are equated (metric invariance for that parameter); others are
#' group-specific (configural).
#'
#' @param model An `MxModel` as returned by `umxRAM()` (or a bare `mxModel`
#'   container holding submodels for multi-group use).  The model must have
#'   raw (row-level) data in `model@data@observed`.
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
run_fastsem <- function(model) {
  if (!requireNamespace("OpenMx", quietly = TRUE))
    stop("The 'OpenMx' package is required for run_fastsem().")
  if (!isTRUE(.fastsem_env$loaded))
    stop("fastsem native library is not loaded. ",
         "Run fastsem_install() or fastsem_load().")

  is_multigroup <- length(model@submodels) > 0

  if (is_multigroup) {
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
    fs_res     <- fastsem_fit(lav_syntax, stacked_data, col_names)

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

  } else {
    data      <- .extract_model_data(model)
    lav_info  <- umx_to_lavaan(model, include_map = TRUE)
    col_names <- colnames(data)
    if (is.null(col_names))
      stop("run_fastsem: data must have column names")
    fs_res <- fastsem_fit(lav_info$syntax, data, col_names)
    model  <- .inject_fastsem_into_model(model, fs_res, lav_info$param_map)
  }

  model
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

  cat(sprintf("  %-3s  %-26s  %12s  %10s  %8s\n",
              "#", "Parameter", "Estimate", "SE", "z"))
  cat("  ", strrep("-", 65), "\n", sep = "")
  for (i in seq_along(res$paramNames)) {
    cat(sprintf("  %-3d  %-26s  %12.6f  %10.6f  %8.3f\n",
                i, res$paramNames[i], res$estimates[i], res$se[i], res$z[i]))
  }

  if (length(res$stdEst) > 0) {
    cat("\n  Standardized Estimates (StdAll)\n")
    cat("  ", strrep("-", 50), "\n", sep = "")
    for (i in seq_along(res$paramNames))
      cat(sprintf("  %-3d  %-26s  %10.6f\n",
                  i, res$paramNames[i], res$stdEst[i]))
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
  model@output$fit          <- neg2logL
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
