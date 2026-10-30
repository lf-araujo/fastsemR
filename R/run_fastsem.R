## run_fastsem.R — fit an OpenMx/umx or lavaan model with fastsem and return the
## same kind of object, so the rest of either workflow runs unchanged.

#' Fit an OpenMx/umx or lavaan model with fastsem.
#'
#' Takes a model specified in OpenMx/umx or lavaan, fits it with the fastsem
#' engine, and returns the same kind of object with the fastsem estimates in
#' place, so downstream tools keep working:
#' * an `MxModel` (`umxRAM()`, `mxModel(type = "RAM")`, or a multi-group
#'   container of RAM models) comes back as an `MxModel` for `summary()`,
#'   `umxSummary()`, `umxCompare()`, `vcov()`, `mxStandardizeRAMpaths()`, plots;
#' * a `lavaan` object built with `do.fit = FALSE` comes back as a fitted
#'   `lavaan` object for `summary()`, `fitMeasures()`, `parameterEstimates()`,
#'   `standardizedSolution()`, `lavTestLRT()`.
#'
#' @section How it works:
#' The model is read into a parameter table (one row per RAM cell and group,
#' like lavaan's `parTable()`), written as fastsem syntax with safe internal
#' names and full-precision numbers, and fitted.  The estimates are then
#' written back and the original package evaluates its own model at them
#' (`mxRun(useOptimizer = FALSE)`; lavaan with `optim.method = "none"`), so
#' fit statistics, algebras and standardised solutions are computed natively.
#' For OpenMx, fastsem's standard errors and parameter covariance are attached
#' (`@output$standardErrors`, `$vcov`, `$hessian`); lavaan computes its own
#' (including robust ones for MLR and friends).  The fastsem result list is kept
#' in `model@output$fastsem` / `fit@external$fastsem`.
#'
#' @section Self-check:
#' The −2lnL that OpenMx/lavaan compute at the fastsem estimates must equal the
#' engine's own.  A difference means the translation lost something; `check`
#' decides whether that warns, errors, or is ignored.  Bounds are enforced by
#' the engine (an estimate sits exactly on a bound when it binds, as in
#' OpenMx); inequalities between two parameters cannot be sent and are only
#' checked, the same way.
#'
#' @section Supported models:
#' RAM models with any mix of labelled / unlabelled, fixed and free paths,
#' (co)variances, means and thresholds; equality through shared labels (within
#' and across groups); multi-group models with shared and group-specific
#' parameters; definition variables on paths from latent variables; `mxConstraint`
#' equalities between labels and lavaan `==` / `:=`.  Not supported (clear
#' error): umx twin supermodels built from algebras (`umxACE`, `umxCP`, ...),
#' cells bound to algebras, definition variables on variances or means,
#' lavaan models with ordered indicators.  Means fixed at a non-zero value
#' need at least one free mean in the model.
#'
#' @param model An `MxModel` or an unfitted `lavaan` object.
#' @param control Optional named list of run-control settings (`maxIter`, `tol`,
#'   `ftol`, `tryHard`, `quadrature`); see *Model directives* in [fastsem_fit()].
#' @param check What to do when the self-check fails: `"warn"` (default),
#'   `"error"` or `"off"`.  Defaults to `getOption("fastsemR.check", "warn")`.
#' @param engine Deprecated and ignored.
#'
#' @return An object of the same class as `model`, fitted by fastsem.
#'
#' @seealso [fastsem_fit()] to fit fastsem syntax directly.
#'
#' @export
#' @examples
#' \dontrun{
#' # ── umx / OpenMx ─────────────────────────────────────────────────────────
#' library(umx)
#' df <- as.data.frame(scale(iris[, c(1, 3, 4)]))
#' names(df) <- c("sl", "pl", "pw")
#' m <- umxRAM("CFA",
#'   umxPath(from = "g", to = c("sl", "pl", "pw")),
#'   umxPath(v1m0 = "g"),
#'   umxPath(var  = c("sl", "pl", "pw")),
#'   data = df, autoRun = FALSE)
#' m_fit <- run_fastsem(m)
#' umxSummary(m_fit)
#'
#' # Multi-group: same labels are equal across groups, as in OpenMx
#' hs <- lavaan::HolzingerSwineford1939
#' grp <- function(nm, d) umxRAM(nm,
#'   umxPath(from = "g", to = c("x1", "x2", "x3")),
#'   umxPath(v1m0 = "g"), umxPath(var = c("x1", "x2", "x3")),
#'   umxPath(means = c("x1", "x2", "x3")),
#'   data = d[, c("x1", "x2", "x3")], autoRun = FALSE)
#' mg <- umxSuperModel("HS", grp("pasteur", hs[hs$school == "Pasteur", ]),
#'                           grp("grant",   hs[hs$school == "Grant-White", ]),
#'                     autoRun = FALSE)
#' summary(run_fastsem(mg))
#'
#' # ── lavaan ───────────────────────────────────────────────────────────────
#' library(lavaan)
#' fit0 <- cfa("visual =~ x1 + x2 + x3", HolzingerSwineford1939,
#'             meanstructure = TRUE, do.fit = FALSE)
#' fit <- run_fastsem(fit0)
#' summary(fit, fit.measures = TRUE)
#' }
run_fastsem <- function(model, control = list(),
                        check = getOption("fastsemR.check", "warn"),
                        engine = NULL) {
  if (!is.null(engine))
    warning("run_fastsem(engine = ) is deprecated and ignored.", call. = FALSE)
  check <- match.arg(check, c("warn", "error", "off"))
  if (!isTRUE(.fastsem_env$loaded))
    stop("fastsem native library is not loaded. ",
         "Run fastsem_install() or fastsem_load().")
  .require_engine_caps(.bridge_caps)
  .check_table_schema()
  if (methods::is(model, "MxModel")) {
    pt <- .ptab_from_mx(model)
    fs <- .ptab_fit(pt, control)
    return(.mx_from_fastsem(model, pt, fs, check))
  }
  if (methods::is(model, "lavaan")) {
    pt <- .ptab_from_lavaan(model)
    fs <- .ptab_fit(pt, control)
    return(.lav_from_fastsem(model, pt, fs, check))
  }
  stop("run_fastsem(): 'model' must be an OpenMx MxModel or a lavaan object, ",
       "not ", class(model)[1], ".", call. = FALSE)
}

.bridge_caps <- c("likelihood", "intercepts", "mg-instances", "bounds",
                  "fixed-means", "table")

# The bridge relies on engine fixes that older builds lack (labelled
# intercepts were silently dropped, multi-group labels ignored); refuse to run
# rather than return a wrong fit.
.require_engine_caps <- function(caps) {
  have <- tryCatch({
    sym <- getNativeSymbolInfo("engineCapsR", .fastsem_env$dll_handle)
    .Call(sym)
  }, error = function(e) character(0))
  miss <- setdiff(caps, have)
  if (length(miss))
    stop("run_fastsem() needs a newer fastsem engine (missing: ",
         paste(miss, collapse = ", "), "). Run fastsem_update(), or ",
         "fastsem_load() a current build.", call. = FALSE)
  invisible(TRUE)
}

# Send the table to the engine, fit, and attach the per-parameter results
# (the engine returns them per table row; rows sharing a plabel agree).
.ptab_fit <- function(pt, control) {
  a   <- .ptab_engine_args(pt, control)
  sym <- getNativeSymbolInfo("fitSemTableR", .fastsem_env$dll_handle)
  fs  <- .Call(sym, a$rows, a$cons, a$latent, a$manifest, a$ngroups,
               a$directives, a$data, colnames(a$data))
  if (is.null(fs))
    stop("run_fastsem(): the engine returned no result", call. = FALSE)
  if (!is.null(fs$error)) stop(fs$error, call. = FALSE)
  nf <- length(fs$freeNames)
  if (nf > 0 && length(fs$vcov) == nf * nf)
    fs$vcov <- matrix(fs$vcov, nf, nf, byrow = TRUE,
                      dimnames = list(fs$freeNames, fs$freeNames))

  first <- which(pt$free & !duplicated(pt$plabel))
  est <- setNames(fs$rowEstimates[first], pt$plabel[first])
  se  <- setNames(fs$rowSEs[first], pt$plabel[first])
  if (anyNA(est))
    stop("run_fastsem(): the engine did not return parameter(s) ",
         paste(names(est)[is.na(est)], collapse = ", "),
         ". Generated syntax:\n", fs$syntax, call. = FALSE)
  fs$fs_est <- est
  fs$fs_se  <- se
  fs
}

.self_check <- function(check, what, pkg_m2ll, fs, pt, tol = NULL) {
  if (check == "off") return(invisible())
  ordinal <- any(pt$mat == "T")
  tol <- tol %||% (if (ordinal) 0.5 else 1e-3 + 1e-8 * abs(fs$neg2logL))
  msgs <- character(0)
  d <- pkg_m2ll - fs$neg2logL
  if (is.na(pkg_m2ll)) d <- 0   # not comparable (see .lav_exo_m2ll)
  if (!is.finite(d) || abs(d) > tol)
    msgs <- c(msgs, sprintf(
      "%s's -2lnL at the fastsem estimates (%.6f) differs from fastsem's (%.6f) by %.3g%s",
      what, pkg_m2ll, fs$neg2logL, d,
      if (ordinal) " (ordinal: raise control = list(quadrature = 9))" else ""))
  msgs <- c(msgs, .ptab_check_bounds(pt, fs$fs_est))
  if (!length(msgs)) return(invisible())
  txt <- paste0("run_fastsem() self-check:\n  ", paste(msgs, collapse = "\n  "))
  if (check == "error") stop(txt, call. = FALSE) else warning(txt, call. = FALSE)
}

# ── OpenMx back end ───────────────────────────────────────────────────────────

.mx_from_fastsem <- function(model, pt, fs, check) {
  groups <- .mx_groups(model)
  is_mg  <- length(model@submodels) > 0L
  est    <- fs$fs_est
  devmap <- attr(pt, "devmap")
  for (g in seq_along(groups)) {
    sub <- groups[[g]]
    ex  <- sub$expectation
    dev_g <- if (is.null(devmap)) NULL else devmap[devmap$group == g, , drop = FALSE]
    rows <- pt[pt$group == g & pt$free, , drop = FALSE]
    for (r in seq_len(nrow(rows))) {
      v  <- est[[rows$plabel[r]]]
      rv <- rows$row[r]; cv <- rows$col[r]
      switch(rows$mat[r],
        A = { sub@matrices[[ex$A]]@values[rv, cv] <- v },
        S = { sub@matrices[[ex$S]]@values[rv, cv] <- v
              sub@matrices[[ex$S]]@values[cv, rv] <- v },
        M = { sub@matrices[[ex$M]]@values[1, rv] <- v },
        T = if (is.null(dev_g) || !nrow(dev_g))
              sub@matrices[[ex$thresholds]]@values[as.integer(sub("^t", "", rv)), cv] <- v)
    }
    # Deviation thresholds: d_1 = t_1, d_i = t_i - t_{i-1}, into the free cells.
    if (!is.null(dev_g) && nrow(dev_g)) {
      th <- pt[pt$group == g & pt$mat == "T", , drop = FALSE]
      for (v in unique(dev_g$var)) {
        tv <- th[th$col == v, , drop = FALSE]
        tv <- tv[order(as.integer(sub("^t", "", tv$row))), , drop = FALSE]
        t_hat <- ifelse(tv$free, est[tv$plabel], tv$value)
        d_hat <- c(t_hat[1], diff(t_hat))
        dv <- dev_g[dev_g$var == v & dev_g$free, , drop = FALSE]
        Ynm <- dv$Y[1]
        for (r in seq_len(nrow(dv)))
          sub@matrices[[Ynm]]@values[dv$i[r], v] <- d_hat[dv$i[r]]
      }
    }
    if (is_mg) model@submodels[[sub@name]] <- sub else model <- sub
  }
  if (is_mg && is.null(model@fitfunction))
    model <- OpenMx::mxModel(model,
               OpenMx::mxFitFunctionMultigroup(names(model@submodels)))

  # Let OpenMx evaluate its own model at the fastsem solution: fit, algebras
  # and expected moments come out exactly as after mxRun().
  model <- OpenMx::mxRun(model, useOptimizer = FALSE, silent = TRUE,
                         suppressWarnings = TRUE)

  pn <- names(OpenMx::omxGetParameters(model))
  V  <- .mx_param_vcov(pn, pt, fs)
  se <- if (is.null(V)) rep(NA_real_, length(pn)) else sqrt(pmax(diag(V), 0))
  model@output$standardErrors <- matrix(se, ncol = 1L, dimnames = list(pn, "SE"))
  if (!is.null(V) && !anyNA(V)) {
    model@output$vcov <- V
    H <- tryCatch(2 * solve(V), error = function(e) NULL)
    if (!is.null(H)) model@output$hessian <- H
  }
  model@output$fastsem <- fs
  .self_check(check, "OpenMx", model@output$fit, fs, pt)
  model
}

# Covariance of OpenMx's free parameters from the engine's covariance over
# (structural parameters, thresholds) — `vcovAug` for ordinal models, `vcov`
# otherwise — through V = G Va Gᵀ.  A structural or plain-threshold parameter
# is one column; a umx deviation d_i = t_i − t_{i−1} is +1/−1 on two threshold
# columns.  Rows that cannot be mapped are NA.
.mx_param_vcov <- function(pn, pt, fs) {
  an <- fs$vcovAugNames
  Va <- if (length(an) && length(fs$vcovAug) == length(an)^2)
          matrix(fs$vcovAug, length(an), length(an), byrow = TRUE,
                 dimnames = list(an, an))
        else if (is.matrix(fs$vcov)) fs$vcov
  if (is.null(Va)) return(NULL)
  vmap   <- setNames(fs$varTo, fs$varFrom)
  devmap <- attr(pt, "devmap")
  G <- matrix(0, length(pn), ncol(Va), dimnames = list(pn, colnames(Va)))
  ok <- logical(length(pn))
  set1 <- function(p, nm, w) {
    if (!nm %in% colnames(Va)) return(FALSE)
    G[p, nm] <<- G[p, nm] + w; TRUE
  }
  for (p in pn) {
    r <- which(pt$free & pt$omx == p)
    if (length(r)) {
      r <- r[1]
      ok[match(p, pn)] <- if (pt$mat[r] == "T")
        set1(p, paste0(vmap[[pt$col[r]]], "|", pt$row[r]), 1)
      else set1(p, fs$rowPlabels[r], 1)
      next
    }
    d <- if (is.null(devmap)) integer(0) else which(devmap$free & devmap$omx == p)
    if (length(d)) {
      d  <- devmap[d[1], ]
      vt <- vmap[[d$var]]
      a  <- set1(p, paste0(vt, "|t", d$i), 1)
      b  <- d$i == 1 || set1(p, paste0(vt, "|t", d$i - 1), -1)
      ok[match(p, pn)] <- a && b
    }
  }
  V <- G %*% Va %*% t(G)
  V[!ok, ] <- NA; V[, !ok] <- NA
  V
}

# ── lavaan back end ───────────────────────────────────────────────────────────

.lav_from_fastsem <- function(obj, pt, fs, check) {
  P   <- lavaan::parTable(obj)
  est <- fs$fs_est
  free_rows <- pt[pt$free & !is.na(pt$lav_id), , drop = FALSE]
  i <- match(free_rows$lav_id, P$id)
  P$start[i] <- est[free_rows$plabel]
  P$est[i]   <- est[free_rows$plabel]
  # Let lavaan evaluate its own model at the fastsem solution (no optimiser):
  # test statistics, SEs (incl. robust), implied moments come out natively.
  fit <- .lav_refit(obj, P)
  fit@external$fastsem <- fs
  # With fixed.x lavaan reports the likelihood conditional on the exogenous x;
  # fastsem fitted their moments freely (same estimates), so add the x
  # marginal at its MLE before comparing.
  m2ll <- -2 * as.numeric(lavaan::logLik(fit)) + .lav_exo_m2ll(obj, pt)
  .self_check(check, "lavaan", m2ll, fs, pt)
  fit
}

# Re-run lavaan on the same parameter table, data and options with
# `optim.method = "none"` at the given estimates.  Built from the object, not
# from its call: update() re-evaluates the call in the wrong frame (data held
# in a local variable are not found) and re-derives fixed.x; slot-based reuse
# trips over the `se = "none"` that do.fit = FALSE stores.
.lav_refit <- function(obj, P) {
  opt <- obj@Options
  ng  <- lavaan::lavInspect(obj, "ngroups")
  X   <- lavaan::lavInspect(obj, "data")
  if (ng == 1L) X <- list(X)
  glab <- lavaan::lavInspect(obj, "group.label")
  df <- do.call(rbind, lapply(seq_len(ng), function(g) {
    d <- as.data.frame(X[[g]])
    if (ng > 1L) d$.fastsem_group <- glab[g]
    d
  }))
  # do.fit = FALSE stores se = "none" (and test = "none" for plain ML), and
  # passing a resolved test list back changes lavaan's information settings;
  # take the user's literal choice from the call, else the estimator default.
  lit  <- function(x) if (is.character(x)) x else "default"
  se   <- lit(obj@call$se)
  test <- lit(obj@call$test)
  args <- list(model = as.data.frame(obj@ParTable), data = df, start = P,
               optim.method = "none",
               estimator = opt$estimator.orig %||% opt$estimator,
               se = se, test = test, missing = opt$missing,
               fixed.x = opt$fixed.x, conditional.x = opt$conditional.x,
               likelihood = opt$likelihood)
  if (ng > 1L) { args$group <- ".fastsem_group"; args$group.label <- glab }
  do.call(lavaan::lavaan, args)
}

# −2lnL of the exogenous observed covariates at their sample moments (N
# denominator), summed over groups; 0 without fixed.x, NA when x has missing
# values (the marginal is then not closed-form and the check is skipped).
.lav_exo_m2ll <- function(obj, pt) {
  if (!isTRUE(obj@Options$fixed.x)) return(0)
  data <- attr(pt, "data")
  total <- 0
  for (g in seq_along(data)) {
    xv <- lavaan::lavNames(obj, "ov.x", block = g)
    if (!length(xv)) next
    X <- as.matrix(data[[g]][, xv, drop = FALSE])
    if (anyNA(X)) return(NA_real_)
    n <- nrow(X); k <- ncol(X)
    S <- stats::cov(X) * (n - 1) / n
    total <- total + n * (k * log(2 * pi) +
                            as.numeric(determinant(S, logarithm = TRUE)$modulus) + k)
  }
  total
}
