## partable.R — the internal parameter table shared by the OpenMx and lavaan
## front ends.
##
## One row per non-zero or free RAM cell per group:
##   id      row identifier
##   group   1-based group number
##   mat     "A" (directed), "S" (symmetric, upper triangle), "M" (means),
##           "T" (thresholds: row = threshold index, col = variable)
##   row,col variable names (A[row, col] is the path col -> row)
##   free    free parameter
##   value   fixed value, or start value when free
##   label   user label (OpenMx / lavaan), NA if none
##   defvar  definition-variable data column, NA if none
##   lbound, ubound   bounds (NA = none)
##   key     parameter identity: rows sharing a key are one parameter
##   plabel  safe internal label sent to the engine ("p<k>"), NA when fixed
##   omx     OpenMx free-parameter name (OpenMx front end only)
##
## Attributes carry what is not per-cell: "manifest", "latent" (per group),
## "ordered" (variables), "constraints" (data.frame lhs/op/rhs, labels already
## mapped to plabels), "data" (list of per-group data frames), "source".

.ptab_cols <- c("group", "mat", "row", "col", "free", "value", "label",
                "defvar", "lbound", "ubound", "omx")

.ptab_row <- function(group, mat, row, col, free, value, label = NA_character_,
                      defvar = NA_character_, lbound = NA_real_,
                      ubound = NA_real_, omx = NA_character_) {
  data.frame(group = group, mat = mat, row = row, col = col, free = free,
             value = value, label = label, defvar = defvar, lbound = lbound,
             ubound = ubound, omx = omx, stringsAsFactors = FALSE)
}

.lab_ok <- function(x) !is.na(x) & nzchar(x)

# ── OpenMx front end ──────────────────────────────────────────────────────────

# The RAM groups of an MxModel: the model itself, or its submodels when it is a
# multi-group container.  Anything that is not a RAM expectation is rejected
# here, so umx twin supermodels (umxACE, umxCP, ...: MxExpectationNormal over
# algebras in a 'top' submodel) fail with an explicit message.
.mx_groups <- function(model) {
  is_ram <- function(m) inherits(m$expectation, "MxExpectationRAM")
  if (length(model@submodels) == 0L) {
    if (!is_ram(model))
      stop("run_fastsem(): '", model@name, "' is not a RAM model (expectation: ",
           class(model$expectation)[1], "). fastsem fits RAM models ",
           "(umxRAM / mxModel(type = 'RAM')).", call. = FALSE)
    return(list(model))
  }
  subs <- model@submodels
  ram  <- vapply(subs, is_ram, logical(1))
  if (!all(ram))
    stop("run_fastsem(): submodel(s) ",
         paste(names(subs)[!ram], collapse = ", "), " of '", model@name,
         "' are not RAM models. fastsem fits RAM models; umx twin supermodels ",
         "(umxACE, umxCP, umxIP, ...) build their covariance from algebras ",
         "and are not supported.", call. = FALSE)
  unname(subs)
}

.mx_ptab_group <- function(sub, g, locs) {
  ex   <- sub$expectation
  A    <- sub[[ex$A]];  S <- sub[[ex$S]]
  M    <- if (!is.na(ex$M)) sub[[ex$M]] else NULL
  vars <- rownames(A@values)
  omx_name <- function(mat, i, j) {
    hit <- locs[locs$model == sub@name & locs$matrix == mat &
                locs$row == i & locs$col == j, , drop = FALSE]
    if (nrow(hit) == 0L) return(NA_character_)
    if (.lab_ok(hit$label[1])) hit$label[1] else
      sprintf("%s.%s[%d,%d]", sub@name, mat, i, j)
  }
  rows <- list()
  devs <- list()
  add  <- function(...) rows[[length(rows) + 1L]] <<- .ptab_row(...)
  cell <- function(mobj, mname, i, j, matcode, rowv, colv) {
    fr  <- isTRUE(mobj@free[i, j])
    val <- mobj@values[i, j]
    lab <- mobj@labels[i, j]
    if (.lab_ok(lab) && grepl("[", lab, fixed = TRUE))
      stop("run_fastsem(): ", sub@name, ".", mname, "[", i, ",", j,
           "] is bound to '", lab, "' (an algebra or matrix element); fastsem ",
           "cannot evaluate OpenMx algebras inside a RAM model.", call. = FALSE)
    dv <- NA_character_
    if (.lab_ok(lab) && startsWith(lab, "data.")) {
      dv <- substring(lab, 6L); fr <- FALSE; lab <- NA_character_
    }
    if (!fr && is.na(dv) && (is.na(val) || val == 0)) return(invisible())
    lb <- if (fr && !is.null(mobj@lbound) && length(mobj@lbound)) mobj@lbound[i, j] else NA_real_
    ub <- if (fr && !is.null(mobj@ubound) && length(mobj@ubound)) mobj@ubound[i, j] else NA_real_
    add(g, matcode, rowv, colv, fr, if (is.na(val)) 0 else val,
        if (fr) lab else NA_character_, dv, lb, ub,
        if (fr) omx_name(mname, i, j) else NA_character_)
  }
  for (j in seq_along(vars)) for (i in seq_along(vars)) {
    if (i == j) next
    cell(A, ex$A, i, j, "A", vars[i], vars[j])
  }
  for (i in seq_along(vars)) for (j in i:length(vars)) {
    # S is symmetric; OpenMx keeps free/labels on both triangles.
    cell(S, ex$S, i, j, "S", vars[i], vars[j])
  }
  if (!is.null(M)) {
    mv <- colnames(M@values)
    for (j in seq_along(mv)) cell(M, ex$M, 1L, j, "M", mv[j], mv[j])
  }
  thr <- ex$thresholds
  if (!is.na(thr)) {
    ti  <- .mx_threshold_source(sub, thr)
    obs <- sub@data@observed
    for (j in seq_along(ti$vars)) {
      v    <- ti$vars[j]
      nlev <- length(levels(obs[[v]])) - 1L
      if (ti$type == "matrix") {
        Tm <- ti$T
        for (i in seq_len(nlev)) {
          fr  <- isTRUE(Tm@free[i, j])
          add(g, "T", paste0("t", i), v, fr, Tm@values[i, j],
              if (fr) Tm@labels[i, j] else NA_character_, NA_character_,
              NA_real_, NA_real_, if (fr) omx_name(thr, i, j) else NA_character_)
        }
      } else {
        # Deviation form (umxThresholdMatrix): t_i = d_1 + ... + d_i.  t_i is
        # free when any of its deviations is; its identity is the sequence of
        # deviation identities (label, cell, or fixed value), so thresholds
        # built from shared deviation labels (twins, groups) stay one parameter.
        Y <- ti$Y
        devkey <- vapply(seq_len(nlev), function(r) {
          if (isTRUE(Y@free[r, j])) {
            if (.lab_ok(Y@labels[r, j])) Y@labels[r, j]
            else sprintf("@%d%s[%d,%d]", g, ti$Yname, r, j)
          } else paste0("=", formatC(Y@values[r, j], digits = 17, format = "g"))
        }, character(1))
        for (i in seq_len(nlev)) {
          fr <- any(Y@free[seq_len(i), j])
          add(g, "T", paste0("t", i), v, fr, sum(Y@values[seq_len(i), j]),
              if (fr) paste0("dev:", paste(devkey[seq_len(i)], collapse = "+"))
              else NA_character_, NA_character_, NA_real_, NA_real_, NA_character_)
          devs[[length(devs) + 1L]] <- data.frame(
            group = g, Y = ti$Yname, var = v, i = i, free = isTRUE(Y@free[i, j]),
            omx = if (isTRUE(Y@free[i, j])) omx_name(ti$Yname, i, j) else NA_character_,
            stringsAsFactors = FALSE)
        }
      }
    }
  }
  out <- do.call(rbind, rows)
  attr(out, "devmap") <- if (length(devs)) do.call(rbind, devs) else NULL
  out
}

# Where a RAM expectation's thresholds come from: a plain matrix, or umx's
# deviation algebra `X %*% Y` (umxThresholdMatrix: X a fixed lower-triangular
# matrix of ones, Y the deviations).  Recognised by structure, not by name.
.mx_threshold_source <- function(sub, thr) {
  get <- function(nm) {
    o <- tryCatch(sub[[nm]], error = function(e) NULL)
    if (is.null(o)) o <- tryCatch(sub[[sub("^[^.]+\\.", "", nm)]], error = function(e) NULL)
    o
  }
  obj <- get(thr)
  if (methods::is(obj, "MxMatrix"))
    return(list(type = "matrix", T = obj, vars = colnames(obj@values)))
  if (methods::is(obj, "MxAlgebra")) {
    f <- obj@formula
    if (is.call(f) && identical(f[[1]], as.name("%*%")) &&
        is.name(f[[2]]) && is.name(f[[3]])) {
      X <- get(as.character(f[[2]])); Y <- get(as.character(f[[3]]))
      if (methods::is(X, "MxMatrix") && methods::is(Y, "MxMatrix") &&
          !any(X@free) && nrow(X@values) == ncol(X@values) &&
          nrow(X@values) == nrow(Y@values) &&
          all(X@values == 1 * lower.tri(X@values, diag = TRUE)))
        return(list(type = "dev", Y = Y, Yname = Y@name, vars = colnames(Y@values)))
    }
  }
  stop("run_fastsem(): thresholds of '", sub@name, "' come from '", thr,
       "', which is neither a threshold matrix nor umx's deviation form ",
       "(lowerOnes %*% deviations).", call. = FALSE)
}

.ptab_from_mx <- function(model) {
  if (!requireNamespace("OpenMx", quietly = TRUE))
    stop("The 'OpenMx' package is required for run_fastsem() on an MxModel.")
  groups <- .mx_groups(model)
  locs   <- OpenMx::omxLocateParameters(model)
  locs$model <- sub("\\..*$", "", locs$model)
  parts <- lapply(seq_along(groups), function(g) .mx_ptab_group(groups[[g]], g, locs))
  pt <- do.call(rbind, parts)
  devmap <- do.call(rbind, lapply(parts, attr, "devmap"))
  # Parameter identity: OpenMx equates every free cell carrying the same label,
  # across matrices and groups; an unlabelled free cell is its own parameter.
  pt$key <- ifelse(pt$free & .lab_ok(pt$label), pt$label,
                   ifelse(pt$free, paste0("@", pt$group, pt$mat, ":", pt$row,
                                          ":", pt$col), NA_character_))
  data <- lapply(groups, function(s) .mx_group_data(s))
  attr(pt, "manifest") <- lapply(groups, function(s) s@manifestVars)
  attr(pt, "latent")   <- lapply(groups, function(s) s@latentVars)
  attr(pt, "ordered")  <- unique(unlist(lapply(groups, function(s)
    names(Filter(is.ordered, s@data@observed)))))
  attr(pt, "constraints") <- .mx_constraints(model)
  attr(pt, "data")     <- data
  attr(pt, "source")   <- "OpenMx"
  attr(pt, "devmap")   <- devmap
  .ptab_check_ordinal_id(pt)
  .ptab_finish(pt)
}

.mx_group_data <- function(sub) {
  dat <- sub@data
  if (is.null(dat))
    stop("run_fastsem(): model '", sub@name, "' has no data.", call. = FALSE)
  if (!is.null(dat@type) && dat@type != "raw")
    stop("run_fastsem(): only raw (row-level) data is supported; '", sub@name,
         "' has type = '", dat@type, "'.", call. = FALSE)
  if (is.null(dat@observed))
    stop("run_fastsem(): '", sub@name, "'@data@observed is NULL.", call. = FALSE)
  as.data.frame(dat@observed)
}

# mxConstraints between labelled parameters (a == b, a < b, a > 0) map onto
# fastsem's label constraints; anything else cannot be expressed.
.mx_constraints <- function(model) {
  cons <- c(model@constraints,
            unlist(lapply(model@submodels, function(s) s@constraints)))
  if (length(cons) == 0L)
    return(data.frame(lhs = character(0), op = character(0), rhs = character(0)))
  do.call(rbind, lapply(cons, function(cn) {
    lhs <- deparse(cn@formula[[2]]); rhs <- deparse(cn@formula[[3]])
    op  <- as.character(cn@formula[[1]])
    if (!grepl("^[A-Za-z_][A-Za-z0-9_.]*$", lhs) ||
        !(grepl("^[A-Za-z_][A-Za-z0-9_.]*$", rhs) ||
          !is.na(suppressWarnings(as.numeric(rhs)))))
      stop("run_fastsem(): mxConstraint '", cn@name, "' (", lhs, " ", op, " ",
           rhs, ") is not a relation between parameter labels; fastsem ",
           "supports only label == label, label < value and label > value.",
           call. = FALSE)
    data.frame(lhs = lhs, op = op, rhs = rhs, stringsAsFactors = FALSE)
  }))
}

# ── lavaan front end ──────────────────────────────────────────────────────────

.ptab_from_lavaan <- function(obj) {
  if (!requireNamespace("lavaan", quietly = TRUE))
    stop("The 'lavaan' package is required for run_fastsem() on a lavaan object.")
  opt <- obj@Options
  if (length(lavaan::lavNames(obj, "ov.ord")))
    stop("run_fastsem(): lavaan models with ordered indicators are not supported ",
         "yet (lavaan's DWLS/WLSMV and fastsem's ordinal FIML are different ",
         "estimators). Fit the model through OpenMx/umx instead.", call. = FALSE)
  if (!toupper(opt$estimator) %in% c("ML", "MLR", "MLM", "MLMV", "MLMVS", "MLF"))
    stop("run_fastsem(): lavaan estimator '", opt$estimator, "' is not ML-based; ",
         "fastsem returns ML/FIML point estimates.", call. = FALSE)
  clus <- tryCatch(lavaan::lavInspect(obj, "cluster"), error = function(e) NULL)
  if (length(clus) || isTRUE(nzchar(opt$sampling.weights)))
    stop("run_fastsem(): clustered or weighted lavaan models are not supported ",
         "through run_fastsem() yet.", call. = FALSE)

  P <- lavaan::parTable(obj)
  ng <- lavaan::lavInspect(obj, "ngroups")
  lv <- lavaan::lavNames(obj, "lv")
  rows <- list()
  for (r in seq_len(nrow(P))) {
    op <- P$op[r]
    if (!op %in% c("=~", "~", "~~", "~1", "|")) next
    g   <- if (!is.null(P$group)) P$group[r] else 1L
    fr  <- P$free[r] > 0L || (isTRUE(P$exo[r] == 1L))
    val <- if (!is.na(P$start[r])) P$start[r] else 0
    if (!fr && val == 0) next
    lab <- if (nzchar(P$label[r])) P$label[r] else NA_character_
    cellrc <- switch(op,
      "=~" = c("A", P$rhs[r], P$lhs[r]),
      "~"  = c("A", P$lhs[r], P$rhs[r]),
      "~~" = c("S", P$lhs[r], P$rhs[r]),
      "~1" = c("M", P$lhs[r], P$lhs[r]),
      "|"  = c("T", P$rhs[r], P$lhs[r]))
    rows[[length(rows) + 1L]] <- .ptab_row(
      g, cellrc[1], cellrc[2], cellrc[3], fr, val,
      if (fr) lab else NA_character_, NA_character_,
      if (!is.null(P$lower) && is.finite(P$lower[r])) P$lower[r] else NA_real_,
      if (!is.null(P$upper) && is.finite(P$upper[r])) P$upper[r] else NA_real_)
    rows[[length(rows)]]$lav_id <- P$id[r]
  }
  # Without a mean structure lavaan's likelihood uses the saturated (sample)
  # means; fastsem gets free intercepts on every observed variable to match.
  ov <- lavaan::lavNames(obj, "ov")
  if (!isTRUE(opt$meanstructure))
    for (g in seq_len(ng)) for (v in ov) {
      rows[[length(rows) + 1L]] <- .ptab_row(g, "M", v, v, TRUE, 0)
      rows[[length(rows)]]$lav_id <- NA_integer_
    }
  pt <- do.call(rbind, rows)
  # Identity: a user label is one parameter everywhere; otherwise lavaan's
  # plabel.  `==` rows between labels/plabels (group.equal and friends) are
  # folded in with union-find so equal parameters share one key.
  plab <- setNames(P$plabel, P$id)
  pt$key <- ifelse(pt$free, ifelse(.lab_ok(pt$label), pt$label,
                                   ifelse(is.na(pt$lav_id),
                                          paste0("@mean", pt$group, ":", pt$row),
                                          plab[as.character(pt$lav_id)])),
                   NA_character_)
  eqs  <- P[P$op == "==", c("lhs", "rhs")]
  key_of <- function(x) {   # a label or plabel → the key it names
    hit <- which(P$label == x | P$plabel == x)
    if (!length(hit))
      stop("run_fastsem(): lavaan equality constraint on '", x, "' is not ",
           "between parameters; fastsem supports label == label only.",
           call. = FALSE)
    i <- which(pt$lav_id %in% P$id[hit])
    if (length(i)) pt$key[i[1]] else x
  }
  if (nrow(eqs)) {
    parent <- list()
    find <- function(x) { while (!is.null(parent[[x]])) x <- parent[[x]]; x }
    for (e in seq_len(nrow(eqs))) {
      a <- find(key_of(eqs$lhs[e])); b <- find(key_of(eqs$rhs[e]))
      if (a != b) parent[[b]] <- a
    }
    pt$key <- vapply(pt$key, function(k) if (is.na(k)) k else find(k), "")
  }
  # Thresholds: lavaan names them t1, t2, ... per variable (rhs of `|`).
  data <- lapply(seq_len(ng), function(g) {
    X <- if (ng > 1L) lavaan::lavInspect(obj, "data")[[g]] else lavaan::lavInspect(obj, "data")
    as.data.frame(X)
  })
  attr(pt, "manifest") <- rep(list(ov), ng)
  attr(pt, "latent")   <- rep(list(lv), ng)
  attr(pt, "ordered")  <- character(0)
  attr(pt, "constraints") <- P[P$op %in% c(":=", "<", ">"), c("lhs", "op", "rhs")]
  attr(pt, "data")     <- data
  attr(pt, "source")   <- "lavaan"
  attr(pt, "likelihood") <- if (identical(opt$likelihood, "wishart")) "wishart" else "normal"
  .ptab_finish(pt)
}

# Ordinal FIML needs each ordinal variable's location and scale pinned: two of
# {mean fixed, residual variance fixed, a fixed threshold (each)} — theta
# (mean 0, residual 1), Mehta (two thresholds fixed), binary conventions.  The
# engine applies the same rule; checking here gives the message in OpenMx terms.
.ptab_check_ordinal_id <- function(pt) {
  th <- pt[pt$mat == "T", , drop = FALSE]
  if (!nrow(th)) return(invisible())
  for (g in unique(th$group)) for (v in unique(th$col[th$group == g])) {
    in_g <- pt$group == g
    # (parenthesised: `!` binds more loosely than `+` in R)
    pins <- (!any(in_g & pt$mat == "M" & pt$row == v & pt$free)) +
            (!any(in_g & pt$mat == "S" & pt$row == v & pt$col == v & pt$free)) +
            sum(th$group == g & th$col == v & !th$free)
    if (pins < 2)
      stop("run_fastsem(): ordinal variable '", v, "' (group ", g, ") is not ",
           "identified: fix two of its mean, residual variance and thresholds ",
           "(theta: mean 0 and residual variance 1; Mehta: fix two thresholds). ",
           "A total-variance constraint (delta parameterisation via ",
           "mxConstraint/algebra) is not translated.", call. = FALSE)
  }
  invisible()
}

# ── Shared finishing + validation ─────────────────────────────────────────────

.ptab_finish <- function(pt) {
  pt$id <- seq_len(nrow(pt))
  keys  <- unique(pt$key[pt$free])
  pt$plabel <- ifelse(pt$free, paste0("p", match(pt$key, keys)), NA_character_)
  .ptab_validate(pt)
  pt
}

.ptab_validate <- function(pt) {
  where <- function(i) sprintf("%s[%s, %s] (group %d)", pt$mat[i], pt$row[i],
                               pt$col[i], pt$group[i])
  bad <- which(!is.na(pt$defvar) & pt$mat != "A")
  if (length(bad))
    stop("run_fastsem(): definition variable '", pt$defvar[bad[1]], "' on ",
         where(bad[1]), "; fastsem supports definition variables on paths and ",
         "loadings (A) only.", call. = FALSE)
  latent_all <- unique(unlist(attr(pt, "latent")))
  bad <- which(!is.na(pt$defvar) & !(pt$col %in% latent_all))
  if (length(bad))
    stop("run_fastsem(): definition variable '", pt$defvar[bad[1]], "' on the ",
         "path ", pt$col[bad[1]], " -> ", pt$row[bad[1]], "; fastsem supports ",
         "definition variables on paths from latent variables only.", call. = FALSE)
  # Definition-variable cells must be the same in every group (the engine's
  # def-var mappings are not group-specific).
  dv <- pt[!is.na(pt$defvar), , drop = FALSE]
  if (nrow(dv) && length(unique(pt$group)) > 1L) {
    sig <- table(paste(dv$row, dv$col, dv$defvar))
    if (any(sig != length(unique(pt$group))))
      stop("run_fastsem(): definition variables differ between groups; fastsem ",
           "needs the same definition-variable paths in every group.", call. = FALSE)
  }
  invisible(pt)
}
