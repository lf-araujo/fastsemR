## serialize.R — parameter table → the engine's table input.
##
## The table travels as typed columns (fitSemTableR); the engine lowers it to
## its own syntax with internal names and round-trip numbers, fits it, and
## returns estimates and SEs per row, so nothing here writes syntax.  The
## column set is the engine's row type (`ParRow` / `ConRow` in partable.nim);
## .table_schema() is checked against the engine's tableSchemaR().

.table_schema <- function() list(
  rows = c(group = "integer", mat = "character", row = "character",
           col = "character", free = "logical", value = "double",
           label = "character", defvar = "character", lbound = "double",
           ubound = "double"),
  constraints = c(lhs = "character", op = "character", rhs = "character"))

.check_table_schema <- function() {
  sym <- getNativeSymbolInfo("tableSchemaR", .fastsem_env$dll_handle)
  eng <- .Call(sym)
  if (!identical(eng, .table_schema()))
    stop("run_fastsem(): the loaded fastsem engine expects a different ",
         "parameter-table layout than this fastsemR version sends. Run ",
         "fastsem_update() or update fastsemR.", call. = FALSE)
  invisible(TRUE)
}

# The arguments of fitSemTableR / tableSyntaxR for a parameter table.
.ptab_engine_args <- function(pt, control = list()) {
  chr <- function(x) ifelse(is.na(x), "", as.character(x))
  rows <- list(
    group  = as.integer(pt$group), mat = as.character(pt$mat),
    row    = as.character(pt$row), col = as.character(pt$col),
    free   = as.logical(pt$free),  value = as.numeric(pt$value),
    label  = chr(ifelse(pt$free, pt$key, NA)), defvar = chr(pt$defvar),
    lbound = as.numeric(pt$lbound), ubound = as.numeric(pt$ubound))
  cons <- attr(pt, "constraints")
  cons <- list(lhs = as.character(cons$lhs %||% character(0)),
               op  = as.character(cons$op  %||% character(0)),
               rhs = as.character(cons$rhs %||% character(0)))
  manifest <- unique(unlist(attr(pt, "manifest")))
  directives <- c(paste0("likelihood: ", attr(pt, "likelihood") %||% "normal"),
                  .fastsem_control_directives(control))
  list(rows = rows, cons = cons,
       latent = as.character(unique(unlist(attr(pt, "latent")))),
       manifest = as.character(manifest),
       ngroups = length(attr(pt, "data")),
       directives = paste(directives, collapse = "\n"),
       data = .ptab_engine_data(pt, manifest))
}

# Stacked numeric data under the model's own column names: manifest
# variables (ordered factors as 0-based codes), definition columns, and the
# group column `fastsem_grp` (1..n in data order) for multi-group models.
.ptab_engine_data <- function(pt, manifest) {
  data <- attr(pt, "data")
  dvs  <- unique(pt$defvar[!is.na(pt$defvar)])
  blocks <- lapply(seq_along(data), function(g) {
    d    <- data[[g]]
    miss <- setdiff(c(manifest, dvs), names(d))
    if (length(miss))
      stop("run_fastsem(): column(s) ", paste(miss, collapse = ", "),
           " not found in the data of group ", g, ".", call. = FALSE)
    out <- lapply(c(manifest, dvs), function(v) {
      x <- d[[v]]
      if (is.factor(x)) {
        if (!is.ordered(x))
          stop("run_fastsem(): variable '", v, "' is an unordered factor.",
               call. = FALSE)
        as.numeric(as.integer(x) - 1L)
      } else as.numeric(x)
    })
    m <- do.call(cbind, out)
    colnames(m) <- c(manifest, dvs)
    if (length(data) > 1L) m <- cbind(m, fastsem_grp = g)
    m
  })
  do.call(rbind, blocks)
}

# The fastsem syntax the engine lowers a table to (internal names), for
# inspection: list(syntax, var_map = model name -> internal, plabel per row).
.ptab_syntax <- function(pt) {
  a   <- .ptab_engine_args(pt)
  sym <- getNativeSymbolInfo("tableSyntaxR", .fastsem_env$dll_handle)
  out <- .Call(sym, a$rows, a$cons, a$latent, a$manifest, a$ngroups, a$directives)
  if (!is.null(out$error)) stop(out$error, call. = FALSE)
  list(syntax = out$syntax, var_map = setNames(out$varTo, out$varFrom),
       plabel = out$plabel)
}

# Free parameters whose estimate falls outside the bound declared in OpenMx /
# lavaan.  The engine enforces the bounds it is sent, so this is a safety net,
# plus the only check for inequalities between two parameters.
.ptab_check_bounds <- function(pt, est) {
  fr <- pt[pt$free & (!is.na(pt$lbound) | !is.na(pt$ubound)), , drop = FALSE]
  msgs <- character(0)
  if (nrow(fr)) {
    v   <- est[fr$plabel]
    out <- (!is.na(fr$lbound) & v < fr$lbound - 1e-8) |
           (!is.na(fr$ubound) & v > fr$ubound + 1e-8)
    out[is.na(out)] <- FALSE
    if (any(out)) {
      fr <- fr[out, , drop = FALSE]
      msgs <- unique(sprintf("%s = %.4g outside [%s, %s]",
        ifelse(is.na(fr$omx), ifelse(is.na(fr$label), fr$key, fr$label), fr$omx),
        v[out], ifelse(is.na(fr$lbound), "-Inf", format(fr$lbound)),
        ifelse(is.na(fr$ubound), "Inf", format(fr$ubound))))
    }
  }
  # Inequality constraints between labels / against values (mxConstraint,
  # lavaan `<` / `>` lines).
  cons <- attr(pt, "constraints")
  if (!is.null(cons) && nrow(cons)) {
    lab_val <- function(x) {
      n <- suppressWarnings(as.numeric(x))
      if (!is.na(n)) return(n)
      i <- which(pt$label == x | pt$key == x)
      if (length(i)) est[[pt$plabel[i[1]]]] else NA_real_
    }
    for (i in which(cons$op %in% c("<", ">"))) {
      a <- lab_val(cons$lhs[i]); b <- lab_val(cons$rhs[i])
      ok <- if (cons$op[i] == "<") a < b + 1e-8 else a > b - 1e-8
      if (isFALSE(ok))
        msgs <- c(msgs, sprintf("%s %s %s violated (%.4g vs %.4g)", cons$lhs[i],
                                cons$op[i], cons$rhs[i], a, b))
    }
  }
  msgs
}
