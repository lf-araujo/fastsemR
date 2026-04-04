## test-umx-bridge.R — umx_to_lavaan() and run_fastsem() tests

skip_if_no_openmx <- function() {
  skip_if_not_installed("OpenMx")
  skip_if_not_installed("umx")
}

skip_if_not_loaded <- function() {
  skip_if_not(isTRUE(fastsemR:::.fastsem_env$loaded),
               "fastsem native library not loaded")
}

# ── umx_to_lavaan: syntax generation ─────────────────────────────────────────

test_that("umx_to_lavaan produces non-empty syntax for a CFA model", {
  skip_if_no_openmx()

  df <- as.data.frame(scale(iris[, c(1, 3, 4)]))
  names(df) <- c("sl", "pl", "pw")

  m <- umx::umxRAM("cfa_test",
    umx::umxPath(from = "g", to = c("sl", "pl", "pw")),
    umx::umxPath(v1m0 = "g"),
    umx::umxPath(var  = c("sl", "pl", "pw")),
    data = df, autoRun = FALSE
  )

  syn <- umx_to_lavaan(m)
  expect_type(syn, "character")
  expect_true(nzchar(syn))
  expect_true(grepl("g =~", syn))
  expect_true(grepl("sl ~~ 1\\*sl|g ~~ 1\\*g", syn))
})

test_that("umx_to_lavaan include_map returns param_map", {
  skip_if_no_openmx()

  df <- as.data.frame(scale(iris[, c(1, 3, 4)]))
  names(df) <- c("sl", "pl", "pw")

  m <- umx::umxRAM("cfa_map",
    umx::umxPath(from = "g", to = c("sl", "pl", "pw")),
    umx::umxPath(v1m0 = "g"),
    umx::umxPath(var  = c("sl", "pl", "pw")),
    data = df, autoRun = FALSE
  )

  out <- umx_to_lavaan(m, include_map = TRUE)
  expect_type(out, "list")
  expect_true("syntax"    %in% names(out))
  expect_true("param_map" %in% names(out))
})

# ── run_fastsem: single-group CFA ─────────────────────────────────────────────

test_that("run_fastsem CFA returns an MxModel with non-zero loadings", {
  skip_if_no_openmx()
  skip_if_not_loaded()

  df <- as.data.frame(scale(iris[, c(1, 3, 4)]))
  names(df) <- c("sl", "pl", "pw")

  m <- umx::umxRAM("cfa_run",
    umx::umxPath(from = "g", to = c("sl", "pl", "pw")),
    umx::umxPath(v1m0 = "g"),
    umx::umxPath(var  = c("sl", "pl", "pw")),
    data = df, autoRun = FALSE
  )

  m_fit <- run_fastsem(m)

  expect_s4_class(m_fit, "MxModel")
  expect_true(isTRUE(m_fit@.wasRun))

  pars <- OpenMx::omxGetParameters(m_fit)
  loadings <- pars[grep("g_to_", names(pars))]
  expect_true(all(loadings > 0.5))
})

# ── run_fastsem: path model matches OpenMx estimates ─────────────────────────

test_that("run_fastsem path model estimates match OpenMx to 2 decimal places", {
  skip_if_no_openmx()
  skip_if_not_loaded()

  df <- as.data.frame(scale(mtcars[, c("mpg", "hp", "wt")]))

  m_omx <- umx::umxRAM("path_omx",
    umx::umxPath(from = c("hp", "wt"), to = "mpg"),
    umx::umxPath(cov  = c("hp", "wt")),
    umx::umxPath(var  = c("hp", "wt", "mpg")),
    umx::umxPath(means = c("hp", "wt", "mpg")),
    data = df
  )

  m_fs <- umx::umxRAM("path_fs",
    umx::umxPath(from = c("hp", "wt"), to = "mpg"),
    umx::umxPath(cov  = c("hp", "wt")),
    umx::umxPath(var  = c("hp", "wt", "mpg")),
    umx::umxPath(means = c("hp", "wt", "mpg")),
    data = df, autoRun = FALSE
  )
  m_fs <- run_fastsem(m_fs)

  omx_est <- OpenMx::omxGetParameters(m_omx)
  fs_est  <- OpenMx::omxGetParameters(m_fs)

  shared <- intersect(names(omx_est), names(fs_est))
  expect_true(length(shared) >= 2)
  expect_true(all(abs(omx_est[shared] - fs_est[shared]) < 0.05))
})

# ── run_fastsem: multi-group CFA ──────────────────────────────────────────────

test_that("run_fastsem multi-group CFA returns reasonable loadings per group", {
  skip_if_no_openmx()
  skip_if_not_loaded()
  skip_if_not_installed("lavaan")

  hs <- lavaan::HolzingerSwineford1939
  mg_p <- umx::umxRAM("mg_p",
    umx::umxPath(from = "g", to = c("x1","x2","x3")),
    umx::umxPath(v1m0 = "g"),
    umx::umxPath(var  = c("x1","x2","x3")),
    umx::umxPath(means = c("x1","x2","x3")),
    data = hs[hs$school == "Pasteur",     c("x1","x2","x3")]
  )
  mg_gw <- umx::umxRAM("mg_gw",
    umx::umxPath(from = "g", to = c("x1","x2","x3")),
    umx::umxPath(v1m0 = "g"),
    umx::umxPath(var  = c("x1","x2","x3")),
    umx::umxPath(means = c("x1","x2","x3")),
    data = hs[hs$school == "Grant-White", c("x1","x2","x3")]
  )

  combined <- OpenMx::mxModel("hs_mg", mg_p, mg_gw)
  r        <- run_fastsem(combined)

  p_pars  <- OpenMx::omxGetParameters(r@submodels$mg_p)
  gw_pars <- OpenMx::omxGetParameters(r@submodels$mg_gw)

  # All loadings should be between 0.3 and 1.2
  for (nm in c("g_to_x1","g_to_x2","g_to_x3")) {
    expect_gt(p_pars[[nm]],  0.3)
    expect_lt(p_pars[[nm]],  1.2)
    expect_gt(gw_pars[[nm]], 0.3)
    expect_lt(gw_pars[[nm]], 1.2)
  }
})
