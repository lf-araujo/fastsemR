## test-umx-bridge.R — run_fastsem() on OpenMx / umx RAM models.
##
## Each model is fitted by OpenMx and by fastsem; the fastsem MxModel must
## carry the same -2lnL, estimates and standard errors, and work with the
## OpenMx/umx tools downstream.

skip_if_no_openmx <- function() {
  skip_if_not_installed("OpenMx")
  skip_if_not_installed("umx")
}

expect_same_fit <- function(omx, fs, tol_est = 1e-4, tol_se = 1e-4) {
  expect_s4_class(fs, "MxModel")
  expect_equal(fs$output$fit, omx$output$fit, tolerance = 1e-6)
  p_o <- OpenMx::omxGetParameters(omx)
  p_f <- OpenMx::omxGetParameters(fs)
  expect_setequal(names(p_f), names(p_o))
  expect_equal(p_f[names(p_o)], p_o, tolerance = tol_est)
  se_o <- omx$output$standardErrors[, 1]
  se_f <- fs$output$standardErrors[names(se_o), 1]
  expect_equal(se_f, se_o, tolerance = tol_se)
}

hs3 <- function() lavaan::HolzingerSwineford1939
v3  <- c("x1", "x2", "x3")

cfa_mx <- function(name, d, var_labels = NA) {
  OpenMx::mxModel(name, type = "RAM", manifestVars = v3, latentVars = "g",
    OpenMx::mxPath("g", v3, values = .5, labels = c("l1", "l2", "l3")),
    OpenMx::mxPath("g", arrows = 2, free = FALSE, values = 1),
    OpenMx::mxPath(v3, arrows = 2, values = 1, labels = var_labels),
    OpenMx::mxPath("one", v3),
    OpenMx::mxData(d[, v3], type = "raw"))
}

test_that("single-group umx CFA matches OpenMx and works downstream", {
  skip_if_no_openmx(); skip_if_no_bridge_engine(); skip_if_not_installed("lavaan")
  hs <- hs3()
  m <- umx::umxRAM("cfa",
    umx::umxPath(from = "g", to = v3), umx::umxPath(v1m0 = "g"),
    umx::umxPath(var = v3), umx::umxPath(means = v3),
    data = hs[, v3], autoRun = FALSE)
  omx <- OpenMx::mxRun(m, silent = TRUE)
  fs  <- run_fastsem(m, check = "error")
  expect_same_fit(omx, fs)
  expect_true(fs@.wasRun)
  expect_s3_class(summary(fs), "summary.mxmodel")
  expect_equal(dim(vcov(fs)), rep(length(OpenMx::omxGetParameters(fs)), 2))
  expect_silent(OpenMx::mxStandardizeRAMpaths(fs, SE = FALSE))
  expect_false(is.null(fs$output$fastsem$neg2logL))
})

test_that("labelled, shared means fit like OpenMx (engine used to drop them)", {
  skip_if_no_openmx(); skip_if_no_bridge_engine()
  set.seed(11); n <- 300; f <- rnorm(n)
  d <- data.frame(x1 = 2 + .8 * f + rnorm(n, sd = .6), x2 = 2 + .7 * f + rnorm(n, sd = .6),
                  x3 = 5 + .6 * f + rnorm(n, sd = .6))
  m <- cfa_mx("mlab", d)
  m$M$labels[1, c("x1", "x2")] <- "mu"
  expect_same_fit(OpenMx::mxRun(m, silent = TRUE), run_fastsem(m, check = "error"))
})

test_that("multi-group: shared labels equal, unlabelled group-specific", {
  skip_if_no_openmx(); skip_if_no_bridge_engine(); skip_if_not_installed("lavaan")
  hs <- hs3()
  mg <- OpenMx::mxModel("mg",
    cfa_mx("pasteur", hs[hs$school == "Pasteur", ]),
    cfa_mx("grant",   hs[hs$school == "Grant-White", ]),
    OpenMx::mxFitFunctionMultigroup(c("pasteur", "grant")))
  omx <- OpenMx::mxRun(mg, silent = TRUE)
  fs  <- run_fastsem(mg, check = "error")
  expect_equal(length(OpenMx::omxGetParameters(fs)), 15L)
  expect_same_fit(omx, fs)
})

test_that("umxSuperModel with identical labels equates everything", {
  skip_if_no_openmx(); skip_if_no_bridge_engine(); skip_if_not_installed("lavaan")
  hs <- hs3()
  grp <- function(nm, d) umx::umxRAM(nm,
    umx::umxPath(from = "g", to = v3), umx::umxPath(v1m0 = "g"),
    umx::umxPath(var = v3), umx::umxPath(means = v3),
    data = d[, v3], autoRun = FALSE)
  mg <- umx::umxSuperModel("HS", grp("pasteur", hs[hs$school == "Pasteur", ]),
                                 grp("grant", hs[hs$school == "Grant-White", ]),
                           autoRun = FALSE)
  expect_same_fit(OpenMx::mxRun(mg, silent = TRUE), run_fastsem(mg, check = "error"))
})

test_that("a container without a fit function gets a multigroup one", {
  skip_if_no_openmx(); skip_if_no_bridge_engine(); skip_if_not_installed("lavaan")
  hs <- hs3()
  mg <- OpenMx::mxModel("mg2", cfa_mx("p", hs[hs$school == "Pasteur", ]),
                               cfa_mx("g2", hs[hs$school == "Grant-White", ]))
  fs <- run_fastsem(mg, check = "error")
  expect_true(is.finite(fs$output$fit))
})

test_that("FIML with missing data and definition variables match OpenMx", {
  skip_if_no_openmx(); skip_if_no_bridge_engine(); skip_if_not_installed("lavaan")
  hm <- hs3(); hm$x1[1:40] <- NA
  m <- cfa_mx("miss", hm)
  expect_same_fit(OpenMx::mxRun(m, silent = TRUE), run_fastsem(m, check = "error"))

  set.seed(3); n <- 400; age <- runif(n, -1, 1); f <- rnorm(n); h <- rnorm(n, sd = .7)
  dd <- data.frame(y1 = .8 * f + age * h + rnorm(n, sd = .6),
                   y2 = .7 * f + rnorm(n, sd = .6), y3 = .6 * f + rnorm(n, sd = .6),
                   age = age)
  y <- c("y1", "y2", "y3")
  dm <- OpenMx::mxModel("dv", type = "RAM", manifestVars = y, latentVars = c("f", "h"),
    OpenMx::mxPath("f", y, values = .5, labels = c("b1", "b2", "b3")),
    OpenMx::mxPath("h", "y1", free = FALSE, labels = "data.age"),
    OpenMx::mxPath("f", arrows = 2, free = FALSE, values = 1),
    OpenMx::mxPath("h", arrows = 2, values = .5, labels = "vh"),
    OpenMx::mxPath(y, arrows = 2, values = 1), OpenMx::mxPath("one", y),
    OpenMx::mxData(dd, type = "raw"))
  # def-var FIML SEs come from finite-difference Hessians on both sides
  expect_same_fit(OpenMx::mxRun(dm, silent = TRUE), run_fastsem(dm, check = "error"),
                  tol_se = 1e-3)
})

test_that("mxConstraint equality fits like a shared label", {
  skip_if_no_openmx(); skip_if_no_bridge_engine(); skip_if_not_installed("lavaan")
  mc <- OpenMx::mxModel(cfa_mx("con", hs3()), OpenMx::mxConstraint(l1 == l2, name = "eq"))
  fs <- run_fastsem(mc, check = "error")
  shared <- OpenMx::mxRun(OpenMx::omxSetParameters(cfa_mx("con", hs3()), labels = "l2",
                                                  newlabels = "l1"), silent = TRUE)
  expect_equal(fs$output$fit, shared$output$fit, tolerance = 1e-6)
  expect_equal(OpenMx::omxGetParameters(fs)[["l1"]],
               OpenMx::omxGetParameters(shared)[["l1"]], tolerance = 1e-4)
})

test_that("ordinal model with a plain threshold matrix matches OpenMx", {
  skip_if_no_openmx(); skip_if_no_bridge_engine()
  set.seed(7); n <- 800; f <- rnorm(n)
  d <- data.frame(y1 = .8 * f + rnorm(n), y2 = .6 * f + rnorm(n), y3 = .5 * f + rnorm(n))
  d[] <- lapply(d, function(x) OpenMx::mxFactor(as.integer(x > 0.2), levels = 0:1))
  y <- names(d)
  m <- OpenMx::mxModel("ordm", type = "RAM", manifestVars = y, latentVars = "g",
    OpenMx::mxPath("g", y, values = .5, labels = c("l1", "l2", "l3")),
    OpenMx::mxPath("g", arrows = 2, free = FALSE, values = 1),
    OpenMx::mxPath(y, arrows = 2, free = FALSE, values = 1),
    OpenMx::mxPath("one", y, free = FALSE, values = 0),
    OpenMx::mxThreshold(y, nThresh = 1, free = TRUE, values = .1,
                        labels = c("t_y1", "t_y2", "t_y3")),
    OpenMx::mxData(d, type = "raw"))
  omx <- OpenMx::mxRun(m, silent = TRUE)
  fs  <- run_fastsem(m, control = list(quadrature = 9), check = "error")
  expect_equal(fs$output$fit, omx$output$fit, tolerance = 1e-5)
  expect_equal(OpenMx::omxGetParameters(fs), OpenMx::omxGetParameters(omx), tolerance = 1e-2)
})

test_that("bounds bind exactly as in OpenMx (Heywood case)", {
  skip_if_no_openmx(); skip_if_no_bridge_engine()
  df <- as.data.frame(scale(iris[, c(1, 3, 4)])); names(df) <- c("sl", "pl", "pw")
  m <- umx::umxRAM("hey", umx::umxPath(from = "g", to = names(df)),
                   umx::umxPath(v1m0 = "g"), umx::umxPath(var = names(df)),
                   data = df, autoRun = FALSE)
  omx <- OpenMx::mxRun(m, silent = TRUE)
  fs  <- run_fastsem(m, check = "error")
  expect_equal(fs$output$fit, omx$output$fit, tolerance = 1e-6)
  expect_identical(OpenMx::omxGetParameters(fs)[["pl_with_pl"]], 0)
  p_o <- OpenMx::omxGetParameters(omx)
  expect_equal(OpenMx::omxGetParameters(fs)[names(p_o)], p_o, tolerance = 1e-4)
})

test_that("unsupported OpenMx constructs fail loudly", {
  skip_if_no_openmx(); skip_if_no_bridge_engine(); skip_if_not_installed("lavaan")
  set.seed(1)
  sim <- function(rA, n) {
    A1 <- rnorm(n); A2 <- rA * A1 + sqrt(1 - rA^2) * rnorm(n)
    data.frame(T_1 = A1 + rnorm(n), T_2 = A2 + rnorm(n))
  }
  ace <- umx::umxACE(selDVs = "T", sep = "_", dzData = sim(.5, 300),
                     mzData = sim(1, 300), autoRun = FALSE)
  expect_error(run_fastsem(ace), "not RAM models")
})

test_that("a mean fixed at a non-zero value matches OpenMx", {
  skip_if_no_openmx(); skip_if_no_bridge_engine(); skip_if_not_installed("lavaan")
  m <- cfa_mx("fixmean", hs3())
  m$M$free[1, "x1"] <- FALSE; m$M$values[1, "x1"] <- 5
  expect_same_fit(OpenMx::mxRun(m, silent = TRUE), run_fastsem(m, check = "error"))
})

test_that("fixed values keep full precision", {
  skip_if_no_openmx(); skip_if_no_bridge_engine(); skip_if_not_installed("lavaan")
  m <- cfa_mx("prec", hs3())
  m$S$values["g", "g"] <- 0.123456789012
  syn <- fastsemR:::.ptab_syntax(fastsemR:::.ptab_from_mx(m))$syntax
  expect_match(syn, "0.123456789012", fixed = TRUE)
})

test_that("deprecated entry points still work and warn", {
  skip_if_no_openmx(); skip_if_no_bridge_engine(); skip_if_not_installed("lavaan")
  m <- cfa_mx("dep", hs3())
  expect_warning(syn <- umx_to_lavaan(m), "deprecated")
  expect_match(syn, "=~")
  expect_warning(res <- umx_to_fastsem(m), "deprecated")
  expect_true(is.finite(res$neg2logL))
  expect_warning(fs <- run_fastsem_ram(m), "deprecated")
  expect_s4_class(fs, "MxModel")
})
