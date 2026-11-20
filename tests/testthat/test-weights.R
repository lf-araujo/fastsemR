## test-weights.R — row weights through run_fastsem().
##
## OpenMx mxData(weight =) and lavaan sampling.weights must reach the engine:
## both used to be dropped silently (the fit came back unweighted).

test_that("OpenMx row weights match mxRun(), single- and multi-group", {
  skip_if_not_installed("OpenMx"); skip_if_no_bridge_engine()
  set.seed(42)
  n <- 400; f <- rnorm(n)
  X <- sapply(c(.8, .7, .6, .75), function(l) l * f + rnorm(n, 0, .6))
  colnames(X) <- paste0("x", 1:4)
  X[matrix(runif(n * 4) < .1, n, 4)] <- NA
  D <- data.frame(X, wt = runif(n, .5, 2))
  mk <- function(name, dat) OpenMx::mxModel(name, type = "RAM",
    manifestVars = paste0("x", 1:4), latentVars = "F",
    OpenMx::mxPath("F", paste0("x", 1:4), values = .7, free = c(FALSE, TRUE, TRUE, TRUE),
                   labels = paste0("l", 1:4)),
    OpenMx::mxPath("F", arrows = 2, values = 1, labels = "vF"),
    OpenMx::mxPath(paste0("x", 1:4), arrows = 2, values = .5, labels = paste0("e", 1:4)),
    OpenMx::mxPath("one", paste0("x", 1:4), values = 0, labels = paste0("m", 1:4)), dat)
  m  <- mk("w", OpenMx::mxData(D, "raw", weight = "wt"))
  m0 <- OpenMx::mxRun(m, silent = TRUE)
  m1 <- run_fastsem(m, check = "error")
  expect_equal(m1$output$fit, m0$output$fit, tolerance = 1e-6)
  expect_equal(OpenMx::omxGetParameters(m1), OpenMx::omxGetParameters(m0), tolerance = 1e-4)
  # weights in one group only
  mg <- OpenMx::mxModel("mg",
    mk("g1", OpenMx::mxData(D[1:200, ], "raw", weight = "wt")),
    mk("g2", OpenMx::mxData(D[201:400, ], "raw")),
    OpenMx::mxFitFunctionMultigroup(c("g1", "g2")))
  g0 <- OpenMx::mxRun(mg, silent = TRUE)
  g1 <- run_fastsem(mg, check = "error")
  expect_equal(g1$output$fit, g0$output$fit, tolerance = 1e-6)
  expect_equal(OpenMx::omxGetParameters(g1), OpenMx::omxGetParameters(g0), tolerance = 1e-4)
})

test_that("OpenMx frequency weights equal duplicated rows", {
  skip_if_not_installed("OpenMx"); skip_if_no_bridge_engine()
  set.seed(7)
  n <- 300; f <- rnorm(n)
  X <- sapply(c(.8, .7, .6), function(l) l * f + rnorm(n, 0, .6))
  colnames(X) <- paste0("x", 1:3)
  D <- data.frame(X, fq = sample(1:3, n, TRUE))
  mk <- function(dat) OpenMx::mxModel("f", type = "RAM", manifestVars = paste0("x", 1:3),
    latentVars = "F",
    OpenMx::mxPath("F", paste0("x", 1:3), values = .7, free = c(FALSE, TRUE, TRUE)),
    OpenMx::mxPath("F", arrows = 2, values = 1),
    OpenMx::mxPath(paste0("x", 1:3), arrows = 2, values = .5),
    OpenMx::mxPath("one", paste0("x", 1:3), values = 0), dat)
  dup <- OpenMx::mxRun(mk(OpenMx::mxData(D[rep(seq_len(n), D$fq), ], "raw")), silent = TRUE)
  fs  <- suppressWarnings(run_fastsem(mk(OpenMx::mxData(D, "raw", frequency = "fq")),
                                      check = "off"))
  expect_equal(OpenMx::omxGetParameters(fs), OpenMx::omxGetParameters(dup), tolerance = 1e-4)
})

test_that("lavaan sampling weights match lavaan, robust SEs included", {
  skip_if_not_installed("lavaan"); skip_if_no_bridge_engine()
  set.seed(1)
  hs <- lavaan::HolzingerSwineford1939
  hs$w <- runif(nrow(hs), .5, 2)
  mod <- "visual =~ x1 + x2 + x3\ntextual =~ x4 + x5 + x6"
  for (grp in list(NULL, "school")) {
    ref <- lavaan::cfa(mod, hs, sampling.weights = "w", group = grp)
    fs  <- run_fastsem(lavaan::cfa(mod, hs, sampling.weights = "w", group = grp,
                                   do.fit = FALSE), check = "error")
    expect_equal(lavaan::coef(fs), lavaan::coef(ref), tolerance = 1e-4)
    expect_equal(sqrt(diag(lavaan::vcov(fs))), sqrt(diag(lavaan::vcov(ref))),
                 tolerance = 1e-3)
    expect_identical(fs@Options$se, ref@Options$se)
  }
  ref <- lavaan::sem("x9 ~ x1 + x4 + x7", hs, sampling.weights = "w")
  fs  <- run_fastsem(lavaan::sem("x9 ~ x1 + x4 + x7", hs, sampling.weights = "w",
                                 do.fit = FALSE), check = "error")
  expect_equal(lavaan::coef(fs), lavaan::coef(ref), tolerance = 1e-4)
})
