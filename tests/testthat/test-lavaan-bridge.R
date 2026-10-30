## test-lavaan-bridge.R — run_fastsem() on unfitted lavaan objects.
##
## The returned object must be a fitted lavaan object equal to lavaan's own fit
## (logLik, estimates, SEs) and work with lavaan's tools downstream.

skip_if_no_lavaan <- function() skip_if_not_installed("lavaan")

expect_same_lav <- function(ref, fs) {
  expect_s4_class(fs, "lavaan")
  expect_equal(as.numeric(lavaan::logLik(fs)), as.numeric(lavaan::logLik(ref)),
               tolerance = 1e-7)
  expect_equal(lavaan::coef(fs), lavaan::coef(ref), tolerance = 1e-4)
  expect_equal(sqrt(diag(lavaan::vcov(fs))), sqrt(diag(lavaan::vcov(ref))),
               tolerance = 1e-4)
}

hs_model <- "visual =~ x1 + x2 + x3\ntextual =~ x4 + x5 + x6"

test_that("single-group CFA matches lavaan and works downstream", {
  skip_if_no_lavaan(); skip_if_no_bridge_engine()
  hs  <- lavaan::HolzingerSwineford1939
  ref <- lavaan::cfa(hs_model, hs, meanstructure = TRUE)
  fs  <- run_fastsem(lavaan::cfa(hs_model, hs, meanstructure = TRUE, do.fit = FALSE),
                     check = "error")
  expect_same_lav(ref, fs)
  expect_equal(lavaan::fitMeasures(fs, c("chisq", "df", "cfi", "rmsea")),
               lavaan::fitMeasures(ref, c("chisq", "df", "cfi", "rmsea")),
               tolerance = 1e-4)
  expect_equal(lavaan::standardizedSolution(fs)$est.std,
               lavaan::standardizedSolution(ref)$est.std, tolerance = 1e-4)
  expect_false(is.null(fs@external$fastsem))
})

test_that("multi-group with group.equal matches lavaan; lavTestLRT works", {
  skip_if_no_lavaan(); skip_if_no_bridge_engine()
  hs <- lavaan::HolzingerSwineford1939
  m  <- paste(hs_model, "textual ~ visual", sep = "\n")
  ref_eq <- lavaan::sem(m, hs, group = "school", group.equal = "loadings")
  fs_eq  <- run_fastsem(lavaan::sem(m, hs, group = "school", group.equal = "loadings",
                                    do.fit = FALSE), check = "error")
  expect_same_lav(ref_eq, fs_eq)
  fs_free <- run_fastsem(lavaan::sem(m, hs, group = "school", do.fit = FALSE),
                         check = "error")
  lrt <- lavaan::lavTestLRT(fs_free, fs_eq)
  expect_equal(lrt[["Df diff"]][2], 4)
})

test_that("FIML, fixed.x covariates and := match lavaan", {
  skip_if_no_lavaan(); skip_if_no_bridge_engine()
  d <- as.data.frame(scale(mtcars[, c("mpg", "hp", "wt")])); d$wt[1:3] <- NA
  m <- "wt ~ a*hp\nmpg ~ b*wt + c*hp\nab := a*b"
  ref <- lavaan::sem(m, d, missing = "ml")
  fs  <- run_fastsem(lavaan::sem(m, d, missing = "ml", do.fit = FALSE), check = "error")
  expect_same_lav(ref, fs)
  pe_f <- lavaan::parameterEstimates(fs); pe_r <- lavaan::parameterEstimates(ref)
  expect_equal(pe_f$est[pe_f$op == ":="], pe_r$est[pe_r$op == ":="], tolerance = 1e-4)
})

test_that("robust estimators keep lavaan's own robust SEs and tests", {
  skip_if_no_lavaan(); skip_if_no_bridge_engine()
  hs  <- lavaan::HolzingerSwineford1939
  ref <- lavaan::cfa(hs_model, hs, estimator = "MLR")
  fs  <- run_fastsem(lavaan::cfa(hs_model, hs, estimator = "MLR", do.fit = FALSE),
                     check = "error")
  expect_same_lav(ref, fs)
  expect_equal(lavaan::fitMeasures(fs, "chisq.scaled"),
               lavaan::fitMeasures(ref, "chisq.scaled"), tolerance = 1e-4)
})

test_that("unsupported lavaan models fail loudly", {
  skip_if_no_lavaan(); skip_if_no_bridge_engine()
  hs <- lavaan::HolzingerSwineford1939
  hs$x1o <- ordered(cut(hs$x1, 3))
  expect_error(run_fastsem(lavaan::cfa("f =~ x1o + x2 + x3", hs, ordered = "x1o",
                                       do.fit = FALSE)), "ordered")
  expect_error(run_fastsem(lavaan::cfa(hs_model, hs, estimator = "ULS", do.fit = FALSE)),
               "not ML-based")
  expect_error(run_fastsem(list()), "MxModel or a lavaan object")
})
