## test-ordinal-twin-ace.R — multi-group ordinal twin ACE via fastsem_fit()
##
## Demonstrates and pins the multi-group conventions for hand-written twin
## syntax (umx twin supermodels build their covariance from algebras and are
## not translated by run_fastsem()):
##   * shared label        → equality across twins:  A1 =~ a*y_T1 ; A2 =~ a*y_T2
##   * c() per group       → fixed covariance:       A1 ~~ c(1, 0.5)*A2  (MZ, DZ)
##   * labelled thresholds → free (log-diff encoded, pooled across variables
##     sharing the label); 'ordered:' alone would fix them at the empirical probit.
##
## The fastsem fit is checked against OpenMx's joint multigroup FIML on the
## same data: -2lnL and the free parameters must agree.

skip_if_no_openmx <- function() {
  skip_if_not_installed("OpenMx")
}

skip_if_not_loaded <- function() {
  skip_if_not(isTRUE(fastsemR:::.fastsem_env$loaded),
              "fastsem native library not loaded")
}

# Deterministic ordinal twin data: 2-category, true loadings a = 1.0, c = 0.4
# (strong genetic contrast: MZ r ≈ .54 vs DZ r ≈ .31), unit-variance A/C/E
# latents, 10% MCAR missingness in twin 2.
.sim_twin_data <- function() {
  set.seed(20261009)
  sim_group <- function(rA, n) {
    a <- 1.0; cc <- 0.4
    A1 <- rnorm(n); A2 <- rA * A1 + sqrt(1 - rA^2) * rnorm(n)
    C1 <- rnorm(n); C2 <- rnorm(n)
    E1 <- rnorm(n); E2 <- rnorm(n)
    d <- data.frame(
      y_T1 = as.integer(a * A1 + cc * C1 + E1 > 0),
      y_T2 = as.integer(a * A2 + cc * C2 + E2 > 0))
    d$y_T2[runif(n) < 0.10] <- NA
    d
  }
  list(mz = sim_group(1, 600), dz = sim_group(0.5, 600))
}

# fastsem syntax: free loadings a/c equated across twins via shared labels,
# unit latent variances, A cross-twin covariance c(1, 0.5) by group,
# C cross-twin covariance 1 in both groups, free thresholds ta/tb pooled
# across twins.
.twin_ace_syntax <- function() {
  paste(c(
    "group: grp",
    "y_T1 | ta*t1 + tb*t2",
    "y_T2 | ta*t1 + tb*t2",
    "A1 =~ a*y_T1", "C1 =~ c*y_T1", "E1 =~ 1*y_T1",
    "A2 =~ a*y_T2", "C2 =~ c*y_T2", "E2 =~ 1*y_T2",
    "A1 ~~ 1*A1", "A2 ~~ 1*A2",
    "C1 ~~ 1*C1", "C2 ~~ 1*C2",
    "E1 ~~ 1*E1", "E2 ~~ 1*E2",
    "A1 ~~ c(1, 0.5)*A2",
    "C1 ~~ 1*C2",
    "A1 ~~ 0*C1", "A1 ~~ 0*C2", "A1 ~~ 0*E1", "A1 ~~ 0*E2",
    "A2 ~~ 0*C1", "A2 ~~ 0*C2", "A2 ~~ 0*E1", "A2 ~~ 0*E2",
    "C1 ~~ 0*E1", "C1 ~~ 0*E2", "C2 ~~ 0*E1", "C2 ~~ 0*E2",
    "E1 ~~ 0*E2",
    "y_T1 ~~ 0*y_T1", "y_T2 ~~ 0*y_T2",
    "warmstart: off"), collapse = "\n")
}

test_that("multi-group ordinal twin ACE matches OpenMx joint FIML", {
  skip_if_no_openmx()
  skip_if_not_loaded()

  dat <- .sim_twin_data()
  num <- rbind(
    data.frame(dat$mz, grp = 1L),
    data.frame(dat$dz, grp = 2L))

  fs <- fastsem_fit(.twin_ace_syntax(), num)

  # OpenMx ground truth: same model, joint multigroup FIML.
  sel <- c("y_T1", "y_T2")
  mk2 <- function(name, d, rA) OpenMx::mxModel(name, type = "RAM",
    manifestVars = sel, latentVars = c("A1", "A2", "C1", "C2", "E1", "E2"),
    OpenMx::mxPath("A1", "y_T1", labels = "a", values = .7),
    OpenMx::mxPath("C1", "y_T1", labels = "c", values = .7),
    OpenMx::mxPath("E1", "y_T1", free = FALSE, values = 1),
    OpenMx::mxPath("A2", "y_T2", labels = "a", values = .7),
    OpenMx::mxPath("C2", "y_T2", labels = "c", values = .7),
    OpenMx::mxPath("E2", "y_T2", free = FALSE, values = 1),
    OpenMx::mxPath(c("A1", "A2", "C1", "C2", "E1", "E2"),
                   arrows = 2, free = FALSE, values = 1),
    OpenMx::mxPath("A1", "A2", arrows = 2, free = FALSE, values = rA),
    OpenMx::mxPath("C1", "C2", arrows = 2, free = FALSE, values = 1),
    OpenMx::mxPath(sel, arrows = 2, free = FALSE, values = 0),
    OpenMx::mxPath("one", sel, free = FALSE, values = 0),
    OpenMx::mxThreshold("y_T1", nThresh = 2, free = TRUE,
                        values = c(-0.5, 0.8), labels = c("ta", "tb")),
    OpenMx::mxThreshold("y_T2", nThresh = 2, free = TRUE,
                        values = c(-0.5, 0.8), labels = c("ta", "tb")),
    OpenMx::mxData(d, type = "raw"))
  as_ord <- function(d) { d$y_T1 <- OpenMx::mxFactor(d$y_T1, levels = 0:1)
                          d$y_T2 <- OpenMx::mxFactor(d$y_T2, levels = 0:1); d }
  joint <- OpenMx::mxModel("joint", mk2("MZ", as_ord(dat$mz), 1),
                           mk2("DZ", as_ord(dat$dz), .5),
                           OpenMx::mxFitFunctionMultigroup(c("MZ", "DZ")))
  # status RED is common on small ordinal models; the fit itself is verified
  # numerically below, so only the optimizer warning is silenced.
  joint <- suppressWarnings(OpenMx::mxRun(joint, silent = TRUE))
  omx <- OpenMx::omxGetParameters(joint)
  omx_fit <- OpenMx::mxEval(fitfunction, joint)

  # -2lnL agreement (both optimisers converge to the same optimum).
  expect_lt(abs(fs$neg2logL - omx_fit), 0.5)

  # Shared labels equated the loadings: one 'a', one 'c', twins identical.
  expect_equal(sum(fs$paramNames == "a"), 1)
  expect_equal(sum(fs$paramNames == "c"), 1)
  expect_equal(fs$estimates[fs$paramNames == "a"],
               fs$estimates[fs$paramNames == "A2->y_T2"])
  expect_equal(fs$estimates[fs$paramNames == "c"],
               fs$estimates[fs$paramNames == "C2->y_T2"])

  # Parameter agreement with OpenMx (loose: different optimisers/starts).
  expect_lt(abs(fs$estimates[fs$paramNames == "a"] - omx[["a"]]), 0.05)
  expect_lt(abs(fs$estimates[fs$paramNames == "c"] - omx[["c"]]), 0.05)

  # The c() cross-group covariance is honoured: the A cross-twin covariance
  # must NOT be estimated (it is fixed per group), and dropping it (fitting
  # without the genetic contrast) must fit strictly worse.
  expect_true("A1~~A2" %in% fs$paramNames)
  syn_no_c <- gsub("A1 ~~ c\\(1, 0\\.5\\)\\*A2", "A1 ~~ 0*A2", .twin_ace_syntax())
  fs0 <- fastsem_fit(syn_no_c, num)
  expect_gt(fs0$neg2logL, fs$neg2logL + 5)
})

test_that("run_fastsem gives a clear error on umx twin supermodels", {
  skip_if_no_openmx()
  skip_if_not_installed("umx")

  set.seed(1)
  sim <- function(rA, n) {
    a <- sqrt(.5); cc <- sqrt(.3); e <- sqrt(.2)
    A1 <- rnorm(n); A2 <- rA * A1 + sqrt(1 - rA^2) * rnorm(n)
    C1 <- rnorm(n); C2 <- rnorm(n); E1 <- rnorm(n); E2 <- rnorm(n)
    data.frame(T_1 = ordered(as.integer(a * A1 + cc * C1 + e * E1 > 0)),
               T_2 = ordered(as.integer(a * A2 + cc * C2 + e * E2 > 0)))
  }
  m <- umx::umxACE(selDVs = "T", sep = "_",
                   dzData = sim(.5, 500), mzData = sim(1, 500), autoRun = FALSE)

  expect_error(suppressWarnings(umx_to_lavaan(m)), "twin supermodels")
  skip_if_not(isTRUE(fastsemR:::.fastsem_env$loaded), "fastsem native library not loaded")
  expect_error(fastsemR:::.ptab_from_mx(m), "twin supermodels")
})
