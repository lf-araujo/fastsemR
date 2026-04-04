## test-fit.R — basic fastsem_fit() tests using built-in R datasets

skip_if_not_loaded <- function() {
  skip_if_not(isTRUE(fastsemR:::.fastsem_env$loaded),
               "fastsem native library not loaded")
}

# ── fastsem_fit: path model on scaled mtcars ───────────────────────────────────

test_that("path model converges and returns expected structure", {
  skip_if_not_loaded()

  df <- as.data.frame(scale(mtcars[, c("mpg", "hp", "wt")]))
  syntax <- "
    mpg ~ hp + wt
    hp ~~ wt
    hp ~~ hp
    wt ~~ wt
    mpg ~~ mpg
    mpg ~ 1
    hp  ~ 1
    wt  ~ 1
  "
  res <- fastsem_fit(syntax, df)

  expect_type(res, "list")
  expect_true(!is.null(res$estimates))
  expect_true(!is.null(res$paramNames))
  expect_true(!is.null(res$se))
  expect_true(length(res$estimates) == length(res$paramNames))
})

test_that("path model hp->mpg coefficient is negative", {
  skip_if_not_loaded()

  df <- as.data.frame(scale(mtcars[, c("mpg", "hp", "wt")]))
  syntax <- "
    mpg ~ hp + wt
    hp ~~ wt
    hp ~~ hp;  wt ~~ wt;  mpg ~~ mpg
    mpg ~ 1;   hp ~ 1;    wt ~ 1
  "
  res  <- fastsem_fit(syntax, df)
  idx  <- which(res$paramNames == "hp->mpg")
  expect_length(idx, 1L)
  expect_lt(res$estimates[[idx]], 0)
})

# ── fastsem_fit: single-factor CFA on iris ─────────────────────────────────────

test_that("one-factor CFA on iris converges", {
  skip_if_not_loaded()

  df <- as.data.frame(scale(iris[, c(1, 3, 4)]))
  names(df) <- c("sl", "pl", "pw")
  syntax <- "
    g =~ 1*sl + pl + pw
    g ~~ g
    sl ~~ sl;  pl ~~ pl;  pw ~~ pw
  "
  res <- fastsem_fit(syntax, df)

  expect_true(all(is.finite(res$estimates)))
  # All loadings should be positive (positively-correlated iris variables)
  for (nm in c("g->pl", "g->pw")) {
    idx <- which(res$paramNames == nm)
    if (length(idx) > 0)
      expect_gt(res$estimates[[idx]], 0)
  }
})

# ── fastsem_fit: FIML with missing data ────────────────────────────────────────

test_that("FIML runs when data has NAs", {
  skip_if_not_loaded()

  set.seed(42)
  df <- as.data.frame(scale(mtcars[, c("mpg", "hp", "wt")]))
  df$mpg[sample(32, 4)] <- NA

  syntax <- "
    mpg ~ hp + wt
    hp ~~ wt
    hp ~~ hp;  wt ~~ wt;  mpg ~~ mpg
    mpg ~ 1;   hp ~ 1;    wt ~ 1
  "
  res <- fastsem_fit(syntax, df)

  expect_true(!is.null(res$neg2logL))
  expect_true(is.finite(res$neg2logL))
  expect_equal(res$estimator, "FIML")
})

# ── fastsem_fit: standardized estimates returned ───────────────────────────────

test_that("standardized estimates have same length as paramNames", {
  skip_if_not_loaded()

  df <- as.data.frame(scale(mtcars[, c("mpg", "hp", "wt")]))
  syntax <- "
    mpg ~ hp + wt
    hp ~~ wt
    hp ~~ hp;  wt ~~ wt;  mpg ~~ mpg
    mpg ~ 1;   hp ~ 1;    wt ~ 1
  "
  res <- fastsem_fit(syntax, df)

  expect_equal(length(res$stdEst), length(res$paramNames))
})

# ── fastsem_sem_vars ───────────────────────────────────────────────────────────

test_that("fastsem_sem_vars extracts manifest variable names", {
  skip_if_not_loaded()

  syntax <- "f =~ y1 + y2 + y3\nf ~~ f\ny1 ~~ y1\ny2 ~~ y2\ny3 ~~ y3"
  cols   <- c("id", "y1", "y2", "y3", "y4")
  vars   <- fastsem_sem_vars(syntax, cols)

  expect_true(all(c("y1", "y2", "y3") %in% vars))
  expect_false("y4" %in% vars)
  expect_false("id" %in% vars)
})
