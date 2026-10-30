## test-partable.R — the internal parameter table and its fastsem syntax
## (no engine needed).

test_that("OpenMx RAM model → table: cells, labels, keys, plabels", {
  skip_if_not_installed("OpenMx"); skip_if_not_installed("lavaan")
  hs <- lavaan::HolzingerSwineford1939
  v  <- c("x1", "x2", "x3")
  m <- OpenMx::mxModel("t", type = "RAM", manifestVars = v, latentVars = "g",
    OpenMx::mxPath("g", v, values = .5, labels = c("l", "l", NA)),
    OpenMx::mxPath("g", arrows = 2, free = FALSE, values = 1),
    OpenMx::mxPath(v, arrows = 2, values = 1),
    OpenMx::mxPath("one", v),
    OpenMx::mxData(hs[, v], type = "raw"))
  pt <- fastsemR:::.ptab_from_mx(m)
  expect_setequal(unique(pt$mat), c("A", "S", "M"))
  a <- pt[pt$mat == "A", ]
  expect_equal(nrow(a), 3L)
  # shared label → one parameter; unlabelled → its own
  expect_equal(length(unique(a$plabel)), 2L)
  expect_equal(a$plabel[a$row == "x1"], a$plabel[a$row == "x2"])
  expect_equal(a$omx[a$row == "x3"], "t.A[3,4]")
  # fixed latent variance kept as a fixed cell
  expect_true(any(pt$mat == "S" & pt$row == "g" & !pt$free & pt$value == 1))
  expect_equal(length(unique(pt$plabel[pt$free])), 8L)
})

test_that("engine syntax uses internal names only and full-precision numbers", {
  skip_if_not_installed("OpenMx"); skip_if_not_installed("lavaan")
  skip_if_no_bridge_engine()
  hs <- lavaan::HolzingerSwineford1939
  d  <- hs[, c("x1", "x2", "x3")]; names(d) <- c("a_r1c1", "start", "data_x")
  v  <- names(d)
  m <- OpenMx::mxModel("odd", type = "RAM", manifestVars = v, latentVars = "g",
    OpenMx::mxPath("g", v, values = .5, labels = c("a_r1c1", "t", "c(1)")),
    OpenMx::mxPath("g", arrows = 2, free = FALSE, values = 1 / 3),
    OpenMx::mxPath(v, arrows = 2, values = 1), OpenMx::mxPath("one", v),
    OpenMx::mxData(d, type = "raw"))
  pt  <- fastsemR:::.ptab_from_mx(m)
  syn <- fastsemR:::.ptab_syntax(pt)$syntax
  expect_false(any(grepl("a_r1c1|start\\*|data_x|c\\(1\\)", syn)))
  expect_match(syn, "0.3333333333333333", fixed = TRUE)
  expect_match(syn, "likelihood: normal", fixed = TRUE)
  # the data keep the model's own column names; the engine renames them
  expect_equal(colnames(fastsemR:::.ptab_engine_args(pt)$data), v)
})

test_that("multi-group table → c() only where groups differ", {
  skip_if_not_installed("OpenMx"); skip_if_not_installed("lavaan")
  skip_if_no_bridge_engine()
  hs <- lavaan::HolzingerSwineford1939
  v  <- c("x1", "x2", "x3")
  mk <- function(nm, d) OpenMx::mxModel(nm, type = "RAM", manifestVars = v,
    latentVars = "g",
    OpenMx::mxPath("g", v, values = .5, labels = c("l1", "l2", "l3")),
    OpenMx::mxPath("g", arrows = 2, free = FALSE, values = 1),
    OpenMx::mxPath(v, arrows = 2, values = 1), OpenMx::mxPath("one", v),
    OpenMx::mxData(d[, v], type = "raw"))
  mg <- OpenMx::mxModel("mg", mk("p", hs[hs$school == "Pasteur", ]),
                              mk("q", hs[hs$school == "Grant-White", ]))
  pt  <- fastsemR:::.ptab_from_mx(mg)
  syn <- fastsemR:::.ptab_syntax(pt)$syntax
  dat <- fastsemR:::.ptab_engine_args(pt)$data
  expect_match(syn, "group: fastsem_grp", fixed = TRUE)
  lines <- strsplit(syn, "\n")[[1]]
  # shared loadings: one label, no c(); unlabelled variances: c(pA, pB)
  expect_false(any(grepl("=~.*c\\(", lines)))
  expect_true(any(grepl("~~ c\\(p[0-9]+, p[0-9]+\\)\\*", lines)))
  expect_equal(sort(unique(dat[, "fastsem_grp"])), c(1, 2))
  expect_equal(unname(dat[1, "fastsem_grp"]), 1)
})

test_that("validation rejects what fastsem cannot express", {
  skip_if_not_installed("OpenMx"); skip_if_not_installed("lavaan")
  hs <- lavaan::HolzingerSwineford1939
  v  <- c("x1", "x2", "x3")
  base <- function() OpenMx::mxModel("b", type = "RAM", manifestVars = v,
    latentVars = "g", OpenMx::mxPath("g", v, values = .5),
    OpenMx::mxPath("g", arrows = 2, free = FALSE, values = 1),
    OpenMx::mxPath(v, arrows = 2, values = 1), OpenMx::mxPath("one", v),
    OpenMx::mxData(hs[, c(v, "ageyr")], type = "raw"))
  m <- base(); m$S$labels["x1", "x1"] <- "data.ageyr"; m$S$free["x1", "x1"] <- FALSE
  expect_error(fastsemR:::.ptab_from_mx(m), "paths and loadings \\(A\\) only")
  m <- base(); m$A$labels["x2", "x1"] <- "data.ageyr"
  expect_error(fastsemR:::.ptab_from_mx(m), "from latent variables only")
  m <- base(); m$A$labels["x1", "g"] <- "alg[1,1]"
  expect_error(fastsemR:::.ptab_from_mx(m), "algebra or matrix element")
})

test_that("the engine's table schema matches what fastsemR sends", {
  skip_if_no_bridge_engine()
  expect_true(fastsemR:::.check_table_schema())
})

test_that("umx deviation thresholds become threshold rows (Mehta, 4 levels)", {
  skip_if_not_installed("OpenMx"); skip_if_not_installed("umx")
  set.seed(4); n <- 200; f <- rnorm(n)
  y <- paste0("y", 1:3)
  d <- as.data.frame(sapply(1:3, function(i)
    as.integer(cut(.7 * f + rnorm(n, sd = .7), c(-Inf, -.6, .1, .8, Inf))) - 1L))
  names(d) <- y
  for (v in y) d[[v]] <- OpenMx::mxFactor(d[[v]], levels = 0:3)
  m <- suppressMessages(umx::umxRAM("dev", umx::umxPath("g", to = y),
         umx::umxPath(v1m0 = "g"), umx::umxPath(v.m. = y), data = d, autoRun = FALSE))
  pt <- fastsemR:::.ptab_from_mx(m)
  th <- pt[pt$mat == "T" & pt$col == "y1", ]
  th <- th[order(th$row), ]
  expect_equal(th$row, c("t1", "t2", "t3"))
  # umx's Mehta default: first two deviations fixed, the third free
  expect_equal(th$free, c(FALSE, FALSE, TRUE))
  Y <- m$deviations_for_thresh
  expect_equal(th$value, unname(cumsum(Y$values[1:3, "y1"])))
  dm <- attr(pt, "devmap")
  expect_true(all(c("y1", "y2", "y3") %in% dm$var))
  expect_equal(sum(dm$free & dm$var == "y1"), 1L)
})

test_that("an ordinal variable without two pins is rejected in OpenMx terms", {
  skip_if_not_installed("OpenMx")
  set.seed(5); n <- 200; f <- rnorm(n)
  d <- data.frame(y1 = OpenMx::mxFactor(as.integer(f + rnorm(n) > 0), levels = 0:1),
                  y2 = rnorm(n), y3 = rnorm(n))
  m <- OpenMx::mxModel("bad", type = "RAM", manifestVars = names(d), latentVars = "g",
    OpenMx::mxPath("g", names(d), values = .5), OpenMx::mxPath("g", arrows = 2, free = FALSE, values = 1),
    OpenMx::mxPath(names(d), arrows = 2, values = 1), OpenMx::mxPath("one", names(d)),
    OpenMx::mxThreshold("y1", nThresh = 1, free = TRUE, values = 0),
    OpenMx::mxData(d, type = "raw"))
  expect_error(fastsemR:::.ptab_from_mx(m), "'y1' \\(group 1\\) is not identified")
})
