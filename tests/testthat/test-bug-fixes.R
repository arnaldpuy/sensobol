library(testthat)

# Regression tests for bugs fixed in 1.2.0.

# ---- Bug 1: third-order convergence plot used the wrong label ("Sijk") ----

test_that("sobol_convergence third-order plot is populated (Sijl, not Sijk)", {
  N <- 2^9
  params <- paste0("X", 1:3)
  mat <- sobol_matrices(N = N, params = params, order = "third")
  Y <- ishigami_Fun(mat)
  out <- sobol_convergence(matrices = c("A", "B", "AB"), Y = Y, N = N,
                           sub.sample = c(100, 200, 400), params = params,
                           first = "saltelli", total = "jansen",
                           order = "third", plot.order = "third")
  expect_true("plot.third" %in% names(out))
  expect_s3_class(out$plot.third, "ggplot")
  # The data underlying the third-order plot must not be empty.
  expect_gt(nrow(out$plot.third$data), 0)
  expect_true(all(out$plot.third$data$sensitivity == "Sijl"))
})

# ---- Bug 2: single parameter / single group dropped the matrix to a vector ----

test_that("sobol_indices works with a single parameter", {
  N <- 256
  m <- sobol_matrices(N = N, params = "X1")
  Y <- m[, 1]^2 + m[, 1]
  ind <- sobol_indices(Y = Y, N = N, params = "X1")
  expect_s3_class(ind$results, "data.table")
  expect_equal(nrow(ind$results[sensitivity == "Si"]), 1)
  expect_equal(nrow(ind$results[sensitivity == "Ti"]), 1)
})

test_that("sobol_indices works when groups collapse to a single group", {
  N <- 256
  params <- paste0("X", 1:3)
  g <- list(all = params)
  m <- sobol_matrices(N = N, params = params, groups = g)
  Y <- ishigami_Fun(m)
  ind <- sobol_indices(Y = Y, N = N, params = params, groups = g)
  expect_s3_class(ind$results, "data.table")
  expect_equal(nrow(ind$results[sensitivity == "Si"]), 1)
  expect_equal(ind$results[sensitivity == "Si", parameters], "all")
})

test_that("single-parameter sobol_indices works on the four-matrix design", {
  N <- 256
  m <- sobol_matrices(N = N, params = "X1", matrices = c("A", "B", "AB", "BA"))
  Y <- m[, 1]^2
  ind <- sobol_indices(matrices = c("A", "B", "AB", "BA"),
                       Y = Y, N = N, params = "X1",
                       first = "owen", total = "owen")
  expect_s3_class(ind$results, "data.table")
  expect_equal(nrow(ind$results[sensitivity == "Si"]), 1)
})

# ---- Bug 3: higher-order plots were silently empty without bootstrap ----

test_that("higher-order plot is populated without bootstrap", {
  N <- 2^9
  params <- paste0("X", 1:3)
  mat <- sobol_matrices(N = N, params = params, order = "second")
  Y <- ishigami_Fun(mat)
  ind <- sobol_indices(Y = Y, N = N, params = params, order = "second")
  gg <- plot(ind, order = "second")
  expect_s3_class(gg, "ggplot")
  # Without bootstrap, all Sij point estimates should be plotted.
  expect_equal(nrow(gg$data), nrow(ind$results[sensitivity == "Sij"]))
})

test_that("higher-order plot still works with bootstrap", {
  N <- 2^9
  params <- paste0("X", 1:3)
  mat <- sobol_matrices(N = N, params = params, order = "second")
  Y <- ishigami_Fun(mat)
  ind <- sobol_indices(Y = Y, N = N, params = params, order = "second",
                       boot = TRUE, R = 50)
  gg <- plot(ind, order = "second")
  expect_s3_class(gg, "ggplot")
})

# ---- Bug 4: discrepancy_ersatz dropped points with an exact-0 x coordinate ----

test_that("discrepancy_ersatz handles an exact 0 in the first column", {
  set.seed(1)
  N <- 100
  x <- c(0, stats::runif(N - 1))   # an exact 0 in the input column
  Y <- stats::rnorm(N)
  val <- sensobol:::s_ersatz(cbind(x, Y))
  expect_true(is.numeric(val))
  expect_false(is.na(val))
  # The point at x = 0 must be counted, not silently dropped: compare with
  # the same configuration where the 0 is nudged to a tiny positive value.
  x2 <- x; x2[1] <- 1e-12
  val_nudged <- sensobol:::s_ersatz(cbind(x2, Y))
  expect_equal(val, val_nudged)
})
