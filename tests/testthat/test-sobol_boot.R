library(testthat)

# ---- sobol_indices: basic workflow ----

test_that("sobol_indices returns a data.table with expected columns", {
  N <- 100
  params <- paste("X", 1:3, sep = "")
  mat <- sobol_matrices(N = N, params = params)
  Y <- ishigami_Fun(mat)
  ind <- sobol_indices(Y = Y, N = N, params = params)
  expect_s3_class(ind$results, "data.table")
  expect_true(all(c("sensitivity", "parameters", "original") %in% colnames(ind$results)))
})

test_that("sobol_indices first-order indices are between -1 and 1", {
  N <- 200
  params <- paste("X", 1:3, sep = "")
  mat <- sobol_matrices(N = N, params = params)
  Y <- ishigami_Fun(mat)
  ind <- sobol_indices(Y = Y, N = N, params = params)
  si <- ind$results[sensitivity == "Si", original]
  expect_true(all(si >= -1 & si <= 1))
})

test_that("sobol_indices total-order indices are between 0 and 1 (approx)", {
  N <- 200
  params <- paste("X", 1:3, sep = "")
  mat <- sobol_matrices(N = N, params = params)
  Y <- ishigami_Fun(mat)
  ind <- sobol_indices(Y = Y, N = N, params = params)
  ti <- ind$results[sensitivity == "Ti", original]
  expect_true(all(ti >= -0.1))  # allow small numerical negatives
})

# ---- sobol_indices: bootstrap ----

test_that("sobol_indices bootstrap adds confidence interval columns", {
  N <- 100
  params <- paste("X", 1:3, sep = "")
  mat <- sobol_matrices(N = N, params = params)
  Y <- ishigami_Fun(mat)
  ind <- sobol_indices(Y = Y, N = N, params = params, boot = TRUE, R = 10)
  expect_true(all(c("low.ci", "high.ci") %in% colnames(ind$results)))
})

# ---- sobol_indices: error paths ----

test_that("sobol_indices stops when Y contains NA", {
  N <- 100
  params <- paste("X", 1:3, sep = "")
  mat <- sobol_matrices(N = N, params = params)
  Y <- ishigami_Fun(mat)
  Y[5] <- NA
  expect_error(sobol_indices(Y = Y, N = N, params = params),
               "NA or NaN")
})

test_that("sobol_indices stops when boot = TRUE but R is NULL", {
  N <- 100
  params <- paste("X", 1:3, sep = "")
  mat <- sobol_matrices(N = N, params = params)
  Y <- ishigami_Fun(mat)
  expect_error(sobol_indices(Y = Y, N = N, params = params, boot = TRUE),
               regexp = NULL)
})

test_that("sobol_indices stops with invalid estimator for AB matrices", {
  N <- 100
  params <- paste("X", 1:3, sep = "")
  mat <- sobol_matrices(N = N, params = params)
  Y <- ishigami_Fun(mat)
  expect_error(
    sobol_indices(Y = Y, N = N, params = params, first = "azzini"),
    "Revise the correspondence"
  )
})

# ---- sobol_dummy ----

test_that("sobol_dummy returns a data.table with Si and Ti", {
  N <- 100
  params <- paste("X", 1:3, sep = "")
  mat <- sobol_matrices(N = N, params = params)
  Y <- ishigami_Fun(mat)
  out <- sobol_dummy(Y = Y, N = N, params = params)
  expect_s3_class(out, "data.table")
  expect_true(all(c("Si", "Ti") %in% out$sensitivity))
})

test_that("sobol_dummy stops when Y contains NA", {
  N <- 100
  params <- paste("X", 1:3, sep = "")
  mat <- sobol_matrices(N = N, params = params)
  Y <- ishigami_Fun(mat)
  Y[1] <- NaN
  expect_error(sobol_dummy(Y = Y, N = N, params = params), "NA or NaN")
})

# ---- sobol_matrices: input validation ----

test_that("sobol_matrices stops when N is not a positive integer", {
  params <- paste("X", 1:3, sep = "")
  expect_error(sobol_matrices(N = -1, params = params), "'N' must be a single positive integer")
  expect_error(sobol_matrices(N = 1.5, params = params), "'N' must be a single positive integer")
  expect_error(sobol_matrices(N = "10", params = params), "'N' must be a single positive integer")
})

# ---- vars_matrices: input validation ----

test_that("vars_matrices stops when star.centers is not a positive integer", {
  params <- paste("X", 1:3, sep = "")
  expect_error(vars_matrices(star.centers = -1, params = params), "'star.centers' must be a single positive integer")
  expect_error(vars_matrices(star.centers = 2.5, params = params), "'star.centers' must be a single positive integer")
})

# ---- load_packages ----

test_that("load_packages stops with informative error for missing package", {
  expect_error(load_packages("_package_that_does_not_exist_xyz_"),
               "required but not installed")
})

# ---- test functions (ishigami, sobol_Fun) ----

test_that("ishigami_Fun returns a numeric vector of correct length", {
  N <- 50
  params <- paste("X", 1:3, sep = "")
  mat <- sobol_matrices(N = N, params = params)
  Y <- ishigami_Fun(mat)
  expect_type(Y, "double")
  expect_length(Y, nrow(mat))
})

test_that("sobol_Fun returns a numeric vector of correct length", {
  N <- 50
  params <- paste("X", 1:8, sep = "")
  mat <- sobol_matrices(N = N, params = params)
  Y <- sobol_Fun(mat)
  expect_type(Y, "double")
  expect_length(Y, nrow(mat))
})

# ---- high-order indices ----

test_that("sobol_indices computes second-order indices", {
  N <- 100
  params <- paste("X", 1:3, sep = "")
  mat <- sobol_matrices(N = N, params = params, order = "second")
  Y <- ishigami_Fun(mat)
  ind <- sobol_indices(Y = Y, N = N, params = params, order = "second")
  expect_true("Sij" %in% ind$results$sensitivity)
  expect_equal(nrow(ind$results[sensitivity == "Sij"]),
               ncol(utils::combn(params, 2)))
})

test_that("sobol_indices computes third-order indices", {
  N <- 100
  params <- paste("X", 1:3, sep = "")
  mat <- sobol_matrices(N = N, params = params, order = "third")
  Y <- ishigami_Fun(mat)
  ind <- sobol_indices(Y = Y, N = N, params = params, order = "third")
  expect_true("Sijl" %in% ind$results$sensitivity)
})

# ---- New estimators: Owen, Martinez, Mauntz ----

# helper: produce Ishigami output on the appropriate sampling design
.ish_data <- function(N, params, matrices = c("A", "B", "AB"), order = "first") {
  mat <- sobol_matrices(N = N, params = params, matrices = matrices, order = order)
  list(Y = ishigami_Fun(mat), mat = mat)
}

# -- mauntz (requires A, B, AB) --

test_that("mauntz Si matches saltelli Si in expectation on Ishigami", {
  N <- 2^13
  params <- paste0("X", 1:3)
  d <- .ish_data(N, params)
  ind_s <- sobol_indices(Y = d$Y, N = N, params = params,
                         first = "saltelli", total = "jansen")
  ind_m <- sobol_indices(Y = d$Y, N = N, params = params,
                         first = "mauntz", total = "jansen")
  si_s <- ind_s$results[sensitivity == "Si", original]
  si_m <- ind_m$results[sensitivity == "Si", original]
  expect_equal(si_s, si_m, tolerance = 0.01)
})

test_that("mauntz Si runs with higher-order designs", {
  N <- 200; params <- paste0("X", 1:3)
  d <- .ish_data(N, params, order = "third")
  ind <- sobol_indices(Y = d$Y, N = N, params = params,
                       first = "mauntz", total = "jansen", order = "third")
  expect_true("Sij"  %in% ind$results$sensitivity)
  expect_true("Sijl" %in% ind$results$sensitivity)
})

test_that("mauntz Si runs with groups and bootstrap", {
  N <- 200; params <- paste0("X", 1:3)
  groups <- list(g1 = "X1", g23 = c("X2", "X3"))
  mat <- sobol_matrices(N = N, params = params, groups = groups)
  Y <- ishigami_Fun(mat)
  ind <- sobol_indices(Y = Y, N = N, params = params, groups = groups,
                       first = "mauntz", total = "jansen",
                       boot = TRUE, R = 25)
  expect_equal(nrow(ind$results[sensitivity == "Si"]), length(groups))
  expect_true(all(c("low.ci", "high.ci") %in% colnames(ind$results)))
})

# -- martinez (requires A, B, AB) --

test_that("martinez Si matches saltelli Si in expectation on Ishigami", {
  N <- 2^13; params <- paste0("X", 1:3)
  d <- .ish_data(N, params)
  ind_s <- sobol_indices(Y = d$Y, N = N, params = params,
                         first = "saltelli", total = "jansen")
  ind_m <- sobol_indices(Y = d$Y, N = N, params = params,
                         first = "martinez", total = "jansen")
  expect_equal(ind_s$results[sensitivity == "Si", original],
               ind_m$results[sensitivity == "Si", original],
               tolerance = 0.02)
})

test_that("martinez Si runs with higher-order designs", {
  N <- 200; params <- paste0("X", 1:3)
  d <- .ish_data(N, params, order = "second")
  ind <- sobol_indices(Y = d$Y, N = N, params = params,
                       first = "martinez", total = "jansen", order = "second")
  expect_true("Sij" %in% ind$results$sensitivity)
})

test_that("martinez Si runs with groups and bootstrap", {
  N <- 200; params <- paste0("X", 1:3)
  groups <- list(g1 = "X1", g23 = c("X2", "X3"))
  mat <- sobol_matrices(N = N, params = params, groups = groups)
  Y <- ishigami_Fun(mat)
  ind <- sobol_indices(Y = Y, N = N, params = params, groups = groups,
                       first = "martinez", total = "jansen",
                       boot = TRUE, R = 25)
  expect_equal(nrow(ind$results[sensitivity == "Si"]), length(groups))
  expect_true(all(c("low.ci", "high.ci") %in% colnames(ind$results)))
})

# -- owen Si and Ti (require A, B, AB, BA) --

test_that("owen Si matches saltelli Si in expectation on Ishigami", {
  N <- 2^13; params <- paste0("X", 1:3)
  d <- .ish_data(N, params, matrices = c("A", "B", "AB", "BA"))
  ind_s <- sobol_indices(matrices = c("A","B","AB","BA"),
                         Y = d$Y, N = N, params = params,
                         first = "saltelli", total = "jansen")
  ind_o <- sobol_indices(matrices = c("A","B","AB","BA"),
                         Y = d$Y, N = N, params = params,
                         first = "owen", total = "owen")
  expect_equal(ind_s$results[sensitivity == "Si", original],
               ind_o$results[sensitivity == "Si", original],
               tolerance = 0.02)
  expect_equal(ind_s$results[sensitivity == "Ti", original],
               ind_o$results[sensitivity == "Ti", original],
               tolerance = 0.02)
})

test_that("owen Si and Ti run with higher-order designs", {
  N <- 200; params <- paste0("X", 1:3)
  d <- .ish_data(N, params, matrices = c("A","B","AB","BA"), order = "second")
  ind <- sobol_indices(matrices = c("A","B","AB","BA"),
                       Y = d$Y, N = N, params = params,
                       first = "owen", total = "owen", order = "second")
  expect_true("Sij" %in% ind$results$sensitivity)
})

test_that("owen Si and Ti run with groups and bootstrap", {
  N <- 200; params <- paste0("X", 1:3)
  groups <- list(g1 = "X1", g23 = c("X2", "X3"))
  mat <- sobol_matrices(N = N, params = params,
                        matrices = c("A","B","AB","BA"), groups = groups)
  Y <- ishigami_Fun(mat)
  ind <- sobol_indices(matrices = c("A","B","AB","BA"),
                       Y = Y, N = N, params = params, groups = groups,
                       first = "owen", total = "owen",
                       boot = TRUE, R = 25)
  expect_equal(nrow(ind$results[sensitivity == "Si"]), length(groups))
  expect_equal(nrow(ind$results[sensitivity == "Ti"]), length(groups))
  expect_true(all(c("low.ci", "high.ci") %in% colnames(ind$results)))
})

# -- error paths: wrong matrix design --

test_that("owen first errors with c('A','B','AB') design", {
  N <- 100; params <- paste0("X", 1:3)
  d <- .ish_data(N, params)
  expect_error(
    sobol_indices(Y = d$Y, N = N, params = params,
                  first = "owen", total = "jansen"),
    "Revise the correspondence"
  )
})

test_that("owen total errors with c('A','B','AB') design", {
  N <- 100; params <- paste0("X", 1:3)
  d <- .ish_data(N, params)
  expect_error(
    sobol_indices(Y = d$Y, N = N, params = params,
                  first = "saltelli", total = "owen"),
    "Revise the correspondence"
  )
})

test_that("martinez and mauntz error with c('A','B','BA') design", {
  N <- 100; params <- paste0("X", 1:3)
  mat <- sobol_matrices(N = N, params = params, matrices = c("A","B","BA"))
  Y <- ishigami_Fun(mat)
  expect_error(
    sobol_indices(matrices = c("A","B","BA"), Y = Y, N = N, params = params,
                  first = "martinez", total = "saltelli"),
    "Revise the correspondence"
  )
  expect_error(
    sobol_indices(matrices = c("A","B","BA"), Y = Y, N = N, params = params,
                  first = "mauntz", total = "saltelli"),
    "Revise the correspondence"
  )
})

# ---- sobol_indices with groups ----

test_that("sobol_indices with groups returns one Si/Ti per group", {
  N <- 200
  params <- paste("X", 1:3, sep = "")
  groups <- list(g1 = "X1", g23 = c("X2", "X3"))
  mat <- sobol_matrices(N = N, params = params, groups = groups)
  Y <- ishigami_Fun(mat)
  ind <- sobol_indices(Y = Y, N = N, params = params, groups = groups)
  expect_equal(nrow(ind$results[sensitivity == "Si"]), length(groups))
  expect_equal(nrow(ind$results[sensitivity == "Ti"]), length(groups))
  expect_setequal(ind$results[sensitivity == "Si", parameters], names(groups))
})

test_that("sobol_indices with trivial groups equals ungrouped result", {
  N <- 200
  params <- paste("X", 1:3, sep = "")
  trivial <- as.list(params); names(trivial) <- params
  set.seed(42)
  mat1 <- sobol_matrices(N = N, params = params, type = "R")
  set.seed(42)
  mat2 <- sobol_matrices(N = N, params = params, type = "R", groups = trivial)
  Y1 <- ishigami_Fun(mat1)
  Y2 <- ishigami_Fun(mat2)
  ind1 <- sobol_indices(Y = Y1, N = N, params = params)
  ind2 <- sobol_indices(Y = Y2, N = N, params = params, groups = trivial)
  expect_equal(ind1$results$original, ind2$results$original)
})

test_that("sobol_indices with groups + second order gives one Sij per pair", {
  N <- 200
  params <- paste("X", 1:3, sep = "")
  groups <- list(g1 = "X1", g23 = c("X2", "X3"))
  mat <- sobol_matrices(N = N, params = params, groups = groups, order = "second")
  Y <- ishigami_Fun(mat)
  ind <- sobol_indices(Y = Y, N = N, params = params, groups = groups,
                       order = "second")
  expect_equal(nrow(ind$results[sensitivity == "Sij"]),
               choose(length(groups), 2))
  expect_equal(ind$results[sensitivity == "Sij", parameters],
               unlist(lapply(utils::combn(names(groups), 2, simplify = FALSE),
                             function(x) paste0(x, collapse = "."))))
})

test_that("sobol_indices with groups bootstraps without error", {
  N <- 100
  params <- paste("X", 1:3, sep = "")
  groups <- list(g1 = "X1", g23 = c("X2", "X3"))
  mat <- sobol_matrices(N = N, params = params, groups = groups)
  Y <- ishigami_Fun(mat)
  ind <- sobol_indices(Y = Y, N = N, params = params, groups = groups,
                       boot = TRUE, R = 10)
  expect_true(all(c("low.ci", "high.ci") %in% colnames(ind$results)))
  expect_equal(nrow(ind$results[sensitivity == "Si"]), length(groups))
})

# ---- sobol_convergence ----

test_that("sobol_convergence returns a list with convergence data and a plot", {
  N <- 500
  params <- paste("X", 1:3, sep = "")
  mat <- sobol_matrices(N = N, params = params)
  Y <- ishigami_Fun(mat)
  sub.sample <- c(100, 200, 300)
  out <- sobol_convergence(matrices = c("A", "B", "AB"), Y = Y, N = N,
                           sub.sample = sub.sample, params = params,
                           first = "saltelli", total = "jansen",
                           order = "first", plot.order = "first")
  expect_true(all(c("convergence", "plot.first") %in% names(out)))
  expect_s3_class(out$convergence, "data.table")
  expect_s3_class(out$plot.first, "ggplot")
})

test_that("sobol_convergence stops when sub.sample is unsorted", {
  N <- 500
  params <- paste("X", 1:3, sep = "")
  mat <- sobol_matrices(N = N, params = params)
  Y <- ishigami_Fun(mat)
  expect_error(
    sobol_convergence(matrices = c("A", "B", "AB"), Y = Y, N = N,
                      sub.sample = c(300, 100, 200), params = params,
                      first = "saltelli", total = "jansen",
                      order = "first", plot.order = "first"),
    "sorted"
  )
})

test_that("sobol_convergence stops when sub.sample exceeds N", {
  N <- 200
  params <- paste("X", 1:3, sep = "")
  mat <- sobol_matrices(N = N, params = params)
  Y <- ishigami_Fun(mat)
  expect_error(
    sobol_convergence(matrices = c("A", "B", "AB"), Y = Y, N = N,
                      sub.sample = c(100, 300), params = params,
                      first = "saltelli", total = "jansen",
                      order = "first", plot.order = "first"),
    regexp = NULL
  )
})

# ---- vars_matrices and vars_to ----

test_that("vars_matrices returns a matrix with correct number of columns", {
  params <- paste("X", 1:3, sep = "")
  mat <- vars_matrices(star.centers = 10, params = params, h = 0.1)
  expect_true(is.matrix(mat))
  expect_equal(ncol(mat), length(params))
})

test_that("vars_to returns a list with results of class vars", {
  params <- paste("X", 1:3, sep = "")
  mat <- vars_matrices(star.centers = 10, params = params, h = 0.1)
  Y <- ishigami_Fun(mat)
  ind <- vars_to(Y = Y, star.centers = 10, params = params, h = 0.1)
  expect_s3_class(ind, "vars")
  expect_s3_class(ind$results, "data.table")
  expect_equal(nrow(ind$results), length(params))
})

# ---- discrepancy_ersatz ----

test_that("discrepancy_ersatz returns a data.table with one row per parameter", {
  N <- 2^7
  params <- paste("X", 1:8, sep = "")
  mat <- sobol_matrices(N = N, params = params, matrices = "A")
  Y <- sobol_Fun(mat)
  ind <- discrepancy_ersatz(mat = mat, Y = Y, params = params)
  expect_s3_class(ind, "data.table")
  expect_equal(nrow(ind), length(params))
  expect_true("value" %in% colnames(ind))
})

# ---- metafunction ----

test_that("metafunction returns a numeric vector and is reproducible with epsilon", {
  N <- 50
  params <- paste("X", 1:10, sep = "")
  mat <- sobol_matrices(N = N, params = params)
  Y1 <- metafunction(mat, epsilon = 42)
  Y2 <- metafunction(mat, epsilon = 42)
  expect_type(Y1, "double")
  expect_length(Y1, nrow(mat))
  expect_equal(Y1, Y2)
})

# ---- plotting ----

test_that("plot.sensobol returns a ggplot object", {
  N <- 100
  params <- paste("X", 1:3, sep = "")
  mat <- sobol_matrices(N = N, params = params)
  Y <- ishigami_Fun(mat)
  ind <- sobol_indices(Y = Y, N = N, params = params)
  gg <- plot(ind)
  expect_s3_class(gg, "ggplot")
})
