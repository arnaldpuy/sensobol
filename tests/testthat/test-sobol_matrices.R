context("test-sobol_matrices")
library(testthat)

test_that("Output is a matrix", {
  N <- 10
  params <- paste("X", 1:3, sep = "")
  mat <- sobol_matrices(N = N, params = params)
  expect_is(mat, "matrix")
})

test_that("Sample size for first and total-order indices 1", {
  N <- 10
  params <- paste("X", 1:3, sep = "")
  matrices <- c("A", "B", "AB")
  mat <- sobol_matrices(N = N, params = params, matrices = matrices)
  expect_equal(nrow(mat), N * (length(params) + 2))
})

test_that("Sample size for first and total-order indices 2", {
  N <- 10
  params <- paste("X", 1:3, sep = "")
  matrices <- c("A", "B", "AB", "BA")
  mat <- sobol_matrices(N = N, params = params, matrices = matrices)
  expect_equal(nrow(mat), N * (2 * length(params) + 2))
})

test_that("Sample size for first and total-order indices 3", {
  N <- 10
  params <- paste("X", 1:3, sep = "")
  matrices <- c("A", "AB", "BA")
  mat <- sobol_matrices(N = N, params = params, matrices = matrices)
  expect_equal(nrow(mat), N * (2 * length(params) + 1))
})


test_that("Sample size for first, total and second-order indices", {
  N <- 10
  params <- paste("X", 1:3, sep = "")
  matrices <- c("A", "B", "AB")
  order <- "second"
  mat <- sobol_matrices(N = N, params = params, matrices = matrices, order = order)
  expect_equal(nrow(mat), N * (length(params) + 2) +
                 N * ncol(utils::combn(params, 2)))
})

test_that("Sample size for first, total, second and third-order indices", {
  N <- 10
  params <- paste("X", 1:3, sep = "")
  matrices <- c("A", "B", "AB")
  order <- "third"
  mat <- sobol_matrices(N = N, params = params, matrices = matrices, order = order)
  expect_equal(nrow(mat), N * (length(params) + 2) +
                 N * ncol(utils::combn(params, 2)) +
                 N * ncol(utils::combn(params, 3)))
})
