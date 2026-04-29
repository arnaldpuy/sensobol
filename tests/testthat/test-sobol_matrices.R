context("test-sobol_matrices")
library(testthat)

test_that("Output is a matrix", {
  N <- 10
  params <- paste("X", 1:3, sep = "")
  mat <- sobol_matrices(N = N, params = params)
  expect_true(is.matrix(mat))
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

# ---- groups: regression (groups = NULL is identical to old behaviour) ----

test_that("groups = NULL produces the same matrix as the ungrouped default", {
  N <- 50
  params <- paste("X", 1:4, sep = "")
  set.seed(1)
  m1 <- sobol_matrices(N = N, params = params, type = "R")
  set.seed(1)
  m2 <- sobol_matrices(N = N, params = params, type = "R", groups = NULL)
  expect_identical(m1, m2)
})

test_that("trivial partition (one parameter per group) matches ungrouped", {
  N <- 50
  params <- paste("X", 1:4, sep = "")
  trivial <- as.list(params); names(trivial) <- params
  set.seed(2)
  m1 <- sobol_matrices(N = N, params = params, type = "R")
  set.seed(2)
  m2 <- sobol_matrices(N = N, params = params, type = "R", groups = trivial)
  expect_equal(m1, m2)
})

# ---- groups: row counts ----

test_that("groups reduces the number of conditional blocks (AB)", {
  N <- 10
  params <- paste("X", 1:4, sep = "")
  groups <- list(g1 = "X1", g2 = c("X2", "X3"), g3 = "X4")
  mat <- sobol_matrices(N = N, params = params, groups = groups,
                        matrices = c("A", "B", "AB"))
  expect_equal(nrow(mat), N * (length(groups) + 2))
})

test_that("groups + BA produces 2*G conditional blocks", {
  N <- 10
  params <- paste("X", 1:4, sep = "")
  groups <- list(g1 = "X1", g2 = c("X2", "X3"), g3 = "X4")
  mat <- sobol_matrices(N = N, params = params, groups = groups,
                        matrices = c("A", "B", "AB", "BA"))
  expect_equal(nrow(mat), N * (2 * length(groups) + 2))
})

test_that("groups + second-order adds choose(G, 2) blocks", {
  N <- 10
  params <- paste("X", 1:4, sep = "")
  groups <- list(g1 = "X1", g2 = c("X2", "X3"), g3 = "X4")
  G <- length(groups)
  mat <- sobol_matrices(N = N, params = params, groups = groups,
                        matrices = c("A", "B", "AB"), order = "second")
  expect_equal(nrow(mat), N * (G + 2) + N * choose(G, 2))
})

# ---- groups: AB block correctness ----

test_that("AB conditional block swaps every column of the group", {
  N <- 5
  params <- paste("X", 1:4, sep = "")
  groups <- list(g1 = "X1", g23 = c("X2", "X3"), g4 = "X4")
  mat <- sobol_matrices(N = N, params = params, groups = groups, type = "R",
                        matrices = c("A", "B", "AB"))
  A <- mat[1:N, ]
  B <- mat[(N + 1):(2 * N), ]
  # AB blocks start at row 2N+1; one block of N rows per group, in list order.
  AB_g1  <- mat[(2 * N + 1):(3 * N), ]
  AB_g23 <- mat[(3 * N + 1):(4 * N), ]
  AB_g4  <- mat[(4 * N + 1):(5 * N), ]

  expect_equal(AB_g1[,  "X1"], B[, "X1"])
  expect_equal(AB_g1[,  c("X2", "X3", "X4")], A[, c("X2", "X3", "X4")])

  expect_equal(AB_g23[, c("X2", "X3")], B[, c("X2", "X3")])
  expect_equal(AB_g23[, c("X1", "X4")], A[, c("X1", "X4")])

  expect_equal(AB_g4[,  "X4"], B[, "X4"])
  expect_equal(AB_g4[,  c("X1", "X2", "X3")], A[, c("X1", "X2", "X3")])
})

# ---- groups: vector form ----

test_that("vector form of groups is equivalent to the list form", {
  N <- 20
  params <- paste("X", 1:4, sep = "")
  groups_list <- list(g1 = "X1", g2 = c("X2", "X3"), g3 = "X4")
  groups_vec  <- c("g1", "g2", "g2", "g3")
  set.seed(3)
  m_list <- sobol_matrices(N = N, params = params, type = "R", groups = groups_list)
  set.seed(3)
  m_vec  <- sobol_matrices(N = N, params = params, type = "R", groups = groups_vec)
  expect_equal(m_list, m_vec)
})

# ---- groups: validation errors ----

test_that("groups must be a strict partition (missing parameter)", {
  expect_error(
    sobol_matrices(N = 10, params = paste0("X", 1:4),
                   groups = list(g1 = "X1", g2 = c("X2", "X3"))),
    "strict partition"
  )
})

test_that("groups must be a strict partition (duplicated parameter)", {
  expect_error(
    sobol_matrices(N = 10, params = paste0("X", 1:4),
                   groups = list(g1 = "X1", g2 = c("X1", "X2"), g3 = c("X3", "X4"))),
    "strict partition"
  )
})

test_that("groups list must be named", {
  expect_error(
    sobol_matrices(N = 10, params = paste0("X", 1:3),
                   groups = list("X1", c("X2", "X3"))),
    "named list"
  )
})

test_that("groups vector must have the same length as params", {
  expect_error(
    sobol_matrices(N = 10, params = paste0("X", 1:3),
                   groups = c("g1", "g2")),
    "same length as 'params'"
  )
})
