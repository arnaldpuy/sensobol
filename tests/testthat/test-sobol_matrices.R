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

# ---- scrambling: regression (default == "none") ----

test_that("scrambling = 'none' matches the pre-existing default", {
  N <- 50; params <- paste0("X", 1:3)
  m_old  <- sobol_matrices(N = N, params = params)
  m_none <- sobol_matrices(N = N, params = params, scrambling = "none")
  expect_identical(m_old, m_none)
})

# ---- scrambling: Cranley-Patterson shift ----

test_that("scrambling = 'shift' with same seed is reproducible", {
  N <- 100; params <- paste0("X", 1:3)
  m1 <- sobol_matrices(N = N, params = params, scrambling = "shift", seed = 7)
  m2 <- sobol_matrices(N = N, params = params, scrambling = "shift", seed = 7)
  expect_equal(m1, m2)
})

test_that("scrambling = 'shift' with different seeds differs", {
  N <- 100; params <- paste0("X", 1:3)
  m1 <- sobol_matrices(N = N, params = params, scrambling = "shift", seed = 7)
  m2 <- sobol_matrices(N = N, params = params, scrambling = "shift", seed = 99)
  expect_false(isTRUE(all.equal(m1, m2)))
})

test_that("scrambling = 'shift' output stays in [0, 1]", {
  N <- 100; params <- paste0("X", 1:3)
  m <- sobol_matrices(N = N, params = params, scrambling = "shift", seed = 1)
  expect_true(all(m >= 0 & m <= 1))
})

test_that("shift column means converge to 0.5 and variances to 1/12", {
  N <- 2^13; params <- paste0("X", 1:3)
  m <- sobol_matrices(N = N, params = params, scrambling = "shift", seed = 1)
  A <- m[1:N, ]
  expect_equal(unname(colMeans(A)), rep(0.5, 3), tolerance = 0.02)
  expect_equal(unname(apply(A, 2, var)), rep(1 / 12, 3), tolerance = 0.02)
})

# ---- scrambling: in-house Sobol' + Owen ----

test_that("scrambling = 'owen' with same seed is reproducible", {
  N <- 100; params <- paste0("X", 1:3)
  m1 <- sobol_matrices(N = N, params = params, scrambling = "owen", seed = 7)
  m2 <- sobol_matrices(N = N, params = params, scrambling = "owen", seed = 7)
  expect_equal(m1, m2)
})

test_that("scrambling = 'owen' with different seeds differs", {
  N <- 100; params <- paste0("X", 1:3)
  m1 <- sobol_matrices(N = N, params = params, scrambling = "owen", seed = 7)
  m2 <- sobol_matrices(N = N, params = params, scrambling = "owen", seed = 99)
  expect_false(isTRUE(all.equal(m1, m2)))
})

test_that("scrambling = 'owen' output stays in [0, 1]", {
  N <- 100; params <- paste0("X", 1:3)
  m <- sobol_matrices(N = N, params = params, scrambling = "owen", seed = 1)
  expect_true(all(m >= 0 & m <= 1))
})

test_that("owen column means converge to 0.5 and variances to 1/12", {
  N <- 2^13; params <- paste0("X", 1:3)
  m <- sobol_matrices(N = N, params = params, scrambling = "owen", seed = 1)
  A <- m[1:N, ]
  expect_equal(unname(colMeans(A)), rep(0.5, 3), tolerance = 0.02)
  expect_equal(unname(apply(A, 2, var)), rep(1 / 12, 3), tolerance = 0.02)
})

test_that("owen QMC integral converges faster than plain MC", {
  # Integral of sum(x_i^2) over [0,1]^d is d/3.
  N <- 2^12; d <- 5
  mat <- sobol_matrices(N = N, params = paste0("X", 1:d),
                        matrices = "A", scrambling = "owen", seed = 1)
  qmc_err <- abs(mean(rowSums(mat ^ 2)) - d / 3)
  set.seed(1)
  mc_err <- abs(mean(replicate(N, sum(stats::runif(d) ^ 2))) - d / 3)
  expect_lt(qmc_err, mc_err)
})

# ---- scrambling: integration with other features ----

test_that("scrambling integrates with groups and order = 'second'", {
  N <- 100; params <- paste0("X", 1:4)
  g <- list(g1 = "X1", g23 = c("X2", "X3"), g4 = "X4")

  for (sc in c("shift", "owen")) {
    m <- sobol_matrices(N = N, params = params, scrambling = sc, seed = 1,
                        groups = g, order = "second")
    G <- length(g)
    expect_equal(nrow(m), N * (G + 2 + choose(G, 2)),
                 info = paste("scrambling =", sc))
  }
})

test_that("Ishigami Sobol' indices converge under all scrambling modes", {
  N <- 2^13; params <- paste0("X", 1:3)
  for (sc in c("none", "shift", "owen")) {
    mat <- sobol_matrices(N = N, params = params, scrambling = sc, seed = 1)
    Y <- ishigami_Fun(mat)
    ind <- sobol_indices(Y = Y, N = N, params = params)
    si <- ind$results[sensitivity == "Si", original]
    # S1 ~ 0.38, S2 ~ 0.001, S3 ~ 0 for sensobol's Ishigami parametrization
    expect_lt(abs(si[1] - 0.383), 0.05, label = paste(sc, "S1"))
    expect_lt(abs(si[2] - 0.001), 0.05, label = paste(sc, "S2"))
    expect_lt(abs(si[3] - 0.000), 0.05, label = paste(sc, "S3"))
  }
})

# ---- scrambling: error and validation paths ----

test_that("scrambling = 'owen' errors when dim > 250", {
  expect_error(
    sobol_matrices(N = 10, params = paste0("X", 1:200),
                   scrambling = "owen", seed = 1),
    "up to 250 dimensions"
  )
})

test_that("scrambling with type != 'QRN' emits a warning and falls back", {
  expect_warning(
    sobol_matrices(N = 50, params = paste0("X", 1:3),
                   type = "R", scrambling = "shift", seed = 1),
    "only used when type"
  )
})

test_that("invalid 'seed' is rejected", {
  expect_error(
    sobol_matrices(N = 50, params = paste0("X", 1:3),
                   scrambling = "shift", seed = "abc"),
    "'seed' must be a single integer"
  )
})

test_that("invalid 'scrambling' value is rejected", {
  expect_error(
    sobol_matrices(N = 50, params = paste0("X", 1:3),
                   scrambling = "bogus", seed = 1),
    "'arg'"
  )
})

test_that("scrambling does not clobber the user's global RNG state", {
  set.seed(123); a <- stats::runif(1)
  set.seed(123); .junk <- sobol_matrices(N = 50, params = paste0("X", 1:3),
                                         scrambling = "shift", seed = 7)
  b <- stats::runif(1)
  expect_equal(a, b)

  set.seed(123); .junk <- sobol_matrices(N = 50, params = paste0("X", 1:3),
                                         scrambling = "owen", seed = 7)
  d <- stats::runif(1)
  expect_equal(a, d)
})
