test_that("simdesign_factor_copula generates correct shape and names", {
  skip_if_not_installed("simdata")
  skip_if_not_installed("copula")

  loadings <- matrix(c(0.9, 0.6, 0.7), nrow = 3, ncol = 1)
  uniq <- c(0.4, 0.7, 0.6)

  dist <- list(
    function(u) qnorm(u, mean = 0, sd = 1),
    function(u) qexp(u, rate = 0.5),
    function(u) qbeta(u, shape1 = 2, shape2 = 5)
  )

  dsgn <- simdesign_factor_copula(
    loadings = loadings,
    uniq_var = uniq,
    dist = dist,
    names_final = c("x", "t", "b")
  )

  x <- simdata::simulate_data(dsgn, n_obs = 200, seed = 1)
  expect_true(is.data.frame(x) || is.matrix(x))
  expect_equal(nrow(x), 200)
  expect_equal(colnames(x), c("x", "t", "b"))
})

test_that("simulate_data is reproducible with seed (factor)", {
  skip_if_not_installed("simdata")
  skip_if_not_installed("copula")

  loadings <- matrix(c(0.8, 0.5, 0.4), nrow = 3, ncol = 1)
  uniq <- c(0.6, 0.7, 0.8)

  dist <- list(function(u) qnorm(u), function(u) qnorm(u), function(u) qnorm(u))

  dsgn <- simdesign_factor_copula(
    loadings = loadings,
    uniq_var = uniq,
    dist = dist,
    names_final = c("a", "b", "c")
  )

  x1 <- simdata::simulate_data(dsgn, n_obs = 300, seed = 999)
  x2 <- simdata::simulate_data(dsgn, n_obs = 300, seed = 999)
  expect_equal(x1, x2)
})

test_that("simdesign_factor_copula_from_data preserves rank correlation approximately", {
  skip_if_not_installed("simdata")
  skip_if_not_installed("copula")

  loadings0 <- matrix(c(0.9, 0.7, 0.6), nrow = 3, ncol = 1)
  uniq0 <- c(0.5, 0.6, 0.7)

  dist0 <- list(
    function(u) qexp(u, rate = 0.2),
    function(u) qnorm(u, mean = 10, sd = 2),
    function(u) qlnorm(u, meanlog = 0, sdlog = 0.5)
  )

  dsgn0 <- simdesign_factor_copula(
    loadings = loadings0,
    uniq_var = uniq0,
    dist = dist0,
    names_final = c("x1", "x2", "x3")
  )

  pilot <- simdata::simulate_data(dsgn0, n_obs = 2000, seed = 202601)

  dsgn_fit <- simdesign_factor_copula_from_data(
    data = pilot,
    vars = c("x1", "x2", "x3"),
    n_factors = 1
  )

  sim <- simdata::simulate_data(dsgn_fit, n_obs = 4000, seed = 42)

  r_pilot <- cor(copula::pobs(as.matrix(pilot)), method = "spearman")
  r_sim <- cor(copula::pobs(as.matrix(sim)), method = "spearman")

  expect_true(max(abs(r_pilot - r_sim)) < 0.12)
})

test_that("simdesign_latent_threshold_gaussian produces requested marginals approximately", {
  skip_if_not_installed("simdata")

  loadings <- matrix(c(0.9, 0.6, 0.7), nrow = 3, ncol = 1)
  uniq <- c(0.4, 0.7, 0.6)

  p_y1 <- 0.30
  p_y2 <- 0.55
  p_y3 <- c(0.20, 0.50, 0.30)

  thresholds <- list(
    qnorm(p_y1),
    qnorm(p_y2),
    qnorm(c(p_y3[1], p_y3[1] + p_y3[2]))
  )

  dsgn <- simdesign_latent_threshold_gaussian(
    loadings = loadings,
    uniq_var = uniq,
    thresholds = thresholds,
    names_final = c("y1", "y2", "y3")
  )

  x <- simdata::simulate_data(dsgn, n_obs = 5000, seed = 202602)

  prop_y1 <- mean(x$y1 == 1L)
  prop_y2 <- mean(x$y2 == 1L)
  prop_y3 <- as.numeric(table(factor(x$y3, levels = 1:3))) / nrow(x)

  expect_true(abs(prop_y1 - p_y1) < 0.03)
  expect_true(abs(prop_y2 - p_y2) < 0.03)
  expect_true(max(abs(prop_y3 - p_y3)) < 0.03)
})

test_that("simdesign_latent_threshold_gaussian_from_data is reproducible with seed", {
  skip_if_not_installed("simdata")
  skip_if_not_installed("copula")

  loadings <- matrix(c(0.9, 0.6, 0.7), nrow = 3, ncol = 1)
  uniq <- c(0.4, 0.7, 0.6)

  thresholds <- list(qnorm(0.3), qnorm(0.55), qnorm(c(0.2, 0.7)))

  dsgn0 <- simdesign_latent_threshold_gaussian(
    loadings = loadings,
    uniq_var = uniq,
    thresholds = thresholds,
    names_final = c("y1", "y2", "y3")
  )

  dat <- simdata::simulate_data(dsgn0, n_obs = 1200, seed = 202603)

  set.seed(999)
  d1 <- simdesign_latent_threshold_gaussian_from_data(dat, vars = c("y1", "y2", "y3"), n_factors = 1)
  set.seed(999)
  d2 <- simdesign_latent_threshold_gaussian_from_data(dat, vars = c("y1", "y2", "y3"), n_factors = 1)

  x1 <- simdata::simulate_data(d1, n_obs = 800, seed = 123)
  x2 <- simdata::simulate_data(d2, n_obs = 800, seed = 123)

  expect_equal(x1, x2)
})