test_that("simdesign_vine_copula_from_data works (rvine) and preserves rank correlation approximately", {
  skip_if_not_installed("simdata")
  skip_if_not_installed("copula")
  skip_if_not_installed("VineCopula")

  set.seed(202503)
  n_pilot <- 1500

  cop_true <- copula::normalCopula(param = 0.6, dim = 4, dispstr = "ex")
  u <- copula::rCopula(n_pilot, cop_true)

  pilot <- data.frame(
    x1 = qexp(u[, 1], rate = 0.3),
    x2 = qnorm(u[, 2], mean = 10, sd = 2),
    x3 = qlnorm(u[, 3], meanlog = 0, sdlog = 0.4),
    x4 = qbeta(u[, 4], 2, 6)
  )

  dsgn <- simdesign_vine_copula_from_data(
    data = pilot,
    vars = c("x1", "x2", "x3", "x4"),
    vine_type = "rvine",
    familyset = 1:6,
    indeptest = TRUE,
    selectioncrit = "AIC",
    method = "mle"
  )

  sim <- simdata::simulate_data(dsgn, n_obs = 4000, seed = 42)
  expect_true(is.data.frame(sim) || is.matrix(sim))
  expect_equal(nrow(sim), 4000)
  expect_equal(colnames(sim), c("x1", "x2", "x3", "x4"))

  r_pilot <- cor(copula::pobs(as.matrix(pilot)), method = "spearman")
  r_sim <- cor(copula::pobs(as.matrix(sim)), method = "spearman")

  expect_true(max(abs(r_pilot - r_sim)) < 0.12)
})

test_that("simulate_data is reproducible with seed (vine)", {
  skip_if_not_installed("simdata")
  skip_if_not_installed("copula")
  skip_if_not_installed("VineCopula")

  set.seed(202504)
  pilot <- data.frame(
    x1 = rnorm(800),
    x2 = rexp(800),
    x3 = rlnorm(800)
  )

  dsgn <- simdesign_vine_copula_from_data(
    data = pilot,
    vars = c("x1", "x2", "x3"),
    vine_type = "cvine",
    familyset = 1:6,
    indeptest = TRUE,
    method = "mle"
  )

  x1 <- simdata::simulate_data(dsgn, n_obs = 600, seed = 999)
  x2 <- simdata::simulate_data(dsgn, n_obs = 600, seed = 999)
  expect_equal(x1, x2)
})