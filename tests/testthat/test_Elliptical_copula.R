test_that("simdesign_gaussian_copula generates correct shape and names", {
  skip_if_not_installed("simdata")
  skip_if_not_installed("copula")

  dist <- list(
    function(u) qnorm(u, mean = 0, sd = 1),
    function(u) qexp(u, rate = 0.5),
    function(u) qbeta(u, shape1 = 2, shape2 = 5)
  )

  dsgn <- simdesign_gaussian_copula(
    dist = dist,
    rho = 0.3,
    dim = 3,
    structure = "ex",
    names_final = c("x", "t", "b")
  )

  x <- simdata::simulate_data(dsgn, n_obs = 100, seed = 1)
  expect_true(is.data.frame(x) || is.matrix(x))
  expect_equal(nrow(x), 100)
  expect_equal(colnames(x), c("x", "t", "b"))
})

test_that("simulate_data is reproducible with seed", {
  skip_if_not_installed("simdata")
  skip_if_not_installed("copula")

  dist <- list(
    function(u) qnorm(u),
    function(u) qnorm(u),
    function(u) qnorm(u)
  )

  dsgn <- simdesign_gaussian_copula(
    dist = dist,
    rho = 0.5,
    dim = 3,
    structure = "ex",
    names_final = c("a", "b", "c")
  )

  x1 <- simdata::simulate_data(dsgn, n_obs = 200, seed = 999)
  x2 <- simdata::simulate_data(dsgn, n_obs = 200, seed = 999)
  expect_equal(x1, x2)
})

test_that("fit-from-data reproduces rank correlation approximately (Gaussian, ex)", {
  skip_if_not_installed("simdata")
  skip_if_not_installed("copula")

  set.seed(202401)
  n_pilot <- 1500
  rho_true <- 0.55
  cop_true <- copula::normalCopula(param = rho_true, dim = 3, dispstr = "ex")
  u <- copula::rCopula(n_pilot, cop_true)

  pilot <- data.frame(
    x1 = qexp(u[, 1], rate = 0.2),
    x2 = qnorm(u[, 2], mean = 10, sd = 2),
    x3 = qlnorm(u[, 3], meanlog = 0, sdlog = 0.5)
  )

  dsgn <- simdesign_elliptical_copula_from_data(
    data = pilot,
    vars = c("x1", "x2", "x3"),
    family = "gaussian",
    structure = "ex",
    fit_method = "itau"
  )

  sim <- simdata::simulate_data(dsgn, n_obs = 4000, seed = 42)

  r_pilot <- cor(copula::pobs(as.matrix(pilot)), method = "spearman")
  r_sim <- cor(copula::pobs(as.matrix(sim)), method = "spearman")

  expect_true(max(abs(r_pilot - r_sim)) < 0.10)
})