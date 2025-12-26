.assert_is_numeric_matrix <- function(x, name) {
  if (!is.matrix(x) || !is.numeric(x)) {
    stop(name, " must be a numeric matrix.", call. = FALSE)
  }
  if (!all(is.finite(x))) {
    stop(name, " must contain only finite values.", call. = FALSE)
  }
  invisible(TRUE)
}

.assert_is_numeric_vector <- function(x, len, name) {
  if (is.null(x)) {
    stop(name, " must not be NULL.", call. = FALSE)
  }
  if (!is.numeric(x) || length(x) != len) {
    stop(name, " must be a numeric vector of length ", len, ".", call. = FALSE)
  }
  if (!all(is.finite(x))) {
    stop(name, " must contain only finite values.", call. = FALSE)
  }
  invisible(TRUE)
}

.assert_thresholds <- function(thresholds, d) {
  if (!is.list(thresholds) || length(thresholds) != d) {
    stop("thresholds must be a list of length ", d, ".", call. = FALSE)
  }
  for (j in seq_len(d)) {
    thr <- thresholds[[j]]
    if (is.null(thr)) {
      stop("thresholds[[", j, "]] must not be NULL.", call. = FALSE)
    }
    if (!is.numeric(thr) || length(thr) < 1) {
      stop("thresholds[[", j, "]] must be a numeric vector of length >= 1.", call. = FALSE)
    }
    if (!all(is.finite(thr))) {
      stop("thresholds[[", j, "]] contains non-finite values.", call. = FALSE)
    }
    if (is.unsorted(thr, strictly = TRUE)) {
      stop("thresholds[[", j, "]] must be strictly increasing.", call. = FALSE)
    }
  }
  invisible(TRUE)
}

.estimate_thresholds_from_data <- function(y, eps = 1e-6) {
  if (is.factor(y)) {
    y2 <- y
  } else if (is.logical(y)) {
    y2 <- factor(as.integer(y), levels = c(0, 1))
  } else if (is.numeric(y) || is.integer(y)) {
    y2 <- factor(y)
  } else {
    stop("Unsupported variable type for threshold estimation.", call. = FALSE)
  }

  tab <- table(y2, useNA = "no")
  if (any(tab == 0)) {
    stop("Some categories have zero count; cannot estimate thresholds reliably.", call. = FALSE)
  }

  p <- cumsum(as.numeric(tab) / sum(tab))
  if (length(p) < 2) {
    stop("Need at least 2 categories to estimate thresholds.", call. = FALSE)
  }

  cuts <- stats::qnorm(.clip_unit(p[seq_len(length(p) - 1)], eps = eps))
  as.numeric(cuts)
}

.make_ordered_from_codes <- function(code, lev) {
  if (is.null(lev)) {
    return(as.integer(code))
  }
  factor(code, levels = seq_along(lev), labels = lev, ordered = TRUE)
}

.fit_factor_model <- function(z, n_factors = 1, rotation = "none", eps = 1e-6) {
  p <- ncol(z)
  if (p < 2) {
    stop("Need at least 2 variables to fit a factor model.", call. = FALSE)
  }
  if (as.integer(n_factors) >= p) {
    stop("n_factors must be < number of variables.", call. = FALSE)
  }

  if (p < 3) {
    if (as.integer(n_factors) != 1L) {
      stop("For 2 variables, only n_factors = 1 is supported.", call. = FALSE)
    }

    s1 <- stats::sd(z[, 1])
    s2 <- stats::sd(z[, 2])
    if (!is.finite(s1) || !is.finite(s2) || s1 <= 0 || s2 <= 0) {
      stop("Selected vars have zero variance after transformation; cannot fit a 2-variable factor model.", call. = FALSE)
    }

    r <- suppressWarnings(stats::cor(z[, 1], z[, 2], use = "pairwise.complete.obs"))
    if (!is.finite(r)) {
      stop("Cannot compute correlation for 2-variable factor fit.", call. = FALSE)
    }

    r <- max(min(r, 1 - eps), -1 + eps)
    a <- sqrt(abs(r))

    list(
      loadings = matrix(c(a, sign(r) * a), nrow = 2, ncol = 1),
      uniquenesses = rep(1 - abs(r), 2),
      method = "analytic_2var",
      correlation = r
    )
  } else {
    fit <- tryCatch(
      stats::factanal(x = z, factors = as.integer(n_factors), rotation = rotation, scores = "none"),
      error = function(e) stop("factanal failed: ", e$message, call. = FALSE)
    )

    list(
      loadings = as.matrix(fit$loadings),
      uniquenesses = as.numeric(fit$uniquenesses),
      method = "factanal",
      factor_fit = fit
    )
  }
}

simdesign_factor_copula <- function(
  loadings,
  uniq_var = NULL,
  dist,
  names_final = NULL,
  process_final = list(),
  name = "Factor Gaussian copula design",
  eps = 1e-6,
  ...
) {
  .assert_is_numeric_matrix(loadings, "loadings")
  d <- nrow(loadings)
  k <- ncol(loadings)

  if (d < 2 || k < 1) {
    stop("loadings must have nrow >= 2 and ncol >= 1.", call. = FALSE)
  }

  if (is.null(uniq_var)) {
    uniq_var <- rep(1, d)
  }
  .assert_is_numeric_vector(uniq_var, d, "uniq_var")
  if (any(uniq_var <= 0)) {
    stop("uniq_var must be strictly positive.", call. = FALSE)
  }

  if (!is.numeric(eps) || length(eps) != 1 || eps <= 0 || eps >= 0.5) {
    stop("eps must be a single numeric in (0, 0.5).", call. = FALSE)
  }

  .assert_is_function_list(dist, d, "dist")

  sd_j <- sqrt(rowSums(loadings^2) + uniq_var)
  if (!all(is.finite(sd_j)) || any(sd_j <= 0)) {
    stop("Invalid loadings/uniq_var: cannot standardize latent variables.", call. = FALSE)
  }

  generator <- function(n_obs, ...) {
    f <- matrix(stats::rnorm(n_obs * k), nrow = n_obs, ncol = k)
    e <- matrix(stats::rnorm(n_obs * d), nrow = n_obs, ncol = d)
    e <- sweep(e, 2, sqrt(uniq_var), "*")
    z <- f %*% t(loadings) + e
    z <- sweep(z, 2, sd_j, "/")
    u <- stats::pnorm(z)
    .clip_unit(u, eps = eps)
  }

  transform_initial <- function(u) {
    if (!is.matrix(u) && !is.data.frame(u)) {
      stop("Internal error: generator did not return a 2D object.", call. = FALSE)
    }
    u <- as.matrix(u)
    if (ncol(u) != d) {
      stop("Internal error: U has unexpected number of columns.", call. = FALSE)
    }
    x <- lapply(seq_len(d), function(j) dist[[j]](u[, j]))
    as.data.frame(x, optional = TRUE, stringsAsFactors = FALSE)
  }

  simdata::simdesign(
    generator = generator,
    transform_initial = transform_initial,
    names_final = names_final,
    process_final = process_final,
    name = name,
    ...
  )
}

simdesign_factor_copula_from_data <- function(
  data,
  vars,
  n_factors = 1,
  rotation = "none",
  qtype = 8,
  eps = 1e-6,
  name = "Factor Gaussian copula design (fit from data)",
  ...
) {
  if (!is.data.frame(data)) {
    stop("data must be a data.frame.", call. = FALSE)
  }
  if (!is.character(vars) || length(vars) < 2) {
    stop("vars must be a character vector of length >= 2.", call. = FALSE)
  }
  if (!all(vars %in% names(data))) {
    miss <- setdiff(vars, names(data))
    stop("vars missing from data: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  if (!is.numeric(n_factors) || length(n_factors) != 1 || n_factors < 1) {
    stop("n_factors must be a single numeric >= 1.", call. = FALSE)
  }

  df_in <- data[, vars, drop = FALSE]
  for (nm in vars) {
    if (!is.numeric(df_in[[nm]])) {
      stop("All vars must be numeric for MVP. Non-numeric: ", nm, call. = FALSE)
    }
  }

  xmat <- as.matrix(df_in)
  if (!all(is.finite(xmat))) {
    stop("Data contains non-finite values in selected vars.", call. = FALSE)
  }

  dist <- lapply(df_in, .make_empirical_q, qtype = qtype, eps = eps)
  u_hat <- copula::pobs(xmat)
  u_hat <- .clip_unit(u_hat, eps = eps)
  z <- stats::qnorm(u_hat)

  fit <- .fit_factor_model(z, n_factors = n_factors, rotation = rotation, eps = eps)

  loadings <- as.matrix(fit$loadings)
  uniq_var <- as.numeric(fit$uniquenesses)
  factor_fit <- if (!is.null(fit$factor_fit)) fit$factor_fit else fit

  simdesign_factor_copula(
    loadings = loadings,
    uniq_var = uniq_var,
    dist = dist,
    names_final = vars,
    name = name,
    eps = eps,
    factor_fit = factor_fit,
    margins = "empirical",
    ...
  )
}

simdesign_latent_threshold_gaussian <- function(
  loadings,
  uniq_var = NULL,
  thresholds,
  levels = NULL,
  names_final = NULL,
  process_final = list(),
  name = "Latent Gaussian threshold design",
  eps = 1e-6,
  ...
) {
  .assert_is_numeric_matrix(loadings, "loadings")
  d <- nrow(loadings)
  k <- ncol(loadings)

  if (d < 2 || k < 1) {
    stop("loadings must have nrow >= 2 and ncol >= 1.", call. = FALSE)
  }

  if (is.null(uniq_var)) {
    uniq_var <- rep(1, d)
  }
  .assert_is_numeric_vector(uniq_var, d, "uniq_var")
  if (any(uniq_var <= 0)) {
    stop("uniq_var must be strictly positive.", call. = FALSE)
  }

  if (!is.numeric(eps) || length(eps) != 1 || eps <= 0 || eps >= 0.5) {
    stop("eps must be a single numeric in (0, 0.5).", call. = FALSE)
  }

  .assert_thresholds(thresholds, d)

  if (!is.null(levels)) {
    if (!is.list(levels) || length(levels) != d) {
      stop("levels must be NULL or a list of length ", d, ".", call. = FALSE)
    }
    for (j in seq_len(d)) {
      lev <- levels[[j]]
      if (!is.null(lev)) {
        if (!is.character(lev)) {
          stop("levels[[", j, "]] must be NULL or a character vector.", call. = FALSE)
        }
        if (length(lev) != length(thresholds[[j]]) + 1) {
          stop("levels[[", j, "]] must have length length(thresholds[[j]]) + 1.", call. = FALSE)
        }
      }
    }
  }

  sd_j <- sqrt(rowSums(loadings^2) + uniq_var)
  if (!all(is.finite(sd_j)) || any(sd_j <= 0)) {
    stop("Invalid loadings/uniq_var: cannot standardize latent variables.", call. = FALSE)
  }

  generator <- function(n_obs, ...) {
    f <- matrix(stats::rnorm(n_obs * k), nrow = n_obs, ncol = k)
    e <- matrix(stats::rnorm(n_obs * d), nrow = n_obs, ncol = d)
    e <- sweep(e, 2, sqrt(uniq_var), "*")
    z <- f %*% t(loadings) + e
    sweep(z, 2, sd_j, "/")
  }

  transform_initial <- function(z) {
    if (!is.matrix(z) && !is.data.frame(z)) {
      stop("Internal error: generator did not return a 2D object.", call. = FALSE)
    }
    z <- as.matrix(z)
    if (ncol(z) != d) {
      stop("Internal error: latent Z has unexpected number of columns.", call. = FALSE)
    }

    out <- vector("list", d)
    for (j in seq_len(d)) {
      thr <- thresholds[[j]]
      code <- findInterval(z[, j], thr) + 1L
      lev <- if (is.null(levels)) NULL else levels[[j]]
      out[[j]] <- .make_ordered_from_codes(code, lev)
    }

    as.data.frame(out, optional = TRUE, stringsAsFactors = FALSE)
  }

  simdata::simdesign(
    generator = generator,
    transform_initial = transform_initial,
    names_final = names_final,
    process_final = process_final,
    name = name,
    ...
  )
}

simdesign_latent_threshold_gaussian_from_data <- function(
  data,
  vars,
  n_factors = 1,
  rotation = "none",
  tie_jitter = 1e-8,
  eps = 1e-6,
  name = "Latent Gaussian threshold design (fit from data)",
  ...
) {
  if (!is.data.frame(data)) {
    stop("data must be a data.frame.", call. = FALSE)
  }
  if (!is.character(vars) || length(vars) < 2) {
    stop("vars must be a character vector of length >= 2.", call. = FALSE)
  }
  if (!all(vars %in% names(data))) {
    miss <- setdiff(vars, names(data))
    stop("vars missing from data: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  if (!is.numeric(n_factors) || length(n_factors) != 1 || n_factors < 1) {
    stop("n_factors must be a single numeric >= 1.", call. = FALSE)
  }
  if (!is.numeric(tie_jitter) || length(tie_jitter) != 1 || tie_jitter < 0) {
    stop("tie_jitter must be a single non-negative numeric.", call. = FALSE)
  }

  df_in <- data[, vars, drop = FALSE]
  n <- nrow(df_in)
  if (n < 5) {
    stop("Not enough rows in data to fit a latent threshold model.", call. = FALSE)
  }

  thresholds <- vector("list", length(vars))
  levels_out <- vector("list", length(vars))

  x_num <- vector("list", length(vars))
  for (j in seq_along(vars)) {
    y <- df_in[[j]]

    thresholds[[j]] <- .estimate_thresholds_from_data(y, eps = eps)

    if (is.factor(y)) {
      levels_out[[j]] <- levels(y)
      x_num[[j]] <- as.integer(y)
    } else if (is.logical(y)) {
      levels_out[[j]] <- NULL
      x_num[[j]] <- as.integer(y)
    } else if (is.numeric(y) || is.integer(y)) {
      levels_out[[j]] <- NULL
      x_num[[j]] <- as.numeric(y)
    } else {
      stop("Unsupported variable type: ", vars[[j]], call. = FALSE)
    }
  }

  xmat <- do.call(cbind, x_num)
  if (!all(is.finite(xmat))) {
    stop("Selected vars contain non-finite values.", call. = FALSE)
  }

  if (tie_jitter > 0) {
    xmat <- xmat + matrix(stats::runif(n * ncol(xmat), min = -tie_jitter, max = tie_jitter), nrow = n)
  }

  u_hat <- copula::pobs(xmat)
  u_hat <- .clip_unit(u_hat, eps = eps)
  z <- stats::qnorm(u_hat)

  fit <- .fit_factor_model(z, n_factors = n_factors, rotation = rotation, eps = eps)

  loadings <- as.matrix(fit$loadings)
  uniq_var <- as.numeric(fit$uniquenesses)
  factor_fit <- if (!is.null(fit$factor_fit)) fit$factor_fit else fit

  simdesign_latent_threshold_gaussian(
    loadings = loadings,
    uniq_var = uniq_var,
    thresholds = thresholds,
    levels = levels_out,
    names_final = vars,
    name = name,
    eps = eps,
    factor_fit = factor_fit,
    margins = "threshold",
    ...
  )
}