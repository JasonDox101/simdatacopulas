simdesign_empirical_copula <- function(
  u_data,
  dist,
  names_final = NULL,
  process_final = list(),
  replace = TRUE,
  jitter = 0,
  name = "Empirical copula design",
  eps = 1e-6,
  ...
) {
  if (!is.matrix(u_data) && !is.data.frame(u_data)) {
    stop("u_data must be a matrix or data.frame.", call. = FALSE)
  }
  u_data <- as.matrix(u_data)
  if (ncol(u_data) < 2) {
    stop("u_data must have at least 2 columns.", call. = FALSE)
  }
  if (nrow(u_data) < 2) {
    stop("u_data must have at least 2 rows.", call. = FALSE)
  }
  if (!all(is.finite(u_data))) {
    stop("u_data contains non-finite values.", call. = FALSE)
  }
  if (!is.logical(replace) || length(replace) != 1) {
    stop("replace must be a single logical.", call. = FALSE)
  }
  if (!is.numeric(jitter) || length(jitter) != 1 || jitter < 0) {
    stop("jitter must be a single non-negative numeric.", call. = FALSE)
  }
  if (!is.numeric(eps) || length(eps) != 1 || eps <= 0 || eps >= 0.5) {
    stop("eps must be a single numeric in (0, 0.5).", call. = FALSE)
  }

  dim <- ncol(u_data)
  .assert_is_function_list(dist, dim, "dist")

  u_data <- .clip_unit(u_data, eps = eps)
  n_u <- nrow(u_data)

  generator <- function(n_obs, ...) {
    if (!replace && n_obs > n_u) {
      stop("When replace=FALSE, n_obs cannot exceed nrow(u_data).", call. = FALSE)
    }
    idx <- sample.int(n_u, size = n_obs, replace = replace)
    u <- u_data[idx, , drop = FALSE]

    if (jitter > 0) {
      u <- u + matrix(
        stats::runif(n_obs * dim, min = -jitter, max = jitter),
        nrow = n_obs,
        ncol = dim
      )
    }

    .clip_unit(u, eps = eps)
  }

  transform_initial <- function(u) {
    if (!is.matrix(u) && !is.data.frame(u)) {
      stop("Internal error: generator did not return a 2D object.", call. = FALSE)
    }
    u <- as.matrix(u)
    if (ncol(u) != dim) {
      stop("Internal error: U has unexpected number of columns.", call. = FALSE)
    }

    x <- lapply(seq_len(dim), function(j) dist[[j]](u[, j]))
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

simdesign_empirical_copula_from_data <- function(
  data,
  vars,
  qtype = 8,
  replace = TRUE,
  jitter = 0,
  eps = 1e-6,
  name = "Empirical copula design (fit from data)",
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

  simdesign_empirical_copula(
    u_data = u_hat,
    dist = dist,
    names_final = vars,
    replace = replace,
    jitter = jitter,
    name = name,
    eps = eps,
    margins = "empirical",
    ...
  )
}