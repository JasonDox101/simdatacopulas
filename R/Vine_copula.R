.assert_vinecopula_installed <- function() {
  if (!requireNamespace("VineCopula", quietly = TRUE)) {
    stop("Package 'VineCopula' is required for vine copulas. Install it via install.packages('VineCopula').", call. = FALSE)
  }
  invisible(TRUE)
}

.get_vine_dim <- function(vine) {
  m <- tryCatch(vine$Matrix, error = function(e) NULL)
  if (!is.null(m) && is.matrix(m) && nrow(m) == ncol(m) && nrow(m) >= 2) {
    return(as.integer(nrow(m)))
  }
  stop("Unable to determine vine dimension (expecting vine$Matrix to be a square matrix).", call. = FALSE)
}

#' Create a simdata design using a vine copula (C-/D-/R-vine) as joint generator
#'
#' @param vine A vine model object from VineCopula (typically an RVineMatrix), e.g. returned by RVineStructureSelect().
#' @param dist List of marginal quantile functions. Each maps u in (0,1) to x.
#' @param names_final Optional character vector of final variable names.
#' @param process_final Optional simdata post-processing list passed to simdata::simdesign.
#' @param name Optional design name.
#' @param eps Numeric in (0,1). Clipping applied to U to avoid 0/1 in quantiles.
#' @param ... Further arguments stored in the simdata design object.
#'
#' @return A 'simdesign' object usable with simdata::simulate_data().
#' @export
simdesign_vine_copula <- function(
  vine,
  dist,
  names_final = NULL,
  process_final = list(),
  name = "Vine copula design",
  eps = 1e-6,
  ...
) {
  .assert_vinecopula_installed()

  dim <- .get_vine_dim(vine)
  .assert_is_function_list(dist, dim, "dist")

  generator <- function(n_obs, ...) {
    u <- VineCopula::RVineSim(n_obs, vine)
    u <- .clip_unit(u, eps = eps)
    u
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
    x <- as.data.frame(x, optional = TRUE, stringsAsFactors = FALSE)
    x
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

#' Fit a vine copula design from data (empirical margins)
#'
#' Notes:
#' - R-vine and C-vine are fit via VineCopula::RVineStructureSelect(type = 0/1).
#' - D-vine is fit with a fixed order via VineCopula::RVineCopSelect() using a D-vine structure matrix derived from VineCopula::D2RVine().
#'
#' @param data data.frame containing variables.
#' @param vars character vector of variable names to model.
#' @param vine_type One of "rvine", "cvine", "dvine".
#' @param order Optional integer vector specifying node order (used for fixed C-/D-vine). If NULL:
#'   - cvine uses automatic root selection (RVineStructureSelect(type=1))
#'   - dvine defaults to 1:d
#' @param familyset Integer vector of pair-copula families (VineCopula codes), or NA for default (all families).
#' @param selectioncrit "AIC", "BIC", or "logLik".
#' @param indeptest Logical; whether to perform independence tests.
#' @param level Significance level for independence tests.
#' @param trunclevel Optional truncation level.
#' @param method "mle" or "itau".
#' @param rotations Logical; whether to include rotations of families.
#' @param qtype Quantile type for empirical margins.
#' @param eps Clipping for pseudo-observations and simulation U.
#' @param name Optional simdesign name.
#' @param ... Stored in the returned simdesign (e.g., for metadata).
#'
#' @return A 'simdesign' object usable with simdata::simulate_data().
#' @export
simdesign_vine_copula_from_data <- function(
  data,
  vars,
  vine_type = c("rvine", "cvine", "dvine"),
  order = NULL,
  familyset = NA,
  selectioncrit = "AIC",
  indeptest = FALSE,
  level = 0.05,
  trunclevel = NA,
  method = c("mle", "itau"),
  rotations = TRUE,
  qtype = 8,
  eps = 1e-6,
  name = "Vine copula design (fit from data)",
  ...
) {
  .assert_vinecopula_installed()

  vine_type <- match.arg(vine_type)
  method <- match.arg(method)

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

  d <- ncol(u_hat)

  if (!is.null(order)) {
    if (!is.numeric(order) || length(order) != d) {
      stop("order must be a numeric vector of length ", d, ".", call. = FALSE)
    }
    order <- as.integer(order)
  }

  vine_fit <- NULL

  if (vine_type == "rvine") {
    vine_fit <- VineCopula::RVineStructureSelect(
      data = u_hat,
      familyset = familyset,
      type = 0,
      selectioncrit = selectioncrit,
      indeptest = indeptest,
      level = level,
      trunclevel = trunclevel,
      progress = FALSE,
      rotations = rotations,
      method = method
    )
  } else if (vine_type == "cvine") {
    if (is.null(order)) {
      vine_fit <- VineCopula::RVineStructureSelect(
        data = u_hat,
        familyset = familyset,
        type = 1,
        selectioncrit = selectioncrit,
        indeptest = indeptest,
        level = level,
        trunclevel = trunclevel,
        progress = FALSE,
        rotations = rotations,
        method = method
      )
    } else {
      dd <- d * (d - 1) / 2
      rvm0 <- VineCopula::C2RVine(
        order = order,
        family = rep(0, dd),
        par = rep(0, dd),
        par2 = rep(0, dd)
      )

      vine_fit <- VineCopula::RVineCopSelect(
        data = u_hat,
        familyset = familyset,
        Matrix = rvm0$Matrix,
        selectioncrit = selectioncrit,
        indeptest = indeptest,
        level = level,
        trunclevel = trunclevel,
        rotations = rotations,
        method = method
      )
    }
  } else {
    if (is.null(order)) {
      order <- seq_len(d)
    }

    dd <- d * (d - 1) / 2
    rvm0 <- VineCopula::D2RVine(
      order = order,
      family = rep(0, dd),
      par = rep(0, dd),
      par2 = rep(0, dd)
    )

    vine_fit <- VineCopula::RVineCopSelect(
      data = u_hat,
      familyset = familyset,
      Matrix = rvm0$Matrix,
      selectioncrit = selectioncrit,
      indeptest = indeptest,
      level = level,
      trunclevel = trunclevel,
      rotations = rotations,
      method = method
    )
  }

  simdesign_vine_copula(
    vine = vine_fit,
    dist = dist,
    names_final = vars,
    name = name,
    eps = eps,
    vine_fit = vine_fit,
    vine_type = vine_type,
    margins = "empirical",
    ...
  )
}