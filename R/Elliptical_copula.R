.make_empirical_q <- function(x, qtype = 8, eps = 1e-6) { # Build an "empirical quantile function": map u∈(0,1) to sample quantiles of x
  force(x) # Force-capture x to avoid lazy evaluation / external shadowing inside the closure
  function(u) { # Return a function: input u (probability) -> output the corresponding quantile
    u <- pmin(pmax(u, eps), 1 - eps) # Clip u to (eps, 1-eps) to avoid boundary issues at 0/1 in quantile()
    as.numeric(stats::quantile(x, probs = u, type = qtype, names = FALSE)) # Compute empirical quantiles and return a plain numeric vector
  } # End inner function
} # End .make_empirical_q

.assert_is_function_list <- function(x, expected_len, name) { # Assert: x must be a function list of expected length (used for marginal quantiles dist)
  if (!is.list(x) || length(x) != expected_len) { # Check that x is a list and has the expected length
    stop(name, " must be a list of length ", expected_len, ".", call. = FALSE) # Throw an error (without call stack)
  }
  ok <- vapply(x, is.function, logical(1)) # Check each element is a function, returning a logical vector
  if (!all(ok)) { # Fail if any element is not a function
    stop(name, " must be a list of functions.", call. = FALSE) # Throw an error if not a function list
  }
  invisible(TRUE) # Return an invisible TRUE so it can be used inside pipelines/functions
} # End .assert_is_function_list

.clip_unit <- function(u, eps = 1e-6) { # Clip numeric values to (eps, 1-eps), used for pseudo-observations or simulated U(0,1)
  pmin(pmax(u, eps), 1 - eps) # Vectorized clipping: lower bound then upper bound
} # End .clip_unit

.get_copula_dim <- function(copula) { # Get the copula object's dimension (workaround for copula::dim not being exported)
  d <- tryCatch(base::dim(copula), error = function(e) NULL) # Try base::dim first (works if the class implements a dim method)
  if (!is.null(d) && length(d) == 1 && is.finite(d)) { # If a single finite value is returned
    return(as.integer(d)) # Return as integer
  }

  d <- tryCatch(methods::slot(copula, "dimension"), error = function(e) NULL) # Then try reading S4 slot "dimension"
  if (!is.null(d) && length(d) == 1 && is.finite(d)) { # Validate it is a single finite value
    return(as.integer(d)) # Return as integer
  }

  stop("Unable to determine copula dimension.", call. = FALSE) # Both methods failed
} # End .get_copula_dim

#' Create a simdata design using an elliptical copula as joint generator
#'
#' This is the core building block: a copula generates dependent U(0,1)
#' samples, then marginal quantile functions map each column to the final scale.
#'
#' @param copula An object inheriting from class 'copula' (e.g., normalCopula, tCopula).
#' @param dist List of marginal quantile functions. Each function maps u in (0,1) to x.
#' @param names_final Optional character vector of final variable names.
#' @param process_final Optional simdata post-processing list passed to simdata::simdesign.
#' @param name Optional design name.
#' @param eps Numeric in (0,1). Clipping applied to U to avoid 0/1 in quantiles.
#' @param ... Further arguments stored in the simdata design object.
#'
#' @return A 'simdesign' object usable with simdata::simulate_data().
#' @export
simdesign_elliptical_copula <- function( # Build a simdata design: copula generates dependent U(0,1), then marginal quantiles map to the target scales
  copula, # Copula object (e.g., normalCopula / tCopula) defining the dependence structure
  dist, # List of marginal quantile functions: dist[[j]](u) maps u∈(0,1) to variable j
  names_final = NULL, # Optional output column names
  process_final = list(), # Optional simdata post-processing steps, executed after transform_initial
  name = "Elliptical copula design", # Design name
  eps = 1e-6, # Numerical-stability clipping for U to avoid 0/1 in quantiles
  ... # Further arguments passed to simdata::simdesign (e.g., metadata)
) { # Begin function body
  if (!inherits(copula, "copula")) { # Validate copula inherits from class "copula"
    stop("copula must inherit from class 'copula'.", call. = FALSE) # Error if invalid
  }

  dim <- .get_copula_dim(copula) # Infer dimension d (number of variables) from the copula object
  .assert_is_function_list(dist, dim, "dist") # Validate dist is a function list of length d

  generator <- function(n_obs, ...) { # simdata generator: produces the "initial state" random draws (here: a U matrix)
    u <- copula::rCopula(n_obs, copula) # Sample an n_obs×d dependent U(0,1) matrix from the copula
    u <- .clip_unit(u, eps = eps) # Clip to (eps, 1-eps) to avoid extreme quantiles (±Inf)
    u # Return U as input to transform_initial
  } # End generator

  transform_initial <- function(u) { # simdata transform_initial: map U to the final variable scales
    if (!is.matrix(u) && !is.data.frame(u)) { # Defensive check: generator should return a 2D object
      stop("Internal error: generator did not return a 2D object.", call. = FALSE) # Internal error if not 2D
    }
    u <- as.matrix(u) # Coerce to matrix for column-wise processing
    if (ncol(u) != dim) { # Check number of columns matches the copula dimension
      stop("Internal error: U has unexpected number of columns.", call. = FALSE) # Internal error if mismatch
    }

    x <- lapply(seq_len(dim), function(j) dist[[j]](u[, j])) # Apply each marginal quantile function to its U column
    x <- as.data.frame(x, optional = TRUE, stringsAsFactors = FALSE) # Assemble as a data.frame (do not force factors)

    x # Return mapped data (simdata handles naming and post-processing)
  } # End transform_initial

  simdata::simdesign( # Build and return the simdata design object
    generator = generator, # Initial generator
    transform_initial = transform_initial, # Mapping from initial state to final variables
    names_final = names_final, # Final column names
    process_final = process_final, # Final post-processing pipeline
    name = name, # Design name
    ... # Pass through remaining arguments
  ) # Return simdesign
} # End simdesign_elliptical_copula

#' Convenience constructor for Gaussian copula simdata designs
#'
#' @param dist List of marginal quantile functions.
#' @param rho Scalar correlation used when structure = "ex".
#' @param Sigma Correlation matrix used when structure = "un".
#' @param dim Dimension (required for structure = "ex"; inferred from Sigma for "un").
#' @param structure Correlation structure: "ex" (exchangeable) or "un" (unstructured).
#' @param names_final Optional names for final variables.
#' @param ... Passed to simdesign_elliptical_copula().
#'
#' @return A 'simdesign' object.
#' @export
simdesign_gaussian_copula <- function( # Convenience constructor: quickly build a simdata design using a Gaussian (normal) copula
  dist, # List of marginal quantile functions (length equals dimension)
  rho = 0, # Scalar correlation used when structure = "ex" (exchangeable)
  Sigma = NULL, # Correlation matrix used when structure = "un" (unstructured)
  dim = NULL, # Required for "ex"; inferred from Sigma for "un"
  structure = c("ex", "un"), # Correlation structure: "ex"=exchangeable; "un"=unstructured
  names_final = NULL, # Optional final column names
  ... # Passed through to simdesign_elliptical_copula
) { # Begin function body
  structure <- match.arg(structure) # Normalize structure

  if (structure == "ex") { # Exchangeable correlation: one rho for all off-diagonal correlations
    if (is.null(dim) || !is.numeric(dim) || length(dim) != 1 || dim < 2) { # Validate dim
      stop("For structure='ex', dim must be a single numeric >= 2.", call. = FALSE) # Error if invalid
    }
    cop <- copula::normalCopula(param = rho, dim = as.integer(dim), dispstr = "ex") # Build exchangeable Gaussian copula
  } else { # Unstructured correlation: determined by full correlation matrix Sigma
    if (is.null(Sigma) || !is.matrix(Sigma) || nrow(Sigma) != ncol(Sigma)) { # Validate Sigma is square
      stop("For structure='un', Sigma must be a square matrix.", call. = FALSE) # Error if not square
    }
    dim <- ncol(Sigma) # Dimension is determined by Sigma's size
    if (dim < 2) { # Dimension must be at least 2
      stop("Sigma must have dimension >= 2.", call. = FALSE) # Error otherwise
    }
    if (!all(is.finite(Sigma))) { # Check for non-finite values
      stop("Sigma must be finite.", call. = FALSE) # Error if non-finite
    }
    if (max(abs(diag(Sigma) - 1)) > 1e-10) { # Check diag is 1 (correlation matrix requirement)
      stop("Sigma must be a correlation matrix with diag = 1.", call. = FALSE) # Error if violated
    }
    cop <- copula::normalCopula(param = copula::P2p(Sigma), dim = dim, dispstr = "un") # Convert Sigma to parameter vector and build unstructured Gaussian copula
  }

  simdesign_elliptical_copula( # Reuse the generic elliptical copula design constructor
    copula = cop, # Built Gaussian copula
    dist = dist, # Marginals
    names_final = names_final, # Final names
    ... # Pass through remaining arguments
  ) # Return simdesign
} # End simdesign_gaussian_copula

#' Fit an elliptical copula design from data (empirical margins)
#'
#' MVP supports Gaussian copula + exchangeable/unstructured structure with robust defaults.
#'
#' @param data data.frame containing variables.
#' @param vars character vector of variable names to model.
#' @param family "gaussian" or "t" (t currently requires df provided).
#' @param structure "ex" or "un".
#' @param fit_method "itau" or "itau_mpl". "itau_mpl" tries MPL after ITAU init.
#' @param df Degrees of freedom if family="t".
#' @param qtype Quantile type for empirical margins.
#' @param eps Clipping for pseudo-observations and simulation U.
#' @param name Optional simdesign name.
#' @param ... Stored in the returned simdesign (e.g., for metadata).
#'
#' @return A 'simdesign' object usable with simdata::simulate_data().
#' @export
simdesign_elliptical_copula_from_data <- function( # Fit an elliptical copula from data and build a simdata design with empirical margins (quantile functions)
  data, # Input data (data.frame)
  vars, # Variable names to model (length >= 2)
  family = c("gaussian", "t"), # Copula family: Gaussian or t; t requires df
  structure = c("ex", "un"), # Correlation structure: exchangeable or unstructured
  fit_method = c("itau", "itau_mpl"), # Fit method: itau only (robust) or itau then try mpl
  df = NULL, # Degrees of freedom for t copula (required for family="t", must be > 2)
  qtype = 8, # Quantile type for empirical margins (stats::quantile algorithm)
  eps = 1e-6, # Numerical stability clipping for pobs and simulated U
  name = "Elliptical copula design (fit from data)", # Design name
  ... # Additional fields stored as metadata in the simdesign
) { # Begin function body
  family <- match.arg(family) # Normalize family
  structure <- match.arg(structure) # Normalize structure
  fit_method <- match.arg(fit_method) # Normalize fit_method

  if (!is.data.frame(data)) { # Validate input data type
    stop("data must be a data.frame.", call. = FALSE) # Error if not a data.frame
  }
  if (!is.character(vars) || length(vars) < 2) { # Validate vars
    stop("vars must be a character vector of length >= 2.", call. = FALSE) # vars must be a character vector of length >= 2
  }
  if (!all(vars %in% names(data))) { # Check vars exist in data
    miss <- setdiff(vars, names(data)) # Identify missing columns
    stop("vars missing from data: ", paste(miss, collapse = ", "), call. = FALSE) # Error with missing columns
  }

  df_in <- data[, vars, drop = FALSE] # Subset to modeled variables (keep as data.frame)
  for (nm in vars) { # Minimal viable checks for each variable
    if (!is.numeric(df_in[[nm]])) { # Current MVP supports numeric variables only (for pobs and correlation)
      stop("All vars must be numeric for MVP. Non-numeric: ", nm, call. = FALSE) # Error with variable name
    }
  }
  xmat <- as.matrix(df_in) # Convert to numeric matrix for downstream calculations
  if (!all(is.finite(xmat))) { # Check for NA/NaN/Inf
    stop("Data contains non-finite values in selected vars.", call. = FALSE) # Error if non-finite values exist
  }

  dist <- lapply(df_in, .make_empirical_q, qtype = qtype, eps = eps) # Build empirical marginal quantile functions (map simulated U back to original scale)

  u_hat <- copula::pobs(xmat) # Convert samples to pseudo-observations, approximating U(0,1) per column
  u_hat <- .clip_unit(u_hat, eps = eps) # Clip to (eps, 1-eps) to avoid 0/1 boundaries affecting fit and quantiles

  dim <- ncol(u_hat) # Copula dimension d (number of variables)

  if (family == "gaussian") { # Fit a Gaussian copula
    if (structure == "ex") { # Exchangeable: single parameter
      cop0 <- copula::normalCopula(param = 0, dim = dim, dispstr = "ex") # Initialize at zero correlation
    } else { # Unstructured: d(d-1)/2 parameters
      cop0 <- copula::normalCopula(param = rep(0, dim * (dim - 1) / 2), dim = dim, dispstr = "un") # Initialize all correlations at zero
    }
  } else { # Fit a t copula (stronger tail dependence)
    if (is.null(df) || !is.numeric(df) || length(df) != 1 || df <= 2) { # df constraints for t copula
      stop("For family='t', df must be a single numeric > 2.", call. = FALSE) # Error if invalid
    }
    if (structure == "ex") { # Exchangeable t copula
      cop0 <- copula::tCopula(param = 0, dim = dim, dispstr = "ex", df = df, df.fixed = TRUE) # df.fixed=TRUE means df is not estimated
    } else { # Unstructured t copula
      cop0 <- copula::tCopula(param = rep(0, dim * (dim - 1) / 2), dim = dim, dispstr = "un", df = df, df.fixed = TRUE) # Initialize correlation parameters at zero
    }
  }

  fit_itau <- copula::fitCopula(cop0, data = u_hat, method = "itau") # Robust fit/initialization via inversion of Kendall's tau (ITAU)
  fit_final <- fit_itau # Use itau as the default final fit

  if (fit_method == "itau_mpl") { # If selected, try MPL (maximum pseudo-likelihood) after itau
    fit_final <- tryCatch( # MPL can fail for some data/starts; fall back to itau on error
      copula::fitCopula(fit_itau@copula, data = u_hat, method = "mpl"), # Use itau-fitted copula as MPL starting point
      error = function(e) fit_itau # Fall back to itau fit if MPL errors
    )
  }

  simdesign_elliptical_copula( # Build a simdesign from fitted dependence + empirical margins for direct simulate_data()
    copula = fit_final@copula, # Fitted copula object (parameters/structure)
    dist = dist, # Empirical marginal quantile functions (preserve original marginal shapes)
    names_final = vars, # Output names match modeled variables
    name = name, # Design name
    eps = eps, # U clipping (also used during simulation)
    copula_fit = fit_final, # Store fit object for inspection (parameters/convergence)
    copula_family = family, # Metadata: copula family
    copula_structure = structure, # Metadata: correlation structure
    margins = "empirical", # Metadata: empirical margins
    ... # Pass through additional metadata fields
  ) # Return simdesign
} # End simdesign_elliptical_copula_from_data
