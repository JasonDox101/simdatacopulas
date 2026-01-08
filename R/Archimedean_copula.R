#' Create a simdata design using an Archimedean copula as joint generator
#'
#' @param copula An object inheriting from class 'copula' (e.g., claytonCopula, gumbelCopula, frankCopula, joeCopula).
#' @param dist List of marginal quantile functions. Each function maps u in (0,1) to x.
#' @param names_final Optional character vector of final variable names.
#' @param process_final Optional simdata post-processing list passed to simdata::simdesign.
#' @param name Optional design name.
#' @param eps Numeric in (0,1). Clipping applied to U to avoid 0/1 in quantiles.
#' @param ... Further arguments stored in the simdata design object.
#'
#' @return A 'simdesign' object usable with simdata::simulate_data().
#' @export
simdesign_archimedean_copula <- function( # Entry point: build a simdata design using an Archimedean copula (reuses the generic copula design)
  copula, # Copula object (e.g., claytonCopula / gumbelCopula / frankCopula / joeCopula) defining dependence and tail behavior
  dist, # List of marginal quantile functions: dist[[j]](u) maps u∈(0,1) to variable j's scale
  names_final = NULL, # Optional output column names
  process_final = list(), # Optional simdata post-processing steps, executed after transform_initial
  name = "Archimedean copula design", # Design name
  eps = 1e-6, # Numerical stability clipping for U to prevent 0/1 from producing ±Inf in quantiles
  ... # Passed through to simdata::simdesign (via simdesign_elliptical_copula)
) { # Begin function body
  simdesign_elliptical_copula( # Reuse the generic copula -> U -> marginal-quantile mapping (works for any copula-class object)
    copula = copula, # Dependence structure
    dist = dist, # Marginal quantiles
    names_final = names_final, # Final names
    process_final = process_final, # Final post-processing
    name = name, # Design name
    eps = eps, # U clipping
    ... # Pass through additional metadata/arguments
  ) # Return simdesign
} # End simdesign_archimedean_copula

#' Convenience constructor for Clayton copula simdata designs
#'
#' @param dist List of marginal quantile functions.
#' @param theta Copula parameter.
#' @param dim Dimension (>= 2).
#' @param names_final Optional names for final variables.
#' @param ... Passed to simdesign_archimedean_copula().
#'
#' @return A 'simdesign' object.
#' @export
simdesign_clayton_copula <- function( # Convenience constructor: build a simdata design using a Clayton copula (common lower-tail dependence)
  dist, # List of marginal quantile functions
  theta, # Clayton copula parameter controlling dependence strength (range defined by the copula package)
  dim, # Dimension d (number of variables, >= 2)
  names_final = NULL, # Optional final column names
  ... # Passed through to simdesign_archimedean_copula
) { # Begin function body
  if (is.null(dim) || !is.numeric(dim) || length(dim) != 1 || dim < 2) { # Validate dim
    stop("dim must be a single numeric >= 2.", call. = FALSE) # Error if invalid
  }
  cop <- copula::claytonCopula(param = theta, dim = as.integer(dim)) # Build Clayton copula object (dependence structure)

  simdesign_archimedean_copula( # Wrap into a simdesign via the generic Archimedean entry point
    copula = cop, # Clayton copula
    dist = dist, # Marginals
    names_final = names_final, # Final names
    ... # Pass through remaining arguments
  ) # Return simdesign
} # End simdesign_clayton_copula

#' Convenience constructor for Gumbel copula simdata designs
#'
#' @param dist List of marginal quantile functions.
#' @param theta Copula parameter.
#' @param dim Dimension (>= 2).
#' @param names_final Optional names for final variables.
#' @param ... Passed to simdesign_archimedean_copula().
#'
#' @return A 'simdesign' object.
#' @export
simdesign_gumbel_copula <- function( # Convenience constructor: build a simdata design using a Gumbel copula (common upper-tail dependence)
  dist, # List of marginal quantile functions
  theta, # Gumbel copula parameter controlling dependence strength (often theta >= 1)
  dim, # Dimension d (number of variables, >= 2)
  names_final = NULL, # Optional final column names
  ... # Passed through to simdesign_archimedean_copula
) { # Begin function body
  if (is.null(dim) || !is.numeric(dim) || length(dim) != 1 || dim < 2) { # Validate dim
    stop("dim must be a single numeric >= 2.", call. = FALSE) # Error if invalid
  }
  cop <- copula::gumbelCopula(param = theta, dim = as.integer(dim)) # Build Gumbel copula object

  simdesign_archimedean_copula( # Wrap into simdesign
    copula = cop, # Gumbel copula
    dist = dist, # Marginals
    names_final = names_final, # Final names
    ... # Pass through remaining arguments
  ) # Return simdesign
} # End simdesign_gumbel_copula

#' Convenience constructor for Frank copula simdata designs
#'
#' @param dist List of marginal quantile functions.
#' @param theta Copula parameter (cannot be 0).
#' @param dim Dimension (>= 2).
#' @param names_final Optional names for final variables.
#' @param ... Passed to simdesign_archimedean_copula().
#'
#' @return A 'simdesign' object.
#' @export
simdesign_frank_copula <- function( # Convenience constructor: build a simdata design using a Frank copula (symmetric dependence, no tail dependence)
  dist, # List of marginal quantile functions
  theta, # Frank copula parameter (theta=0 degenerates to independence; range defined by the copula package)
  dim, # Dimension d (number of variables, >= 2)
  names_final = NULL, # Optional final column names
  ... # Passed through to simdesign_archimedean_copula
) { # Begin function body
  if (is.null(dim) || !is.numeric(dim) || length(dim) != 1 || dim < 2) { # Validate dim
    stop("dim must be a single numeric >= 2.", call. = FALSE) # Error if invalid
  }
  cop <- copula::frankCopula(param = theta, dim = as.integer(dim)) # Build Frank copula object

  simdesign_archimedean_copula( # Wrap into simdesign
    copula = cop, # Frank copula
    dist = dist, # Marginals
    names_final = names_final, # Final names
    ... # Pass through remaining arguments
  ) # Return simdesign
} # End simdesign_frank_copula

#' Convenience constructor for Joe copula simdata designs
#'
#' @param dist List of marginal quantile functions.
#' @param theta Copula parameter.
#' @param dim Dimension (>= 2).
#' @param names_final Optional names for final variables.
#' @param ... Passed to simdesign_archimedean_copula().
#'
#' @return A 'simdesign' object.
#' @export
simdesign_joe_copula <- function( # Convenience constructor: build a simdata design using a Joe copula (upper-tail dependence with sharper behavior)
  dist, # List of marginal quantile functions
  theta, # Joe copula parameter controlling dependence strength (range defined by the copula package)
  dim, # Dimension d (number of variables, >= 2)
  names_final = NULL, # Optional final column names
  ... # Passed through to simdesign_archimedean_copula
) { # Begin function body
  if (is.null(dim) || !is.numeric(dim) || length(dim) != 1 || dim < 2) { # Validate dim
    stop("dim must be a single numeric >= 2.", call. = FALSE) # Error if invalid
  }
  cop <- copula::joeCopula(param = theta, dim = as.integer(dim)) # Build Joe copula object

  simdesign_archimedean_copula( # Wrap into simdesign
    copula = cop, # Joe copula
    dist = dist, # Marginals
    names_final = names_final, # Final names
    ... # Pass through remaining arguments
  ) # Return simdesign
} # End simdesign_joe_copula

#' Fit an Archimedean copula design from data (empirical margins)
#'
#' @param data data.frame containing variables.
#' @param vars character vector of variable names to model.
#' @param family One of "clayton", "gumbel", "frank", "joe".
#' @param fit_method "itau" or "itau_mpl". "itau_mpl" tries MPL after ITAU init.
#' @param qtype Quantile type for empirical margins.
#' @param eps Clipping for pseudo-observations and simulation U.
#' @param name Optional simdesign name.
#' @param ... Stored in the returned simdesign (e.g., for metadata).
#'
#' @return A 'simdesign' object usable with simdata::simulate_data().
#' @export
simdesign_archimedean_copula_from_data <- function( # Fit an Archimedean copula dependence structure from data and build a simdesign with empirical margins
  data, # Input data (data.frame)
  vars, # Variable names to model (length >= 2)
  family = c("clayton", "gumbel", "frank", "joe"), # Archimedean family (different tail dependence behaviors)
  fit_method = c("itau", "itau_mpl"), # Fit method: itau (robust) or itau then try mpl (fallback on failure)
  qtype = 8, # Empirical quantile type (stats::quantile algorithm)
  eps = 1e-6, # Numerical stability clipping for pobs and simulated U
  name = "Archimedean copula design (fit from data)", # Design name
  ... # Additional fields stored as metadata in the simdesign
) { # Begin function body
  family <- match.arg(family) # Normalize family
  fit_method <- match.arg(fit_method) # Normalize fit_method

  if (!is.data.frame(data)) { # Validate data type
    stop("data must be a data.frame.", call. = FALSE) # Error if not a data.frame
  }
  if (!is.character(vars) || length(vars) < 2) { # Validate vars
    stop("vars must be a character vector of length >= 2.", call. = FALSE) # Must be a character vector of length >= 2
  }
  if (!all(vars %in% names(data))) { # Check vars exist in data
    miss <- setdiff(vars, names(data)) # Identify missing columns
    stop("vars missing from data: ", paste(miss, collapse = ", "), call. = FALSE) # Error with missing columns
  }

  df_in <- data[, vars, drop = FALSE] # Subset to modeled variables (keep as data.frame)
  for (nm in vars) { # MVP: numeric variables only
    if (!is.numeric(df_in[[nm]])) { # If a column is non-numeric
      stop("All vars must be numeric for MVP. Non-numeric: ", nm, call. = FALSE) # Error with variable name
    }
  }
  xmat <- as.matrix(df_in) # Convert to matrix for downstream calculations
  if (!all(is.finite(xmat))) { # Check for NA/NaN/Inf
    stop("Data contains non-finite values in selected vars.", call. = FALSE) # Error if non-finite values exist
  }

  dist <- lapply(df_in, .make_empirical_q, qtype = qtype, eps = eps) # Build empirical marginal quantile functions (preserve original marginal shapes in simulation)

  u_hat <- copula::pobs(xmat) # Convert samples to pseudo-observations U (approx. U(0,1) per column) for copula fitting
  u_hat <- .clip_unit(u_hat, eps = eps) # Clip to (eps, 1-eps) to avoid 0/1 boundaries affecting fit

  dim <- ncol(u_hat) # Variable dimension d

  cop0 <- switch( # Choose an initial copula object by family (used to initialize fitCopula)
    family, # Branch selector
    clayton = copula::claytonCopula(param = 0.5, dim = dim), # Clayton init (example: moderate dependence)
    gumbel  = copula::gumbelCopula(param = 1.1, dim = dim), # Gumbel init (slightly above independence: theta=1 is independence)
    frank   = copula::frankCopula(param = 1, dim = dim), # Frank init (theta=0 is independence; 1 is a mild start)
    joe     = copula::joeCopula(param = 1.1, dim = dim), # Joe init (theta=1 is independence; 1.1 is a mild start)
    stop("Unsupported family.", call. = FALSE) # Error for unsupported family
  )

  fit_itau <- copula::fitCopula(cop0, data = u_hat, method = "itau") # Robust fit/initialization via ITAU (inversion of Kendall's tau)
  fit_final <- fit_itau # Use itau as the default final fit

  if (fit_method == "itau_mpl") { # If selected, try MPL (maximum pseudo-likelihood) after itau
    fit_final <- tryCatch( # MPL can be unstable for some data/starts; fall back on failure
      copula::fitCopula(fit_itau@copula, data = u_hat, method = "mpl"), # Use itau-fitted copula as MPL starting point
      error = function(e) fit_itau # Fall back to itau if MPL errors
    )
  }

  simdesign_archimedean_copula( # Wrap "fitted dependence + empirical margins" into a simdesign ready for simulation
    copula = fit_final@copula, # Fitted copula object (estimated parameters)
    dist = dist, # Empirical marginal quantile function list
    names_final = vars, # Output names match modeled variables
    name = name, # Design name
    eps = eps, # U clipping
    copula_fit = fit_final, # Store fit object for inspecting parameters/fit information
    copula_family = family, # Metadata: copula family
    margins = "empirical", # Metadata: empirical margins
    ... # Pass through remaining metadata
  ) # Return simdesign
} # End simdesign_archimedean_copula_from_data
