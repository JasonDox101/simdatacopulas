\name{simdesign_factor_copula}
\alias{simdesign_factor_copula}
\alias{simdesign_factor_copula_from_data}
\alias{simdesign_latent_threshold_gaussian}
\alias{simdesign_latent_threshold_gaussian_from_data}
\title{Factor and latent-threshold copula simdesigns}
\usage{
simdesign_factor_copula(loadings, uniq_var = NULL, dist, names_final = NULL, process_final = list(), name = "Factor Gaussian copula design", eps = 1e-6, ...)
simdesign_factor_copula_from_data(data, vars, n_factors = 1, rotation = "none", qtype = 8, eps = 1e-6, name = "Factor Gaussian copula design (fit from data)", ...)
simdesign_latent_threshold_gaussian(loadings, uniq_var = NULL, thresholds, levels = NULL, names_final = NULL, process_final = list(), name = "Latent Gaussian threshold design", eps = 1e-6, ...)
simdesign_latent_threshold_gaussian_from_data(data, vars, n_factors = 1, rotation = "none", tie_jitter = 1e-8, eps = 1e-6, name = "Latent Gaussian threshold design (fit from data)", ...)
}
\arguments{
\item{loadings}{Numeric loading matrix.}
\item{uniq_var}{Numeric vector of uniqueness variances.}
\item{dist}{List of marginal quantile functions mapping \eqn{u \in (0,1)} to the target scale.}
\item{names_final}{Optional character vector of output names.}
\item{process_final}{Optional simdata post-processing list.}
\item{name}{Design name.}
\item{eps}{Numeric in (0, 0.5) used for clipping.}
\item{...}{Additional fields stored in the simdesign.}
\item{data}{data.frame containing variables.}
\item{vars}{Character vector of variable names.}
\item{n_factors}{Number of latent factors to fit.}
\item{rotation}{Factor rotation passed to \code{factanal}.}
\item{qtype}{Quantile type for empirical margins.}
\item{thresholds}{List of threshold vectors for each variable.}
\item{levels}{Optional list of character vectors defining ordered levels.}
\item{tie_jitter}{Non-negative jitter applied before pseudo-observations.}
}
\value{A \code{simdesign} object usable with \code{simdata::simulate_data()}.}
\description{Construct simdata designs using factor-based Gaussian copulas and latent threshold Gaussian models.}
