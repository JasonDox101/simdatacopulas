\name{simdesign_empirical_copula}
\alias{simdesign_empirical_copula}
\alias{simdesign_empirical_copula_from_data}
\title{Empirical copula simdesigns}
\usage{
simdesign_empirical_copula(u_data, dist, names_final = NULL, process_final = list(), replace = TRUE, jitter = 0, name = "Empirical copula design", eps = 1e-6, ...)
simdesign_empirical_copula_from_data(data, vars, qtype = 8, replace = TRUE, jitter = 0, eps = 1e-6, name = "Empirical copula design (fit from data)", ...)
}
\arguments{
\item{u_data}{Matrix or data.frame of pseudo-observations in (0,1).}
\item{dist}{List of marginal quantile functions mapping \eqn{u \in (0,1)} to the target scale.}
\item{names_final}{Optional character vector of output names.}
\item{process_final}{Optional simdata post-processing list.}
\item{replace}{Logical; sample with replacement when resampling U.}
\item{jitter}{Non-negative jitter applied to U.}
\item{name}{Design name.}
\item{eps}{Numeric in (0, 0.5) used for clipping.}
\item{...}{Additional fields stored in the simdesign.}
\item{data}{data.frame containing variables.}
\item{vars}{Character vector of variable names.}
\item{qtype}{Quantile type for empirical margins.}
}
\value{A \code{simdesign} object usable with \code{simdata::simulate_data()}.}
\description{Construct simdata designs using an empirical copula with empirical margins.}
