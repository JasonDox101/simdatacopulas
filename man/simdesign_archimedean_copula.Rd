\name{simdesign_archimedean_copula}
\alias{simdesign_archimedean_copula}
\alias{simdesign_clayton_copula}
\alias{simdesign_gumbel_copula}
\alias{simdesign_frank_copula}
\alias{simdesign_joe_copula}
\alias{simdesign_archimedean_copula_from_data}
\title{Archimedean copula simdesigns}
\usage{
simdesign_archimedean_copula(copula, dist, names_final = NULL, process_final = list(), name = "Archimedean copula design", eps = 1e-6, ...)
simdesign_clayton_copula(dist, theta, dim, names_final = NULL, ...)
simdesign_gumbel_copula(dist, theta, dim, names_final = NULL, ...)
simdesign_frank_copula(dist, theta, dim, names_final = NULL, ...)
simdesign_joe_copula(dist, theta, dim, names_final = NULL, ...)
simdesign_archimedean_copula_from_data(data, vars, family = c("clayton", "gumbel", "frank", "joe"), fit_method = c("itau", "itau_mpl"), qtype = 8, eps = 1e-6, name = "Archimedean copula design (fit from data)", ...)
}
\arguments{
\item{copula}{An object inheriting from class \code{"copula"}.}
\item{dist}{List of marginal quantile functions mapping \eqn{u \in (0,1)} to the target scale.}
\item{names_final}{Optional character vector of output names.}
\item{process_final}{Optional simdata post-processing list.}
\item{name}{Design name.}
\item{eps}{Numeric in (0, 0.5) used for clipping.}
\item{...}{Additional fields stored in the simdesign.}
\item{theta}{Copula parameter.}
\item{dim}{Dimension (>= 2).}
\item{data}{data.frame containing variables.}
\item{vars}{Character vector of variable names.}
\item{family}{Archimedean family: \code{"clayton"}, \code{"gumbel"}, \code{"frank"}, or \code{"joe"}.}
\item{fit_method}{Fit method: \code{"itau"} or \code{"itau_mpl"}.}
\item{qtype}{Quantile type for empirical margins.}
}
\value{A \code{simdesign} object usable with \code{simdata::simulate_data()}.}
\description{Construct simdata designs based on Archimedean copulas with parametric or data-fitted dependence.}
