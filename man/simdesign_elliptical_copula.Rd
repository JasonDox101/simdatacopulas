\name{simdesign_elliptical_copula}
\alias{simdesign_elliptical_copula}
\alias{simdesign_gaussian_copula}
\alias{simdesign_elliptical_copula_from_data}
\title{Elliptical copula simdesigns}
\usage{
simdesign_elliptical_copula(copula, dist, names_final = NULL, process_final = list(), name = "Elliptical copula design", eps = 1e-6, ...)
simdesign_gaussian_copula(dist, rho = 0, Sigma = NULL, dim = NULL, structure = c("ex", "un"), names_final = NULL, ...)
simdesign_elliptical_copula_from_data(data, vars, family = c("gaussian", "t"), structure = c("ex", "un"), fit_method = c("itau", "itau_mpl"), df = NULL, qtype = 8, eps = 1e-6, name = "Elliptical copula design (fit from data)", ...)
}
\arguments{
\item{copula}{An object inheriting from class \code{"copula"}.}
\item{dist}{List of marginal quantile functions mapping \eqn{u \in (0,1)} to the target scale.}
\item{names_final}{Optional character vector of output names.}
\item{process_final}{Optional simdata post-processing list.}
\item{name}{Design name.}
\item{eps}{Numeric in (0, 0.5) used for clipping.}
\item{...}{Additional fields stored in the simdesign.}
\item{rho}{Scalar correlation for exchangeable structure.}
\item{Sigma}{Correlation matrix for unstructured structure.}
\item{dim}{Dimension for exchangeable structure.}
\item{structure}{Correlation structure: \code{"ex"} or \code{"un"}.}
\item{data}{data.frame containing variables.}
\item{vars}{Character vector of variable names.}
\item{family}{Copula family: \code{"gaussian"} or \code{"t"}.}
\item{fit_method}{Fit method: \code{"itau"} or \code{"itau_mpl"}.}
\item{df}{Degrees of freedom for t copula when \code{family="t"}.}
\item{qtype}{Quantile type for empirical margins.}
}
\value{A \code{simdesign} object usable with \code{simdata::simulate_data()}.}
\description{Construct simdata designs based on elliptical copulas with parametric or data-fitted dependence.}
