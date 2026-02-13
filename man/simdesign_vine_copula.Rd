\name{simdesign_vine_copula}
\alias{simdesign_vine_copula}
\alias{simdesign_vine_copula_from_data}
\title{Vine copula simdesigns}
\usage{
simdesign_vine_copula(vine, dist, names_final = NULL, process_final = list(), name = "Vine copula design", eps = 1e-6, ...)
simdesign_vine_copula_from_data(data, vars, vine_type = c("rvine", "cvine", "dvine"), order = NULL, familyset = NA, selectioncrit = "AIC", indeptest = FALSE, level = 0.05, trunclevel = NA, method = c("mle", "itau"), rotations = TRUE, qtype = 8, eps = 1e-6, name = "Vine copula design (fit from data)", ...)
}
\arguments{
\item{vine}{A VineCopula vine object.}
\item{dist}{List of marginal quantile functions mapping \eqn{u \in (0,1)} to the target scale.}
\item{names_final}{Optional character vector of output names.}
\item{process_final}{Optional simdata post-processing list.}
\item{name}{Design name.}
\item{eps}{Numeric in (0, 0.5) used for clipping.}
\item{...}{Additional fields stored in the simdesign.}
\item{data}{data.frame containing variables.}
\item{vars}{Character vector of variable names.}
\item{vine_type}{One of \code{"rvine"}, \code{"cvine"}, \code{"dvine"}.}
\item{order}{Optional integer vector specifying node order.}
\item{familyset}{VineCopula family set; use \code{NA} for default families.}
\item{selectioncrit}{Selection criterion.}
\item{indeptest}{Logical; perform independence tests.}
\item{level}{Significance level for independence tests.}
\item{trunclevel}{Truncation level.}
\item{method}{Estimation method: \code{"mle"} or \code{"itau"}.}
\item{rotations}{Logical; allow rotated families.}
\item{qtype}{Quantile type for empirical margins.}
}
\value{A \code{simdesign} object usable with \code{simdata::simulate_data()}.}
\description{Construct simdata designs using vine copulas with empirical margins or data-fitted dependence.}
