\name{beat.1cv}
\alias{beat.1cv}
\title{The function returns a dataframe with planned and actual coefficients of variation (CV) in a multivariate multi-domain allocation problem.}
\usage{
beat.1cv(stratif, errors, alloc, minnumstrat = 2)
}
\arguments{
\item{stratif}{A dataframe with strata information.}

\item{errors:}{A dataframe with planned CV.}

\item{alloc:}{A vector with an allocation into strata.}

\item{minnumstrat:}{Optional. A number indicating the lower bound for the number of allocated units in each stratum (default is 2).}
}
\value{
A dataframe with planned and actual CV for each combination of domain, domain category and auxiliary variable.
}
\description{
The function returns a dataframe with planned and actual coefficients of variation (CV) in a multivariate multi-domain allocation problem.
}
