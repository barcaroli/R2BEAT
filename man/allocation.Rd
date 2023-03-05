\name{allocation}
\Rdversion{1.0}
\alias{allocation}
\docType{data}
\title{
Sample sizes for each stratum
}
\description{
Example data frame containing a given allocation.
}
\usage{
data(beat.example)
}
\format{
  The Strata data frame contains a row per each stratum with the following variables:
  \describe{
    \item{SIZE}{
    Stratum sample size (numeric)
    }
  }
 }
\details{
Note: the names of the variables must be the ones indicated above.
}
\examples{
# Load example data
data(beat.example)
allocation
str(allocation)
}
\keyword{datasets}

