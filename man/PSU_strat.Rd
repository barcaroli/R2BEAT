\name{PSU_strat}
\Rdversion{1.0}
\alias{PSU_strat}
\docType{data}
\title{
Information on Primary Stage Units (PSUs) stratification   
}
\description{
Example data frame containing information on Primary Stage Units (PSUs) stratification.
}
\usage{
data(beat.example)
}
\format{
The PSU_strat data frame contains a row for each Primary Stage Units (PSUs) with the following variables: 
 \describe{
    \item{STRATUM}{
    Identifier of the stratum (numeric)
    } 
\item{PSU_MOS}{
    Measure of size of the primary stage unit (numeric)
    } 
\item{PSU_ID}{
    Identifier of the primary stage unit (numeric)
    } 
  }
}
\details{
Note: the names of the variables must be the ones indicated above.
}

\examples{
\dontrun{
# Load example data
data(beat.example)
PSU_strat
str(PSU_strat)
}
}
\keyword{datasets}
