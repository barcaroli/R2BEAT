\name{design}
\Rdversion{1.0}
\alias{design}
\docType{data}
\title{
Sampling design variables 
}
\description{
Example data frame containing variables for describing the sampling design. 
}
\usage{
data(beat.example)
}
\format{
  The design data frame contains a row per each stratum with the following variables: 
  \describe{
    \item{STRATUM}{
    Identifier of the stratum (numeric)
    } 
    \item{STRAT_MOS}{
    Measure of size of the stratum (numeric)
    }     
    \item{DELTA}{
    The average size of Secondary Stage Units (SSU) in the strata. With respect to the sample on which we are interested in, it could be equal or greater than 1 (numeric). See details for a depth explanation.
   } 
   \item{MINIMUM}{
    the minimum number of SSU to be selected in each PSU. It could be different in each stratum (numeric)
   } 
}
}
\details{
Note: the names of the variables must be the ones indicated above.

The sample design can be defined through a measure of size of the stratum, the average size of each SSU (>=1) and the minimum number of SSU to be selected in each PSU.
In particular, if SSU are not cluster DELTA=1 and the sample size determined will be given in term of SSU.
Instead, when SSUs are, in turn, clusters (for instance, households composed by individuals), defining DELTA equal to the average size of SSUs, enables to derive a sample in term of individuals.

Furthermore, modifying the MINIMUM it is possible to tune the number of PSU in the sample (see the example in \link{beat.1st}). 
In fact, considering the same sample size, increasing the MINIMUM, less PSU will be involved in the sample, but worst estimates in term of expected coefficient of variations will be provided.
On the contrary, decreasing the MINIMUM, more PSU will be involved in the sample and better estimates will be obtained.
Instead, increasing the MINIMUM for obtaining the same expected errors, requests less PSU, but much more SSU. The contrary occurs decreasing the MINIMUM.
}

\examples{
\dontrun{
# Load example data
data(beat.example)
design
str(design)
}
}
\keyword{datasets}