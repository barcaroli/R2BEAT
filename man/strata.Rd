\name{strata}
\Rdversion{1.0}
\alias{strata}
\docType{data}
\title{
Strata characteristics
}
\description{
Example data frame containing information on strata characteristics.
}
\usage{
data(beat.example)
}
\format{
  The Strata data frame contains a row per each stratum with the following variables: 
  \describe{
    \item{STRATUM}{
    Identifier of the stratum (numeric).
    } 
    \item{N}{
    Stratum population size (numeric).
    } 
    \item{M1}{
    Mean in the stratum of the first variable (numeric).
    }
    \item{Mj}{
    Mean in the stratum of the j-th variable (numeric).
    }
    \item{Mn}{
    Mean in the stratum of the last variable (numeric).
    }
    \item{S1}{
    Standard deviation in the stratum of the first variable (numeric).
    }
    \item{Sj}{
    Standard deviation in the stratum of the j-th variable (numeric).
    }
    \item{Sn}{
    Standard deviation in the stratum of the last variable (numeric).
    }
    \item{CENS}{flag (1 indicates a take all stratum, 0 a sampling stratum, usually 0 ) (numeric).
    }
    \item{COST}{Cost per interview in each stratum, usually 0  (numeric).
    }
    \item{DOM1}{
    Domain value to which the stratum belongs for the first type of domain (factor or numeric).
    }
    \item{DOMa}{
    Domain value to which the stratum belongs for the a-th type of domain  (factor or numeric).
    }
    \item{DOMk}{
    Domain value to which the stratum belongs for the k-th type of domain  (factor or numeric).
    }
  }
 }
\details{
Note: the names of the variables must be the ones indicated above.
}
\keyword{datasets}

\examples{
\dontrun{
# Load example data
data(beat.example)
strata
str(strata)
}
}

