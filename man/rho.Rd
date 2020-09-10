\name{rho}
\Rdversion{1.0}
\alias{rho}
\docType{data}
\title{
Intraclass correlation coefficients for self and non self representative in the strata
}
\description{
Example data frame containing intraclass correlation \eqn{\rho} in Self Representative (SR) and Non Self Representative (NSR) strata.
}
\usage{
data(beat.example)
}
\format{
  The intraclass correlation coefficienta (\eqn{\rho}) data frame contains a row per each stratum with the following variables: 
  \describe{
    \item{STRATUM}{
    Identifier of the stratum (numeric)
    } 
    \item{RHO_AR1}{
    intraclass correlation of the elementary units for each primary stage unit of the self representing area belonging to the stratum for the first variable. 
    }
    \item{RHO_ARj}{
    intraclass correlation of the elementary units for each primary stage unit of the self representing area belonging to the stratum for the j-th variable. 
    }
    \item{RHO_ARn}{
    intraclass correlation of the elementary units for each primary stage unit of the self representing area belonging to the stratum for the n-th variable.
    }
    \item{RHO_NAR1}{
    intraclass correlation of the elementary units for each primary stage unit of the non self representing area belonging to the stratum for the first variable. 
    }
    \item{RHO_NARj}{
    intraclass correlation of the elementary units for each primary stage unit of the non self representing area belonging to the stratum for the j-th variable. 
    }
    \item{RHO_NARn}{
    intraclass correlation of the elementary units for each primary stage unit of the non self representing area belonging to the stratum for the n-th variable. 
    }
  }
 }
\details{
Note: the names of the variables must be the ones indicated above.

Intraclass correlation, \eqn{\rho}, provide a measure of the cluster heterogeneity and they have a direct impact on the design effect (\link{design}).
It can be indirectly computed from the design effect and the average minimum number of interviews in the Primary Stage Units (PSUs). 
%\out{b<sup>*</sup>}:
%\deqn{\rho= \frac{deff - 1}{b^{*}-1}}.
%\out{
%<i>&#961=<div class="fraction"><span class="fup">(deff - 1)</span><span class="bar">/</span><span class="fdn">(<i>b<sup>*</sup>-1)</span></div></i>
%}
%<p><span class="math display">\[\deqn{\rho= \frac{deff - 1}{b^{*}-1}}\]</span></p>
The ideal situation is when all the clusters in which the population is divided are more heterogeneous possible within them. 
At the limit, if each cluster were a reduced copy of the population then it would be sufficient to extract one just to have the same information that would be obtained from a complete survey.
Then, more similar the units in the cluster are, higher the sample size must be (Cochran, 1977, Chapter 8).

By definition, in SR strata \eqn{\rho}, is equal to 1, because there is just a single PSU in SR strata.
In NSR strata usually, \eqn{\rho} is usual higher than 1, because a double stage of selection is needed. 
}

\keyword{datasets}

\references{
  Cochran, W. (1977) 
  \emph{Sampling Techniques.} 
  John Wiley & Sons, Inc., New York.
}

\examples{
# Load example data
data(beat.example)
rho
str(rho)
}