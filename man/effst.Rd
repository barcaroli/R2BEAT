\name{effst}
\Rdversion{1.0}
\alias{effst}
\docType{data}
\title{
Estimator effect
}
\description{
Example data frame containing estimator effect, (\eqn{effst}), in each stratum for each variable.
}
\usage{
data(beat.example)
}
\format{
  The estimator effect data frame contains a row per each stratum with the following variables: 
  \describe{
    \item{STRATUM}{
    Identifier of the stratum (numeric).
    } 
    \item{EFFST1}{
    Estimator effect in the stratum of the first variable (numeric).
    }
    \item{EFFSTj}{
    Estimator effect in the stratum of the j-th variable (numeric).
    }
    \item{EFFSTn}{
    Estimator effect in the stratum of the last variable (numeric).
    }
  }
 }
\details{
Note: the names of the variables must be the ones indicated above.

The estimator effect, (\eqn{effst}), provides a measure of the variance inflaction or reduction due to the use of a different estimator from the HT (Horvitz and Thompson, 1952).
It is equal to the ratio between the sampling variance of the estimator planned to be used and the sampling variance  of the HT. 
%In formula
%\deqn{effst = \frac{\text{var} \left(\hat{Y}_{\cdot}\right)}{\text{var} \left(\hat{Y}_{HT}\right)}}.

Then, when the HT is used, (\eqn{effst}) is equal to 1. 
However, always more often, different estimators, such as calibration estimator (Deville and Särndal, 1992) or generalized regression estimator GREG (Fuller, 2002 and references therein), are used.
Usually this kind of estimators take into account auxiliary variables that enables to increase the accuracy of the estimates, that is, they reduce their errors (CV). Then, their \eqn{effst} is usually lower than 1.

Therefore, taking into account the estimator effect when planning the survey can help in saving sample size or at least to more properly evatuate the allocation.  
}
\keyword{datasets}

\references{
  Deville, J.C., Särndal, C.E. (1992). 
  Calibration estimators in survey sampling. 
  \emph{Journal of the American statistical Association}, 87(418): 376-38.
  
  Fuller, W.A.. (2002). 
  Regression estimation for survey samples.
  Survey Methodology 28(1): 5-23.
  
  Horvitz, D.G., Thompson, D.J. (1952) 
  A generalization of sampling without replacement from a finite universe. 
  \emph{Journal of the American statistical Association}, 47(260): 663-685.
}


\examples{
# Load example data
data(beat.example)
effst
str(effst)
}

