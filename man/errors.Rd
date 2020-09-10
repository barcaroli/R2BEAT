\name{errors}
\Rdversion{1.0}
\alias{errors}
\docType{data}
\title{
Precision constraints (maximum CVs) as input for Bethel allocation
}
\description{
Example data frame containing precision levels (expressed in terms of acceptable CV's).
}
\usage{
data(beat.example)
}
\format{
  The constraint data frame (errors) contains a row per each type of domain with the following variables:
  \describe{
 \item{DOM}{Type of domain code (factor).} 
 \item{CV1}{Planned coefficient of variation for first variable (numeric).} 
 \item{CVj}{Planned coefficient of variation for j-th variable (numeric).}
 \item{CVn}{Planned coefficient of variation for last variable (numeric).} 
  }
}
\details{
Note: the names of the variables must be the ones indicated above.

The coefficient of variation (CV) is a standardized measure of variance. 
It is often expressed as a percentage and is defined as the ratio between the standard deviation of the estimate and the estimate (or its absolute value). 
%In formula:
%\deqn{CV(\hat{Y}) = \frac{\sqrt{\text{var} \left(\hat{Y}\right)}}{\hat{Y}} * 100}.
}

\examples{
# Load example data
data(beat.example)
errors
str(errors)
}
\keyword{datasets}
