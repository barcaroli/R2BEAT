\name{beat.1st}
\alias{beat.1st}
\title{
Compute one stage multivariate optimal allocation.
}
\description{
Compute multivariate optimal allocation for different domains in one stage stratified sample design
}
\usage{
     beat.1st(stratif, errors, minnumstrat=2, maxiter=200, maxiter1=25, epsilon=10^(-11))
}

\arguments{
 \item{stratif}{Data frame of survey strata, for more details see, e.g.,\link{strata}.}
 \item{errors}{Data frame of expected coefficients of variation (CV) for each domain, for more details see, e.g.,\link{errors}.}
 \item{minnumstrat}{Minimum number of elementary units per strata (default=2).}
 \item{maxiter}{ Maximum number of iterations (default=200) of the general procedure. This kind of iteration may be required by the fact that when in a stratum the number of allocated units is greater or equal to its population, that stratum is set as "census stratum", and the whole procedure is re-initialised.}
 \item{maxiter1}{Maximum number of iterations in Chromy algorithm (default=25).}
 \item{epsilon}{Tollerance for the maximum absolute differences between the expected CV and the realised CV with the allocation obtained in the last iteraction for all domains. The default is 10^(-11).}
}

\details{The methodology is a generalization of Bethel multivariate allocation (1989) that extended the Neyman (1959) - Tchuprov (1923) allocation for multi-purpose and multi-domains surveys.
The generalized Bethel’s algorithm allows to determine the optimal sample size for each stratum in a stratified sample design. The overall sample size and the allocation among the different strata is determined starting from the accuracy constraints imposed in the survey on interest estimates.}

\value{Object of class \code{list}. The list contains 4 objects:
\item{n}{Vector with the optimal sample size for each stratum.}
\item{file_strata}{Data frame corresponding to the input data.frame \code{stratif} with the \code{n} optimal sample size column added.}
\item{alloc}{Data frame with optimal (\code{ALLOC}), proportional (\code{PROP}), equal (\code{EQUAL}) sample size allocation.}
\item{sensitivity}{Data frame with a summary of expected coefficients of variation (\code{Planned CV}), realized coefficients of variation with the given optimal allocation (\code{Actual CV}) and the sensitivity at 10\% for each domain and each variable. Sensitivity can be a useful tool to help in finding the best allocation, because it provides a hint of the expected sample size variation for a 10\% change in planned CVs.}
}
\references{

  Bethel, J. (1989)
  \emph{Sample allocation in multivariate surveys.}
  Survey methodology, 15.1: 47-57.

  Cochran, W. (1977) 
  \emph{Sampling Techniques.} 
  John Wiley & Sons, Inc., New York
  
  Neyman, J. (1934). 
  On the two different aspects of the representative method: the method of stratified sampling and the method of purposive selection. 
  \emph{Journal of the Royal Statistical Society}, 97(4): 558-625.
  
  Tschuprow, A. A. (1923). On the mathematical expectation of the moments of frequency distributions in the case of correlated observation. 
  (Chapters 4-6). \emph{Metron}, 2: 646-683.
}
\examples{
# Load example data
data(beat.example)

## Example 1
# Allocate the sample
allocation_1 <- beat.1st(stratif=strata, errors=errors)

# The total sample size is
sum(allocation_1$n)

## Example 2
# Assume 5700 units is the maximum sample size to stick to our budget.
# Looking at allocation_1$sensitivity we can see that most of the 
# sensitivity is in DOM1 for REG1 and REG2 due to V1.
allocation_1$sensitivity
# We can relax the constraints increasing the expected coefficients of variation for X1 by 10%
errors1 <- errors 
errors1[1,2] <- errors[1,2]+errors[1,2]*0.1

# Try the new allocation 
allocation_2 <- beat.1st(stratif=strata, errors=errors1)
sum(allocation_2$n)

## Example 3
# On the contrary, if we tighten the constraints decreasing the expected coefficients of variation 
# for X1 by 10%
errors2 <- errors 
errors2[1,2] <- errors[1,2]-errors[1,2]*0.1

# The new allocation leads to a larger sample than the first example 
allocation_3 <- beat.1st(stratif=strata, errors=errors2)
sum(allocation_3$n)
}
\author{Developed by Stefano Falorsi, Andrea Fasulo, Alessio Guandalini, Daniela Pagliuca, Marco D. Terribili.}
