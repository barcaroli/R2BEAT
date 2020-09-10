\name{beat.2st}
\alias{beat.2st}
\title{
Multivariate optimal allocation for different domains in two stage statified sample design
}
\description{
Compute multivariate optimal allocation for different domains corrected considering stratified two stages design
}
\usage{
	beat.2st(stratif, errors, des_file, psu_file, rho, deft_start = NULL, 
    effst = NULL, epsilon1 = 5, mmdiff_deft = 1,maxi = 20, 
    epsilon = 10^(-11), minnumstrat = 2, maxiter = 200, maxiter1 = 25)
}
\arguments{
 \item{stratif}{Data frame of survey strata, for more details see, e.g.,\link{strata}.}
 \item{errors}{Data frame of coefficients of variation for each domain, for more details see, e.g.,\link{errors}.}
 \item{des_file}{Data frame containing information on sampling design variables, for more details see, e.g.,\link{design}.}
 \item{psu_file}{Data frame containing information on primary stage units stratification, for more details see, e.g.,\link{PSU_strat}.}
 \item{rho}{Data frame of survey strata, for more details see, e.g.,\link{rho}.}
 \item{deft_start}{Data frame of survey strata, for taking into account the initial design effect on each variable, for more details see, e.g.,\link{deft_start}.}
 \item{effst}{Data frame of survey strata, for taking into account the estimator effect on each variable, for more details see, e.g.,\link{effst}.}
 \item{epsilon1}{First stop condition: sample sizes differences beetween two iterations; iteration continues until the maximum of sample sizes differences is greater than the default value. The default is 5.} 
 \item{mmdiff_deft}{Second stop condition: defts differences beetween two iterations; iteration continues until the maximum of defts largest differences is greater than the default value. The default is 0.06.}
 \item{maxi}{Third stop condition: maximum number of allowed iterations. The default is 20.}
 \item{epsilon}{The same as in function \link{beat.1st}.}
 \item{minnumstrat}{The same as in function \link{beat.1st}.}
 \item{maxiter}{The same as in function \link{beat.1st}.}
 \item{maxiter1}{The same as in function \link{beat.1st}.}
}
\details{
The methodology is a generalization of Bethel multivariate allocation (1989) that extended the Neyman (1959) - Tchuprov (1923) allocation for multi-purpose and multi-domains surveys.
The generalized Bethel’s algorithm allows to determine the optimal sample size for each stratum in a stratified sample design. The overall sample size and the allocation among the different strata is determined starting from the accuracy constraints imposed in the survey on interest estimates.
The optimal allocation is obtained throught a procedure that converge in few iteractions:

The first iteration is a computation of an initial allocation with the multivariate optimal allocation for different domains in one stages statified sample design (the methodology is a generalization for multidomains and multistages designs of Bethel multivariate allocation, 1989).

The correction of the initial allocation is based on an iterative method calculating new allocations and is based on an inflaction of strata variances  using the design effect (Ganninger, 2010).

}

\value{
Object of class \code{list}. The list contains 8 objects:
\item{iteractions}{Data frame that for each iteraction provides a summary with the number of Primary Stage Units (\code{PSU_Total}) distinguish between Self-Representative (\code{PSU_SR}) from Non-Self-Representative (\code{PSU_NSR}) and the number of Secondary Stage Units (\code{SSU}). This output is also printed to the screen.}
\item{file_strata}{Input data frame in \code{stratif} with the design effect for each variables in each stratum (\code{DEFT1 - DEFTn}) and the optimal sample size columns.}
\item{alloc}{Data frame with optimal (\code{ALLOC}), proportional (\code{PROP}), equal (\code{EQUAL}) sample size allocation.}
\item{planned}{Data frame with a summary of expected coefficients of variation for each variable in each domain.}
\item{expected}{Data frame with a summary of realized coefficients of variation with the given optimal allocation for each variable in each domain.}
\item{sensitivity}{Data frame with a summary of the sensitivity at 10\% for each domain and each variable. Sensitivity can be a useful tool to help in finding the best allocation, because it provides a hint of the expected sample size variation for a 10\% change in planned CVs.}
\item{deft_c}{Data frame with the design effect for each variable in each domain in each iteraction. Note that \code{DEFT1_0 - DEFTn_0} is always equal to 1 if \code{deft_start} is \code{NULL}. Instead is equal to \code{deft_start}. While \code{DEFT1 - DEFTn} are the final design effect related to the given allocation.}
\item{param_alloc}{A vector with a resume of all the parameter given for the allocation.}
}
\references{
  Cochran, W. (1977) 
  \emph{Sampling Techniques.} 
  John Wiley & Sons, Inc., New York
  
  Ganninger, M. (2010). 
  \emph{Design effects: model-based versus design-based approach.} 
  Vol. 3, p. 174. DEU.
  
  Neyman, J. (1934). 
  On the two different aspects of the representative method: the method of stratified sampling and the method of purposive selection. 
  \emph{Journal of the Royal Statistical Society}, 97(4), 558-625.
  
  Tschuprow, A. A. (1923). On the mathematical expectation of the moments of frequency distributions in the case of correlated observation. 
  (Chapters 4-6). \emph{Metron}, 1923, 2: 646-683.
}
\examples{
# Load example data
data(beat.example)

## Example 1
# Allocate the sample
allocation2st_1 <- beat.2st(stratif=strata, errors=errors,
des_file=design, psu_file=PSU_strat,rho=rho)
# The total ammount of sample size is 191 PSU (36 SR + 155 NSR) and 15147 SSU. 

## Example 2
# Assume 13000 SSUs is the maximum sample size to stick to our budget.
# Look at the sensitivity is in DOM1 for REG1 and REG2 due to V1.
allocation2st_1$sensitivity
# We can relax the constraints increasing the expected coefficients of variation for X1 by 10%
errors1 <- errors 
errors1[1,2] <- errors[1,2]+errors[1,2]*0.1

# Try the new allocation 
allocation2st_2 <- beat.2st(stratif=strata, errors=errors1,
des_file=design, psu_file=PSU_strat,rho=rho)

## Example 3
# On the contrary, if we tighten the constraints decreasing the expected coefficients of variation 
# for X1 by 10%
errors2 <- errors 
errors2[1,2] <- errors[1,2]-errors[1,2]*0.1

# The new allocation leads to a larger sample than the first example (around 18000)
allocation2st_3 <- beat.2st(stratif=strata, errors=errors2,
des_file=design, psu_file=PSU_strat,rho=rho)

## Example 4
# Sometimes some budget constraints concern the number of PSU involved in the survey.
# Tuning the PSUs number is possible modyfing the MINIMUM in des_file. 
# Assume to increase the MINIMUM from 48 to 60
design1 <- design 
design1[,4] <- 60
allocation2st_4 <- beat.2st(stratif=strata, errors=errors2,
des_file=design1, psu_file=PSU_strat, rho=rho)

# The PSUs numer is decreased, while the SSUs number increased 
# due to cluster intra-correlation effect.
# Under the same expected errors, to offset a slight reduction of PSUs (from 221 to 207) 
# an increase of SSUs involved is observed.
allocation2st_3$expected
allocation2st_4$expected

## Example 5
# On the contrary, assume to decrease the MINIMUM from 48 to 24.
# The SSUs number strongly decrease in the face of an increase of PSUs,
# always under the same expected errors.
design2 <- design 
design2[,4] <- 24
allocation2st_5 <- beat.2st(stratif=strata, errors=errors2,
des_file=design2, psu_file=PSU_strat, rho=rho)
allocation2st_4$expected
allocation2st_5$expected

# Note that MINIMUM can be different for each stratum.

## Example 6
# Assume that the SSUs are in turn clusters, for instance households composed by individuals.
# In the previous examples we always derived optimal allocations
# for sample of SSUs (i.e. households, because
# DELTA = 1).
design
design1
design2
# For obtaining a sample in terms of the elements composing SSUs
# (i.e., individuals) is just sufficient to 
# modify the DELTA in des_file.
design3 <- design 
design3$DELTA <- 2.31
# DELTA_IND=2.31, the average size of household in Italy.
allocation2st_6 <- beat.2st(stratif=strata, errors=errors,
des_file=design3, psu_file=PSU_strat, rho=rho)
}
\author{
Developed by Stefano Falorsi, Andrea Fasulo, Alessio Guandalini, Daniela Pagliuca, Marco D. Terribili.
}
