\name{beat.cv}
\alias{beat.cv}
\title{
Computation of coefficient of variation (CV) for a given multivariate multiple allocation
}
\description{
Compute the coefficients of variation considering a given multivariate optimal allocation.
}
\usage{ 
beat.cv(n_file, stratif, errors, des_file, psu_file, rho, epsilon)
}
 
\arguments{
 \item{n_file}{Data frame containing the sample size allocated in each stratum, for more details, e.g.,\link{allocation}.}
  \item{stratif}{Data frame of survey strata, for more details see, e.g.,\link{strata}.}
 \item{errors}{Data frame of expected coefficients of variation (CV) for each domain, for more details see, e.g.,\link{errors}.}
 \item{des_file}{Data frame containing information on sampling design variables, for more details see, e.g.,\link{design}.}
 \item{psu_file}{Data frame containing information on primary stage units stratification, for more details see, e.g.,\link{PSU_strat}.}
 \item{rho}{Data frame of survey strata, for more details see, e.g.,\link{rho}.}
 \item{epsilon}{The same as in function \link{beat.1st}.}
}
\details{This function enables to derive the expected coefficient of variation (CV) from a given allocation.
The function \code{beat.cv} returns the estimates expected accuracy in terms of coefficient of variation, for several variables in different domains, given a certain allocation among the different strata.
}
\value{
Object of class \code{list}. The list contains a set of \code{data.frame}, as many of the cross product between domain and interest variables, containing total estimates, population, variance and expected coefficient of variation for every domain modality.
For each domain and each variable is defined
\item{Tot1}{Total estimate}
\item{N}{Measure of size}
\item{Varfin}{The sample variance of the total estimate}
\item{CV}{The coefficient of variation of the}
}
\examples{
\dontrun{

# Load example data
data(beat.example)

## Example 1
# Calculate coefficients of variation, for two variables in two domains,
# given an allocation among the different strata.
allocation
cv1<-beat.cv( n_file=allocation, stratif=strata, errors=errors,
des_file=design, psu_file=PSU_strat, rho=rho)

## Example 2
# Take the example 1 in beat.2st.
allocation2st_1 <- beat.2st(stratif=strata, errors=errors,
des_file=design, psu_file=PSU_strat,rho=rho)
# The allocation obtained is 
allocation2st_1$alloc
# with these precision constraints 
errors
# and these expected coefficient of variation
allocation2st_1$expected

# Now, fit the output of beat.2st to allocation, that is
SIZE <- allocation2st_1$alloc[-18,c(2)]
allocation1 <- data.frame(SIZE)
# If apply beat.cv the same error in allocation2st_1$expected should be obtained.
# In fact
cv2<-beat.cv( n_file=allocation1, stratif=strata, errors=errors,
des_file=design, psu_file=PSU_strat, rho=rho)
cv2

# Please, note that some very slightly differences may occur.
}
}
\author{Developed by Stefano Falorsi, Andrea Fasulo, Alessio Guandalini, Daniela Pagliuca, Marco D. Terribili.}
