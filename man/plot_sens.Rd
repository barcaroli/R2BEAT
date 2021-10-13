\name{plot_sens}
\alias{plot_sens}
\title{
Plot of the sensitivity analysis for some parameters by means of grid search 
}
\description{
This function allows to plot the results of the simulation carried out by using the function 'sensitivity'.

}
\usage{
plot_sens ( 
  x,
  min,
  max) 
}

\arguments{
 \item{x}{The result of the 'sensitivity' function.}
 \item{min}{minimum value of the parameter.}
 \item{max}{maximum value of the parameter.}
}

\value{
A list containing the (i) vector of allocated PSUs in the iterations and (ii) the vector of allocated SSUs in the iterations
}

\examples{
\dontrun{
library(R2BEAT)
load("pop.RData")
cv <- as.data.frame(list(DOM=c("DOM1","DOM2"),
                         CV1=c(0.02,0.03),
                         CV2=c(0.03,0.05),
                         CV3=c(0.03,0.05),
                         CV4=c(0.04,0.08)))
cv
# parameters
samp_frame <- pop
errors <- cv
id_PSU <- "municipality"  
id_SSU <- "id_ind"        
strata_var <- "stratum"  
target_vars <- c("income_hh","active","inactive","unemployed")   # more than one
deff_var <- "stratum"        
domain_var <- "region"        
minimum <- 50 # minimum number of SSUs to be interviewed in each selected PSU
# average dimension of the SSU in terms of elementary survey units
#delta =  nrow(pop) /length(unique(pop$id_hh))     
delta =  1  # average dimension of the SSU in terms of elementary survey units
deff_sugg <- 1.5
min <- 30
max <- 80
sens <- sensitivity_min_SSU (
  samp_frame,
  errors,
  id_PSU,
  id_SSU,
  strata_var,
  target_vars,
  deff_var,
  domain_var,
  minimum,
  delta,
  deff_sugg,
  min,
  max)
plot_sens(sens,min,max)
}
}
\author{Giulio Barcaroli}
