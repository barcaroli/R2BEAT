\name{sensitivity_min_SSU}
\alias{sensitivity_min_SSU}
\title{
Sensitivity analysis for choosing minimum number of SSUs per PSU
}
\description{
This function allows to analyse the different results in terms 
of first stage size (number of PSUs) and second stage size (number of SSUs),
when varying the values of the minimum number of SSU per single PSU)
The name of the parameter has to be given, together with the minimum and
maximum value. On the basis of these minimum and maximum values, 10 different values will be used for carrying out the allocation.
The output will be a graphical one.
To be used only in the scenario when no previous rounds of the survey are 
available, and a frame complete with values of target variables is available.

}
\usage{
sensitivity_min_SSU ( 
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
  min,
  max,
  plot=TRUE) 
}

\arguments{
 \item{samp_frame}{The dataframe containing sampling units in the reference population.}
 \item{errors}{Precision constraints.}
 \item{id_PSU}{variables used as identifiers in sampling frame.}
 \item{id_SSU}{variables used as identifiers in sampling frame.}
 \item{strata_var}{stratification variable used in sampling frame.}
 \item{target_vars}{target variables.}
 \item{deff_var}{stratification variable to be used when calculating deff.}
 \item{domain_var}{the variable used to identify the domain of interest.}
 \item{delta}{average number of analysis units per sampling unit.}
 \item{deff_sugg}{suggestion for deff value.}
 \item{min}{minimum value of the parameter.}
 \item{max}{maximum value of the parameter.}
  \item{plot}{plot (TRUE/FALSE) the final result.}
}

\value{
A list containing the (i) vector of allocated PSUs in the iterations and (ii) the vector of allocated SSUs in the iterations
}

\examples{
\dontrun{
library(R2BEAT)
data(pop)
cv <- as.data.frame(list(DOM=c("DOM1","DOM2"),
                         CV1=c(0.03,0.04),
                         CV2=c(0.06,0.08),
                         CV3=c(0.06,0.08),
                         CV4=c(0.06,0.08)))
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
# average dimension of the SSU in terms of elementary survey units
#delta =  nrow(pop) /length(unique(pop$id_hh))     
delta =  1  # average dimension of the SSU in terms of elementary survey units
deff_sugg <- 1.5  # deff (suggested)
sensitivity_min_SSU(
            samp_frame,
            errors,
            id_PSU,
            id_SSU,
            strata_var,
            target_vars,
            deff_var,
            domain_var,
            delta,
            deff_sugg,
            min=1,
            max=2)
}
}
\author{Giulio Barcaroli}
