\name{sensitivity}
\alias{sensitivity}
\title{
Sensitivity analysis for some parameters by means of grid search 
}
\description{
This function allows to analyse the different results in terms 
of first stage size (number of PSUs) and second stage size (number of SSUs),
when varying the values of the following parameters: 

- deff_sugg (suggested value of deff) 
- minimum (minimum number of SSU per single PSU)
- f (suggested sample fraction)

The name of the parameter has to be given, together with the minimum and
maximum value. On the basis of these minimum and maximum values, 10 different values will be used for carrying out the allocation.
The output will be a graphical one.
It can be time consuming.

To be used only in the scenario when no previous rounds of the survey are 
available, and a frame complete with values of target variables is available.

}
\usage{
sensitivity ( 
  samp_frame,
  id_PSU,
  id_SSU,
  strata_var,
  target_vars,
  deff_var,
  domain_var,
  minimum,
  delta,
  f,
  deff_sugg,
  search=c("deff","min_SSU","sample_fraction"),
  min,
  max) 
}

\arguments{
 \item{samp_frame}{The dataframe containing sampling units in the reference population.}
 \item{id_PSU}{variables used as identifiers in sampling frame.}
 \item{id_SSU}{variables used as identifiers in sampling frame.}
 \item{strata_var}{stratification variable used in sampling frame.}
 \item{target_vars}{target variables.}
 \item{deff_var}{stratification variable to be used when calculating deff.}
 \item{domain_var}{the variable used to identify the domain of interest.}
 \item{minimum}{minimum number of SSU to be selected from a PSU.}
 \item{delta}{average number of analysis units per sampling unit.}
 \item{f}{expected (initial) sampling rate.}
 \item{deff_sugg}{suggested value of the deff.}
 \item{search}{parameter for which the analysis must be done.}
 \item{min}{minimum value of the parameter.}
 \item{max}{maximum value of the parameter.}
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
id_PSU <- "municipality"  
id_SSU <- "id_ind"        
strata_var <- "stratum"  
target_vars <- c("income_hh","active","inactive","unemployed")   # more than one
deff_var <- "stratum"        
domain_var <- "region"        
minimum <- 50 # minimum number of SSUs to be interviewed in each selected PSU
delta =  nrow(pop) /length(unique(pop$id_hh))     # average dimension of the SSU in terms of elementary survey units
delta =  1     # average dimension of the SSU in terms of elementary survey units
f = 0.05 # sampling fraction (suggested)
deff_sugg <- 1.5
sensitivity(samp_frame,
            id_PSU,
            id_SSU,
            strata_var,
            target_vars,
            deff_var,
            domain_var,
            minimum,
            delta,
            f,
            deff_sugg,
            search=c("deff"),
            min=1,
            max=2)
}
}
\author{Giulio Barcaroli}
