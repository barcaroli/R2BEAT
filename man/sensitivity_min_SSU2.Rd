\name{sensitivity_min_SSU2}
\alias{sensitivity_min_SSU2}
\title{
Sensitivity analysis for choosing minimum number of SSUs per PSU (no sampling frame available)
}
\description{
This function is similar to the function sensitivity_min_SSU, as it allows to analyse the different results in terms 
of first stage size (number of PSUs) and second stage size (number of SSUs),
when varying the values of the minimum number of SSU per single PSU)
The name of the parameter has to be given, together with the minimum and
maximum value. Differently from sensitivity_min_SSU, it requires in input all the outputs already prepared by prepareInputToAllocation1 or prepareInputToAllocation2.

}
\usage{
sensitivity_min_SSU2 (strata, 
                      des_file,
                      psu_file,
                      rho,
                      effst,
                      errors, 
                      min, 
                      max, 
                      plot) 
}

\arguments{
 \item{strata}{Data frame of survey strata, for more details see, e.g.,\link{strata}.}
 \item{des_file}{Data frame containing information on sampling design variables, for more details see, e.g.,\link{design}.}
 \item{psu_file}{Data frame containing information on primary stage units stratification, for more details see, e.g.,\link{PSU_strat}.}
 \item{rho}{Data frame of survey strata, for more details see, e.g.,\link{rho}.}
 \item{effst}{Data frame of survey strata, for taking into account the estimator effect on each variable, for more details see, e.g.,\link{effst}.}
 \item{errors}{Data frame of coefficients of variation for each domain, for more details see, e.g.,\link{errors}.}
 \item{min}{starting value for the minimum value of SSUs per PSU.}
 \item{max}{ending value for the minimum value of SSUs per PSU.}
 \item{plot}{plot (TRUE/FALSE) the final result.}
}

\value{
A list containing the (i) vector of allocated PSUs in the iterations and (ii) the vector of allocated SSUs in the iterations
}

\examples{
\dontrun{
library(readr)
samp <- read_rds("https://github.com/barcaroli/R2BEAT_workflows/blob/master/sample.rds?raw=true")
pop <- read_rds("https://github.com/barcaroli/R2BEAT_workflows/blob/master/pop.RDS?raw=true")
library(R2BEAT)
str(samp)
## Sample design description
library(ReGenesees)
samp$stratum_2 <- as.factor(samp$stratum_2)
sample.des <- e.svydesign(samp, 
                          ids= ~ municipality + id_hh, 
                          strata = ~ stratum_2, 
                          weights = ~ weight,
                          self.rep.str = ~ SR,
                          check.data = TRUE)
## Find and collapse lonely strata
ls <- find.lon.strata(sample.des)
if (!is.null(ls)) sample.des <- collapse.strata(sample.des)
## Calibration with known totals
totals <- pop.template(sample.des,
             calmodel = ~ sex : cl_age, 
             partition = ~ region)
totals <- fill.template(pop, totals, mem.frac = 10)
sample.cal <- e.calibrate(sample.des, 
                          totals,
                          calmodel = ~ sex : cl_age, 
                          partition = ~ region,
                          calfun = "logit",
                          bounds = c(0.3, 2.6), 
                          aggregate.stage = 2,
                          force = FALSE)
samp_frame <- pop
RGdes <- sample.des
RGcal <- sample.cal
strata_var <- c("stratum")      
target_vars <- c("income_hh",
                 "active",
                 "inactive",
                 "unemployed")   
weight_var <- "weight"
deff_var <- "stratum"            
id_PSU <- c("municipality")      
id_SSU <- c("id_hh")             
domain_var <- c("region") 
delta <- 1                   
minimum <- 50     
inp <- prepareInputToAllocation2 ( 
  samp_frame,
  RGdes, 
  RGcal, 
  id_PSU, 
  id_SSU, 
  strata_var, 
  target_vars, 
  deff_var, 
  domain_var,
  delta,
  minimum) 
cv <- as.data.frame(list(DOM=c("DOM1","DOM2"),
                         CV1=c(0.02,0.03),
                         CV2=c(0.03,0.06),
                         CV3=c(0.03,0.06),
                         CV4=c(0.05,0.08)))
sens <- sensitivity_min_SSU2 (strata = inp$strata, 
                         des_file = inp$des_file,
                         psu_file = inp$psu_file,
                         rho = inp$rho,
                         effst = inp$effst,
                         errors = cv, 
                         min = 15, 
                         max = 25, 
                         plot = TRUE) 
}
}
\author{Giulio Barcaroli}
