\name{prepareInputToAllocation2}
\alias{prepareInputToAllocation2}
\title{
Input dataframes for R2BEAT two-stages sample design when sampling frame is available
}
\description{
In case of scenario 2 (at least one previous round of the survey available), prepares the following input dataframes for R2BEAT two-stages sample design starting from the sampling frame:
1. strata
2. deff
3. effst
4. rho
5. PSU_file
6. des_file
}
\usage{
prepareInputToAllocation2 ( 
  samp_frame,
  RGdes, 
  RGcal, 
  id_PSU, 
  id_SSU, 
  strata_vars, 
  target_vars, 
  deff_vars, 
  domain_vars,
  delta,
  minimum) 
}

\arguments{
 \item{samp_frame}{The dataframe containing sampling units in the reference population.}
 \item{RGdes}{The 'design' ReGenesees object.}
 \item{RGcal}{The 'calibration' ReGenesees object.}
 \item{id_PSU}{variables used as identifiers in sampling frame.}
 \item{id_SSU}{variables used as identifiers in sampling frame.}
 \item{strata_vars}{stratification variable used in sampling frame.}
 \item{target_vars}{target variables.}
 \item{deff_vars}{stratification variable to be used when calculating deff.}
 \item{domain_var}{the variable used to identify the domain of interest.}
 \item{minimum}{minimum number of SSU to be selected from a PSU.}
 \item{delta}{average number of analysis units per sampling unit.}
 \item{deff_sugg}{suggested value of the deff.}
}

\value{
A list containing: (1) strata, (2) deff, (3) effst, (4) rho, (5) PSU_file, (6) des_file
}

\examples{
\dontrun{
library(R2BEAT)
# the following two RData are available at the link:
# https://github.com/barcaroli/R2BEAT/tree/master/data
load("pop.RData")                 # sampling frame  
load("R2BEAT_ReGenesees.RData")   # ReGenesees design and calibration objects plus PSU data

samp_frame <- pop
RGdes <- des                           # ReGenesees design object
RGcal <- cal                           # ReGenesees calibrated object

strata_vars <- c("stratum")            # variables of stratification
target_vars <- c("income_hh",
                 "active",
                 "inactive",
                 "unemployed")         # target variables
deff_vars <- "stratum"    # stratification variables for calculating deff and effst 
                          # (n.b: must coincide or be a subset of variables of stratification)
id_PSU <- c("municipality")            # identification variable of PSUs
id_SSU <- c("id_hh")                   # identification variable of SSUs
domain_vars <- c("region")             # domain variables
delta=1                      # Average number of SSUs for each selection unit
minimum <- 50                # Minimum number of SSUs to be selected in each PSU
inp <- prepareInputToAllocation2 ( 
  samp_frame,
  RGdes, 
  RGcal, 
  id_PSU, 
  id_SSU, 
  strata_vars, 
  target_vars, 
  deff_vars, 
  domain_vars,
  delta,
  minimum) 
}
}
\author{Giulio Barcaroli}
