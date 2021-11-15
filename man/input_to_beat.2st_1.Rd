\name{input_to_beat.2st_1}
\alias{input_to_beat.2st_1}
\title{
Input dataframes for R2BEAT two-stages sample design (when a previous round of the survey is available, but no sampling frame)
}
\description{
Prepares the following input dataframes for R2BEAT two-stages sample design starting from ReGenesees design and/or calibrated objects:
1. strata
2. deff
3. effst
4. rho
}
\usage{
     input_to_beat.2st_1(RGdes,
                        RGcal,
                        id_PSU,
                        id_SSU,
                        strata_vars,
                        target_vars,
                        deff_vars,
                        domain_vars)
}

\arguments{
 \item{RGdes}{ReGenesees design object.}
 \item{RGcal}{ReGenesees calibrated object.}
 \item{id_PSU}{variables used as identifiers in ReGenesees objects.}
 \item{id_SSU}{variables used as identifiers in ReGenesees objects.}
 \item{strata_vars}{stratification variables used in ReGenesees objects.}
 \item{target_vars}{target variables.}
 \item{deff_vars}{stratification variables to be used when calculating deff.}
 \item{domain_vars}{the variables used to identify the domain of interest.}
}


\examples{
\dontrun{
library(R2BEAT)

load("R2BEAT_ReGenesees.RData")   # ReGenesees design and calibration objects plus PSU data

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
inp1 <- input_to_beat.2st_1(RGdes,
                            RGcal,
                            id_PSU,
                            id_SSU,
                            strata_vars,
                            target_vars,
                            deff_vars,
                            domain_vars)							
inp1$strata
inp1$deff
inp1$effst
inp1$rho
}
}
\author{Giulio Barcaroli}
