\name{prepareInputToAllocation1}
\alias{prepareInputToAllocation1}
\title{
Input dataframes for R2BEAT two-stages sample design when sampling frame is available
}
\description{
In case of scenario 1 (no previous round of the survey available), prepares the following input dataframes for R2BEAT two-stages sample design starting from the sampling frame:
1. strata
2. deff
3. effst
4. rho
5. PSU_file
6. des_file
}
\usage{
prepareInputToAllocation1 ( 
  samp_frame,
  id_PSU,
  id_SSU,
  strata_var,
  target_vars,
  deff_var,
  domain_var,
  minimum,
  delta,
  deff_sugg) 
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
 \item{deff_sugg}{suggested value of the deff.}
}

\value{
A list containing  (i) the vector of allocated PSUs in the iterations and (ii) the vector of allocated SSUs in the iterations
}

\examples{
\dontrun{
library(R2BEAT)
library(SamplingStrata)
library(sampling)
load("pop.RData")
pop$active <- ifelse(pop$work==1,1,0)
pop$inactive <- ifelse(pop$work==2,1,0)
pop$one <- 1
# parameters
samp_frame <- pop
id_PSU <- "municipality"  
id_SSU <- "id_ind"        
strata_var <- "stratum"  
target_vars <- c("income_hh","active","inactive","unemployed") 
deff_var <- "stratum"        
domain_var <- "region"        
minimum <- 50 # minimum number of SSUs to be interviewed in each selected PSU
# average dimension of the SSU in terms of elementary survey units
# delta =  nrow(pop) /length(unique(pop$id_hh))     
delta =  1     
deff_sugg <- 1.5
# prepare inputs
inp <- prepareInputToAllocation1(samp_frame,
                                id_PSU,
                                id_SSU,
                                strata_var,
                                target_vars,
                                deff_var,
                                domain_var,
                                minimum,
                                delta,
                                deff_sugg)
}
}
\author{Giulio Barcaroli}
