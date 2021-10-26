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
  strata_var, 
  target_vars, 
  deff_var, 
  domain_var,
  delta,
  minimum) 
}

\arguments{
 \item{samp_frame}{The dataframe containing sampling units in the reference population.}
 \item{RGdes}{The 'design' ReGenesees object.}
 \item{RGcal}{The 'calibration' ReGenesees object.}
 \item{id_PSU}{variables used as identifiers in sampling frame.}
 \item{id_SSU}{variables used as identifiers in sampling frame.}
 \item{strata_var}{stratification variable used in sampling frame.}
 \item{target_vars}{target variables.}
 \item{deff_var}{stratification variable to be used when calculating deff.}
 \item{domain_var}{the variable used to identify the domain of interest.}
 \item{delta}{average number of analysis units per sampling unit.}
 \item{minimum}{minimum number of SSU to be selected from a PSU.}
}

\value{
A list containing: (1) strata, (2) deff, (3) effst, (4) rho, (5) PSU_file, (6) des_file
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
head(inp$strata)
head(inp$deff)
head(inp$effst)
head(inp$rho)
head(inp$psu_file)
head(inp$des_file)
}
}
\author{Giulio Barcaroli}
