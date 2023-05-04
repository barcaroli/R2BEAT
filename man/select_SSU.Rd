\name{select_SSU}
\alias{select_SSU}
\title{
Select sample of secondary stage units (SSU)
}
\description{
Select sample of secondary stage units (SSU) from the population frame on the basis of the SSU allocated to each selected PSU
}
\usage{
     select_SSU(df,PSU_code,SSU_code,PSU_sampled)
}

\arguments{
 \item{df}{Dataframe containing sampling units (SSUs)}
 \item{PSU_code}{Identifier of each PSU in dataframe containing sampling units}
 \item{SSU_code}{Identifier of each SSU in dataframe containing sampling units}
 \item{PSU_sampled}{Dataframe containing selected PSUs}
}

\value{
A dataframe containing the units selected in the sample
}

\examples{
\dontrun{
library(readr)
pop <- read_rds("https://github.com/barcaroli/R2BEAT_workflows/blob/master/pop.RDS?raw=true")
library(R2BEAT)
cv <- as.data.frame(list(DOM=c("DOM1","DOM2"),
                         CV1=c(0.02,0.03),
                         CV2=c(0.03,0.06),
                         CV3=c(0.03,0.06),
                         CV4=c(0.05,0.08)))
cv
samp_frame <- pop
samp_frame$one <- 1
id_PSU <- "municipality"  
id_SSU <- "id_ind"        
strata_var <- "stratum"   
target_vars <- c("income_hh","active","inactive","unemployed")   
deff_var <- "stratum"     
domain_var <- "region"  
delta =  1       # households = survey units
minimum <- 50    # minimum number of SSUs to be interviewed in each selected PSU
deff_sugg <- 1.5 # suggestion for the deff value

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
inp$desfile$MINIMUM <- 50
alloc <- beat.2st(stratif = inp$strata, 
                   errors = cv, 
                   des_file = inp$des_file, 
                   psu_file = inp$psu_file, 
                   rho = inp$rho, 
                   deft_start = NULL,
                   effst = inp$effst, 
                   minPSUstrat = 2,
                   minnumstrat = 50
)
sample_1st <- select_PSU(alloc, type="ALLOC", pps=TRUE, plot=TRUE)
samp <- select_SSU(df=pop,
                   PSU_code="municipality",
                   SSU_code="id_ind",
                   PSU_sampled=sample_1st$sample_PSU)
}
}
\author{Giulio Barcaroli}
