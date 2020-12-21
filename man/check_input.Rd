\name{check_input}
\alias{check_input}
\title{
Check of coherence in the inputs for the allocation step
}
\description{
Checks the coherence between the population in the strata dataset and the population calculated by the PSUs dataset
}
\usage{
     check_input(strata,des,strata_var_strata,strata_var_des)
}

\arguments{
 \item{strata}{strata dataset}
 \item{des}{design dataset}
 \item{strata_var_strata}{variable identifying stratum in strata dataset}
 \item{strata_var_des}{variable identifying stratum in design dataset}
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
inp1 <- check_input(RGdes,
                            RGcal,
                            id_PSU,
                            id_SSU,
                            strata_vars,
                            target_vars,
                            deff_vars,
                            domain_vars)							
head(inp1$strata)
head(psu)
psu_id="municipality"        # Identifier of the PSU
stratum_var="stratum"        # Identifier of the stratum
mos_var="ind"                # Variable to be used as 'measure of size'
delta=1                      # Average number of SSUs for each selection unit
minimum <- 50                # Minimum number of SSUs to be selected in each PSU
inp2 <- input_to_beat.2st_2(psu,
                            psu_id,
                            stratum_var,
                            mos_var,
                            delta,
                            minimum)
head(inp2$psu_file)
head(inp2$des_file)
newstrata <- check_input(strata=inp1$strata,
                         des=inp2$des_file,
                         strata_var_strata="STRATUM",
                         strata_var_des="STRATUM")
}
}
\author{Giulio Barcaroli}
