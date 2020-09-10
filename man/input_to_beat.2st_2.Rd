\name{input_to_beat.2st_2}
\alias{input_to_beat.2st_2}
\title{
Prepares the design file for two-stage sample design
}
\description{
Prepares the design file for two-stage sample design on the basis of a dataset containin information on each PSU
}
\usage{
     input_to_beat.2st_2(psu,psu_id,stratum_var,mos_var,delta,minimum)
}

\arguments{
 \item{psu}{Dataframe containing information on each PSU.}
 \item{psu_id}{Identifier of each PSU in PSU dataframe.}
 \item{stratum_var}{Identifier of stratum in PSU dataframe.}
 \item{mos_var}{Variable containing the number of selection units in each PSU.}
 \item{delta}{Average number of final number of SSU per each selection unit.}
 \item{minimum}{Minimum number of selection units to be interviewed in each PSU.}
}


\examples{
\dontrun{
# psu <- read.csv2("psu.csv") # Read the external file containing PSU information
load("R2BEAT_ReGenesees.RData") # ReGenesees design and calibration objects plus PSU data
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
}
}
\author{Giulio Barcaroli}
