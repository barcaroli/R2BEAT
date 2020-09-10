\name{StratSel} 
\alias{StratSel}
\title{Stratification and selection}

\description{
Merge two data frames (whatever PSU population and allocation data frames) and then compute stratification and selection of a fixed number of sample PSUs (Primary Sampling Units) for stratum using the Sampford's method (unequal probabilities, without replacement, fixed sample size), implemented by the UPsampford function of R package \emph{sampling}. About the stratification, the function realizes for each domain estimation, the computation of a size threshold for a given size PSU. PSUs with measure of size that exceeds a calculated threshold are identify like SR (Self Representative) and each constitutes a stratum by itself, so they come into the sample with probability equal to one. The remaining NSR (Non Self Representative) PSUs are ordered for the measure of size and divided into stratum having size approximately constant to the corrected threshold and with PSUs having sizes as homogeneous as possible.\cr
If the number of NSR PSUs in each stratum is greater than the number of sample PSUs for each stratum, then it is indispensable to define the vector of the inclusion probabilities (pik), argument of UPsampford function, vice versa NSR PSUs will become SR PSUs. If some pik values are greater than 1, pik is newly calculated, for those PSUs belonging to the same domain and stratum of the PSUs having pik >1, until pik values will be all less than 1, then applies the UPsampford function.
}

\usage{
StratSel(dataPop, idpsu, dom, final_pop, size, PSUsamplestratum, min_sample, 
		min_sample_index = FALSE, dataAll, domAll, f_sample,
		planned_min_sample = NULL, launch = TRUE)
}

\arguments{
  \item{dataPop}{PSU population data frame.}
  \item{idpsu}{Formula identifying primary sample units identifier.}
  \item{dom}{Formula identifying the variable domain.}  
  \item{final_pop}{Formula identifying secondary stage units.}
  \item{size}{Formula identifying first stage units size.} 
  \item{PSUsamplestratum}{Number of sample PSUs to select to each stratum.}
  \item{min_sample}{Number of final sample units to observe in each PSU.}
  \item{min_sample_index}{Identify the absence of planned_min_sample (default = \code{FALSE})}
  \item{dataAll}{Allocation data frame.}
  \item{domAll}{Formula identifying the domain variable for allocation data frame.}
  \item{f_sample}{Formula identifying final sample units.} 
  \item{planned_min_sample}{Formula identifying planned final sample units to observe in each PSU, variable between domains; \code{NULL} (the default)
  means the existence of a fixed number of them identify from min_sample.}
  \item{launch}{Identify the parameter related to the launch procedure. If default = \code{TRUE} the launch is partial (see 'Details').} 
}

\details{
It is possible to launch the procedure in two separate steps, so that user can see a first output (launch = \code{TRUE}) and decide if modify the input parameters or continue the procedure setting launch = \code{FALSE}.
}

\value{
An object of class \code{list} or \code{data.frame} depending from argument \code{launch}. If launch = \code{TRUE} the only output is, at domain level, a data frame of the Self Representative (SR) PSUs and of the Non Self Representative (NSR) PSUs before stratification. If launch = \code{FALSE} the output is a \code{list} composed by four members of \code{data.frame} class. The first component, at domain level, is the same data frame obtained if launch = \code{TRUE}. The second element
is, at domain level, a data frame of SR and NSR PSUs after stratification, with their totals like the first output, but with additional information such as the final units sample size distinctly for SR and NSR PSUs and their mean. The third component is a data frame that supply, for each stratum, the number of sample PSUs selected and the total number of PSUs. The fourth element provides, for each PSU, some information like the inclusion probability and the sampling fraction.
}

\author{Raffaella Cianchetta
}

\references{ 
S. Falorsi A. Russo (2001), Il disegno di rilevazione per indagini Panel sulle famiglie, \emph{Rivista di Statistica Ufficiale, N. 3}, p. 55-90.\cr
Sampford, M. (1967), On sampling without replacement with unequal probabilities of selection, \emph{Biometrika}, 54:499-513.\cr
Yves Tille' and Alina Matei (2012). sampling: Survey Sampling. R package version 2.5. \cr
http://CRAN.R-project.org/package=sampling.\cr
}

\examples{
# Start StratSel
\dontrun{
load(FS4)
data(population)
data(allocation)

# The function with a fixed number of final sample units (min_sample= 8) 
# to observe in each PSU and partial launch of the procedure with
# only one output data frame(see 'Value')   
Output_list <- StratSel(dataPop= population, idpsu= ~ comune_num, dom= ~ dom,
               final_pop= ~ fam, size= ~ pop, PSUsamplestratum= 6, min_sample= 8,
               min_sample_index= FALSE, dataAll= allocation, domAll= ~ dom,
               f_sample= ~ campdom, planned_min_sample= NULL, launch= TRUE)
			   
# It can also be written as below due to default values:
Output_list <- StratSel(dataPop= population, idpsu= ~ comune_num, dom= ~ dom,
               final_pop= ~ fam, size= ~ pop, PSUsamplestratum= 6, min_sample= 8, 
               dataAll= allocation, domAll= ~ dom, f_sample= ~ campdom)


# The function with a fixed number of final sample units (min_sample= 8)
# to observe in each PSU and full launch of the procedure with 
# the output list composed of four data frames(see 'Value')   
Output_list <- StratSel(dataPop= population, idpsu= ~ comune_num, dom= ~ dom,
               final_pop= ~ fam,  size= ~ pop, PSUsamplestratum= 6, min_sample= 8, 
               dataAll= allocation, domAll= ~ dom, f_sample= ~ campdom, launch= FALSE)
			   

# The function with a variable number of final sample units (planned_min_sample= 
# ~ planned_final_sample) and partial launch of the procedure
Output_list <- StratSel(dataPop= population, idpsu= ~ comune_num, dom= ~ dom, 
               final_pop= ~ fam, size= ~ pop, PSUsamplestratum= 6, min_sample= NULL,
               min_sample_index= TRUE, dataAll= allocation, domAll= ~ dom,
               f_sample= ~ campdom, planned_min_sample= ~ planned_final_sample)

			   
# The function with a variable number of final sample units (planned_min_sample=
# ~ planned_final_sample) and full launch of the procedure with the output list 
#composed of four data frames		
Output_list <- StratSel(dataPop= population, idpsu= ~ comune_num, dom= ~ dom, 
               final_pop= ~ fam, size= ~ pop, PSUsamplestratum= 6, min_sample= NULL,
               min_sample_index= TRUE, dataAll= allocation, domAll= ~ dom,
               f_sample= ~ campdom, planned_min_sample= ~ planned_final_sample,
               launch= FALSE)	   
}
}

\keyword{survey}