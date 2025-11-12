\name{eval_2stage}
\Rdversion{1.4}
\alias{eval_2stage}
\title{
Evaluation of the two-stage sample design optimized solution by simulation
}
\description{
The user can indicate the number of samples that must be selected by the optimized frame. First, the true values of the parameters are calculated from the frame. Then, for each sample the sampling estimates are calculated, together with the differences between them and the true values of the parameters.
At the end, an estimate of the CV is produced for each target variable, in order to compare them
with the precision constraints set at the beginning of the optimization process.
If the flag 'writeFiles' is set to TRUE, boxplots of distribution of the CV's in the different 
domains are produced for each Y variable ('cv.pdf'),
together with boxplot of the distributions of differences between estimates and values of the parameters
in the population ('differences.pdf').

}
\usage{
eval_2stage(df,
            PSU_code,
            SSU_code,
            domain_var,
            target_vars,
            PSU_sampled,
            nsampl = 100, 
            writeFiles = TRUE) 
}
\arguments{
  \item{df}{
The sampling frame.
	}
  \item{PSU_code}{
In the sampling frame, the identifier of the PSU.
	}
	\item{SSU_code}{
In the sampling frame, the identifier of the SSU.
	}
	\item{domain_var}{
In the sampling frame, the identifier of the domain of interest for the estimates.
	}
	\item{target_vars}{
In the sampling frame, the variables used to produce the target estimates.
	}
	\item{PSU_sampled}{
The set of selected PSUs.
	}
  \item{nsampl}{
The number of samples to be drawn from the frame.
	}
  \item{writeFiles}{
A flag to write in the work directory the outputs of the function. Default is TRUE.
	}
}
\value{
A list containing (i) the CV distribution in the domains, (ii) the bias distribution in the domains, (iii) the dataframe containing the sampling
estimates by domain
}

\author{
Giulio Barcaroli
}

\examples{
\dontrun{
library(readr)
pop <- read_rds("https://github.com/barcaroli/R2BEAT_workflows/blob/master/pop.RDS?raw=true")
library(R2BEAT)
cv <- as.data.frame(list(DOM=c("DOM1","DOM2"),
                         CV1=c(0.02,0.03),
                         CV2=c(0.03,0.05),
                         CV3=c(0.03,0.05),
                         CV4=c(0.05,0.08)))
frame <- pop
id_PSU <- "municipality"  
id_SSU <- "id_ind"        
stratum <- "stratum"   
target <- c("income_hh","active","inactive","unemployed")   # more than one
deff_level <- "stratum"     
domain <- "region"  
delta =  1     # households = survey units
minimum <- 50  # minimum number of SSUs to be interviewed in each selected PSU
deff_sugg <- 1.5  # suggestion for the deff value
inp <- prepareInputToAllocation1(frame,
                                 id_PSU,
                                 id_SSU,
                                 stratum,
                                 target,
                                 deff_level,
                                 domain,
                                 minimum,
                                 delta,
                                 deff_sugg)
inp$desfile$MINIMUM <- 50
alloc <- beat.2st(file_strata = inp$strata, 
                  errors = cv, 
                  des_file = inp$des_file, 
                  psu_file = inp$psu_file, 
                  rho = inp$rho, 
                  deft_start = NULL,
                  effst = inp$effst, 
                  minPSUstrat = 2,
                  minnumstrat = 50
)
sample_1st <- select_PSU(alloc, type="OPT", pps=TRUE)

df=pop
df$one <- 1
PSU_code="municipality"
SSU_code="id_ind"
target <- c("income_hh",
                 "active",
                 "inactive",
                 "unemployed")  

PSU_sampled <- sample_1st$sample_PSU
eval <- eval_2stage(df,
                    PSU_code,
                    SSU_code,
                    domain,
                    target,
                    PSU_sampled,
                    nsampl=10, 
                    writeFiles=FALSE) 
eval$coeff_var
eval$rel_bias

}
}
\keyword{ survey }