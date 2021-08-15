\name{eval_2stage}
\Rdversion{1.4}
\alias{eval_2stage}
\title{
Evaluation of the solution produced by the optimization steps of R2BEAT by
selecting a number of samples from the frame with the optimal stratification, and calculating average CV's on the target variables Y's.
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
            writeFiles = TRUE,
            progress = TRUE) 
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
	\item{progress}{
  If set to TRUE, a progress bar is visualised during the execution. Default is TRUE.
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
library(R2BEAT)
load("./data/pop.RData")
cv <- as.data.frame(list(DOM=c("DOM1","DOM2"),
                         CV1=c(0.03,0.04),
                         CV2=c(0.06,0.08),
                         CV3=c(0.06,0.08),
                         CV4=c(0.06,0.08)))
samp_frame <- pop
id_PSU <- "municipality"  
id_SSU <- "id_ind"        
strata_var <- "stratum"   
target_vars <- c("income_hh","active","inactive","unemployed")   # more than one
deff_var <- "stratum"     
domain_var <- "region"  
delta =  1     # households = survey units
minimum <- 50  # minimum number of SSUs to be interviewed in each selected PSU
f = 0.05         # suggestion for the sampling fraction 
deff_sugg <- 1.5  # suggestion for the deff value
inp <- prepareInputToAllocation(samp_frame,
                                id_PSU,
                                id_SSU,
                                strata_var,
                                target_vars,
                                deff_var,
                                domain_var,
                                minimum,
                                delta,
                                f,
                                deff_sugg)
alloc <- beat.2st(stratif = inp$strata, 
                  errors = cv, 
                  des_file = inp$des_file, 
                  psu_file = inp$psu_file, 
                  rho = inp$rho, 
                  deft_start = NULL, 
                  effst = inp$effst,
                  epsilon1 = 5, 
                  mmdiff_deft = 1,maxi = 15, 
                  epsilon = 10^(-11), 
                  minnumstrat = 2, 
                  maxiter = 200, 
                  maxiter1 = 25)
allocat <- alloc$alloc[-nrow(alloc$alloc),]
sample_2st <- StratSel(dataPop= inp$psu_file,
                       idpsu= ~ PSU_ID, 
                       dom= ~ STRATUM, 
                       final_pop= ~ PSU_MOS, 
                       size= ~ PSU_MOS, 
                       PSUsamplestratum= 1, 
                       min_sample= minimum, 
                       min_sample_index= FALSE, 
                       dataAll=allocat,
                       domAll= ~ factor(STRATUM), 
                       f_sample= ~ ALLOC, 
                       planned_min_sample= NULL, 
                       launch= F)
selected_PSU <- sample_2st[[4]]
df=pop
PSU_code="municipality"
SSU_code="id_ind"
PSU_sampled=selected_PSU[selected_PSU$Sampled_PSU==1,]
target_vars <- c("income_hh",
                 "active",
                 "inactive",
                 "unemployed")  
domain_var <- "region"
PSU_sampled <- selected_PSU[selected_PSU$PSU_final_sample_unit > 0,]
eval <- eval_2stage(df,
                    PSU_code,
                    SSU_code,
                    domain_var,
                    target_vars,
                    PSU_sampled,
                    nsampl=10, 
                    writeFiles=TRUE,
                    progress=TRUE) 
eval$coeff_var
eval$rel_bias

}
}
\keyword{ survey }