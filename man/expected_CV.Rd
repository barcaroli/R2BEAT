\name{expected_CV}
\alias{expected_CV}
\title{expected_CV}
\usage{
expected_CV(strata, errors, alloc)
}
\arguments{
\item{strata}{name of the dataframe containing information in the sampling strata.}

\item{errors}{name of the dataframe}

\item{alloc}{vector containing the allocation of sampling units.}
}
\value{
a dataframe containing the maximum expected coefficients of variation in each domain level for each target variable
}
\description{
Function to report the expected coefficients of variation for target variables Ys in a 'strata' dataset 
given an allocation 'alloc' and the current set of precision constraints
}
\examples{
load("./data/sample.RData")
target <- c("active","inactive","unemployed","income_hh")
alloc <- prepareInputToAllocation_beat.1st(frame = samp,
                                                     ID = "id_hh",
                                                     stratum = "stratum_label",
                                                     dom = "region",
                                                     target = target)
strata <- alloc$file_strata

cv <- as.data.frame(list(DOM = c("DOM1","DOM2"),
                         CV1 = c(0.05,0.10),
                         CV2 = c(0.05,0.10),
                         CV3 = c(0.05,0.10),
                         CV4 = c(0.05,0.10)))
allocation <- beat.1st(strata,cv)

alloc <- allocation$alloc$OPT[-nrow(allocation$alloc)]
exp_cv <- expected_CV(strata,cv,alloc)
exp_cv
}
