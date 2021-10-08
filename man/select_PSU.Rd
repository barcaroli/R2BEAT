\name{select_PSU}
\alias{select_PSU}
\title{
Select sample of primary stage units (PSU)
}
\description{
Select sample of primary stage units (PSU) on the basis of the PSU allocated in the allocation step.
}
\usage{
     select_PSU(alloc, type="ALLOC", pps=TRUE, plot=TRUE)
}

\arguments{
 \item{alloc}{Output of the allocation step (beat.st)}
 \item{type}{Type of SSU allocation ("ALLOC" = optimal, "PROP" = proportional to population size, "EQUAL" = equal size in each stratum)}
 \item{pps}{Type of PSU selection in strata (pps = TRUE --> proportional to size, pps= FALSE--> simple random sampling)}
  \item{plot}{If TRUE, a plot of PSUs and SSUs ditribution is produced}
}

\value{
A list including: universe_PSU, sample_PSU
}

\examples{
\dontrun{
PSU_selection <- select_PSU(alloc, type="ALLOC", pps=TRUE, plot=TRUE)
}
}
\author{Alessio Guandalini}
