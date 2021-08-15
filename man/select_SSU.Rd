\name{select_SSU}
\alias{select_SSU}
\title{
Select sample of secondary stage units (SSU)
}
\description{
Select sample of secondary stage units (SSU) from the population frame on the basis of the SSU allocated to each selected PSU
}
\usage{
     select_SSU(df,PSU_code,SSU_code,PSU_sampled,verbose)
}

\arguments{
 \item{df}{Dataframe containing sampling units (SSUs)}
 \item{PSU_code}{Identifier of each PSU in dataframe containing sampling units}
 \item{SSU_code}{Identifier of each SSU in dataframe containing sampling units}
 \item{PSU_sampled}{Dataframe containing selected PSUs}
 \item{verbose}{Flag for producing output information while executing}
}

\value{
A dataframe containing the units selected in the sample
}

\examples{
\dontrun{
samp <- select_SSU(df=pop,
                   PSU_code="municipality",
                   SSU_code="id_ind",
                   PSU_sampled=selected_PSU)
}
}
\author{Giulio Barcaroli}
