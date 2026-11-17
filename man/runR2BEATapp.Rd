\name{runR2BEATapp}
\alias{runR2BEATapp}
\title{Launch the R2BEAT Shiny application}
\usage{
runR2BEATapp(launch.browser = TRUE, ...)
}
\arguments{
\item{launch.browser}{Passed to \code{shiny::runApp}.}
\item{...}{Further arguments passed to \code{shiny::runApp}.}
}
\value{
Called for its side effect (runs the app).
}
\description{
Starts the interactive R2BEAT application, giving access to three workflows:
CV optimisation for one-stage allocation (Pareto and genetic algorithm),
two-stage design from a sampling frame, and two-stage design from a previous
survey (ReGenesees objects). If any of the required packages (shiny, shinyjs,
shinydashboard, DT, ggplot2, genalg) is not installed, it is installed from
CRAN first.
}
\examples{
\dontrun{
runR2BEATapp()
}
}
