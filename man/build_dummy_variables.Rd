\name{build_dummy_variables}
\Rdversion{1.5}
\alias{build_dummy_variables}
\title{
Derives dummy variables from target variables to the sampling frame
}
\description{
This function allow to derive dummy variables to the sampling frame (this is necessary if we want to differentiate precision constraints in the same domain level, for instance in the different regions).

}
\usage{
build_dummy_variables(frame,
                      domain_var,
                      initial_target_vars,
                      cv) 
}
\arguments{
  \item{frame}{
The sampling frame.
	}
  \item{domain_var}{
In the sampling frame, the variable that identifies a given domain level.
	}
	\item{initial_target_vars}{
In the sampling frame, the initial set of target variables (non dummies)
	}
	\item{cv}{
The initial set of precision constraints, related to the initial set of target variables.
	}
}
\value{
A list containing (i) frame: the sampling frame with the derived dummy variables, (ii) new_target_variables: a vector with the names of the total set of target variables (initial ones plus derived variables), (iii) cvnew: the dataframe containing the new set of precision constraints (related to the initial target variables plus the derived dummy variables)
}

\author{
Giulio Barcaroli
}

\examples{
\dontrun{
}
}
\keyword{ survey }