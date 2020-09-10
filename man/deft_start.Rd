\name{deft_start}
\Rdversion{1.0}
\alias{deft_start}
\docType{data}
\title{
Starting values for the Design Effect (\eqn{deft}) 
}
\description{
Example data frame containing the starting values for the Design Effect (\eqn{deft}) .
}
\usage{
data(beat.example)
}
\format{
  The Design Effect data frame contains a row per each stratum with the following variables: 
  \describe{
    \item{STRATUM}{
    Identifier of the stratum (numeric).
    } 
    \item{DEFT1}{
    Starting values for the Design Effect in the stratum of the first variable (numeric).
    }
    \item{DEFTj}{
    Starting values for the Design Effect in the stratum of the j-th variable (numeric).
    }
    \item{DEFTn}{
    Starting values for the Design Effect in the stratum of the last variable (numeric).
    }
  }
 }
\details{
Note: the names of the variables must be the ones indicated above.

This is an optional input. 
The function \code{beat.2st} independently computes and updates the design effect.
However, it is possible to set the starting values of design effect for each variable in each stratum.
The design effect is the square root of the ratio of the actual sampling variance to the variance expected with the simple random sampling (SRS), on equal sample size.
%The \out{<img src="https://latex.codecogs.com/gif.latex?deft&space;=&space;\sqrt{deff}" title="deft = \sqrt{deff}" />}, where 
%\out{<img src="https://latex.codecogs.com/gif.latex?deff=&space;\frac{\text{var}_{design}&space;\left(&space;\hat{Y}&space;\right)}{\text{var}_{SRS}&space;\left(&space;\hat{Y}&space;\right)}" title="deff= \frac{\text{var}_{design} \left( \hat{Y} \right)}{\text{var}_{SRS} \left( \hat{Y} \right)}" />}
%Therefore, 
Under SRS the desing effect is equal to 1. 
Usually, as increasing the stages of selection the design effect increases because it takes into account the "clusterization" of sampling units and the sample size in Self Representative (SR) and Non Self Representative (NSR) strata.
%Expanding the expression of \eqn{deff} this it looks more clear:
%\deqn{deff = \frac{n}{N} \left\lbrace \frac{N_{SR}^2}{n_{SR}^2} \left[1 + \rho_{SR} \left( b_{SR}^* -1 \right) \right] + 
%\frac{N_{NSR}^2}{n_{NSR}^2} \left[1 + \rho_{NSR} \left( b_{NSR}^* -1 \right) \right]\right\rbrace},
%where \out{<img src="https://latex.codecogs.com/gif.latex?N_{\cdot}" title="N_{\cdot}" />} is the number of SSU in the SR or NSR part of the stratum and and \out{<img src="https://latex.codecogs.com/gif.latex?n_{\cdot}" title="n_{\cdot}" />} the sample size. Furthermore, \eqn{<img src="https://latex.codecogs.com/gif.latex?b_{\cdot}^*" title="b_{\cdot}^*" />}
%is the average minimum number of interviews in the Primary Stage Units (PSU) in SR and NSR part of the stratum. Finally, \out{<img src="https://latex.codecogs.com/gif.latex?\rho_{\cdot}" title="\rho_{\cdot}" />} are the intraclass correlation and they meausure the "clusterization" of sampling units in the SR and NSR of the strata (for more details see \link{rho}).
In practice, higher is the intraclass correlation, higher will be the design effect and much more sample size for satisfying the precision constraints is needed with respect to SRS. 
}
\keyword{datasets}

\examples{
# Load example data
data(beat.example)
deft_start
str(deft_start)
}

