#' cvs_hint
#'  
#' @description
#' Function to propose a convenient set of precision constraints,
#' given the different variability of target variables Ys in the strata
#' 
#' @details
#' The function requires in input the 'strata' dataset, plus an initial
#' set of precision constraints. It is suggested to define a 'neutral'
#' set, i.e. equal CVs for all variable, possibly differentiated
#' only with regard to the different domain levels
#'  
#' 
#' @param strata the 'strata' dataset.
#' @param errors the 'errors' dataset containing the current precision constraints.
#' 
#' @return the new 'errors' dataset containing the changed precision constraints
#' 
#' @examples
#' data(beat.example)
#' errors[1,c(2:3)] <- c(0.03,0.03)
#' errors[2,c(2:3)] <- c(0.03,0.03)
#' errors
#' cv <- CVs_hint(strata,errors)
#' cv
CVs_hint <- function(strata,cv) {
  ndom <- length(grep("DOM",colnames(strata)))
  nvar <- ncol(cv) - 1
  S <- paste0("S",c(1:nvar))
  M <- paste0("M",c(1:nvar))
  cv_new <- cv
  for (i in c(1:ndom)) {
    eval(parse(text=paste0("vett <- list(as.numeric(as.factor(as.character(strata$DOM",i,"))))")))
    vett <- unlist(vett)
    a <- aggrStrata(strata,vett,dominio = 1,nvar=nvar,censiti="CENS")
    cvs <- a[,c(S)] / a[,c(M)]
    cv_new[i,c(2:(nvar+1))] <- cv_new[i,c(2:(nvar+1))]*apply(cvs,2,FUN=mean)
  }
  cv_new
}
