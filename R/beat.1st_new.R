#' Compute multivariate optimal allocation for different domains in one-stage stratified sample design
#'
#' @description 
#' This function computes the multivariate optimal allocation for different domains 
#' in a one-stage file_strataied sample design, following the generalized Bethel algorithm.
#' 
#' @usage
#' beat.1st(file_strata, errors, var_label=NULL, minnumstrat = 2, maxiter = 200, maxiter1 = 25, epsilon = 1e-11)
#'
#' @param file_strata: dataframe of survey strata, for more details see, e.g.,\code{\link{strata}}. 
#' @param errors: dataframe of expected coefficients of variation (CV) for each domain, for more details see, e.g.,\code{\link{errors}}.
#' @param var_label: a vector listing the target variable(s) which are used for leading the optimal allocation.
#' @param minnumstrat: Minimum number of elementary units per strata (default=2).
#' @param maxiter: Maximum number of iterations (default=200) of the general procedure. This kind of iteration may be required by the fact that when in a stratum the number of allocated units is greater or equal to its population, that stratum is set as "census stratum", and the whole procedure is re-initialised.
#' @param maxiter1: Maximum number of iterations in Chromy algorithm (default=25).
#' @param epsilon: Tollerance for the maximum absolute differences between the expected CV and the realised CV with the allocation obtained in the last iteraction for all domains. The default is 10^(-11).
#' 
#' @return 
#' It returns a list with the following objects: 
#' 
#' - a dataframe, *file_strata*, with the strata information that is one of the inputs for performing the multivariate and multidomain optimal allocation implemented in the function *beat.1st* in R2BEAT R-package; 
#' 
#' - a vector, *var_label*, listing the target variable(s) which will be used for leading the optimal allocation.
#'
#' @details
#' The methodology is a generalization of Bethel multivariate allocation (1989) that extends 
#' the Neyman (1959) - Tchuprov (1923) allocation for multi-purpose and multi-domain surveys. 
#' The generalized Bethel’s algorithm allows to determine the optimal sample size for each stratum 
#' in a stratified sample design. The overall sample size and the allocation among the different strata 
#' is determined starting from the accuracy constraints imposed in the survey on interest estimates.
#'
#' @return
#' Returns a list with four components:
#' - n: Vector with the optimal sample size for each stratum 
#'  
#' - file_strata: dataframe corresponding to the input data.frame file_strata with the optimal sample size column added 
#' 
#' - alloc: dataframe with optimal (OPT), proportional (PROP), and equal (UNIF) saple size allocation. 
#' 
#' - sensitivity: dataframe with the planned coefficients of variation (PlannedCV), the expected coefficients of variation (ExpectedCV), and the sensitivity at 10% for each variable and each domain category. Sensitivity provides an indication of sample size variation for a 10% change in planned CVs.
#'  
#' - expectedCV: dataframe with the maximum of the expected coefficients of variation (Expected CV), for each variable in each domain. 
#' 
#' - plannedCV: dataframe with the maximum coefficients of variation admissible for each domain and for each variable. It is the input, error, provided by the user.
#' 
#' - param_alloc: a vector summarising all the parameters used for performing the optimal allocation.
#'
#' @examples
#' # Load example data
#'data(beat.example)

## Example 1
# Allocate the sample
#'allocation_1 <- beat.1st(file_strata=strata, errors=errors)
# The total sample size is
#'sum(allocation_1$n)

## Example 2
# Assume 5700 units is the maximum sample size to stick to our budget.
# Looking at allocation_1$sensitivity we can see that most of the 
# sensitivity is in DOM1 for REG1 and REG2 due to V1.
#'allocation_1$sensitivity
# We can relax the constraints increasing the expected coefficients of variation for X1 by 10%
#'errors1 <- errors 
#'errors1[1,2] <- errors[1,2]+errors[1,2]*0.1
# Try the new allocation 
#'allocation_2 <- beat.1st(file_strata=strata, errors=errors1)
#'sum(allocation_2$n)

## Example 3
# On the contrary, if we tighten the constraints decreasing the expected coefficients of variation 
# for X1 by 10%
#'errors2 <- errors 
#'errors2[1,2] <- errors[1,2]-errors[1,2]*0.1
# The new allocation leads to a larger sample than the first example 
#'allocation_3 <- beat.1st(file_strata=strata, errors=errors2)
#'sum(allocation_3$n)

#' @export


beat.1st <- function (file_strata, errors, var_label=NULL, minnumstrat = 2, maxiter = 200, maxiter1 = 25, 
          epsilon = 10^(-11)) 
{
    param_alloc <- as.data.frame(t(c(epsilon, minnumstrat, 
                                     maxiter, maxiter1, 1)))
    names(param_alloc) = c("epsilon", "minnumstrat", 
                           "maxiter", "maxiter1", "num_stages")
  colnames(file_strata) <- toupper(colnames(file_strata))
  colnames(errors) <- toupper(colnames(errors))
  iter1 <- 0
  val = NULL
  m <- NULL
  s <- NULL
  cv = NULL
  nstrat = nrow(file_strata)
  nvar = ncol(errors) - 1
  ndom = nrow(errors)
  # check
  if(sum(grepl("^M[0-9]+$", names(file_strata)))!=nvar){
    msg = "The number of variables in file_strata do not match with the number of variables in the errors dataframe"
    stop(msg)
  }
  if(sum(grepl("^DOM[0-9]+$", names(file_strata)))!=ndom){
    msg = "The number of domains in file_strata do not match with the number of domains in the errors dataframe"
    stop(msg)
  }
  # parameters
  param_alloc$nvar <- nvar
  param_alloc$ndom <- ndom
  param_alloc$nstrat <- nstrat
  varloop <- c(1:nvar)
  strloop <- c(1:nstrat)
  domloop <- c(1:ndom)
  med <- as.matrix(file_strata[, names(file_strata) %in% sapply(1:nvar, 
                                                        function(i) paste("M", i, sep = ""))])
  esse <- as.matrix(file_strata[, names(file_strata) %in% sapply(1:nvar, 
                                                         function(i) paste("S", i, sep = ""))])
  nom_dom <- sapply(1:ndom, function(i) paste("DOM", i, sep = ""))
  dom <- as.vector(file_strata[, names(file_strata) %in% nom_dom])
  N <- as.vector(file_strata$N)
  cens <- as.vector(file_strata$CENS)
  cost <- as.vector(file_strata$COST)
  nocens = 1 - cens
  dom_labels=list()
  if (ndom == 1){
    (nvalues <- nlevels((as.factor(file_strata$DOM1))))
    dom_labels[[1]]=levels(as.factor(file_strata$DOM1))
  }
  if (ndom > 1) {
    nvalues <- NULL
    for (i in c(1:ndom)) {
      eval(parse(text = paste0("nvalues <- c(nvalues,nlevels(as.factor(file_strata$DOM", 
                               i, ")))")))
      eval(parse(text = paste0("dom_labels[[i]] <- levels(as.factor(file_strata$DOM", 
                               i, "))")))
    }
  }
  names(dom_labels)<-paste("DOM", 1:ndom, sep="")
  # allocation algorithm
  crea_disj = function(data, vars) {
    out = NULL
    sapply(vars, function(vari) {
      col = as.factor(data[, vari])
      out <<- cbind(out, outer(col, levels(col), function(y, 
                                                          x) ifelse(y == x, 1, 0)))
    })
    out
  }
  disj <- crea_disj(file_strata, nom_dom)
  nc <- ncol(disj)
  for (i in 1:nc) {
    m <- cbind(m, disj[, i] * med)
    s <- cbind(s, disj[, i] * esse)
  }
  for (k in domloop) {
    cvx <- as.matrix(errors[k, names(errors) %in% sapply(1:nvar, 
                                                         function(i) paste("CV", i, sep = ""))])
    ndomvalues <- c(1:nvalues[k])
    for (k1 in ndomvalues) {
      cv <- cbind(cv, cvx)
    }
  }
  alfa2 <- NULL
  nvar_orig <- nvar
  nvar <- ncol(cv)
  varloop <- c(1:nvar)
  NTOT <- c(rep(0, nvar))
  CVfin <- c(rep(0, nvar))
  varfin <- c(rep(0, nvar))
  totm <- c(rep(0, nvar))
  crea_a = function() {
    numA <- (N^2) * (s^2) * nocens
    denA1 <- colSums(t(t(N * m) * c(cv)))^2
    denA2 <- colSums(N * (s^2) * nocens)
    denA <- denA1 + denA2 + epsilon
    a <- t(t(numA)/denA)
    return(a)
  }
  chromy = function(alfatot, diff, iter, alfa, alfanext, x) {
    while (diff > epsilon && iter < maxiter) {
      iter <- iter + 1
      den1 = sqrt(rowSums(t(t(a) * c(alfa))))
      den2 = sum(sqrt(rowSums(t(t(a * cost) * c(alfa)))))
      x <- sqrt(cost)/(den1 * den2 + epsilon)
      alfatot <- sum(c(alfa) * (t(a) %*% x)^2)
      alfanext <- c(alfa) * (t(a) %*% x)^2/alfatot
      diff <- max(abs(alfanext - alfa))
      alfa <- alfanext
      alfa2 <<- alfanext
    }
    n <- ceiling(1/x)
    return(n)
  }
  a <- crea_a()
  n <- chromy(0, 999, 0, c(rep(1/nvar, nvar)), c(rep(0, nvar)), 
              array(0.1, dim = c(nstrat, 1)))
  contx <- sum(n > N)
  cens[n > N] <- 1
  nocens = 1 - cens
  for (i in strloop) {
    if (n[i] < minnumstrat) {
      n[i] <- min(minnumstrat, N[i])
    }
  }
  while (contx > 0 && iter1 < maxiter1) {
    iter1 = iter1 + 1
    a <- crea_a()
    n <- chromy(0, 999, 0, c(rep(1/nvar, nvar)), c(rep(0, 
                                                       nvar)), array(0.1, dim = c(nstrat, 1)))
    contx <- sum(n > N)
    cens[n > N] <- 1
    nocens = 1 - cens
    for (i in strloop) {
      if (n[i] < minnumstrat) {
        n[i] <- min(minnumstrat, N[i])
      }
    }
  }
  n = (nocens * n) + (cens * N)
  num_strati <- length(N)
  sampleSize <- sum(n)
  popSize <- sum(N)
  sampleSize_nocens <- sum(n[which(file_strata$CENS == 0)])
  popSize_nocens <- sum(N[which(file_strata$CENS == 0)])
  uniform <- ifelse(cens == 0, sampleSize_nocens/sum(file_strata$CENS == 
                                                      0), N)
  proportional = ifelse(cens == 0, sampleSize_nocens * N/popSize_nocens, 
                         N)
  Bethel_sample <- cbind(file_strata, n)
  colnames(Bethel_sample)[length(colnames(Bethel_sample))] <- "n"
  df = NULL
  df <- data.frame(STRATUM=as.character(file_strata$STRATUM), OPT=n, PROP=proportional, 
              UNIF=uniform)
  totdf=data.frame(STRATUM="Total", OPT=sum(n), PROP=sum(proportional), 
                   UNIF=sum(uniform))
  df <- rbind(df, totdf)
  # expected cv
  calcola_cv <- function() {
    NTOT <- c(rep(0, nvar))
    CVfin <- c(rep(0, nvar))
    NTOT <- colSums((m > 0) * N)
    varfin <- rowSums(t((s * N)^2 * (1 - round(n)/N)/round(n))/NTOT^2)
    totm <- rowSums(t(m * N))
    CVfin <- round(sqrt(varfin/(totm/NTOT)^2), digits = 4)
    return(CVfin)
  }
  CVfin <- calcola_cv()
  # sensitivity
  g <- 0
  for (i in strloop) {
    t <- 0
    for (j in varloop) {
      t <- t + alfa2[j] * a[i, j]
    }
    g <- g + sqrt(cost[i] * t)
  }
  g <- g^2
  sens <- 2 * 0.1 * alfa2 * g
  # prepare output
  j <- 0
  outcv <- NULL
  outcv=data.frame(matrix(nrow=sum(nvalues*nvar_orig),ncol=7))
  colnames(outcv) <- c("Type", "Dom", "Dom_label", "Var", "PlannedCV", 
                       "ExpectedCV", "Sensitivity10%")
  for (k in (1:ndom)) {
    for (k1 in 1:(nvalues[k])) {
      for (k2 in 1:nvar_orig) {
        j <- j + 1
        outcv[j,1] <- nom_dom[k]
        outcv[j,2] <- k1
        outcv[j,3] <- dom_labels[[k]][k1]
        outcv[j,4] <- paste("V", k2, sep = "")
        outcv[j,5] <- cv[j]
        outcv[j,6] <- CVfin[j]
        outcv[j,7] <- ceiling(sens[j])
      }
    }
  }
  if(!is.null(var_label)){
    if(length(var_label)==nvar_orig){
      outcv$Var_label <- rep(var_label, times=sum(nvalues))
      outcv=outcv[,c("Type", "Dom", "Dom_label", "Var","Var_label", "PlannedCV", 
      "ExpectedCV", "Sensitivity10%")]
    }
    else{
      warning("The number of elements in var_label does not match with the number of variables used in the allocation", call.=FALSE)
    }
  }
  # expectedCV output
  exp_cv=aggregate(ExpectedCV ~ Type+Var, outcv, max)
  exp_cv=reshape(exp_cv, 
          direction = "wide",
          idvar = "Type",
          timevar = "Var")
  colnames(exp_cv) <- colnames(errors)
  # output
  out <- list(n = n, file_strata = Bethel_sample, 
              alloc = as.data.frame(df), 
              sensitivity = outcv, 
              expectedCV=exp_cv, plannedCV=errors, 
              param_alloc=param_alloc)
  return(out)
}
