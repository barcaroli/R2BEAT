#' The function returns a dataframe with planned and actual coefficients of variation (CV) in a multivariate multi-domain allocation problem.
#'
#' @param stratif A dataframe with strata information.
#' @param errors: A dataframe with planned CV.
#' @param alloc: A vector with an allocation into strata.
#' @param minnumstrat: Optional. A number indicating the lower bound for the number of allocated units in each stratum (default is 2).
#' @return A dataframe with planned and actual CV for each combination of domain, domain category and auxiliary variable.
# ----------------------------------------------------------------------
beat.1cv <- function(stratif, errors,alloc, minnumstrat = 2) 
{

  # Begin body
  # First input data frame
  colnames(stratif) <- toupper(colnames(stratif))
  
  # Second input data frame
  colnames(errors) <- toupper(colnames(errors))
  #
  #check inputs:
  #
  checkData(strata = stratif, errors = errors)
  #
  check_n(allocation=alloc, strata=stratif, minnumstrat)
  #
  ordina_variabili <- function(dati, prefisso, n_var) {
    #order df columns,
    #s.t. we have prefix_1_var_1, prefix_1_var_2, ..., prefix_n_var_1, prefix_n_var_n
    if (!is.data.frame(dati)) 
      stop()
    as.matrix(dati[, paste(prefisso, 1:n_var, sep = ""), 
                   drop = FALSE])
  }
  #
  # --------------------------------------------------------------
  # Initialization of parameters. Attribution of initial
  # values
  # ------------------------------------------------------------
  
  val <- NULL #to store number of categories of each domain
  m <- NULL #to store the means
  s <- NULL #to store the standard deviation
  cv <- NULL #to store the cv
  
  nstrat <- nrow(stratif)
  nvar <- length(grep("CV", names(errors)))
  ndom <- nrow(errors)

  
  # ---------------------------------------------------------
  # Initial Data Structures - selection from input variables
  # ---------------------------------------------------------
  
  # means
  med <- ordina_variabili(stratif, "M", nvar)
  
  # variances estimates
  esse <- ordina_variabili(stratif, "S", nvar)
  
  #check
  if (ncol(med) != ncol(esse)) 
    stop(print("Error: Number of M variables don't match the number of S variables"))
  if (ncol(med) != nvar) 
    stop(print("Error: Number of variables don't match the number of planned CV"))
  
  
  # populations
  N <- as.vector(stratif$N)
  # vector cens
  cens <- as.vector(stratif$CENS)
  # strata with N < minnumstrata are set to take-all strata
  cens[N < minnumstrat] <- 1
  # vector cost
  cost <- as.vector(stratif$COST)
  # check and default for cost and cens
  if (is.null(cost)) 
    cost <- rep(1, nstrat)
  if (is.null(cens)) 
    cens <- rep(0, nstrat)
  nocens <- 1 - cens
  # check variable cens
  if (sum(cens) == length(cens)) {
    warning(print("Warning: Variable CENS always equal 1"))
  }
  
  # domains
  name_dom <- sapply(1:ndom, function(i) paste("DOM", i, sep = ""))
  dom <- ordina_variabili(stratif, "DOM", ndom)
  
  # numbers of different domains of interest (numbers of
  # modalities/categories for each type of domain) if (ndom
  # ==1) (nvalues<-nlevels((as.factor(stratif$DOM1)))) if
  # (ndom>1) {nvalues<-sapply(nom_dom, function(vari)
  # {val<-c(val, nlevels(as.factor(dom[,vari])))})}
  nvalues <- sapply(name_dom, function(vari) {
    val <- c(val, nlevels(as.factor(dom[, vari])))
  })
  
  # -------------------------------------------------
  # disjunctive matrix
  # -------------------------------------------------
  
  crea_disj <- function(data, vars) {
    out <- NULL
    sapply(vars, function(vari) {
      col <- as.factor(data[, vari])
      out <<- cbind(out, outer(col, levels(col), 
                               function(y,x) 
                                 {
                                 ifelse(y == x, 1, 0)
                                 }))
    })
    out
  }
  
  disj <- crea_disj(stratif, name_dom)
  
  # -------------------------------------------------
  # Building means (m) and deviations matrices (s) for
  # different domains of interest
  # -------------------------------------------------
  #(m) and (s)
  nc <- ncol(disj)
  for (i in 1:nc) {
    m <- cbind(m, disj[, i] * med)
    s <- cbind(s, disj[, i] * esse)
  }
  
  #
  # -------------------------------------------------------------
  # computation of the coefficients of variation CVs for
  # different domains of interest
  # -------------------------------------------------------------
 
  cvDom <- NULL  #to be printed
  cvDom2 <- NULL  #to be printed
  for (k in 1:ndom) {
    cvx <- ordina_variabili(errors[k, ], "CV", nvar)
    ndomvalues <- c(1:nvalues[k])
    for (k1 in ndomvalues) {
      cv <- cbind(cv, cvx)
      cvDom <- c(cvDom, rep(name_dom[k], length(cvx)))
      cvDom2 <- c(cvDom2, rep(levels(as.factor(dom[, k]))[k1],
                              length(cvx)))
    }
  }
  
  calcola_cv <- function() {
    
    #
    # -------------------------------------------------------------
    # Populations in strata for different domains of interest
    # NTOTj
    # -------------------------------------------------------------
    NTOT <- colSums((m > 0) * N)
    # NTOT is a row-vector of dimension equal to the number of column of m (or s)
    # with the total number of population units in each domain
    NTOT_all=m
    NTOT_all=t(apply(NTOT_all, 1, function(x) x=t(as.numeric(x>0)*NTOT)))
    
    # NTOT_all is a matrix with the same dimension as m: 
    #1) NTOT_all has been inizialed as equal to m 
    #2) Then in each column j, the non-null value will be filled up with  
    #the j-the element of vector NTOT
    #
    # -------------------------------------------------------------
    # Sampled population in strata for different domains of interest
    # allocTOTj
    # -------------------------------------------------------------
    allocTOT=colSums((m > 0) * round(alloc))
    # allocTOT is a row-vector of dimension equal to the number of column of m (or s)
    # with the total number of sampled units in each domain
    #
    
    # ------------------------------------------------------------
    # Computation of the CVs
    # ------------------------------------------------------------
    #
    # ------------------------------------------------------------
    # Total
    # ------------------------------------------------------------ 
    totm <- rowSums(t(m * N))
    # row-vector of dimension equal to the number of column of m (or s)
    # with the total (mean*number of pop units) in each domain 
    totm_all=m
    totm_all=t(apply(totm_all, 1, function(x) x=t(as.numeric(x>0)*t(totm))))
    # totm_all is a matrix with the same dimension as m: 
    #1) totm_all has been inizialed as equal to m 
    #2) Then in each column j, the non-null value will be filled up with  
    #the j-the element of vector totm
    
    # ------------------------------------------------------------
    # Variance
    # ------------------------------------------------------------ 
    var=s^2*(N-1)
    #correction factor
    diff=N*((1/N)*(m*N) - totm_all/NTOT_all)^2
    diff[is.na(diff)]=0
    var2=sqrt((colSums(var)+colSums(diff))/NTOT) 
    varfin= var2^2*((1 - allocTOT/NTOT)*1/(allocTOT))

    # ------------------------------------------------------------
    # CV
    # ------------------------------------------------------------ 
    
    CVfin <- sqrt(varfin/(totm/NTOT)^2)
    #sd/mean
    return(CVfin)
  }
  
  # -----------------------------------------------------------

  
  #create the dataframe to be returned as output
  out_cv=data.frame(matrix(nrow=nvar*sum(nvalues), ncol=4))
  CVfin <- calcola_cv()
  colnames(out_cv)=c("TYPE", "DOMAIN/VAR", "PLANNED_CV ", 
                     "ACTUAL_CV")
  out_cv[,1]=as.vector(cvDom)
  out_cv[,2]=paste(as.vector(cvDom2), 
                          "/V", c(1:ncol(med)), sep="")
  out_cv[,3]=as.vector(cv)
  out_cv[,4]= as.vector(CVfin)

  return(out_cv)
  
  # End body
}

checkData <- function(strata = NULL, errors = NULL) {
  # controls on strata dataframe
  if (!is.null(strata)) {
    if (sum(grepl("N", colnames(strata))) < 1) 
      stop("In strata dataframe the indication of population (N) is missing")
    
    if (sum(grepl("STRAT", toupper(colnames(strata)), fixed = TRUE)) < 
        1) 
      stop("In strata dataframe the indication of stratum (STRATUM) is missing")
    if (sum(grepl("DOM1", toupper(colnames(strata)), fixed = TRUE)) < 
        1) 
      stop("In strata dataframe the indication of at least Ine domain (DOM1) is missing")
    if(sum(grepl('CENS',toupper(colnames(strata)),fixed=TRUE))
    < 1) stop('In strata dataframe the indication of strata
    to be sampled or censused (CENS) is missing') 
    if(sum(grepl('COST',toupper(colnames(strata)),fixed=TRUE))
    < 1) stop('In strata dataframe the indication of
    interviewing cost in strata (COST) is missing')
    if (sum(grepl("M1", toupper(colnames(strata)), fixed = TRUE)) < 1) 
      stop("In strata dataframe the indication of at least one mean (M1) is missing")
    if (sum(grepl("S1", toupper(colnames(strata)), fixed = TRUE)) < 1) 
      stop("In strata dataframe the indication of at least one standard deviation (S1) is missing")
    if(sum(grepl('^M+[0123456789]',toupper(colnames(strata)),perl=TRUE))!=  sum(grepl('^S+[0123456789]',toupper(colnames(strata)),perl=TRUE)))
    stop('In strata dataframe the number of means (Mx)
    differs from the number of standard deviations (Sx)')
  }
  # controls on errors dataframe
  if (!is.null(errors)) {
    if (sum(grepl("DOM", toupper(colnames(errors)), fixed = TRUE)) < 
        1) 
      stop("In errors dataframe the indication of domain (DOM) is missing")
    if (sum(grepl("CV", toupper(colnames(errors)), fixed = TRUE)) < 
        1) 
      stop("In errors dataframe the indication of at least one constraint (CV) is missing")
  }
  # crossed controls between errors and strata
  if (!is.null(errors) && !is.null(strata)) {
    if
    (sum(grepl('^S+[0123456789]',toupper(colnames(strata)),perl=TRUE)) != sum(grepl('CV',toupper(colnames(errors)),fixed=TRUE)))
    stop('In strata dataframe the number of means and std deviations differs from the number of coefficient of
    variations in errors dataframe')
    if (sum(grepl("DOM", toupper(colnames(strata)), fixed = TRUE)) !=   nrow(errors)) 
      stop("The different domains (DOMx) in strata dataframe are not represented in errors dataframe")
  }
}

# check n<minnumstr
check_n <- function(allocation=NULL, strata=NULL, minnumstrat) {
  if (!is.null(allocation)) {
    if(any(allocation<minnumstrat & strata$N>minnumstrat)){
      stop(paste("The allocated sample size is lower than the required minimum number of sampled unit per stratum, even if the stratum size allows to select a higher size. This happens in the following strata: ", paste(strata$STRATUM[which(n<minnumstrat & strata$N>minnumstrat)], "\n", sep="")))
    }
    if(any(strata$CENS==1 & strata$N!=allocation)){
      stop(paste("The stratum ", strata$STRATUM[which(strata$CENS==1 & strata$N!=n)], " must be censused. \n", sep=""))
    }
    
  }
}
  
