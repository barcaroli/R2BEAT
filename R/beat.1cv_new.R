#' The function returns a dataframe with planned and actual coefficients of variation (CV) in a multivariate multi-domain allocation problem.
#'
#' @param stratif A dataframe with strata information.
#' @param alloc: A vector with an allocation into strata.
#' @param minnumstrat: Optional. A number indicating the lower bound for the number of allocated units in each stratum (default is 2).
#' @return A dataframe with actual CV for each combination of domain, domain category and target variable.
beat.1cv_2 <- function (stratif, alloc, minnumstrat = 2) 
{
  colnames(stratif) <- toupper(colnames(stratif))
  checkData <- function(strata = NULL) {
    if (!is.null(strata)) {
      if (sum(grepl("N", colnames(strata))) < 1) 
        stop("In strata dataframe the indication of population (N) is missing")
      if (sum(grepl("STRAT", toupper(colnames(strata)), 
                    fixed = TRUE)) < 1) 
        stop("In strata dataframe the indication of stratum (STRATUM) is missing")
      if (sum(grepl("DOM1", toupper(colnames(strata)), 
                    fixed = TRUE)) < 1) 
        stop("In strata dataframe the indication of at least Ine domain (DOM1) is missing")
      if (sum(grepl("CENS", toupper(colnames(strata)), 
                    fixed = TRUE)) < 1) 
        stop("In strata dataframe the indication of strata\n    to be sampled or censused (CENS) is missing")
      if (sum(grepl("COST", toupper(colnames(strata)), 
                    fixed = TRUE)) < 1) 
        stop("In strata dataframe the indication of\n    interviewing cost in strata (COST) is missing")
      if (sum(grepl("M1", toupper(colnames(strata)), fixed = TRUE)) < 
          1) 
        stop("In strata dataframe the indication of at least one mean (M1) is missing")
      if (sum(grepl("S1", toupper(colnames(strata)), fixed = TRUE)) < 
          1) 
        stop("In strata dataframe the indication of at least one standard deviation (S1) is missing")
      if (sum(grepl("^M+[0123456789]", toupper(colnames(strata)), 
                    perl = TRUE)) != sum(grepl("^S+[0123456789]", 
                                               toupper(colnames(strata)), perl = TRUE))) 
        stop("In strata dataframe the number of means (Mx)\n    differs from the number of standard deviations (Sx)")
    }
  }
  check_n <- function(allocation = NULL, strata = NULL, minnumstrat) {
    if (!is.null(allocation)) {
      if (any(allocation < minnumstrat & strata$N > minnumstrat)) {
        error_msg = paste(c(paste("The allocated sample size is lower than the required minimum number of sampled unit per stratum (minnumstrat=", 
                                  minnumstrat, "), even if the stratum size allows to select a higher size. This happens in the following strata:", 
                                  sep = ""), paste(strata$STRATUM[which(allocation < 
                                                                          minnumstrat & strata$N > minnumstrat)], collapse = ";")), 
                          collapse = " ")
        stop(error_msg)
      }
      if (any(strata$CENS == 1 & strata$N != allocation)) {
        stop(paste(c(paste("The following stratum(a) must be censused:"), 
                     paste(strata$STRATUM[which(strata$CENS == 
                                                  1 & strata$N != allocation)], collapse = "; ")), 
                   collapse = " "))
      }
    }
  }
  checkData(strata = stratif)
  check_n(allocation = alloc, strata = stratif, minnumstrat)
  nstrat <- nrow(stratif)
  nvar <- sum(grepl("^M[0-9]+$", names(stratif)))
  ndom <-sum(grepl("^DOM[0-9]+$", names(stratif)))
  med <- as.matrix(stratif[, paste("M", 1:nvar, sep = ""), drop = FALSE])
  esse <- as.matrix(stratif[, paste("S", 1:nvar, sep = ""), drop = FALSE])
  if (ncol(med) != ncol(esse)) 
    stop(print("Error: Number of M variables don't match the number of S variables"))
  if (ncol(med) != nvar) 
    stop(print("Error: Number of variables don't match the number of planned CV"))
  N <- as.vector(stratif$N)
  cens <- as.vector(stratif$CENS)
  cens[N < minnumstrat] <- 1
  cost <- as.vector(stratif$COST)
  if (is.null(cost)) 
    cost <- rep(1, nstrat)
  if (is.null(cens)) 
    cens <- rep(0, nstrat)
  nocens <- 1 - cens
  if (sum(cens) == length(cens)) {
    warning(print("Warning: Variable CENS always equal 1"))
  }
  name_dom <- sapply(1:ndom, function(i) paste("DOM", i, sep = ""))
  dom <- as.data.frame(stratif[, paste("DOM", 1:ndom, sep = ""), drop = FALSE])
  val <- NULL
  nvalues <- sapply(name_dom, function(vari) {
    val <- c(val, nlevels(as.factor(dom[, vari])))
  })
  crea_disj <- function(data, vars) {
    out <- NULL
    sapply(vars, function(vari) {
      col <- as.factor(data[, vari])
      out <<- cbind(out, outer(col, levels(col), function(y, 
                                                          x) {
        ifelse(y == x, 1, 0)
      }))
    })
    out
  }
  disj <- crea_disj(stratif, name_dom)
  nc <- ncol(disj)
  m <- NULL
  s <- NULL
  for (i in 1:nc) {
    m <- cbind(m, disj[, i] * med)
    s <- cbind(s, disj[, i] * esse)
  }
  cvDom <- NULL
  cvDom2 <- NULL
  cvDom3 <- NULL
  for (k in 1:ndom) {
    ndomvalues <- c(1:nvalues[k])
    for (k1 in ndomvalues) {
      cvDom <- c(cvDom, rep(name_dom[k], nvar))
      cvDom2 <- c(cvDom2, rep((1:nlevels(as.factor(dom[, k])))[k1], 
                              nvar))
      cvDom3 <- c(cvDom3, rep(levels(as.factor(dom[, k]))[k1], 
                              nvar))
    }
  }
  calcola_cv <- function() {
    NTOT <- c(rep(0, nvar))
    CVfin <- c(rep(0, nvar))
    NTOT <- colSums((m > 0) * N)
    varfin <- rowSums(t((s * N)^2 * (1 - round(alloc)/N)/round(alloc))/NTOT^2)
    totm <- rowSums(t(m * N))
    CVfin <- round(sqrt(varfin/(totm/NTOT)^2), digits = 6)
    return(CVfin)
  }
  out_cv = data.frame(matrix(nrow = nvar * sum(nvalues), ncol = 5))
  CVfin <- calcola_cv()
  colnames(out_cv) = c("Type", 
                       "Dom", "Dom_label", 
                       "Var", "ExpectedCV")
  out_cv[, 1] = as.vector(cvDom)
  out_cv[, 2] = as.vector(cvDom2)
  out_cv[, 3] = as.vector(cvDom3)
  out_cv[, 4] = as.vector(paste("V", c(1:ncol(med)), sep = ""))
  out_cv[, 5] = as.vector(CVfin)
  
  # expectedCV output
  MAXexp_cv=aggregate(ExpectedCV ~ Type+Var, out_cv, max)
  MAXexp_cv=reshape(MAXexp_cv, 
                 direction = "wide",
                 idvar = "Type",
                 timevar = "Var")
  colnames(MAXexp_cv)[-1] = sub("ExpectedCV.V", "", colnames(MAXexp_cv)[-1])
  MAXexp_cv = MAXexp_cv[,c(1, order(as.numeric(colnames(MAXexp_cv)[-1]))+1)]
  colnames(MAXexp_cv) <- c("DOM", paste("V", 1:nvar, sep=""))
  expected_cv=list(expectedCV_spec=out_cv,expectedCV= MAXexp_cv)
  return(expected_cv)
}
