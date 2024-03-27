beat.1st <- function (stratif, errors, minnumstrat = 2, maxiter = 200, maxiter1 = 25, 
          epsilon = 10^(-11)) 
{
  test_stages <- try(test_stages, TRUE)
  if (test_stages == 2) {
    minnumstrat <- param_alloc$p_minnumstrat
    maxiter <- param_alloc$p_maxiter
    maxiter1 <- param_alloc$p_maxiter1
    epsilon <- param_alloc$p_epsilon
  }
  else {
    param_alloc <- as.data.frame(t(c(epsilon, minnumstrat, 
                                     maxiter, maxiter1, 1)))
    names(param_alloc) = c("p_epsilon", "p_minnumstrat", 
                           "p_maxiter", "p_maxiter1", "num_stages")
  }
  colnames(stratif) <- toupper(colnames(stratif))
  colnames(errors) <- toupper(colnames(errors))
  iter1 <- 0
  val = NULL
  m <- NULL
  s <- NULL
  cv = NULL
  nstrat = nrow(stratif)
  nvar = ncol(errors) - 1
  ndom = nrow(errors)
  param_alloc$nvar <- nvar
  param_alloc$ndom <- ndom
  param_alloc$nstrat <- nstrat
  param_alloc <<- param_alloc
  varloop <- c(1:nvar)
  strloop <- c(1:nstrat)
  domloop <- c(1:ndom)
  med <- as.matrix(stratif[, names(stratif) %in% sapply(1:nvar, 
                                                        function(i) paste("M", i, sep = ""))])
  esse <- as.matrix(stratif[, names(stratif) %in% sapply(1:nvar, 
                                                         function(i) paste("S", i, sep = ""))])
  nom_dom <- sapply(1:ndom, function(i) paste("DOM", i, sep = ""))
  dom <- as.vector(stratif[, names(stratif) %in% nom_dom])
  N <- as.vector(stratif$N)
  cens <- as.vector(stratif$CENS)
  cost <- as.vector(stratif$COST)
  nocens = 1 - cens
  if (ndom == 1) 
    (nvalues <- nlevels((as.factor(stratif$DOM1))))
  if (ndom > 1) {
    nvalues <- NULL
    for (i in c(1:ndom)) {
      eval(parse(text = paste0("nvalues <- c(nvalues,nlevels(as.factor(stratif$DOM", 
                               i, ")))")))
    }
  }
  crea_disj = function(data, vars) {
    out = NULL
    sapply(vars, function(vari) {
      col = as.factor(data[, vari])
      out <<- cbind(out, outer(col, levels(col), function(y, 
                                                          x) ifelse(y == x, 1, 0)))
    })
    out
  }
  disj <- crea_disj(stratif, nom_dom)
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
  sampleSize_nocens <- sum(n[which(stratif$CENS==0)])
  popSize_nocens <- sum(N[which(stratif$CENS==0)])
  #originale
  #uguale = rep(sampleSize/num_strati, num_strati)
  uguale <- ifelse(cens==0, sampleSize_nocens/sum(stratif$CENS==0), N)
  proporzionale = ifelse(cens==0, sampleSize_nocens * N/popSize_nocens, N)
  Bethel_sample <- cbind(stratif, n)
  colnames(Bethel_sample)[length(colnames(Bethel_sample))] <- "n"
  nomi <- c("STRATUM", "ALLOC", "PROP", "EQUAL")
  df = NULL
  #originale
  #df <- cbind(df, as.character(stratif$STRATUM), n, sampleSize * 
   #             N/popSize, uguale)
  df <- cbind(df, as.character(stratif$STRATUM), n, proporzionale, uguale)
  tot <- apply(matrix(as.numeric(df[, 2:4]), ncol = 3), 2, 
               sum)
  df <- rbind(df, c("Total", tot))
  colnames(df) <- nomi
  #originale: 
  # for (j in varloop) {
  #   NTOT[j] <- sum((m[, j] > 0) * N)
  # }
  # varfin = rowSums(t((s * N)^2 * (1 - round(n)/N)/round(n))/NTOT^2)
  # totm = rowSums(t(m * N))
  # CVfin <- round(sqrt(varfin/(totm/NTOT)^2), digits = 4)
  #nuovo: 
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
    allocTOT=colSums((m > 0) * round(n))
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
  
  CVfin <- calcola_cv()
  
  
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
  domcard <- c(rep(0, ndom))
  for (k in (1:ndom)) {
    statement <- paste(" domcard[", k, "] <- length(levels(as.factor(stratif$DOM", 
                       k, ")))", sep = "")
    eval(parse(text = statement))
  }
  j <- 0
  outcv <- NULL
  outcv = cbind("Type", "dom", "Var", "Planned CV", "Actual CV", 
                "Sensitivity 10%")
  colnames(outcv) <- c("Type", "Dom", "Var", "Planned CV", 
                       "Actual CV", "Sensitivity 10%")
  for (k in (1:ndom)) {
    valloop <- c(1:(domcard[k]))
    for (k1 in valloop) {
      for (k2 in 1:nvar_orig) {
        j <- j + 1
        outcv <- rbind(outcv, c(nom_dom[k], k1, paste("V", 
                                                      k2, sep = ""), cv[j], CVfin[j], ceiling(sens[j])))
      }
    }
  }
  output_beth <- list(n = n, file_strata = Bethel_sample, 
                      alloc = as.data.frame(df), sensitivity = as.data.frame(outcv)[-1, 
                      ])
  output_beth$alloc[, 2] <- as.numeric(as.character(output_beth$alloc[, 
                                                                      2]))
  output_beth$alloc[, 3] <- as.numeric(as.character(output_beth$alloc[, 
                                                                      3]))
  output_beth$alloc[, 4] <- as.numeric(as.character(output_beth$alloc[, 
                                                                      4]))
  output_beth$sensitivity[, 6] <- as.numeric(as.character(output_beth$sensitivity[, 
                                                                                  6]))
  output_beth$sensitivity[, 5] <- as.numeric(as.character(output_beth$sensitivity[, 
                                                                                  5]))
  output_beth$sensitivity[, 4] <- as.numeric(as.character(output_beth$sensitivity[, 
                                                                                  4]))
  output_beth <<- output_beth
  return(output_beth)
}
