input_to_beat.2st_1 <- function (RGdes, 
                                 RGcal, 
                                 weight_var_des, 
                                 weight_var_cal,
                                 id_PSU, id_SSU, 
                                 strata_vars, 
                                 target_vars, 
                                 deff_vars, 
                                 domain_vars) 
{
  if (requireNamespace("ReGenesees", quietly = TRUE)) {
    svystat <- ReGenesees::svystat
  }
  else {
    stop("The package ReGenesees is needed. \nInstall it by executing the following: \ndevtools::install_github('DiegoZardetto/ReGenesees'")
  }
  for (i in (1:length(target_vars))) {
    tvi <- target_vars[i]
    sw <- FALSE
    st <- paste0("sw <- class(RGdes$variables$", tvi, ") == 'factor'")
    eval(parse(text = st))
    if (sw == TRUE) {
      tab <- NULL
      st <- paste0("tab <- table(RGdes$variables$", tvi, 
                   ")")
      eval(parse(text = st))
      if (dim(tab) > 2) 
        stop("Factor '", tvi, "' is not binary")
      if (dim(tab) == 2) {
        vals <- NULL
        st <- paste0("vals <- c(levels(RGdes$variables$", 
                     tvi, "))")
        eval(parse(text = st))
        if (vals[1] != "0" | vals[2] != "1") 
          stop("Values of factor '", tvi, "' not equal to 0/1")
      }
    }
  }
  options(warn = -1)
  options(scipen = 9999)
  id_vars <- c(id_PSU, id_SSU)
  st <- paste0("RGcal$variables$",weight_var," <- RGcal$variables$",weight_var_cal)
  eval(parse(text = st))
  sv1 <- NULL
  for (i in 1:(length(strata_vars))) {
    if (i < length(strata_vars)) 
      sv1 <- paste0(sv1, strata_vars[i], "+")
    if (i == length(strata_vars)) 
      sv1 <- paste0(sv1, strata_vars[i])
  }
  sv2 <- NULL
  for (i in 1:(length(strata_vars))) {
    if (i < length(strata_vars)) 
      sv2 <- paste0(sv2, strata_vars[i], ":")
    if (i == length(strata_vars)) 
      sv2 <- paste0(sv2, strata_vars[i])
  }
  tv <- NULL
  for (i in 1:(length(target_vars))) {
    if (i < length(target_vars)) 
      tv <- paste0(tv, target_vars[i], "+")
    if (i == length(target_vars)) 
      tv <- paste0(tv, target_vars[i])
  }
  st <- paste("N <- aggregate(", weight_var, " ~ ", sv1, ",RGcal$variables,FUN=sum)", 
              sep = "")
  eval(parse(text = st))
  st <- "N$STRATUM <- paste0("
  for (i in 1:(length(strata_vars))) {
    if (i < length(strata_vars)) 
      st <- paste0(st, "N$", strata_vars[i], ",")
    if (i == length(strata_vars)) 
      st <- paste0(st, "N$", strata_vars[i], ")")
  }
  eval(parse(text = st))
  colnames(N)[length(strata_vars) + 1] <- "N"
  N$N <- round(N$N, 0)
  RGcal$variables$one <- 1
  ssize <- NULL
  st <- paste0("ssize <- aggregate(one ~ ", sv1, ", RGcal$variables, FUN=sum)")
  eval(parse(text = st))
  st <- "ssize$STRATUM <- paste0("
  for (i in 1:(length(strata_vars))) {
    if (i < length(strata_vars)) 
      st <- paste0(st, "ssize$", strata_vars[i], ",")
    if (i == length(strata_vars)) 
      st <- paste0(st, "ssize$", strata_vars[i], ")")
  }
  eval(parse(text = st))
  for (i in 1:(length(strata_vars))) {
    st <- paste0("ssize$", strata_vars[i], "<-NULL")
    eval(parse(text = st))
  }
  M <- NULL
  M$STRATUM <- ssize$STRATUM
  M <- as.data.frame(M)
  for (i in (1:length(target_vars))) {
    tvi <- target_vars[i]
    st <- paste0("sw <- class(RGcal$variables$", tvi, ") == 'factor'")
    eval(parse(text = st))
    if (sw == FALSE) {
      Mi <- NULL
      st <- paste0("Mi <- svystat(RGcal,kind ='TM',estimator='Mean',y= ~", 
                   tvi, ",by= ~", sv2, ",forGVF=FALSE)")
      eval(parse(text = st))
      M <- cbind(M, Mi[, length(strata_vars) + 1])
    }
    if (sw == TRUE) {
      st <- paste0("Mi <- svystat(RGcal,kind ='TM',estimator='Mean',y= ~", 
                   tvi, ",by= ~", sv2, ",forGVF=FALSE)")
      eval(parse(text = st))
      M <- cbind(M, Mi[, length(strata_vars) + 2])
    }
    colnames(M)[1 + i] <- paste0("M", i)
  }
  S <- NULL
  S$STRATUM <- ssize$STRATUM
  S <- as.data.frame(S)
  for (i in (1:length(target_vars))) {
    tvi <- target_vars[i]
    st <- paste0("sw <- class(RGcal$variables$", tvi, ") == 'factor'")
    eval(parse(text = st))
    if (sw == FALSE) {
      d <- RGcal$variables[, c(strata_vars, weight_var)]
      d$x2 <- RGcal$variables[, target_vars[i]]^2
      st <- paste0("m <- aggregate(x2 * ", weight_var, 
                   " ~ ", sv1, ", data=d, FUN=sum) / aggregate(", 
                   weight_var, " ~ ", sv1, ", data=d, FUN=sum)")
      eval(parse(text = st))
      d$x <- RGcal$variables[, target_vars[i]]
      st <- paste0("m2 <- aggregate((x * ", weight_var, 
                   ") ~ ", sv1, ", data=d, FUN=sum) / aggregate(", 
                   weight_var, " ~ ", sv1, ", data=d, FUN=sum)")
      eval(parse(text = st))
      st <- paste0("S$S", i, " <- sqrt(m[,length(strata_vars)+1]-m2[,length(strata_vars)+1]^2)")
      eval(parse(text = st))
    }
    if (sw == TRUE) {
      st <- paste0("S$S", i, " <- sqrt(M$M", i, " * (1-M$M", 
                   i, "))")
      eval(parse(text = st))
    }
  }
  strata <- NULL
  strata <- merge(N, ssize, by = c("STRATUM"))
  strata <- merge(strata, M, by = c("STRATUM"))
  strata <- merge(strata, S, by = c("STRATUM"))
  strata$COST <- 1
  strata$CENS <- 0
  strata$DOM1 <- 1
  RGdes$variables$ones <- 1
  strata_domain <- aggregate(RGdes$variables$ones, RGdes$variables[, 
                                                                   c(paste0(strata_vars), paste0(domain_vars))], sum)
  strata <- merge(strata, strata_domain[, c(paste0(strata_vars), 
                                            paste0(domain_vars))])
  k <- 1
  for (i in (length(domain_vars):1)) {
    k <- k + 1
    colnames(strata)[ncol(strata) - i + 1] <- paste0("DOM", 
                                                     k)
  }
  strata$one <- NULL
  strata
  dv <- NULL
  for (i in 1:(length(deff_vars))) {
    if (i < length(deff_vars)) 
      dv <- paste0(dv, deff_vars[i], "+")
    if (i == length(deff_vars)) 
      dv <- paste0(dv, deff_vars[i])
  }
  deff <- NULL
  for (i in (1:length(target_vars))) {
    tvi <- target_vars[i]
    st <- paste0("sw <- class(RGcal$variables$", tvi, ") == 'factor'")
    eval(parse(text = st))
    if (sw == FALSE) {
      st <- paste0("deffi <- svystat(RGdes,kind ='TM',estimator='Mean',y= ~", 
                   tvi, ",by=~", dv, ",deff=TRUE,forGVF=TRUE)")
      eval(parse(text = st))
      for (j in 1:nrow(deffi)) {
        st <- paste0("deffi$", deff_vars, "[j] <- substr(deffi$name[j],1,(gregexpr(pattern =':',deffi$name[j])[[1]][1])-1)")
        eval(parse(text = st))
      }
      deffi$label <- paste0("DEFF", i)
      deffi <- deffi[, c(deff_vars, "DEFF", "label")]
      deffi$DEFF <- ifelse(is.nan(deffi$DEFF), 1, deffi$DEFF)
      deffi$DEFF <- round(deffi$DEFF, 6)
      deffi$DEFF <- ifelse(deffi$DEFF == 0, 1, deffi$DEFF)
      deff <- rbind(deff, deffi)
    }
    if (sw == TRUE) {
      st <- paste0("deffi <- svystat(RGdes,kind ='TM',estimator='Mean',y= ~", 
                   tvi, ",by=~", dv, ",deff=TRUE,forGVF=TRUE)")
      eval(parse(text = st))
      for (j in 1:nrow(deffi)) {
        st <- paste0("deffi$", deff_vars, "[j] <- substr(deffi$name[j],1,(gregexpr(pattern =':',deffi$name[j])[[1]][1])-1)")
        eval(parse(text = st))
      }
      deffi$label <- paste0("DEFF", i)
      deffi <- deffi[, c(deff_vars, "DEFF", "label")]
      deffi$DEFF <- ifelse(is.nan(deffi$DEFF), 1, deffi$DEFF)
      deffi$DEFF <- round(deffi$DEFF, 6)
      deffi$DEFF <- ifelse(deffi$DEFF == 0, 1, deffi$DEFF)
      deff <- rbind(deff, deffi[c((nrow(deffi)/2 + 1):nrow(deffi)), 
      ])
    }
  }
  ids <- NULL
  for (i in 1:(length(id_vars))) {
    if (i < length(id_vars)) 
      ids <- paste0(ids, id_vars[i], "+")
    if (i == length(id_vars)) 
      ids <- paste0(ids, id_vars[i])
  }
  RGdes$variables$one <- 1
  vars <- paste(dv, id_PSU, sep = "+")
  st <- paste0("b_nar <- aggregate(one ~ ", vars, ", RGdes$variables, FUN=sum)")
  eval(parse(text = st))
  st <- paste0("b_nar <- aggregate(one ~ ", dv, ", b_nar, FUN=mean)")
  eval(parse(text = st))
  colnames(b_nar) <- c(deff_vars, "b_nar")
  st <- "b_nar$STRATUM <- paste0("
  for (i in 1:(length(deff_vars))) {
    if (i < length(deff_vars)) 
      st <- paste0(st, "b_nar$", deff_vars[i], ",")
    if (i == length(deff_vars)) 
      st <- paste0(st, "b_nar$", deff_vars[i], ")")
  }
  eval(parse(text = st))
  deff$progr <- c(1:nrow(deff))
  deff <- merge(deff, b_nar, by = deff_vars, sort = FALSE)
  deff <- deff[order(deff$progr), ]
  deff$progr <- NULL
  for (i in 1:length(deff_vars)) {
    st <- paste0("d=data.frame(", deff_vars[i], "=unique(deff[,paste0(deff_vars)]))")
    eval(parse(text = st))
  }
  for (i in 1:length(target_vars)) {
    st <- paste0("def <- aggregate(DEFF~label+", dv, ",deff[deff$label=='DEFF", 
                 i, "',],FUN=mean)[,2:3]")
    eval(parse(text = st))
    colnames(def)[2] <- paste("DEFF", i, sep = "")
    d <- merge(d, def, by = c(paste0(deff_vars)))
  }
  b <- NULL
  st <- paste0("b <- aggregate(b_nar ~ ", dv, ",deff,FUN=mean)")
  eval(parse(text = st))
  effst <- NULL
  for (i in (1:length(target_vars))) {
    tvi <- target_vars[i]
    st <- paste0("sw <- class(RGcal$variables$", tvi, ") == 'factor'")
    eval(parse(text = st))
    if (sw == FALSE) {
      st <- paste0("effsti <- svystat(RGcal,kind ='TM',estimator='Mean',y= ~", 
                   tvi, ",by=~", dv, ",deff=TRUE,forGVF=TRUE)")
      eval(parse(text = st))
      for (j in 1:nrow(effsti)) {
        st <- paste0("effsti$", deff_vars, "[j] <- substr(effsti$name[j],1,(gregexpr(pattern =':',effsti$name[j])[[1]][1])-1)")
        eval(parse(text = st))
      }
      effsti$label <- paste0("EFFST", i)
      effsti <- effsti[, c(deff_vars, "DEFF", "label")]
      effsti$DEFF <- ifelse(is.nan(effsti$DEFF), 1, effsti$DEFF)
      effsti$DEFF <- round(effsti$DEFF, 6)
      effsti$DEFF <- ifelse(effsti$DEFF == 0, 1, effsti$DEFF)
      effst <- rbind(effst, effsti)
    }
    if (sw == TRUE) {
      st <- paste0("effsti <- svystat(RGdes,kind ='TM',estimator='Mean',y= ~", 
                   tvi, ",by=~", dv, ",deff=TRUE,forGVF=TRUE)")
      eval(parse(text = st))
      for (j in 1:nrow(effsti)) {
        st <- paste0("effsti$", deff_vars, "[j] <- substr(effsti$name[j],1,(gregexpr(pattern =':',effsti$name[j])[[1]][1])-1)")
        eval(parse(text = st))
      }
      effsti$label <- paste0("EFFST", i)
      effsti <- effsti[, c(deff_vars, "DEFF", "label")]
      effsti$DEFF <- ifelse(is.nan(effsti$DEFF), 1, effsti$DEFF)
      effsti$DEFF <- round(effsti$DEFF, 6)
      effsti$DEFF <- ifelse(effsti$DEFF == 0, 1, effsti$DEFF)
      effst <- rbind(effst, effsti[c((nrow(effsti)/2 + 
                                        1):nrow(effsti)), ])
    }
  }
  effst$EFFST <- effst$DEFF/deff$DEFF
  effst$EFFST <- ifelse(effst$EFFST == Inf | effst$EFFST == 
                          -Inf, 1, effst$EFFST)
  for (i in 1:length(deff_vars)) {
    st <- paste0("e=data.frame(", deff_vars[i], "=unique(effst[,paste0(deff_vars)]))")
    eval(parse(text = st))
  }
  for (i in 1:length(target_vars)) {
    st <- paste0("ef <- aggregate(EFFST~label+", dv, ",effst[effst$label=='EFFST", 
                 i, "',],FUN=mean)[,2:3]")
    eval(parse(text = st))
    colnames(ef)[2] <- paste("EFFST", i, sep = "")
    e <- merge(e, ef)
  }
  e <- merge(strata[, c("STRATUM", paste0(strata_vars))], 
             e)
  d <- merge(d, b, by = paste0(deff_vars))
  d <- merge(strata[, c("STRATUM", paste0(strata_vars))], 
             d)
  rho <- d
  for (i in (1:length(target_vars))) {
    st <- paste0("rho$RHO_AR", i, " <- 1")
    eval(parse(text = st))
    st <- paste0("rho$RHO_NAR", i, "<-(rho$DEFF", i, "-1)/(rho$b_nar-1)")
    eval(parse(text = st))
    st <- paste0("rho$RHO_NAR", i, "<-ifelse(rho$RHO_NAR", 
                 i, " == -Inf,0,rho$RHO_NAR", i, ")")
    eval(parse(text = st))
    st <- paste0("rho$RHO_NAR", i, "<-ifelse(rho$RHO_NAR", 
                 i, " == Inf,0,rho$RHO_NAR", i, ")")
    eval(parse(text = st))
  }
  rho$b_nar <- NULL
  for (i in (1:length(target_vars))) {
    st <- paste0("rho$DEFF", i, " <- NULL")
    eval(parse(text = st))
  }
  for (i in (1:length(strata_vars))) {
    st <- paste0("rho$", strata_vars[i], " <- NULL")
    eval(parse(text = st))
  }
  out <- list(strata = strata, deff = d, effst = e, rho = rho)
  return(out)
}
