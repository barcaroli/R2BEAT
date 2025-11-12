prepareInputToAllocation2 <- function (frame,
                                       RGdes,
                                       RGcal,
                                       id_PSU,
                                       id_SSU,
                                       stratum,
                                       target,
                                       deff_level,
                                       domain,
                                       delta,
                                       minimum){
    if (requireNamespace("ReGenesees", quietly = TRUE)) {
      svystat <- ReGenesees::svystat
    }else{
      stop("The package ReGenesees is needed. \nInstall it by executing the following: \ndevtools::install_github('DiegoZardetto/ReGenesees'")
    }
    for (i in (1:length(target))) {
      tvi <- target[i]
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
    if(is.null(RGcal)){
      RGcal=RGdes
      RGcal=des.addvars(RGcal,weight.cal=weights(RGdes))
    }
    wgtcal <- colnames(RGcal$variables)[grep(".cal", colnames(RGcal$variables))]
    if (length(wgtcal) == 0){stop("No calibrated weights in RGcal calibrated object")}
    pos <- unlist(gregexpr(".cal", wgtcal, fixed = TRUE))
    wgtdes <- substr(wgtcal, 1, (pos - 1))
    if (!(wgtdes %in% colnames(RGdes))) {warning("No weights in RGdes design object")}
    options(warn = -1)
    options(scipen = 9999)
    id_vars <- c(id_PSU, id_SSU)
    sv1 <- NULL
    for (i in 1:(length(stratum))) {
      if (i < length(stratum)) 
        sv1 <- paste0(sv1, stratum[i], "+")
      if (i == length(stratum)) 
        sv1 <- paste0(sv1, stratum[i])
    }
    sv2 <- NULL
    for (i in 1:(length(stratum))) {
      if (i < length(stratum)) 
        sv2 <- paste0(sv2, stratum[i], ":")
      if (i == length(stratum)) 
        sv2 <- paste0(sv2, stratum[i])
    }
    tv <- NULL
    for (i in 1:(length(target))) {
      if (i < length(target)) 
        tv <- paste0(tv, target[i], "+")
      if (i == length(target)) 
        tv <- paste0(tv, target[i])
    }
    st <- paste("N <- aggregate(", wgtcal, " ~ ", sv1, ",RGcal$variables,FUN=sum)", 
                sep = "")
    eval(parse(text = st))
    st <- "N$STRATUM <- paste0("
    for (i in 1:(length(stratum))) {
      if (i < length(stratum)) 
        st <- paste0(st, "N$", stratum[i], ",")
      if (i == length(stratum)) 
        st <- paste0(st, "N$", stratum[i], ")")
    }
    eval(parse(text = st))
    colnames(N)[length(stratum) + 1] <- "N"
    N$N <- round(N$N, 0)
    RGcal$variables$one <- 1
    ssize <- NULL
    st <- paste0("ssize <- aggregate(one ~ ", sv1, ", RGcal$variables, FUN=sum)")
    eval(parse(text = st))
    M <- NULL
    eval(parse(text = paste0("M$", sv1, " <- ssize$", sv1, "")))
    M <- as.data.frame(M)
    for (i in (1:length(target))) {
      tvi <- target[i]
      st <- paste0("sw <- class(RGcal$variables$", tvi, ") == 'factor'")
      eval(parse(text = st))
      if (sw == FALSE) {
        Mi <- NULL
        st <- paste0("Mi <- svystat(RGcal,kind ='TM',estimator='Mean',y= ~", 
                     tvi, ",by= ~", sv2, ",forGVF=FALSE)")
        eval(parse(text = st))
        M <- cbind(M, Mi[, length(stratum) + 1])
      }
      if (sw == TRUE) {
        st <- paste0("Mi <- svystat(RGcal,kind ='TM',estimator='Mean',y= ~", 
                     tvi, ",by= ~", sv2, ",forGVF=FALSE)")
        eval(parse(text = st))
        M <- cbind(M, Mi[, length(stratum) + 2])
      }
      colnames(M)[1 + i] <- paste0("M", i)
    }
    S <- NULL
    eval(parse(text = paste0("S$", sv1, " <- ssize$", sv1, "")))
    S <- as.data.frame(S)
    for (i in (1:length(target))) {
      tvi <- target[i]
      st <- paste0("sw <- class(RGcal$variables$", tvi, ") == 'factor'")
      eval(parse(text = st))
      if (sw == FALSE) {
        d <- RGcal$variables[, c(stratum, wgtcal)]
        d$x2 <- RGcal$variables[, target[i]]^2
        st <- paste0("m <- aggregate(x2 * ", wgtcal, " ~ ", 
                     sv1, ", data=d, FUN=sum) / aggregate(", wgtcal, 
                     " ~ ", sv1, ", data=d, FUN=sum)")
        eval(parse(text = st))
        d$x <- RGcal$variables[, target[i]]
        st <- paste0("m2 <- aggregate((x * ", wgtcal, ") ~ ", 
                     sv1, ", data=d, FUN=sum) / aggregate(", wgtcal, 
                     " ~ ", sv1, ", data=d, FUN=sum)")
        eval(parse(text = st))
        st <- paste0("S$S", i, " <- sqrt(m[,length(stratum)+1]-m2[,length(stratum)+1]^2)")
        eval(parse(text = st))
      }
      if (sw == TRUE) {
        st <- paste0("S$S", i, " <- sqrt(M$M", i, " * (1-M$M", 
                     i, "))")
        eval(parse(text = st))
      }
    }
    strata <- NULL
    strata <- merge(N, ssize, by = c(sv1))
    strata <- merge(strata, M, by = c(sv1))
    strata <- merge(strata, S, by = c(sv1))
    strata$COST <- 1
    strata$CENS <- 0
    strata$DOM1 <- 1
    RGdes$variables$ones <- 1
    strata_domain <- aggregate(RGdes$variables$ones, RGdes$variables[, 
                                                                     c(paste0(stratum), paste0(domain))], sum)
    strata <- merge(strata, strata_domain[, c(paste0(stratum), 
                                              paste0(domain))])
    k <- 1
    for (i in (length(domain):1)) {
      k <- k + 1
      colnames(strata)[ncol(strata) - i + 1] <- paste0("DOM",k)
    }
    strata$one <- NULL
    st <- paste("dl <- aggregate(", wgtdes ,"~ ", stratum, "+", deff_level,", RGdes$variables,FUN=sum)",sep = "")
    eval(parse(text = st))       
    st <- paste("dl$", wgtdes, " <- NULL")
    eval(parse(text = st))       
    strata <- merge(strata,dl,by=sv1)
    dv <- NULL
    for (i in 1:(length(deff_level))) {
      if (i < length(deff_level)) 
        dv <- paste0(dv, deff_level[i], "+")
      if (i == length(deff_level)) 
        dv <- paste0(dv, deff_level[i])
    }
    deff <- NULL
    for (i in (1:length(target))) {
      tvi <- target[i]
      st <- paste0("sw <- class(RGcal$variables$", tvi, ") == 'factor'")
      eval(parse(text = st))
      if (sw == FALSE) {
        st <- paste0("deffi <- svystat(RGdes,kind ='TM',estimator='Mean',y= ~", 
                     tvi, ",by=~", dv, ",deff=TRUE,forGVF=TRUE)")
        eval(parse(text = st))
        for (j in 1:nrow(deffi)) {
          st <- paste0("deffi$", deff_level, "[j] <- substr(deffi$name[j],1,(gregexpr(pattern =':',deffi$name[j])[[1]][1])-1)")
          eval(parse(text = st))
        }
        deffi$label <- paste0("DEFF", i)
        deffi <- deffi[, c(deff_level, "DEFF", "label")]
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
          st <- paste0("deffi$", deff_level, "[j] <- substr(deffi$name[j],1,(gregexpr(pattern =':',deffi$name[j])[[1]][1])-1)")
          eval(parse(text = st))
        }
        deffi$label <- paste0("DEFF", i)
        deffi <- deffi[, c(deff_level, "DEFF", "label")]
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
    colnames(b_nar) <- c(deff_level, "b_nar")
    st <- "b_nar$STRATUM <- paste0("
    for (i in 1:(length(deff_level))) {
      if (i < length(deff_level)) 
        st <- paste0(st, "b_nar$", deff_level[i], ",")
      if (i == length(deff_level)) 
        st <- paste0(st, "b_nar$", deff_level[i], ")")
    }
    eval(parse(text = st))
    deff$progr <- c(1:nrow(deff))
    deff <- merge(deff, b_nar, by = deff_level, sort = FALSE)
    deff <- deff[order(deff$progr), ]
    deff$progr <- NULL
    for (i in 1:length(deff_level)) {
      st <- paste0("d=data.frame(", deff_level[i], "=unique(deff[,paste0(deff_level)]))")
      eval(parse(text = st))
    }
    for (i in 1:length(target)) {
      st <- paste0("def <- aggregate(DEFF~label+", dv, ",deff[deff$label=='DEFF",
                   i, "',],FUN=mean)[,2:3]")
      eval(parse(text = st))
      colnames(def)[2] <- paste("DEFF", i, sep = "")
      d <- merge(d, def, by = c(paste0(deff_level)))
    }
    b <- NULL
    st <- paste0("b <- aggregate(b_nar ~ ", dv, ",deff,FUN=mean)")
    eval(parse(text = st))
    effst <- NULL
    for (i in (1:length(target))) {
      tvi <- target[i]
      st <- paste0("sw <- class(RGcal$variables$", tvi, ") == 'factor'")
      eval(parse(text = st))
      if (sw == FALSE) {
        st <- paste0("effsti <- svystat(RGcal,kind ='TM',estimator='Mean',y= ~", 
                     tvi, ",by=~", dv, ",deff=TRUE,forGVF=TRUE)")
        eval(parse(text = st))
        for (j in 1:nrow(effsti)) {
          st <- paste0("effsti$", deff_level, "[j] <- substr(effsti$name[j],1,(gregexpr(pattern =':',effsti$name[j])[[1]][1])-1)")
          eval(parse(text = st))
        }
        effsti$label <- paste0("EFFST", i)
        effsti <- effsti[, c(deff_level, "DEFF", "label")]
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
          st <- paste0("effsti$", deff_level, "[j] <- substr(effsti$name[j],1,(gregexpr(pattern =':',effsti$name[j])[[1]][1])-1)")
          eval(parse(text = st))
        }
        effsti$label <- paste0("EFFST", i)
        effsti <- effsti[, c(deff_level, "DEFF", "label")]
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
    for (i in 1:length(deff_level)) {
      st <- paste0("e=data.frame(", deff_level[i], "=unique(effst[,paste0(deff_level)]))")
      eval(parse(text = st))
    }
    for (i in 1:length(target)) {
      st <- paste0("ef <- aggregate(EFFST~label+", dv, ",effst[effst$label=='EFFST", 
                   i, "',],FUN=mean)[,2:3]")
      eval(parse(text = st))
      colnames(ef)[2] <- paste("EFFST", i, sep = "")
      e <- merge(e, ef)
    }
    e <- merge(strata[,c("STRATUM", paste0(stratum),  paste0(deff_level))], 
               e,by=paste0(deff_level))
    d <- merge(d, b, by = paste0(deff_level))
    d <- merge(strata[, c("STRATUM", paste0(stratum), paste0(deff_level))],d, 
               by=paste0(deff_level))
    rho <- d
    for (i in (1:length(target))) {
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
    for (i in (1:length(target))) {
      st <- paste0("rho$DEFF", i, " <- NULL")
      eval(parse(text = st))
    }
    for (i in (1:length(stratum))) {
      st <- paste0("rho$", stratum[i], " <- NULL")
      eval(parse(text = st))
    }
  one <- list(strata = strata, deff = d, effst = e, rho = rho)
  frame$one <- 1
  psu <- eval(parse(text = paste0("aggregate(one~", id_PSU, 
                                  "+", stratum, ",data=frame,FUN=sum)")))
  mos_var <- "one"
  if (length(stratum) == 1) {
    st <- paste0("frame$stratum_new <- frame$", 
                 stratum[1], "")}
  if (length(stratum) > 1) {
    st <- NULL
    for (i in c(1:length(stratum))) {
      if (i == 1) 
        st1 <- paste0("frame$stratum_new <- cbind(frame[,'", 
                      stratum[i], "']")}
      if (i > 1 & i != length(stratum)){ 
        st1 <- paste0(",frame[,'", stratum[i], 
                      "']")}
      if (i == length(stratum)){ 
        st1 <- paste0(",frame[,'", stratum[i], 
                      "'])")}
      st <- paste0(st, st1)
    }
    eval(parse(text = st))
psu_file <- NULL
psu_file$PSU_ID <- psu[, paste0(id_PSU)]
psu_file$STRATUM <- psu[, paste0(stratum)]
psu_file$PSU_MOS <- psu[, paste0(mos_var)]
psu_file <- as.data.frame(psu_file)
des_file <- aggregate(PSU_MOS ~ STRATUM, data = psu_file, 
                      sum)
des_file$DELTA <- delta
des_file$MINIMUM <- minimum
colnames(des_file)[2] <- "STRAT_MOS"
two <- list(psu_file = psu_file, des_file = des_file)
out <- list(file_strata = one[[1]], deff = one[[2]], effst = one[[3]], 
              rho = one[[4]], psu_file = two[[1]], des_file = two[[2]])
return(out)
}
