#-----------------------------------
# Function to prepare inputs to 
# allocation and PSU selection steps
# (scenario 1)
#-----------------------------------

prepareInputToAllocation1 <- function( 
  samp_frame,
  id_PSU,
  id_SSU,
  strata_var,
  target_vars,
  deff_var,
  domain_var,
  minimum,
  delta,
  f,
  deff_sugg) 
{
  if (is.null(samp_frame$one)) samp_frame$one <- 1

  # strata
  frame <- SamplingStrata::buildFrameDF(df=samp_frame,
                        id=id_SSU,
                        X=strata_var,
                        Y=target_vars,
                        domainvalue=domain_var)
  nvarY <- length(grep("Y",colnames(frame)))
  strata <- SamplingStrata::buildStrataDF(frame,progress=FALSE)
  strata$DOM2 <- strata$DOM1
  strata$DOM1 <- 1
  strata$STRATUM <- as.factor(strata$STRATO)
  strata$STRATO <- strata$X1 <- NULL
  strata <- strata[order(as.numeric(as.character(strata$STRATUM))),]

  # deff
  deff <- strata$STRATUM
  deff <- as.data.frame(deff)
  colnames(deff) <- "STRATUM"
  for (i in c(1:nvarY)) {
    st <- paste0("deff$DEFF",i," <- ",deff_sugg)
    eval(parse(text=st))
  }
  st <- paste0("b_nar <- aggregate(one ~ ",strata_var," + ",id_PSU,", data=samp_frame, FUN=sum)")
  eval(parse(text=st))
  st <- paste0("b_nar <- aggregate(one ~ ",strata_var,", b_nar, FUN=mean)")
  eval(parse(text=st))
  b_nar$one <- b_nar$one * f
  deff <- merge(deff,b_nar,by.x="STRATUM",by.y="stratum")
  colnames(deff)[ncol(deff)] <- "b_nar"
  deff <- deff[order(as.numeric(as.character(deff$STRATUM))),]

  # effst
  effst <- strata$STRATUM
  effst <- as.data.frame(effst)
  colnames(effst)[1] <- "STRATUM"
  for (i in c(1:nvarY)) {
    st <- paste0("effst$EFFST",i," <- 1")
    eval(parse(text=st))
  }
  effst <- effst[order(as.numeric(as.character(effst$STRATUM))),]

  # rho
  rho <- deff
  for (i in c(1:nvarY)) {
    st <- paste0("rho$RHO_AR",i," <- 1")
    eval(parse(text=st))
    st <- paste0("rho$RHO_NAR",i," <- (rho$DEFF", i, "-1)/(rho$b_nar-1)")
    eval(parse(text=st))
  }
  rho$b_nar <- NULL
  rho[,grep("DEFF",colnames(rho))] <- NULL
  rho <- rho[order(as.numeric(as.character(rho$STRATUM))),]

  # psu_file
  strat_mun <- samp_frame[,c(strata_var,id_PSU)]
  strat_mun <- strat_mun[!duplicated(strat_mun),]
  st <- paste0("mun <- aggregate(samp_frame$one,by=list(samp_frame$",id_PSU,"),sum)")
  eval(parse(text=st))
  colnames(mun) <- c(id_PSU,"N")
  mun <- merge(mun,strat_mun)
  psu_file <- mun[,c(id_PSU,strata_var,"N")]
  colnames(psu_file) <- c("PSU_ID","STRATUM","PSU_MOS")

  # des_file
  des_file <- aggregate(psu_file$PSU_MOS,by=list(psu_file$STRATUM),sum)
  colnames(des_file) <- c("STRATUM","STRAT_MOS")
  des_file$DELTA <- delta
  des_file$MINIMUM <- minimum

  out <- list(strata=strata,
              deff=deff,
              effst=effst,
              rho=rho,
              psu_file=psu_file,
              des_file=des_file)
  return(out)
}