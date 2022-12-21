prepareInputToAllocation1 <- function (samp_frame, id_PSU, id_SSU, strata_var, target_vars, 
          deff_var, domain_var, minimum, delta, deff_sugg) 
{
  f = 0.05
  if (is.null(samp_frame$one)) 
    samp_frame$one <- 1
  cat("\nCalculating strata...")
  frame <- buildFrameDF(df = samp_frame, id = id_SSU, X = strata_var, 
                        Y = target_vars, domainvalue = domain_var)
  nvarY <- length(grep("Y", colnames(frame)))
  strata <- buildStrataDF(frame, progress = FALSE)
  strata$DOM2 <- strata$DOM1
  strata$DOM1 <- 1
  strata$STRATUM <- as.factor(strata$STRATO)
  strata$STRATO <- strata$X1 <- NULL
  strata <- strata[order(as.numeric(as.character(strata$STRATUM))), 
  ]
  deff <- strata$STRATUM
  deff <- as.data.frame(deff)
  colnames(deff) <- "STRATUM"
  for (i in c(1:nvarY)) {
    st <- paste0("deff$DEFF", i, " <- ", deff_sugg)
    eval(parse(text = st))
  }
  st <- paste0("b_nar <- aggregate(one ~ ", strata_var, " + ", 
               id_PSU, ", data=samp_frame, FUN=sum)")
  eval(parse(text = st))
  st <- paste0("b_nar <- aggregate(one ~ ", strata_var, ", b_nar, FUN=mean)")
  eval(parse(text = st))
  if (ncol(b_nar) > 1) 
    colnames(b_nar)[1] <- "STRATUM"
  if (ncol(b_nar) == 1) 
    b_nar$STRATUM <- 1
  b_nar$one <- b_nar$one * f
  deff <- merge(deff, b_nar, by = "STRATUM")
  colnames(deff)[ncol(deff)] <- "b_nar"
  deff <- deff[order(as.numeric(as.character(deff$STRATUM))), 
  ]
  effst <- strata$STRATUM
  effst <- as.data.frame(effst)
  colnames(effst)[1] <- "STRATUM"
  for (i in c(1:nvarY)) {
    st <- paste0("effst$EFFST", i, " <- 1")
    eval(parse(text = st))
  }
  effst <- effst[order(as.numeric(as.character(effst$STRATUM))), 
  ]
  cat("\nCalculating rho in strata...")
  rho <- NULL
  rho$STRATUM <- strata$STRATUM
  rho <- as.data.frame(rho)
  for (i in c(1:nvarY)) {
    eval(parse(text = paste0("rho$RHO_AR", i, " <- 1")))
    eval(parse(text = paste0("rho$RHO_NAR", i, " <- NA")))
  }
  L <- NULL
  k <- 0
  for (s in c(rho$STRATUM)) {
    cat("\nStratum ", s)
    k <- k + 1
    for (i in c(1:nvarY)) {
      eval(parse(text = paste0("mu <- mean(samp_frame$", 
                               target_vars[i], "[samp_frame$", strata_var, 
                               " == ", s, "])")))
      eval(parse(text = paste0("D2y <- sum((samp_frame$", 
                               target_vars[i], "[samp_frame$", strata_var, 
                               " == ", s, "]-mu)^2)")))
      eval(parse(text = paste0("L <- unique(samp_frame$", 
                               id_PSU, "[samp_frame$", strata_var, "==s])")))
      D2w <- 0
      for (l in c(L)) {
        eval(parse(text = paste0("mu_L <- mean(samp_frame$", 
                                 target_vars[i], "[samp_frame$", id_PSU, " == ", 
                                 l, "])")))
        eval(parse(text = paste0("D2w <- D2w + sum((samp_frame$", 
                                 target_vars[i], "[samp_frame$", id_PSU, " == ", 
                                 l, "] - mu_L)^2)")))
      }
      eval(parse(text = paste0("rho$RHO_NAR", i, "[", 
                               k, "] <- 1 - (length(L) / (length(L) - 1) * (D2w / D2y))")))
      eval(parse(text = paste0("rho$RHO_NAR", i, "[", 
                               k, "] <- 1 - (D2w / D2y)")))
    }
  }
  strat_mun <- samp_frame[, c(strata_var, id_PSU)]
  strat_mun <- strat_mun[!duplicated(strat_mun), ]
  st <- paste0("mun <- aggregate(samp_frame$one,by=list(samp_frame$", 
               id_PSU, "),sum)")
  eval(parse(text = st))
  colnames(mun) <- c(id_PSU, "N")
  mun <- merge(mun, strat_mun)
  psu_file <- mun[, c(id_PSU, strata_var, "N")]
  colnames(psu_file) <- c("PSU_ID", "STRATUM", "PSU_MOS")
  des_file <- aggregate(psu_file$PSU_MOS, by = list(psu_file$STRATUM), 
                        sum)
  colnames(des_file) <- c("STRATUM", "STRAT_MOS")
  des_file$DELTA <- delta
  des_file$MINIMUM <- minimum
  #--------------------
  # check values in rho
  for (i in c(1:nvarY)) {
    st <- paste0("rho$RHO_NAR",i," <- ifelse(is.nan(rho$RHO_NAR",i,"),0,rho$RHO_NAR",i,")")
    eval(parse(text=st))
    }
  #--------------------
  out <- list(strata = strata, deff = deff, effst = effst, 
              rho = rho, psu_file = psu_file, des_file = des_file)
  return(out)
}
