select_PSU2 <- function (alloc, 
                         type = "ALLOC", 
                         var_ord = NULL, 
                         des_file = des_file,
                         psu_file = psu_file) 
{
  allocation <- alloc$alloc[-nrow(alloc$alloc), c("STRATUM",type)]
  #------------------------------------------------------------------
  strat <- unique(des_file$STRATUM)
  nstrat <- length(strat)
  for (i in c(1:length(unique(des_file$STRATUM)))) {
    minimum <- des_file$MINIMUM[des_file$STRATUM==strat[i]]
    delta <- des_file$DELTA[des_file$STRATUM==strat[i]]
    allocation$PSUs[allocation$STRATUM == strat[i]] <- round(allocation$ALLOC[allocation$STRATUM == strat[i]]/(minimum*delta))
  }
  #-------------------------------------------------------------------
  allocation$PSUs <- ifelse(allocation$PSUs == 0, 1, allocation$PSUs)
  if (is.null(psu_file)) psu_file <- alloc$psu_trs[, c(4, 2, 5)]
  psu_file$ones <- 1
  PSUs <- aggregate(cbind(ones, PSU_MOS) ~ STRATUM, data = psu_file, FUN = sum)
  PSUs <- merge(PSUs, allocation)
  colnames(PSUs) <- c("STRATUM", "PSUs_available", "SSUs_available", "SSUs_allocated", "PSUs_allocated")
  PSUs <- PSUs[, c(1, 2, 5, 3, 4)]
  psu_frame <- merge(psu_file, PSUs, by = c("STRATUM"))
  if (!is.null(var_ord)) {
    psu_frame <- merge(psu_frame, var_ord, by = "PSU_ID")
    vars <- colnames(var_ord)[c(2:ncol(var_ord))]
    stmt <- NULL
    for (v in c(1:length(vars))) {
      if (v == 1) {
        stmt <- paste0("psu_frame <- psu_frame[order(psu_frame$", vars[v], ",")
      }
      if (v > 1 & v < length(vars)) {
        stmt <- paste0(stmt, "psu_frame", vars[v], ",")
      }
      if (v == length(vars)) {
        stmt <- paste0(stmt, "psu_frame$", vars[v], "),]")
      }
    }
    eval(parse(text = stmt))
  }
  PSU_stats <- as.data.frame(list(STRATUM = strat,
                                  PSU=rep(NA,nstrat),
                                  SSU=rep(NA,nstrat),
                                  DELTA=round(des_file$DELTA,2),
                                  Units=rep(NA,nstrat)))
  k = 0
  strats <- c(unique(psu_frame$STRATUM))
  strats <- strats[order(strats)]
  for (i in c(strats)) {
    k <- k+1
    psu_to_be_selected <- PSUs$PSUs_allocated[PSUs$STRATUM == i]
    psu_frame$pik_1st[psu_frame$STRATUM == i] <- inclusionprobabilities(psu_frame$PSU_MOS[psu_frame$STRATUM == i], psu_to_be_selected)
    try(psu_frame$sampled[psu_frame$STRATUM == i] <- UPsystematic(pik = psu_frame$pik_1st[psu_frame$STRATUM == i]))
    PSU_stats$STRATUM[k] <- i
    PSU_stats$PSU[k] <- sum(psu_frame$sampled[psu_frame$STRATUM == i & psu_frame$sampled == 1])
    PSU_stats$SSU[k] <- sum(psu_frame$sampled[psu_frame$STRATUM == i & psu_frame$sampled == 1]) * des_file$MINIMUM[des_file$STRATUM == i]
    PSU_stats$Units[k] <- round(sum(psu_frame$sampled[psu_frame$STRATUM == i & psu_frame$sampled == 1]) * des_file$MINIMUM[des_file$STRATUM == i] * des_file$DELTA[des_file$STRATUM == i])
    cat("\n Stratum: ", i, "  PSUs to be selected: ", psu_to_be_selected, 
        "  PSUs selected: ", sum(psu_frame$sampled[psu_frame$STRATUM == i & psu_frame$sampled == 1]))
  }
  PSU_stats <- rbind(PSU_stats,
                     c("Total",sum(PSU_stats$PSU),
                       sum(PSU_stats$SSU),
                       round(sum(PSU_stats$Units)/sum(PSU_stats$SSU),2),
                       sum(PSU_stats$Units)))
  
  
  sample_PSUs <- psu_frame[psu_frame$sampled == 1, ]
  for (i in c(unique(psu_frame$STRATUM))) {
    sample_PSUs$PSU_final_sample_unit <- des_file$MINIMUM[des_file$STRATUM == i]
    sample_PSUs$pik_2st[sample_PSUs$STRATUM == i] <- des_file$MINIMUM[des_file$STRATUM == i]/sample_PSUs$PSU_MOS[sample_PSUs$STRATUM == i]
  }
  sample_PSUs$weight <- (1/sample_PSUs$pik_1st) * (1/sample_PSUs$pik_2st)
  sample_PSUs$ones <- NULL
  sample_PSUs$PSUs_available <- NULL
  sample_PSUs$PSUs_allocated <- NULL
  sample_PSUs$SSUs_available <- NULL
  sample_PSUs$SSUs_allocated <- NULL
  sample_PSUs$sampled <- NULL
  out <- list(PSU_stats=PSU_stats,
              sample_PSUs=sample_PSUs)
  return(out)
}
