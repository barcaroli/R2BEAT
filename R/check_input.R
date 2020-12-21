check_input <- function(strata,
                        des,
                        strata_var_strata,
                        strata_var_des) {
  a <- merge(strata[,c(strata_var_strata,"N")],des[,c(strata_var_des,"STRAT_MOS")])
  a <- a[order(as.numeric(a$STRATUM)),]
  a$relative_difference <- round((a$STRAT_MOS - a$N) / a$STRAT_MOS,3)
  colnames(a)[2:3] <- c("N_in_strata","N_in_PSUs")
  cat("\n--------------------------------------------------")
  cat("\n Differences between population in strata and PSUs  ")
  cat("\n--------------------------------------------------")
  cat("\n")
  print(a)
  cat("\n--------------------------------------------------")
  cat("\nPopulation of PSUs has been attributed to strata")
  strata <- merge(strata,des[,c(strata_var_des,"STRAT_MOS")])
  strata$N <- strata$STRAT_MOS
  strata$STRAT_MOS <- NULL
  return(strata)
}