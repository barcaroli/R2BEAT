#-----------------------------------
# Function to prepare inputs to 
# allocation and PSU selection steps
# (scenario 2)
#-----------------------------------

# warning("The non-CRAN package, ReGenesees, is needed.")
# Add this line to NAMESPACE when ReGenesees is on CRAN
# importFrom(ReGenesees,svystat)
# and add ReGenesees to Depends in DESCRIPTION

prepareInputToAllocation2 <- function(
  samp_frame,
  RGdes, 
  RGcal, 
  id_PSU, 
  id_SSU, 
  strata_var, 
  target_vars, 
  deff_var, 
  domain_var,
  delta,
  minimum
) 
{
  one <- input_to_beat.2st_1(RGdes, 
                             RGcal, 
                             id_PSU, 
                             id_SSU, 
                             strata_vars=strata_var, 
                             target_vars, 
                             deff_vars=deff_var, 
                             domain_vars=domain_var)
  samp_frame$one <- 1
  psu <- eval(parse(text=paste0("aggregate(one~",id_PSU,"+",strata_var,",data=samp_frame,FUN=sum)")))
  mos_var <- "one"
  if (length(strata_var) == 1) st <- paste0("samp_frame$stratum_new <- samp_frame$",strata_var[1],"")
  if (length(strata_var) > 1) {
    st <- NULL
    for (i in c(1:length(strata_var))) {
      if (i == 1) st1 <- paste0("samp_frame$stratum_new <- cbind(samp_frame[,'",strata_var[i],"']")
      if (i > 1 & i != length(strata_var)) st1 <- paste0(",samp_frame[,'",strata_var[i],"']")
      if (i == length(strata_var)) st1 <- paste0(",samp_frame[,'",strata_var[i],"'])")
      st <- paste0(st,st1)
    }
    eval(parse(text=st))
  }
  two <- input_to_beat.2st_2(psu,
                             psu_id=id_PSU,
                             stratum_var=strata_var,
                             mos_var,
                             delta,
                             minimum)
  out <- list(strata = one[[1]], deff = one[[2]], effst = one[[3]], rho = one[[4]],
              psu_file=two[[1]],des_file=two[[2]])

  return(out)
}
