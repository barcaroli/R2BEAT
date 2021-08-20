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
  strata_vars, 
  target_vars, 
  deff_vars, 
  domain_vars,
  mos_var,
  delta,
  minimum
) 
{
  one <- input_to_beat.2st_1(RGdes, 
                             RGcal, 
                             id_PSU, 
                             id_SSU, 
                             strata_vars, 
                             target_vars, 
                             deff_vars, 
                             domain_vars)
  samp_frame$one <- 1
  psu <- eval(parse(text=paste0("aggregate(one~",id_PSU,"+",strata_vars,",data=samp_frame,FUN=sum)")))
  two <- input_to_beat.2st_2(psu,
                             psu_id=id_PSU,
                             stratum_var=strata_vars,
                             mos_var,
                             delta,
                             minimum)
  out <- list(strata = one[[1]], deff = one[[2]], effst = one[[3]], rho = one[[4]],
              psu_file=two[[1]],des_file=two[[2]])

  return(out)
}