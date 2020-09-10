#####################################
# Prepare inputs for R2BEAT::beat.2st
# Author: Giulio Barcaroli
#####################################
# 1. psu_file
# 2. des_file
#####################################


input_to_beat.2st_2 <- function(psu,psu_id,stratum_var,mos_var,delta,minimum) {
  
  psu_file <- NULL
  psu_file$PSU_ID <- psu[,paste0(psu_id)]
  psu_file$STRATUM <- psu[,paste0(stratum_var)]
  psu_file$PSU_MOS <- psu[,paste0(mos_var)]
  psu_file <- as.data.frame(psu_file)
  
  des_file <- aggregate(PSU_MOS ~ STRATUM, data=psu_file,sum)
  des_file$DELTA <- delta
  des_file$MINIMUM <- minimum
  colnames(des_file)[2] <- "STRAT_MOS"
  out <- list(psu_file=psu_file,des_file=des_file)
  return(out)
}


