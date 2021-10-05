select_PSU <- function(alloc, type="ALLOC", pps=TRUE)
{
  univ <- alloc$psu_trs
  univ <- univ[order(univ$STRATUM,-univ$PSU_MOS),]
  minPSUstr <- alloc$param_alloc$p_minPSUstrat
  minSSUstr <- alloc$param_alloc$p_minnumstrat
  univ$minPSUstr <- alloc$param_alloc$p_minPSUstrat
  univ$minSSUstr <- alloc$param_alloc$p_minnumstrat
  univ <- merge(univ,alloc$file_strata[,c("STRATUM","N")],by="STRATUM")
  univ <- merge(univ,alloc$alloc[-nrow(alloc$alloc),c("STRATUM",type)],by="STRATUM")
  univ$f <- univ[,type]/univ$N
  univ$THRESHOLD_NAR <- univ$THRESHOLD * univ$minPSUstr
  univ$partial <- 0
  univ$SUB[1] <- 1
   for (i in 2:nrow(univ)){
     # SR
     if(univ$STRATUM[i]==univ$STRATUM[i-1] & univ$AR[i]==1){univ$SUB[i]=univ$SUB[i-1]+1}
     #NSR
     ### the first
     if(univ$STRATUM[i]==univ$STRATUM[i-1] & univ$AR[i]==0 & univ$AR[i-1]==1){univ$SUB[i]=univ$SUB[i-1]+1; univ$partial[i] <- univ$PSU_MOS[i]}
     if(univ$STRATUM[i]!=univ$STRATUM[i-1] & univ$AR[i]==0){univ$SUB[i]=1; univ$partial[i] <- univ$PSU_MOS[i]}
     ### from the second  
     if(univ$STRATUM[i]==univ$STRATUM[i-1] & 
        univ$AR[i]==0 & univ$AR[i-1]==0 & 
        univ$partial[i-1] <= univ$THRESHOLD_NAR[i] & 
        (univ$partial[i-1] + univ$PSU_MOS[i]) <= univ$THRESHOLD_NAR[i])
        {univ$SUB[i]=univ$SUB[i-1]; univ$partial[i] <- (univ$partial[i-1] + univ$PSU_MOS[i])}
     if(univ$STRATUM[i]==univ$STRATUM[i-1] & 
        univ$AR[i]==0 & univ$AR[i-1]==0 & 
        univ$partial[i-1] <= univ$THRESHOLD_NAR[i] & 
        (univ$partial[i-1] + univ$PSU_MOS[i]) > univ$THRESHOLD_NAR[i])
     {univ$SUB[i]=univ$SUB[i-1]; univ$partial[i] <- (univ$partial[i-1] + univ$PSU_MOS[i])}
     if(univ$STRATUM[i]==univ$STRATUM[i-1] & 
        univ$AR[i]==0 & univ$AR[i-1]==0 & 
        univ$partial[i-1] > univ$THRESHOLD_NAR[i]) 
        {univ$SUB[i]=univ$SUB[i-1]+1; univ$partial[i] <- univ$PSU_MOS[i]}
   }
   univ$SUBSTRAT <- paste(univ$STRATUM,univ$SUB,sep="")
   psu_strat <- aggregate(PSU_ID~SUBSTRAT,univ,length)
   colnames(psu_strat)[2] <- "PSU_strat"
   univ <- merge(univ[,which(colnames(univ)!="PSU_strat")],psu_strat,by="SUBSTRAT")
   
   # check
   #------------------------------------------------------
   # univ[which(univ$PSU_strat==minPSUstr),]$AR <- 1
   #------------------------------------------------------
   univ$AR <- ifelse(univ$PSU_strat<=minPSUstr,1,univ$AR)
   #------------------------------------------------------
   univ <- univ[order(univ$STRATUM,-univ$PSU_MOS),]
   univ$partial <- 0
   univ$SUB[1] <- 1
   for (i in 2:nrow(univ)){
     # SR
     if(univ$STRATUM[i]==univ$STRATUM[i-1] & univ$AR[i]==1){univ$SUB[i]=univ$SUB[i-1]+1}
     #NSR
     ### the first
     if(univ$STRATUM[i]==univ$STRATUM[i-1] & univ$AR[i]==0 & univ$AR[i-1]==1){univ$SUB[i]=univ$SUB[i-1]+1; univ$partial[i] <- univ$PSU_MOS[i]}
     if(univ$STRATUM[i]!=univ$STRATUM[i-1] & univ$AR[i]==0){univ$SUB[i]=1; univ$partial[i] <- univ$PSU_MOS[i]}
     ### from the second  
     if(univ$STRATUM[i]==univ$STRATUM[i-1] & 
        univ$AR[i]==0 & univ$AR[i-1]==0 & 
        univ$partial[i-1] <= univ$THRESHOLD_NAR[i] & 
        (univ$partial[i-1] + univ$PSU_MOS[i]) <= univ$THRESHOLD_NAR[i])
     {univ$SUB[i]=univ$SUB[i-1]; univ$partial[i] <- (univ$partial[i-1] + univ$PSU_MOS[i])}
     if(univ$STRATUM[i]==univ$STRATUM[i-1] & 
        univ$AR[i]==0 & univ$AR[i-1]==0 & 
        univ$partial[i-1] <= univ$THRESHOLD_NAR[i] & 
        (univ$partial[i-1] + univ$PSU_MOS[i]) > univ$THRESHOLD_NAR[i])
     {univ$SUB[i]=univ$SUB[i-1]; univ$partial[i] <- (univ$partial[i-1] + univ$PSU_MOS[i])}
     if(univ$STRATUM[i]==univ$STRATUM[i-1] & 
        univ$AR[i]==0 & univ$AR[i-1]==0 & 
        univ$partial[i-1] > univ$THRESHOLD_NAR[i]) 
     {univ$SUB[i]=univ$SUB[i-1]+1; univ$partial[i] <- univ$PSU_MOS[i]}
   }
   univ$SUBSTRAT <- paste(univ$STRATUM,univ$SUB,sep="")
   psu_strat <- aggregate(PSU_ID~SUBSTRAT,univ,length)
   colnames(psu_strat)[2] <- "PSU_strat"
   univ <- merge(univ[,which(colnames(univ)!="PSU_strat")],psu_strat,by="SUBSTRAT")
   univ$POPAR = univ$PSU_MOS * univ$AR
   univ$POPNAR = univ$PSU_MOS - univ$POPAR
   # selection
   SUBSTRAT_MOS <- aggregate(PSU_MOS~SUBSTRAT,univ,sum)
   colnames(SUBSTRAT_MOS)[2] <- "SUBSTRAT_MOS"
   univ <- merge(univ, SUBSTRAT_MOS, by="SUBSTRAT")
   univ$pik <- minPSUstr*(univ$PSU_MOS/univ$SUBSTRAT_MOS)
   #-----------------------------------------------------
   univ$pik <- ifelse(univ$AR==1,1,univ$pik)
   univ$pik <- ifelse(univ$pik > 1,1,univ$pik)
   #-----------------------------------------------------
   if (pps==FALSE){univ$pik <- 1/univ$PSU_strat}
   s <- univ[univ$AR==1,]$PSU_ID
   SUBSTRAT <- unique(univ[univ$AR==0,]$SUBSTRAT)
   for (i in 1:length(SUBSTRAT)){
     U <- univ[univ$SUBSTRAT==SUBSTRAT[i],]
     U$I_s <- UPsampford(U$pik)  
     sh <- U[U$I_s==1,]$PSU_ID
     s <- c(s,sh)
   }
   sample <- data.frame(PSU_ID=s,sampled=1)
   univ <- merge(univ,sample,by="PSU_ID",all.x=TRUE)
   univ[is.na(univ$sampled),]$sampled <- 0
   univ$weight_Ist <- 1/univ$pik
   univ$n <- ifelse(univ$sampled==0,min(round(univ$f * univ$PSU_MOS,0)),minSSUstr)
   univ$weight_IIst <- ifelse(univ$sampled==0,0,univ$PSU_MOS/univ$n)
   univ$weight <- univ$weight_Ist * univ$weight_IIst 
   universe_PSU <- univ[,c("PSU_ID","STRATUM","SUB","SUBSTRAT","AR","PSU_strat","ALLOC","SUBSTRAT_MOS","PSU_MOS","pik","weight","sampled","n")]
   sample_PSU <- universe_PSU[universe_PSU$sampled==1,]
   
   #--------------------------------------------------------------------------------
   # Compute SSU to be sampled in each PSU
   population_stratum <- as.data.frame(tapply(univ$PSU_MOS,univ$STRATUM,sum))
   colnames(population_stratum) <- c("STRATUM_MOS")
   population_stratum$STRATUM <- row.names(population_stratum)
   sample_PSU <- merge(sample_PSU,population_stratum)
   
   sample_substratum <- as.data.frame(tapply(sample_PSU$PSU_MOS,sample_PSU$SUBSTRAT,sum))
   colnames(sample_substratum) <- c("SUBSTRAT_MOS2")
   sample_substratum$SUBSTRAT <- row.names(sample_substratum)
   sample_PSU <- merge(sample_PSU,sample_substratum)
   
   PSU_substratum <- as.data.frame(table(sample_PSU$SUBSTRAT))
   colnames(PSU_substratum) <- c("SUBSTRAT","PSU_substrat")
   sample_PSU <- merge(sample_PSU,PSU_substratum)
   sample_PSU$ALLOC_SUBSTR <- round(sample_PSU$ALLOC * sample_PSU$SUBSTRAT_MOS / sample_PSU$STRATUM_MOS)
   sample_PSU$PSU_final_sample_unit <- round(sample_PSU$ALLOC_SUBSTR /  sample_PSU$PSU_substrat)
   sample_PSU$SR <- ifelse(sample_PSU$AR == 1,1,0)
   sample_PSU$nSR <- ifelse(sample_PSU$AR == 0,1,0)
   sample_PSU$stratum <- sample_PSU$SUBSTRAT
   sample_PSU$Pik <- sample_PSU$pik
   sample_PSU$weight_1st <- 1/sample_PSU$Pik
   sample_PSU$weight_2st <- sample_PSU$PSU_MOS/sample_PSU$PSU_final_sample_unit
   sample_PSU$weight <- sample_PSU$weight_1st * sample_PSU$weight_2st
   sample_PSU <- sample_PSU[,c("PSU_ID","STRATUM","stratum","SR","nSR","PSU_final_sample_unit","Pik","weight_1st","weight_2st","weight")] 
   

   out <- list(universe_PSU=universe_PSU,sample_PSU=sample_PSU)
   out
}
