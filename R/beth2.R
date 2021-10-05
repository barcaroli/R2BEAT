beat.2st <- function (stratif, 
                      errors, 
                      des_file, 
                      psu_file, 
                      rho, 
                      deft_start = NULL, 
                      effst = NULL, 
                      epsilon1 = 5, 
                      mmdiff_deft = 1, 
                      maxi = 20, 
                      epsilon = 10^(-11), 
                      minPSUstrat = 2, 
                      minnumstrat = 2, 
                      maxiter = 200, 
                      maxiter1 = 25) 
{
  diffx = 999
  iterx = 0
  stages <- NULL
  test_stages <<- 2
  param_alloc <- as.data.frame(t(c(epsilon, minPSUstrat, minnumstrat, mean(des_file$MINIMUM), maxiter, 
                                   maxiter1, epsilon1, diffx, iterx, mmdiff_deft, maxi, 
                                   2)))
  names(param_alloc) = c("p_epsilon", "p_minPSUstrat", "p_minnumstrat", "minimum",
                         "p_maxiter", "p_maxiter1", "p_epsilon1", 
                         "p_diffx", "p_iterx", "p_mmdiff_deft", 
                         "p_maxi", "num_stages")
  nvar = ncol(errors) - 1
  if (!is.null(effst)) {
    colnames(stratif) <- toupper(colnames(stratif))
    sapply(1:nvar, function(i) {
      stratif[, paste("S_ORIG_no_effst", i, sep = "")] <<- stratif[, 
                                                                   paste("S", i, sep = "")]
    })
    sapply(1:nvar, function(i) {
      stratif[, paste("S_ORIG", i, sep = "")] <<- stratif[, 
                                                          paste("S", i, sep = "")] * effst[, 
                                                                                           paste("EFFST", i, sep = "")]
    })
    sapply(1:nvar, function(i) {
      stratif[, paste("S", i, sep = "")] <<- (stratif[, 
                                                      paste("S", i, sep = "")] * effst[, 
                                                                                       paste("EFFST", i, sep = "")])
    })
  }
  else {
    colnames(stratif) <- toupper(colnames(stratif))
    sapply(1:nvar, function(i) {
      stratif[, paste("S_ORIG", i, sep = "")] <<- stratif[, 
                                                          paste("S", i, sep = "")]
    })
  }
  colnames(errors) <- toupper(colnames(errors))
  colnames(des_file) <- toupper(colnames(des_file))
  colnames(rho) <- toupper(colnames(rho))
  colnames(psu_file) <- toupper(colnames(psu_file))
  deft <- NULL
  deft <- data.frame(stratif$STRATUM)
  names(deft) = "STRATUM"
  if (is.null(deft_start)) {
    sapply(1:nvar, function(i) {
      deft[, paste("DEFT", i, sep = "")] <<- 1
    })
  }
  else {
    deft <- merge(deft, deft_start, by = "STRATUM", 
                  all.x = T)
  }
  nstrat = nrow(stratif)
  strloop <- c(1:nstrat)
  ndom = nrow(errors)
  param_alloc$nvar <- nvar
  param_alloc$ndom <- ndom
  param_alloc$nstrat <- nstrat
  param_alloc <<- param_alloc
  n_1st <- NULL
  ob <- beat.1st(stratif, errors, minnumstrat)
  n <- ob$n
  n_1st <- n
  iterations <- data.frame(cbind(0, 0, 0, 0, sum(n)))
  colnames(iterations) <- c("iterations", "PSU_SR", 
                            "PSU NSR", "PSU Total", "SSU")
  des_file <- des_file[order(des_file$STRATUM), ]
  des_file$CAMP <- n
  infoloop <- merge(deft, stratif[, c("STRATUM", sapply(1:nvar, 
                                                        function(i) (paste("S", i, sep = ""))))], 
                    by = "STRATUM")
  Bethel_sample_n12 <- NULL
  i_arnar_fin <- NULL
  i_diff_DEFT_fin <- NULL
  i_S_DEFT_loop <- NULL
  newDeft <- function(des_file, stratif, errors, rho, psu_file, 
                      deft) {
    while (diffx > epsilon1 && iterx < maxi && mmdiff_deft > 
           0.06) {
      iterx <- iterx + 1
      des_file$F = (des_file$CAMP/des_file$STRAT_MOS)
      des_file$THRESHOLD = ((des_file$MINIMUM/des_file$F) * des_file$DELTA)
      psu_file1 <- merge(des_file[, c("STRATUM", "THRESHOLD")], psu_file, by = "STRATUM")
      psu_file1$AR <- 0
      psu_file1$AR[(psu_file1$PSU_MOS >= psu_file1$THRESHOLD)] <- 1
      table(psu_file1$AR)
      ###
      psu_file1 <- psu_file1[order(psu_file1$STRATUM,-psu_file1$PSU_MOS),]
      psu_file1$SUB <- 0
      psu_file1$partial <- 0
      psu_file1$SUB[1] <- 1
      psu_file1$THRESHOLD_NAR <- psu_file1$THRESHOLD * minPSUstrat
      for (i in 2:nrow(psu_file1)){
        # SR
        if(psu_file1$STRATUM[i]==psu_file1$STRATUM[i-1] & psu_file1$AR[i]==1){psu_file1$SUB[i]=psu_file1$SUB[i-1]+1}
        #NSR
        ### the first
        if(psu_file1$STRATUM[i]==psu_file1$STRATUM[i-1] & psu_file1$AR[i]==0 & psu_file1$AR[i-1]==1){psu_file1$SUB[i]=psu_file1$SUB[i-1]+1; psu_file1$partial[i] <- psu_file1$PSU_MOS[i]}
        if(psu_file1$STRATUM[i]!=psu_file1$STRATUM[i-1] & psu_file1$AR[i]==0){psu_file1$SUB[i]=1; psu_file1$partial[i] <- psu_file1$PSU_MOS[i]}
        ### from the second  
        if(psu_file1$STRATUM[i]==psu_file1$STRATUM[i-1] & 
           psu_file1$AR[i]==0 & psu_file1$AR[i-1]==0 & 
           psu_file1$partial[i-1] <= psu_file1$THRESHOLD_NAR[i] & 
           (psu_file1$partial[i-1] + psu_file1$PSU_MOS[i]) <= psu_file1$THRESHOLD_NAR[i])
        {psu_file1$SUB[i]=psu_file1$SUB[i-1]; psu_file1$partial[i] <- (psu_file1$partial[i-1] + psu_file1$PSU_MOS[i])}
        if(psu_file1$STRATUM[i]==psu_file1$STRATUM[i-1] & 
           psu_file1$AR[i]==0 & psu_file1$AR[i-1]==0 & 
           psu_file1$partial[i-1] <= psu_file1$THRESHOLD_NAR[i] & 
           (psu_file1$partial[i-1] + psu_file1$PSU_MOS[i]) > psu_file1$THRESHOLD_NAR[i])
        {psu_file1$SUB[i]=psu_file1$SUB[i-1]; psu_file1$partial[i] <- (psu_file1$partial[i-1] + psu_file1$PSU_MOS[i])}
        if(psu_file1$STRATUM[i]==psu_file1$STRATUM[i-1] & 
           psu_file1$AR[i]==0 & psu_file1$AR[i-1]==0 & 
           psu_file1$partial[i-1] > psu_file1$THRESHOLD_NAR[i]) 
        {psu_file1$SUB[i]=psu_file1$SUB[i-1]+1; psu_file1$partial[i] <- psu_file1$PSU_MOS[i]}
      }
      psu_file1$SUBSTRAT <- paste(psu_file1$STRATUM,psu_file1$SUB,sep="")
      psu_strat <- aggregate(PSU_ID~SUBSTRAT,psu_file1,length)
      colnames(psu_strat)[2] <- "PSU_strat"
      psu_file1 <- merge(psu_file1,psu_strat,by="SUBSTRAT")
      # check
      psu_file1[which(psu_file1$PSU_strat==minPSUstrat),]$AR <- 1
      psu_file1 <- psu_file1[order(psu_file1$STRATUM,-psu_file1$PSU_MOS),]
      psu_file1$partial <- 0
      psu_file1$SUB[1] <- 1
      ###
      psu_file1$POPAR = psu_file1$PSU_MOS * psu_file1$AR
      psu_file1$POPNAR = psu_file1$PSU_MOS - psu_file1$POPAR
      popolaz <- aggregate(psu_file1[, c("POPAR", "POPNAR")], 
                           by = list(STRATUM = psu_file1$STRATUM), 
                           sum)
      des_file <- merge(des_file, popolaz, by = "STRATUM")
      des_file <- merge(des_file, rho, by = "STRATUM")
      des_file$CAMPAR = round((des_file$F) * (des_file$POPAR))
      des_file$CAMPNAR = des_file$CAMP - des_file$CAMPAR
      des_file$BDIS_AR = 1 * des_file$DELTA
      des_file$BDIS_NAR = des_file$MINIMUM * des_file$DELTA
      sapply(1:nvar, function(i) {
        des_file[, paste("DEFF", i, sep = "")] <<- ((des_file$CAMP/(des_file$STRAT_MOS^2)) * 
                                                      ((((des_file$POPAR^2)/(des_file$CAMPAR + epsilon)) * 
                                                          (1 + ((des_file[, paste("RHO_AR", i, 
                                                                                  sep = "")]) * (des_file$BDIS_AR - 
                                                                                                   1)))) + (((des_file$POPNAR^2)/(des_file$CAMPNAR + 
                                                                                                                                    epsilon)) * (1 + ((des_file[, paste("RHO_NAR", 
                                                                                                                                                                        i, sep = "")]) * (des_file$BDIS_NAR - 
                                                                                                                                                                                            1))))))
      })
      sapply(1:nvar, function(i) {
        des_file[, paste("DEFF", i, sep = "")][des_file[, 
                                                        paste("DEFF", i, sep = "")] < 0] <<- 1
      })
      sapply(1:nvar, function(i) {
        des_file[, paste("DEFT", i, sep = "")] <<- sqrt(des_file[, 
                                                                 paste("DEFF", i, sep = "")])
      })
      stratif <- merge(stratif, des_file[, c("STRATUM", 
                                             sapply(1:nvar, function(i) (paste("DEFT", 
                                                                               i, sep = ""))))], by = "STRATUM")
      loopdeft <- deft
      sapply(1:nvar, function(i) {
        loopdeft[, paste("oldDEFT", i, sep = "")] <<- loopdeft[, 
                                                               paste("DEFT", i, sep = "")]
      })
      sapply(1:nvar, function(i) {
        loopdeft[, paste("DEFT", i, sep = "")] <<- des_file[, 
                                                            paste("DEFT", i, sep = "")]
      })
      sapply(1:nvar, function(i) {
        loopdeft[, paste("diffx_deft", i, sep = "")] <<- (loopdeft[, 
                                                                   paste("DEFT", i, sep = "")] - loopdeft[, 
                                                                                                          paste("oldDEFT", i, sep = "")])
      })
      sapply(1:nvar, function(i) {
        loopdeft[, paste("mdiff_deft", i, sep = "")] <<- max(abs(loopdeft[, 
                                                                          paste("diffx_deft", i, sep = "")]))
      })
      mmdiff_deft <- max(sapply(1:nvar, function(i) {
        loopdeft[, paste("mdiff_deft", i, sep = "")]
      }))
      if (iterx == 1 & !is.null(deft_start)) {
        sapply(1:nvar, function(i) {
          stratif[, paste("S", i, sep = "")] <<- (stratif[, 
                                                          paste("S_ORIG", i, sep = "")] * 
                                                    deft_start[, paste("DEFT", i, sep = "")])
        })
      }
      else {
        sapply(1:nvar, function(i) {
          stratif[, paste("S", i, sep = "")] <<- (stratif[, 
                                                          paste("S_ORIG", i, sep = "")] * 
                                                    des_file[, paste("DEFT", i, sep = "")])
        })
      }
      c <- nvar
      if (iterx == 1) {
        colnames(infoloop)[2:(2 + c - 1)] <- paste("DEFT", 
                                                   1:nvar, "_", iterx - 1, sep = "")
      }
      else {
        colnames(infoloop)[(2 + ((iterx - 1) * c * 2)):(2 + 
                                                          nvar - 1 + ((iterx - 1) * c * 2))] <- paste("DEFT", 
                                                                                                      1:nvar, "_", iterx - 1, sep = "")
      }
      if (iterx == 1) {
        colnames(infoloop)[(2 + c):(2 + 2 * c - 1)] <- paste("S", 
                                                             1:nvar, "_", iterx - 1, sep = "")
      }
      else {
        colnames(infoloop)[((2 + c) + ((iterx - 1) * 
                                         c * 2)):((2 + c) + nvar - 1 + ((iterx - 1) * 
                                                                          c * 2))] <- paste("S", 1:nvar, "_", 
                                                                                            iterx - 1, sep = "")
      }
      infoloop <- merge(infoloop, des_file[, c("STRATUM", 
                                               sapply(1:nvar, function(i) (paste("DEFT", 
                                                                                 i, sep = ""))))], by = "STRATUM", 
                        suffixes = rep("", 2))
      infoloop <- merge(infoloop, stratif[, c("STRATUM", 
                                              sapply(1:nvar, function(i) (paste("S", 
                                                                                i, sep = ""))))], by = "STRATUM", 
                        suffixes = rep("", 2))
      ob <- beat.1st(stratif, errors, minnumstrat)
      n <- ob$n
      #------------------------------------------------------------------------
      # PSU_SR <- sum(psu_file1$AR)
      # PSU_NSR <- round(sum(aggregate(psu_file1$POPNAR/psu_file1$THRESHOLD, 
      #                                by = list(psu_file1$STRATUM), sum)[, 2]))
      psu_file1[which(psu_file1$PSU_strat==minPSUstrat),]$AR <- 1
      PSU_SR <- sum(aggregate(AR~SUBSTRAT,data=psu_file1,sum)[2])
      PSU_NSR <- length(unique(psu_file1$SUBSTRAT[psu_file1$AR==0]))*minPSUstrat
      PSU_Total <- PSU_SR + PSU_NSR
      #------------------------------------------------------------------------
      
      iter_new <- c(iterx, PSU_SR, PSU_NSR, PSU_Total, 
                    sum(n))
      iterations <- rbind(iterations, iter_new)
      diffx <- max(abs(n - des_file$CAMP))
      loopdeft$diff_n <- (abs(n - des_file$CAMP))
      deft <- des_file[, c("STRATUM", sapply(1:nvar, 
                                             function(i) (paste("DEFT", i, sep = ""))))]
      output_beth12 <- list(i_diff_DEFT_fin = loopdeft, 
                            i_S_DEFT_loop = infoloop, i_arnar_fin = des_file, 
                            n_1st = n_1st, psu_trs = psu_file1)
      output_beth12 <<- output_beth12
      des_file <- des_file[, c("STRATUM", "STRAT_MOS", 
                               "MINIMUM", "DELTA")]
      des_file$CAMP <- n
      stratif <- stratif[, c("STRATUM", "N", 
                             sapply(1:nvar, function(i) paste("M", i, 
                                                              sep = "")), sapply(1:nvar, function(i) paste("S", 
                                                                                                           i, sep = "")), sapply(1:nvar, function(i) paste("S_ORIG", 
                                                                                                                                                           i, sep = "")), sapply(1:ndom, function(i) paste("DOM", 
                                                                                                                                                                                                           i, sep = "")), "COST", "CENS")]
    }
    # print(iterations)
    obb <- list(ob = ob, iterations = iterations)
    return(obb)
  }
  ob_fin <- newDeft(des_file, stratif, errors, rho, psu_file,deft)
  iterations <- ob_fin$iteration
  
  ob_fin <- ob_fin$ob
  n_2st <- ob_fin$n
  Bethel_sample_n12 <- cbind(stratif[, c("STRATUM", "N", 
                                         sapply(1:ndom, function(i) paste("DOM", i, sep = "")))], 
                             n_1st, n_2st)
  test_stages <<- 0
  output_beth12$n2_st <- n_2st
  output_beth12$Bethel_sample_n12 <- (cbind(stratif[, c("STRATUM", 
                                                        "N", sapply(1:ndom, function(i) paste("DOM", 
                                                                                              i, sep = "")))], n_1st, n_2st))
  psu_trs <- output_beth12$psu_trs 
  
  #----------------------------------------------------------
  psu_trs[which(psu_trs$PSU_strat<=minPSUstrat),]$AR <- 1
  iterations$`PSU_SR`[nrow(alloc$iterations)] <- sum(aggregate(AR~SUBSTRAT,data=psu_trs,sum)[2])
  iterations$`PSU NSR`[nrow(alloc$iterations)] <- length(unique(psu_trs$SUBSTRAT[psu_trs$AR==0]))*minPSUstrat
  iterations$`PSU Total`[nrow(alloc$iterations)] <- sum(aggregate(AR~SUBSTRAT,data=psu_trs,sum)[2]) + length(unique(psu_trs$SUBSTRAT[psu_trs$AR==0]))*minPSUstrat
  print(iterations)
  #----------------------------------------------------------
  
  
  output_beth12 <<- output_beth12
  file_strata <- ob_fin$file_strata
  alloc <- ob_fin$alloc
  alloc[, -1] <- round(alloc[, -1], digits = 0)
  cv <- ob_fin$sensitivity
  planned <- cv[, c(1:4)]
  varn <- unique(cv[, 3])
  step <- length(varn)
  namesP <- planned[seq(from = 1, to = dim(planned)[1], by = step), 
                    1:2]
  valueP <- matrix(round(planned[, 4], 5), ncol = step, byrow = TRUE)
  planned <- cbind(namesP, valueP)
  colnames(planned) <- c(colnames(namesP), paste(varn))
  expected <- cv[, c(1:3, 5)]
  step <- length(unique(expected[, 3]))
  namesE <- expected[seq(from = 1, to = dim(expected)[1], by = step), 
                     1:2]
  valueE <- matrix(round(expected[, 4], 5), ncol = step, byrow = TRUE)
  expected <- cbind(namesE, valueE)
  colnames(expected) <- c(colnames(namesE), paste(varn))
  sensitivity <- cv[, c(1:3, 6)]
  step <- length(unique(sensitivity[, 3]))
  namesS <- sensitivity[seq(from = 1, to = dim(sensitivity)[1], 
                            by = step), 1:2]
  valueS <- matrix(sensitivity[, 4], ncol = step, byrow = TRUE)
  sensitivity <- cbind(namesS, valueS)
  colnames(sensitivity) <- c(colnames(namesE), paste(varn))
  deft <- output_beth12$i_S_DEFT_loop
  deft_c <- cbind(STRATUM = file_strata$STRATUM, deft[, startsWith(colnames(deft), 
                                                                   "DEFT")])
 
  out <- list(iterations = iterations, file_strata = file_strata, 
              alloc = alloc, psu_trs=psu_trs, planned = planned, expected = expected, 
              sensitivity = sensitivity, deft_c = deft_c, param_alloc = param_alloc)
  return(out)
}

