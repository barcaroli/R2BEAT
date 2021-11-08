select_PSU <- function (alloc, type = "ALLOC", pps = TRUE, plot = TRUE) 
{
  univ <- alloc$psu_trs
  if (length(unique(univ$PSU_ID)) < nrow(univ))
    stop("PSU identifier not unique")
  univ <- univ[order(univ$STRATUM, -univ$PSU_MOS), ]
  minPSUstr <- alloc$param_alloc$p_minPSUstrat
  minSSUstr <- alloc$param_alloc$p_minnumstrat
  minimum <- alloc$param_alloc$minimum
  univ$minPSUstr <- alloc$param_alloc$p_minPSUstrat
  univ$minSSUstr <- alloc$param_alloc$p_minnumstrat
  univ <- merge(univ, alloc$file_strata[, c("STRATUM", "N")], by = "STRATUM")
  univ <- merge(univ, alloc$alloc[-nrow(alloc$alloc), c("STRATUM", type)], by = "STRATUM")
  univ$f <- univ[, type]/univ$N
  univ$THRESHOLD_NAR <- univ$THRESHOLD * univ$minPSUstr
  univ$partial <- 0
  univ$SUB[1] <- 1
  #----------------------------------------------
  univ <- univ[order(univ$STRATUM,-univ$PSU_MOS),]
  #----------------------------------------------
  for (i in 2:nrow(univ)) {
    # if same stratum and is AR, then SUB + 1 
    if (univ$STRATUM[i] == univ$STRATUM[i - 1] & univ$AR[i] == 1) {
      univ$SUB[i] = univ$SUB[i - 1] + 1 
    }
    # if same stratum and not AR and previous is AR, then SUB + 1 and init partial
    if (univ$STRATUM[i] == univ$STRATUM[i - 1] & univ$AR[i] == 0 & univ$AR[i - 1] == 1) {
      univ$SUB[i] = univ$SUB[i - 1] + 1
      univ$partial[i] <- univ$PSU_MOS[i]
    }
    # if different stratum and AR, then SUB + 1 and init partial
    if (univ$STRATUM[i] != univ$STRATUM[i - 1] & univ$AR[i] == 0) {
      univ$SUB[i] = 1
      univ$partial[i] <- univ$PSU_MOS[i]
    }
    # if same stratum and both previous and current not AR and previous partial lt threshold
    #                 and updated partial lt threshold, then same SUB and add to partial
    if (univ$STRATUM[i] == univ$STRATUM[i - 1] 
        & univ$AR[i] == 0 & univ$AR[i - 1] == 0 
        & univ$partial[i - 1] <= univ$THRESHOLD_NAR[i] 
        & (univ$partial[i - 1] + univ$PSU_MOS[i]) <= univ$THRESHOLD_NAR[i]) {
      univ$SUB[i] = univ$SUB[i - 1]
      univ$partial[i] <- (univ$partial[i - 1] + univ$PSU_MOS[i])
    }
    # if same stratum and both previous and current not AR and previous partial lt threshold
    #                 and updated partial gt threshold, then new SUB and init partial
    if (univ$STRATUM[i] == univ$STRATUM[i - 1] 
        & univ$AR[i] == 0 & univ$AR[i - 1] == 0 
        & univ$partial[i - 1] <= univ$THRESHOLD_NAR[i] 
        & (univ$partial[i - 1] + univ$PSU_MOS[i]) > univ$THRESHOLD_NAR[i]) {
      # univ$SUB[i] = univ$SUB[i - 1]
      # univ$partial[i] <- (univ$partial[i - 1] + univ$PSU_MOS[i])
      univ$SUB[i] = univ$SUB[i - 1] + 1
      univ$partial[i] <- univ$PSU_MOS[i]
    }
    # if same stratum and both previous and current not AR and previous partial lt threshold
    #                 and updated partial gt threshold, then same SUB and add to partial
    # if (univ$STRATUM[i] == univ$STRATUM[i - 1] 
    #     & univ$AR[i] == 0 & univ$AR[i - 1] == 0 
    #     & univ$partial[i - 1] > univ$THRESHOLD_NAR[i]) {
    #   univ$SUB[i] = univ$SUB[i - 1] + 1
    #   univ$partial[i] <- univ$PSU_MOS[i]
    # }
  }
  univ$SUBSTRAT <- paste(univ$STRATUM, univ$SUB, sep = "-")
  psu_strat <- aggregate(PSU_ID ~ SUBSTRAT, univ, length)
  colnames(psu_strat)[2] <- "PSU_strat"
  univ <- merge(univ[, which(colnames(univ) != "PSU_strat")], psu_strat, by = "SUBSTRAT")
  #----------------------------------------------
  univ <- univ[order(univ$STRATUM,-univ$PSU_MOS),]
  #----------------------------------------------
  univ$AR <- ifelse(univ$PSU_strat <= minPSUstr, 1, univ$AR)
  # univ$partial <- ifelse(univ$PSU_strat <= minPSUstr, 0, univ$partial)
  univ$SUB[1] <- 1

  for (i in 2:nrow(univ)) {
    # if same stratum and is AR, then SUB + 1 
    if (univ$STRATUM[i] == univ$STRATUM[i - 1] & univ$AR[i] == 1) {
      univ$SUB[i] = univ$SUB[i - 1] + 1
    }
    # if same stratum and not AR and previous is AR, then SUB + 1 and init partial
    if (univ$STRATUM[i] == univ$STRATUM[i - 1] 
        & univ$AR[i] == 0 & univ$AR[i - 1] == 1) {
          univ$SUB[i] = univ$SUB[i - 1] + 1
          univ$partial[i] <- univ$PSU_MOS[i]
    }
    # if different stratum and AR, then SUB + 1 and init partial
    if (univ$STRATUM[i] != univ$STRATUM[i - 1] 
        & univ$AR[i] == 0) {
          univ$SUB[i] = 1
          univ$partial[i] <- univ$PSU_MOS[i]
    }
    # if same stratum and both previous and current not AR and previous partial lt threshold
    #                 and updated partial lt threshold, then same SUB and add to partial
    if (univ$STRATUM[i] == univ$STRATUM[i - 1] 
        & univ$AR[i] == 0 & univ$AR[i - 1] == 0 
        & univ$partial[i - 1] <= univ$THRESHOLD_NAR[i] 
        & (univ$partial[i - 1] + univ$PSU_MOS[i]) <= univ$THRESHOLD_NAR[i]) {
          univ$SUB[i] = univ$SUB[i - 1]
          univ$partial[i] <- (univ$partial[i - 1] + univ$PSU_MOS[i])
    }
    # if same stratum and both previous and current not AR and previous partial lt threshold
    #                 and updated partial gt threshold, then new SUB and init partial
    if (univ$STRATUM[i] == univ$STRATUM[i - 1] 
        & univ$AR[i] == 0 & univ$AR[i - 1] == 0 
        & univ$partial[i - 1] <= univ$THRESHOLD_NAR[i] 
        & (univ$partial[i - 1] + univ$PSU_MOS[i]) > univ$THRESHOLD_NAR[i]) {
      # univ$SUB[i] = univ$SUB[i - 1]
      # univ$partial[i] <- (univ$partial[i - 1] + univ$PSU_MOS[i])
          univ$SUB[i] = univ$SUB[i - 1] + 1
          univ$partial[i] <- univ$PSU_MOS[i]
    }
    # if same stratum and both previous and current not AR and previous partial lt threshold
    #                 and updated partial gt threshold, then same SUB and add to partial
    # if (univ$STRATUM[i] == univ$STRATUM[i - 1] 
    #     & univ$AR[i] == 0 & univ$AR[i - 1] == 0 
    #     & univ$partial[i - 1] > univ$THRESHOLD_NAR[i]) {
    #   univ$SUB[i] = univ$SUB[i - 1] + 1
    #   univ$partial[i] <- univ$PSU_MOS[i]
  }
  #-------------------------------------------------------
  # univ$SUBSTRAT <- paste(univ$STRATUM, univ$SUB, sep = "")
  univ$SUBSTRAT <- paste(as.character(univ$STRATUM), as.character(univ$SUB), sep = "-")
  #-------------------------------------------------------
  psu_strat <- aggregate(PSU_ID ~ SUBSTRAT, univ, length)
  colnames(psu_strat)[2] <- "PSU_strat"
  univ <- merge(univ[, which(colnames(univ) != "PSU_strat")], 
                psu_strat, by = "SUBSTRAT")
  univ$POPAR = univ$PSU_MOS * univ$AR
  univ$POPNAR = univ$PSU_MOS - univ$POPAR
  SUBSTRAT_MOS <- aggregate(PSU_MOS ~ SUBSTRAT, univ, sum)
  colnames(SUBSTRAT_MOS)[2] <- "SUBSTRAT_MOS"
  univ <- merge(univ, SUBSTRAT_MOS, by = "SUBSTRAT")
  univ$pik <- minPSUstr * (univ$PSU_MOS/univ$SUBSTRAT_MOS)
  univ$pik <- ifelse(univ$AR == 1, 1, univ$pik)
  univ$pik <- ifelse(univ$pik > 1, 1, univ$pik)
  if (pps == FALSE) {
    univ$pik <- 1/univ$PSU_strat
  }
  s <- univ[univ$AR == 1, ]$PSU_ID
  SUBSTRAT <- unique(univ[univ$AR == 0, ]$SUBSTRAT)
  for (i in 1:length(SUBSTRAT)) {
    # cat("\nSUBSTR ",SUBSTRAT[i])
    U <- univ[univ$SUBSTRAT == SUBSTRAT[i], ]
    U$I_s <- UPsampford(U$pik)
    sh <- U[U$I_s == 1, ]$PSU_ID
    # cat("\n sh",sh)
    s <- c(s, sh)
  }
  sample <- data.frame(PSU_ID = s, sampled = 1)
  univ <- merge(univ, sample, by = c("PSU_ID"), all.x = TRUE)
  univ[is.na(univ$sampled), ]$sampled <- 0
  univ$weight_Ist <- 1/univ$pik
  univ$n <- ifelse(univ$sampled == 0, min(round(univ$f * univ$PSU_MOS, 
                                                0)), minSSUstr)
  univ$weight_IIst <- ifelse(univ$sampled == 0, 0, univ$PSU_MOS/univ$n)
  univ$weight <- univ$weight_Ist * univ$weight_IIst
  univ <- univ[order(univ$STRATUM,univ$SUB),]
  universe_PSU <- univ[, c("PSU_ID", "STRATUM", "SUB", "SUBSTRAT", 
                           "AR", "PSU_strat", "ALLOC", "SUBSTRAT_MOS", "PSU_MOS", 
                           "pik", "weight", "sampled", "n")]
  sample_PSU <- universe_PSU[universe_PSU$sampled == 1, ]
  population_stratum <- as.data.frame(tapply(univ$PSU_MOS, 
                                             univ$STRATUM, sum))
  colnames(population_stratum) <- c("STRATUM_MOS")
  population_stratum$STRATUM <- row.names(population_stratum)
  sample_PSU <- merge(sample_PSU, population_stratum)
  sample_substratum <- as.data.frame(tapply(sample_PSU$PSU_MOS, 
                                            sample_PSU$SUBSTRAT, sum))
  colnames(sample_substratum) <- c("SUBSTRAT_MOS2")
  sample_substratum$SUBSTRAT <- row.names(sample_substratum)
  sample_PSU <- merge(sample_PSU, sample_substratum)
  PSU_substratum <- as.data.frame(table(sample_PSU$SUBSTRAT))
  colnames(PSU_substratum) <- c("SUBSTRAT", "PSU_substrat")
  sample_PSU <- merge(sample_PSU, PSU_substratum)
  sample_PSU$ALLOC_SUBSTR <- round(sample_PSU$ALLOC * sample_PSU$SUBSTRAT_MOS/sample_PSU$STRATUM_MOS)
  sample_PSU$PSU_final_sample_unit <- round(sample_PSU$ALLOC_SUBSTR/sample_PSU$PSU_substrat)
  sample_PSU$PSU_final_sample_unit <- ifelse(sample_PSU$PSU_final_sample_unit < 
                                               minimum, minimum, sample_PSU$PSU_final_sample_unit)
  sample_PSU$SR <- ifelse(sample_PSU$AR == 1, 1, 0)
  sample_PSU$nSR <- ifelse(sample_PSU$AR == 0, 1, 0)
  sample_PSU$stratum <- sample_PSU$SUBSTRAT
  sample_PSU$Pik <- sample_PSU$pik
  sample_PSU$weight_1st <- 1/sample_PSU$Pik
  sample_PSU$weight_2st <- sample_PSU$PSU_MOS/sample_PSU$PSU_final_sample_unit
  sample_PSU$weight <- sample_PSU$weight_1st * sample_PSU$weight_2st
  sample_PSU <- sample_PSU[, c("PSU_ID", "STRATUM", "stratum", 
                               "SR", "nSR", "PSU_final_sample_unit", "Pik", "weight_1st", 
                               "weight_2st", "weight")]
  PSU_stats <- as.data.frame(table(sample_PSU$STRATUM))
  colnames(PSU_stats) <- c("STRATUM", "PSU")
  PSU_SR <- as.data.frame(table(sample_PSU$STRATUM[sample_PSU$SR == 1]))
  colnames(PSU_SR) <- c("STRATUM", "PSU_SR")
  PSU_stats <- merge(PSU_stats, PSU_SR, all.x=TRUE)
  PSU_stats$PSU_SR <- ifelse(is.na(PSU_stats$PSU_SR),0,PSU_stats$PSU_SR)
  PSU_stats$PSU_NSR <- PSU_stats$PSU - PSU_stats$PSU_SR
  SSU <- aggregate(PSU_final_sample_unit ~ STRATUM, data = sample_PSU, 
                   sum)
  colnames(SSU)[2] <- "SSU"
  SSU_SR <- aggregate(PSU_final_sample_unit ~ STRATUM, data = sample_PSU[sample_PSU$SR == 1, ], sum)
  colnames(SSU_SR)[2] <- "SSU_SR"
  SSU_NSR <- aggregate(PSU_final_sample_unit ~ STRATUM, data = sample_PSU[sample_PSU$nSR == 
                                                                            1, ], sum)
  colnames(SSU_NSR)[2] <- "SSU_NSR"
  PSU_stats <- merge(PSU_stats, SSU, all.x = T)
  PSU_stats <- merge(PSU_stats, SSU_SR, all.x = T)
  PSU_stats <- merge(PSU_stats, SSU_NSR, all.x = T)
  PSU_stats$SSU_SR <- ifelse(is.na(PSU_stats$SSU_SR), 0, PSU_stats$SSU_SR)
  PSU_stats$SSU_NSR <- ifelse(is.na(PSU_stats$SSU_NSR), 0, 
                              PSU_stats$SSU_NSR)
  PSU_stats <- PSU_stats[order(as.numeric(as.character(PSU_stats$STRATUM))), 
  ]
  PSU_stats <- rbind(PSU_stats, rep(NA, 8))
  PSU_stats$STRATUM <- as.character(PSU_stats$STRATUM)
  PSU_stats[nrow(PSU_stats), 1] <- "Total"
  PSU_stats[nrow(PSU_stats), c(2:7)] <- colSums(PSU_stats[-nrow(PSU_stats),2:7])
  # PSU_stats <- PSU_stats[order(as.numeric(PSU_stats$STRATUM)),]
  if (plot == TRUE) {
    des <- PSU_stats[-nrow(PSU_stats), ]
    des2 <- NULL
    des2$STRATUM <- c(des$STRATUM[1:(nrow(PSU_stats) - 1)], 
                      des$STRATUM[1:(nrow(PSU_stats) - 1)])
    des2$SR <- c(rep("SR", (nrow(PSU_stats) - 1)), rep("nSR", 
                                                       (nrow(PSU_stats) - 1)))
    des2$PSU <- rep(NA, (nrow(PSU_stats) - 1))
    des2$SSU <- rep(NA, (nrow(PSU_stats) - 1))
    des2 <- as.data.frame(des2)
    des2$PSU[c(1:(nrow(PSU_stats) - 1))] <- des$PSU_SR
    des2$PSU[c(nrow(PSU_stats):nrow(des2))] <- des$PSU_NSR
    des2$SSU[c(1:(nrow(PSU_stats) - 1))] <- des$SSU_SR
    des2$SSU[c(nrow(PSU_stats):nrow(des2))] <- des$SSU_NSR
    des2$STRATUM <- as.numeric(des2$STRATUM)
    des2
    par(mfrow = c(2, 1))
    barplot(PSU ~ SR + STRATUM, data = des2, main = "PSUs by strata", 
            xlab = "strata", ylab = "PSUs", col = c("black", 
                                                    "grey"), las = 2, cex.names = 0.7)
    legend("topright", legend = c("Non Self Representative", 
                                  "Self Representative"), cex = 0.7, fill = c("black", 
                                                                              "grey"))
    barplot(SSU ~ SR + STRATUM, data = des2, main = "SSUs by strata", 
            xlab = "strata", ylab = "SSUs", col = c("black", 
                                                    "grey"), las = 2, cex.names = 0.7)
    legend("topright", legend = c("Non Self Representative", 
                                  "Self Representative"), cex = 0.7, fill = c("black", 
                                                                              "grey"))
  }
  out <- list(universe_PSU = universe_PSU, sample_PSU = sample_PSU, 
              PSU_stats = PSU_stats[, c("STRATUM", "PSU", "PSU_SR", 
                                        "PSU_NSR", "SSU", "SSU_SR", "SSU_NSR")])
  out
}
