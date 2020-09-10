# Function taken from FS4 package 
# (with changes in lines 427-445)
# Author: Raffaella Cianchetta

StratSel <- function (dataPop, idpsu, dom, final_pop, size, PSUsamplestratum, 
                      min_sample, min_sample_index = FALSE, dataAll, domAll, f_sample, 
                      planned_min_sample = NULL, launch = TRUE) 
{
  if (!inherits(dataPop, "data.frame")) 
    stop("Population data must be supplied as a data frame")
  if (!inherits(idpsu, "formula")) 
    stop("PSU identifier must be supplied as a formula")
  idpsu <- all.vars(idpsu)
  if (length(idpsu) > 1) 
    stop("Idpsu formula must reference only one variable")
  if (!inherits(dom, "formula")) 
    stop("Domain identifier must be supplied as a formula")
  dom <- all.vars(dom)
  if (length(dom) > 1) 
    stop("Population domain formula must reference only one variable")
  if (!inherits(final_pop, "formula")) 
    stop("Final population identifier must be supplied as a formula")
  final_pop <- all.vars(final_pop)
  if (length(final_pop) > 1) 
    stop("Final population formula must reference only one variable")
  if (!inherits(size, "formula")) 
    stop("Size identifier must be supplied as a formula")
  size <- all.vars(size)
  if (length(size) > 1) 
    stop("Size formula must reference only one variable")
  if (!is.numeric(PSUsamplestratum)) 
    stop("PSUsamplestratum variable must be numeric")
  if (PSUsamplestratum < 1) 
    stop("PSUsamplestratum must be >= 1")
  if (!is.logical(min_sample_index)) 
    stop("Parameter 'min_sample_index' must be logical")
  if (min_sample_index == TRUE) {
    if (!is.null(min_sample)) 
      stop("Min_sample must be supplied as NULL")
    if (is.null(planned_min_sample)) 
      stop("Planned_min_sample must be supplied not NULL, but as a formula")
    if (!inherits(planned_min_sample, "formula")) 
      stop("Cluster identifiers must be supplied as a formula")
    planned_min_sample <- all.vars(planned_min_sample)
    if (length(planned_min_sample) > 1) 
      stop("planned_min_sample formula must reference only one variable")
  }
  else {
    if (!is.numeric(min_sample)) 
      stop("min_sample variable must be numeric")
    if (min_sample < 1) 
      stop("min_sample must be >= 1")
    if (!is.null(planned_min_sample)) 
      stop("Planned_min_sample must be supplied as NULL")
  }
  if (!inherits(dataAll, "data.frame")) 
    stop("Allocation data must be supplied as a data frame")
  if (!inherits(domAll, "formula")) 
    stop("Cluster identifiers must be supplied as a formula")
  domAll <- all.vars(domAll)
  if (length(domAll) > 1) 
    stop("Allocation domain formula must reference only one variable")
  if (!inherits(f_sample, "formula")) 
    stop("Cluster identifiers must be supplied as a formula")
  f_sample <- all.vars(f_sample)
  if (length(f_sample) > 1) 
    stop("Final sample formula must reference only one variable")
  if (!is.logical(launch)) 
    stop("Parameter 'launch' must be logical")
  Domain <- as.factor(dataPop[[dom]])
  population_file <- data.frame(dataPop[idpsu], Domain, dataPop[final_pop], 
                                dataPop[size])
  nrow_pop_f <- nrow(population_file)
  domAll <- as.factor(dataAll[[domAll]])
  if (min_sample_index == TRUE) {
    allocation_file <- data.frame(domAll, dataAll[f_sample], 
                                  dataAll[planned_min_sample])
  }
  else {
    allocation_file <- data.frame(domAll, dataAll[f_sample])
  }
  merge_output <- merge(x = population_file, y = allocation_file, 
                        by.x = names(population_file)[2], by.y = names(allocation_file)[1])
  nrow_merge_out <- nrow(merge_output)
  if (nrow_merge_out != nrow_pop_f) {
    stop("Merge result is inconsistent", call. = FALSE)
  }
  else {
    ndom <- nlevels(merge_output[["Domain"]])
    domains <- split(merge_output, merge_output[["Domain"]])
    Domini <- levels((merge_output)[["Domain"]])
    sampling_fraction <- rep(NA, ndom)
    final_sample <- rep(NA, ndom)
    choicedom <- rep(NA, ndom)
    SRdom <- rep(NA, ndom)
    sizeSRdom <- rep(NA, ndom)
    threshold <- rep(NA, ndom)
    nSRdom <- rep(NA, ndom)
    sizenSRdom <- rep(NA, ndom)
    SRdom_nSRdom <- rep(NA, ndom)
    if (min_sample_index == TRUE) {
      planned_min_sample <- rep(NA, ndom)
    }
    sizedom <- tapply(merge_output[[4]], merge_output[[1]], 
                      FUN = sum)
    final_populationdom <- tapply(merge_output[[3]], merge_output[[1]], 
                                  FUN = sum)
    for (i in 1:ndom) {
      final_sample[i] <- domains[[i]][[5]][1]
      if (min_sample_index == TRUE) {
        planned_min_sample[i] <- domains[[i]][[6]][1]
        threshold[i] <- (planned_min_sample[i] * sizedom[i])/final_sample[i]
      }
      else {
        threshold[i] <- (min_sample * sizedom[i])/final_sample[i]
      }
      sampling_fraction[i] <- final_sample[i]/final_populationdom[[i]]
      cSRdom <- 0
      csizeSRdom <- 0
      for (j in 1:nrow(domains[[i]])) {
        domains[[i]]$threshold[j] <- threshold[i]
        domains[[i]]$final_populationdom[j] <- final_populationdom[i]
        domains[[i]]$sampling_fraction[j] <- sampling_fraction[i]
        SR <- 0
        sizeSR <- 0
        domains[[i]]$SR[j] <- SR
        domains[[i]]$SizeSR[j] <- sizeSR
        if (domains[[i]][[size]][j] >= threshold[i]) {
          domains[[i]]$SR[j] <- 1
          cSRdom <- cSRdom + domains[[i]]$SR[j]
          domains[[i]]$SizeSR[j] <- domains[[i]][[size]][j] * 
            domains[[i]]$SR[j]
          csizeSRdom <- csizeSRdom + domains[[i]]$SizeSR[j]
        }
        SRdom[i] <- cSRdom
        sizeSRdom[i] <- csizeSRdom
        nSRdom[i] <- (round(((sizedom[i] - sizeSRdom[i])/threshold[i])/PSUsamplestratum)) * 
          PSUsamplestratum
        sizenSRdom[i] <- sizedom[i] - sizeSRdom[i]
        if (sizenSRdom[i] > 0 && nSRdom[i] == 0) {
          nSRdom[i] <- 1
        }
        SRdom_nSRdom[i] <- SRdom[i] + nSRdom[i]
      }
      df_dominio <- data.frame(sizedom, threshold, SRdom, 
                               sizeSRdom, sizenSRdom, nSRdom, SRdom_nSRdom)
      df_Domainresults <- data.frame(Domini, SRdom, nSRdom, 
                                     SRdom_nSRdom)
      names(df_Domainresults)[4] <- "SRdom+nSRdom"
    }
    names(df_Domainresults)[1] <- dom
    df_Domainresults[[1]] <- as.character(df_Domainresults[[1]])
    df_Domainresults <- rbind(df_Domainresults, c("Total", 
                                                  apply(df_Domainresults[, -1], 2, sum)))
    class(df_Domainresults) <- c("data.frame", "results")
    if (launch) {
      return(df_Domainresults)
    }
    else {
      domain_base <- NULL
      domain_base <- ldply(domains, data.frame)
      domain_base$.id <- NULL
      domain_base_ord <- domain_base[order(domain_base[[1]], 
                                           -domain_base$SR, -domain_base[[4]], decreasing = FALSE), 
      ]
      domains <- split(domain_base_ord, domain_base_ord[["Domain"]])
      ndom <- nlevels(domain_base_ord[["Domain"]])
      overflow <- 0
      remainder <- 0
      for (i in 1:ndom) {
        csizedom <- 0
        cPSUdom <- 0
        cth <- 0
        csizestr <- 0
        stratum <- 0
        first_time <- 0
        NPSUdomsample <- 0
        for (j in 1:nrow(domains[[i]])) {
          if (domains[[i]]$SR[j] == 1) {
            csizedom <- csizedom + domains[[i]][[size]][j]
            cPSUdom <- cPSUdom + 1
            stratum <- stratum + 1
            domains[[i]]$stratum[j] <- paste(domains[[i]]$Domain[j], 
                                             stratum, sep = "")
            domains[[i]]$csizestr[j] <- csizestr
            domains[[i]]$csizedom[j] <- csizedom
            domains[[i]]$cth[j] <- cth
          }
          else {
            NPSUdomsample <- df_dominio$SRdom_nSRdom[i]
            if ((NPSUdomsample - cPSUdom) <= 0) {
              domains[[i]]$stratum[j] <- paste(domains[[i]]$Domain[j], 
                                               stratum, sep = "")
              domains[[i]]$csizestr[j] <- csizestr
              domains[[i]]$csizedom[j] <- csizedom
              domains[[i]]$cth[j] <- cth
              next
            }
            else {
              if ((domains[[i]]$SR[j] == 0) && (first_time == 
                                                0)) {
                stratum <- stratum + 1
                sizedom <- df_dominio$sizedom[i]
                cth <- (((sizedom - csizedom)/(NPSUdomsample - 
                                                 cPSUdom)) * PSUsamplestratum)
                first_time <- 1
              }
              domains[[i]]$stratum[j] <- paste(domains[[i]]$Domain[j], 
                                               stratum, sep = "")
              domains[[i]]$csizestr[j] <- csizestr
              domains[[i]]$csizedom[j] <- csizedom
              domains[[i]]$cth[j] <- cth
              csizestr <- csizestr + domains[[i]][[size]][j]
              if ((csizestr <= cth) || ((NPSUdomsample - 
                                         cPSUdom) == PSUsamplestratum)) {
                domains[[i]]$stratum[j] <- paste(domains[[i]]$Domain[j], 
                                                 stratum, sep = "")
                domains[[i]]$csizestr[j] <- csizestr
                domains[[i]]$csizedom[j] <- csizedom
                domains[[i]]$cth[j] <- cth
                next
              }
              if ((csizestr > cth) && ((NPSUdomsample - 
                                        cPSUdom) > PSUsamplestratum)) {
                remainder <- cth - (csizestr - domains[[i]][[size]][j])
                overflow <- csizestr - cth
                if (remainder >= overflow) {
                  domains[[i]]$stratum[j] <- paste(domains[[i]]$Domain[j], 
                                                   stratum, sep = "")
                  domains[[i]]$csizestr[j] <- csizestr
                  domains[[i]]$csizedom[j] <- csizedom
                  domains[[i]]$cth[j] <- cth
                  csizedom <- csizedom + csizestr
                  cPSUdom <- cPSUdom + PSUsamplestratum
                  csizestr <- 0
                  sizedom <- df_dominio$sizedom[i]
                  cth <- (((sizedom - csizedom)/(NPSUdomsample - 
                                                   cPSUdom)) * PSUsamplestratum)
                  stratum <- stratum + 1
                }
                else {
                  csizedom <- csizedom + (csizestr - 
                                            domains[[i]][[size]][j])
                  cPSUdom <- cPSUdom + PSUsamplestratum
                  sizedom <- df_dominio$sizedom[i]
                  cth <- (((sizedom - csizedom)/(NPSUdomsample - 
                                                   cPSUdom)) * PSUsamplestratum)
                  csizestr <- domains[[i]][[size]][j]
                  stratum <- stratum + 1
                  domains[[i]]$cth[j] <- cth
                  domains[[i]]$stratum[j] <- paste(domains[[i]]$Domain[j], 
                                                   stratum, sep = "")
                  domains[[i]]$csizedom[j] <- csizedom
                  domains[[i]]$csizestr[j] <- csizestr
                }
              }
            }
          }
        }
      }
      DF_stratum <- NULL
      DF_stratum <- ldply(domains, data.frame)
      DF_stratum$.id <- NULL
      DF_stratum$cth <- NULL
      DF_stratum$csizedom <- NULL
      DF_stratum$csizestr <- NULL
      if (min_sample_index == TRUE) {
        DF_stratum[6] <- NULL
      }
      piece <- NULL
      Countstr <- ddply(DF_stratum, .(DF_stratum$Domain, 
                                      DF_stratum$stratum), transform, N_PSU_Stratum = NROW(piece))
      columnname <- names(DF_stratum)[4]
      names(DF_stratum)[4] <- "size"
      SomSizePSUstr <- ddply(DF_stratum, .(DF_stratum$Domain, 
                                           DF_stratum$stratum), transform, somma = sum(size))
      names(SomSizePSUstr)[6] <- columnname
      columnname <- names(DF_stratum)[3]
      final_population <- NULL
      names(DF_stratum)[3] <- "final_population"
      SumFinalPopulationstr <- ddply(DF_stratum, .(DF_stratum$Domain, 
                                                   DF_stratum$stratum), transform, sumfinalPopstr = sum(final_population))
      names(SumFinalPopulationstr)[5] <- columnname
      Countstr <- cbind(SomSizePSUstr$somma, SumFinalPopulationstr$sumfinalPopstr, 
                        Countstr)
      PSUsamplestratum <- rep(PSUsamplestratum, nrow(Countstr))
      Countstr <- cbind(Countstr, PSUsamplestratum)
      Countstr[4] <- NULL
      Countstr[3] <- NULL
      CountstrE <- Countstr[Countstr$N_PSU_Stratum > Countstr$PSUsamplestratum[1], 
      ]
      CountstrNEs <- Countstr[Countstr$N_PSU_Stratum <= 
                                Countstr$PSUsamplestratum[1], ]
      numrowNEs <- nrow(CountstrNEs)
      if (numrowNEs != 0) {
        vUps <- rep(1, numrowNEs)
        CountstrNEs$SR <- 1
        vPIK <- rep(1, numrowNEs)
        DFAR <- cbind(vPIK, CountstrNEs)
        names(DFAR)[1] <- "PIK"
        DFAR <- cbind(vUps, DFAR)
        names(DFAR)[4] <- "SumFinalPopulationstr"
        names(DFAR)[3] <- "Size_Stratum"
      }
      else {
        DFAR <- NULL
      }
      PIK <- (CountstrE[[6]]/CountstrE[[1]]) * CountstrE$PSUsamplestratum[1]
      CountstrE <- cbind(PIK, CountstrE)
      Countstr_Over_1 <- CountstrE[CountstrE$PIK > 1, 
      ]
      end_PIK <- FALSE
      CountstrSTRATI_Add <- NULL
      while (end_PIK == FALSE) {
        if (nrow(Countstr_Over_1) == 0) {
          ListCountstrE <- split(CountstrE, CountstrE[[4]], 
                                 drop = TRUE)
          ndomList <- length(ListCountstrE)
          newStratum <- 0
          first_time <- 0
          vPIK <- c()
          vUps <- c()
          for (i in 1:ndomList) {
            newStratum <- 0
            x <- 1
            vPIK <- c()
            UPs <- c()
            for (j in 1:nrow(ListCountstrE[[i]])) {
              if (ListCountstrE[[i]]$stratum[j] != newStratum) {
                if (newStratum > 0) {
                  Pik <- vPIK
                  UPs = UPsampford(Pik)
                  vPIK <- c()
                  x <- 1
                  vUps <- c(vUps, UPs)
                }
                newStratum <- ListCountstrE[[i]]$stratum[j]
                vPIK[x] <- ListCountstrE[[i]][[1]][j]
                x <- x + 1
              }
              else {
                newStratum <- ListCountstrE[[i]]$stratum[j]
                vPIK[x] <- ListCountstrE[[i]][[1]][j]
                x <- x + 1
              }
            }
            Pik <- vPIK
            UPs = UPsampford(Pik)
            vUps <- c(vUps, UPs)
          }
          DFUps <- NULL
          DFUps <- ldply(ListCountstrE, data.frame)
          DFUps <- cbind(vUps, DFUps)
          DFUps[2] <- NULL
          names(DFUps)[4] <- "SumFinalPopulationstr"
          names(DFUps)[3] <- "Size_Stratum"
          end_PIK <- TRUE
        }
        else {
          Countstr_Over_1 <- ddply(Countstr_Over_1, 
                                   .(Countstr_Over_1$Domain, Countstr_Over_1$stratum), 
                                   transform, N_Stratum = NROW(piece))
          CountstrE <- CountstrE[CountstrE$PIK < 1, 
          ]
          CountstrSTRATI_Add <- NULL
          for (i in 1:nrow(Countstr_Over_1)) {
            domain <- 0
            stratum <- 0
            CountstrSTRATI <- NULL
            domain <- Countstr_Over_1$Domain[i]
            stratum <- Countstr_Over_1$stratum[i]
            CountstrSTRATI <- CountstrE[((CountstrE$Domain == 
                                            domain) & (CountstrE$stratum == stratum)), 
            ]
            if ((nrow(CountstrSTRATI)) != 0) {
              CountstrSTRATI$PSUsamplestratum <- CountstrSTRATI$PSUsamplestratum - 
                Countstr_Over_1$N_Stratum[i]
              vSomSizePSUstr <- 0
              vSomSizePSUstr <- sum(CountstrSTRATI[[7]])
              CountstrSTRATI[[2]] <- vSomSizePSUstr
              CountstrE <- CountstrE[!((CountstrE$Domain == 
                                          domain) & (CountstrE$stratum == stratum)), 
              ]
              CountstrSTRATI_Add <- rbind(CountstrSTRATI_Add, 
                                          CountstrSTRATI)
            }
            else {
              next
            }
          }
          CountstrE <- rbind(CountstrE, CountstrSTRATI_Add)
          CountstrE <- CountstrE[order(CountstrE$Domain, 
                                       -CountstrE$SR, -CountstrE[[7]], decreasing = FALSE), 
          ]
          CountstrE$PIK <- (CountstrE[[7]]/CountstrE[[2]]) * 
            CountstrE$PSUsamplestratum
          Countstr_Over_1$N_Stratum <- NULL
          Countstr_Over_1[3] <- 1
          Countstr_Over_1[2] <- NULL
          Countstr_Over_1[1] <- NULL
          vUps <- rep(1, nrow(Countstr_Over_1))
          Countstr_Over_1$SR <- 1
          Countstr_Over_1 <- cbind(vUps, Countstr_Over_1)
          names(Countstr_Over_1)[4] <- "SumFinalPopulationstr"
          names(Countstr_Over_1)[3] <- "Size_Stratum"
          DFAR <- rbind(DFAR, Countstr_Over_1)
          Countstr_Over_1 <- NULL
          Countstr_Over_1 <- CountstrE[CountstrE$PIK > 
                                         1, ]
        }
      }
      if (length(DFAR) == 0) {
        DFGlobal <- rbind(DFUps)
      }
      else {
        DFGlobal <- rbind(DFAR, DFUps)
      }
      DFGlobal_ord <- DFGlobal[order(DFGlobal[[5]], -DFGlobal$SR, 
                                     -DFGlobal$SizeSR, decreasing = FALSE), ]
      drops <- c("DF_stratum$Domain", "DF_stratum$stratum")
      DFGlobal_ord <- DFGlobal_ord[, !(names(DFGlobal_ord) %in% 
                                         drops)]
      colnames(DFGlobal_ord)[1:3] <- c("Sampled_PSU", 
                                       "Pik", "Size_Stratum")
      PSU_final_sample_unit <- rep(NA, nrow(DFGlobal_ord))
      nSR <- rep(NA, nrow(DFGlobal_ord))
      # From here: changes to let this function work also with R 4.0 -----------------------
      # change [[i]] in [i]
      for (i in 1:nrow(DFGlobal_ord)) {
        if (DFGlobal_ord$N_PSU_Stratum[i] <= DFGlobal_ord$PSUsamplestratum[i]) 
          DFGlobal_ord$PSU_final_sample_unit[i] <- round((DFGlobal_ord$sampling_fraction[i] * 
                                                              DFGlobal_ord[[4]][i] * DFGlobal_ord$Sampled_PSU[i])/DFGlobal_ord$N_PSU_Stratum[i])
        else {
          DFGlobal_ord$PSU_final_sample_unit[i] <- round((DFGlobal_ord$sampling_fraction[i] * 
                                                              DFGlobal_ord[[4]][i] * DFGlobal_ord$Sampled_PSU[i])/Countstr$PSUsamplestr[1])
        }
        if ((DFGlobal_ord$Sampled_PSU[i] == 1) && 
            (DFGlobal_ord$Pik[i] < 1)) {
          DFGlobal_ord$nSR[i] <- 1
        }
        else {
          DFGlobal_ord$nSR[i] <- 0
        }
      }
      # ------ until here ---------------------------------------------------------------
      DFGlobal_ord[4] <- NULL
      a <- aggregate(DFGlobal_ord$SR, list(Domain = DFGlobal_ord$Domain), 
                     FUN = sum)
      b <- aggregate(DFGlobal_ord$nSR, list(Domain = DFGlobal_ord$Domain), 
                     FUN = sum)
      names(a)[2] <- "SRdom"
      names(b)[2] <- "nSRdom"
      b <- cbind(a, b)
      b[3] <- NULL
      for (i in 1:ndom) {
        b$SRdom_nSRdom[i] <- b$SRdom[i] + b$nSRdom[i]
      }
      aggnSR <- aggregate(DFGlobal_ord$PSU_final_sample_unit, 
                          list(Domain = DFGlobal_ord$Domain, nSR = DFGlobal_ord$nSR), 
                          FUN = sum)
      AGGnSR <- aggnSR[aggnSR$nSR == 1, ]
      AGGSR <- aggnSR[aggnSR$nSR == 0, ]
      MERGEAGG <- merge(AGGSR, AGGnSR, by = "Domain", 
                        all = TRUE, sort = TRUE)
      MERGEAGG[is.na(MERGEAGG)] <- 0
      MERGEAGG[4] <- NULL
      MERGEAGG[2] <- NULL
      MERGEAGG[1] <- NULL
      names(MERGEAGG)[1] <- "SR_PSU_final_sample_unit"
      names(MERGEAGG)[2] <- "NSR_PSU_final_sample_unit"
      b <- cbind(b, MERGEAGG)
      b[[1]] <- as.character(b[[1]])
      Total <- c("Total", round(apply(b[, -1], 2, sum)))
      b[[2]] <- as.character(b[[2]])
      b[[3]] <- as.character(b[[3]])
      b[[4]] <- as.character(b[[4]])
      names(b)[4] <- "SRdom+nSRdom"
      Mean <- c("Mean", " ", " ", " ", round(apply(b[-(1:4)], 
                                                   2, mean)))
      b <- rbind(b, Total, Mean)
      aggPSUsampleStr <- aggregate(DFGlobal_ord$Sampled_PSU, 
                                   list(Dominio = DFGlobal_ord$Domain, Stratum = as.numeric(DFGlobal_ord$stratum)), 
                                   FUN = sum)
      aggPSUsampleStr <- aggPSUsampleStr[order(aggPSUsampleStr$Dominio), 
      ]
      not_dupl_N_PSU_Stratum <- DFGlobal_ord[!duplicated(DFGlobal_ord$stratum, 
                                                         DFGlobal_ord$N_PSU_Stratum), ]
      aggPSUsampleStr <- cbind(aggPSUsampleStr, not_dupl_N_PSU_Stratum$N_PSU_Stratum)
      aggPSUsampleStr[1] <- NULL
      names(aggPSUsampleStr)[3] <- "N_PSU_Stratum"
      names(aggPSUsampleStr)[2] <- "PSU_Sample_Stratum"
      class(b) <- c("data.frame", "results")
      class(aggPSUsampleStr) <- c("data.frame", "results")
      names(DFGlobal_ord)[4] <- dom
      DFGlobal_ord$PSUsamplestratum <- NULL
      class(DFGlobal_ord) <- c("data.frame", "results")
      output <- list(df_Domainresults, b, aggPSUsampleStr, 
                     DFGlobal_ord)
      return(output)
      gc()
    }
  }
}
