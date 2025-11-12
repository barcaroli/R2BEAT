#' Compute multivariate optimal allocation for different domains in two-stage stratified sample design
#'
#' @description 
#' This function computes the multivariate optimal allocation for different domains 
#' in a one-stage file_strataied sample design, following the generalized Bethel algorithm.
#' 
#' @usage
#' beat.1st(file_strata,errors,des_file,psu_file,rho,deft_start,effst,epsilon1,mmdiff_deft,maxi,epsilon,minPSUstrat,minnumstrat,maxiter,maxiter1)
#'
#' @param file_strata: dataframe of survey strata, for more details see, e.g.,\code{\link{strata}}. 
#' @param errors: dataframe of coefficients of variation for each domain, for more details see, e.g.,\code{\link {errors}}.
#' @param des_file: dataframe containing information on sampling design variables, for more details see, e.g.,\code{\link {design}}.
#' @param psu_file: dataframe containing information on primary stage units stratification, for more details see, e.g.,\code{\link {PSU_strat}}.
#' @param rho: dataframe of survey strata, for more details see, e.g.,\code{\link rho}.
#' @param deft_start: dataframe of survey strata, for taking into account the initial design effect on each variable, for more details see, e.g.,\code{\link {deft_start}}.
#' @param effst: dataframe of survey strata, for taking into account the estimator effect on each variable, for more details see, e.g.,\code{\link {effst}}.
#' @param epsilon1: first stop condition: sample sizes differences beetween two iterations; iteration continues until the maximum of sample sizes differences is greater than the default value. The default is 5. 
#' @param mmdiff_deft: second stop condition: defts differences beetween two iterations; iteration continues until the maximum of defts largest differences is greater than the default value. The default is 0.06.
#' @param maxi: third stop condition: maximum number of allowed iterations. The default is 20.
#' @param epsilon: the same as in function \code{\link {beat.1st}}.
#' @param minPSUstrat: minimum number of non-self-represenative PSUs to be selected in each stratum.
#' @param minnumstrat: the same as in function \code{\link {beat.1st}}.
#' @param maxiter: the same as in function \code{\link {beat.1st}}.
#' @param maxiter1: the same as in function \code{\link{beat.1st}}.
#' 
#' @return
#' Returns a list with four components:
#' 
#' - iteractions: dataframe that for each iteraction provides a summary with the number of Primary Stage Units (\code{PSU_Total}) distinguish between Self-Representative (\code{PSU_SR}) from Non-Self-Representative (\code{PSU_NSR}) and the number of Secondary Stage Units (\code{SSU}). This output is also printed to the screen.
#' 
#' - file_strata: Input dataframe in \code{file_strata} with the design effect for each variables in each stratum (\code{DEFT1 - DEFTn}) and the optimal sample size columns.
#' 
#' - alloc: dataframe with optimal (\code{OPT}), proportional (\code{PROP}), equal (\code{UNIF}) sample size allocation.
#' 
#' - planned: Data frame with a summary of expected coefficients of variation for each variable in each domain.
#' 
#' - expected: Data frame with a summary of realized coefficients of variation with the given optimal allocation for each variable in each domain.
#' 
#' - sensitivity: Data frame with a summary of the sensitivity at 10\% for each domain and each variable. Sensitivity can be a useful tool to help in finding the best allocation, because it provides a hint of the expected sample size variation for a 10\% change in planned CVs.
#' 
#' - deft_c: Data frame with the design effect for each variable in each domain in each iteraction. Note that \code{DEFT1_0 - DEFTn_0} is always equal to 1 if \code{deft_start} is \code{NULL}. Instead is equal to \code{deft_start}. While \code{DEFT1 - DEFTn} are the final design effect related to the given allocation.
#' 
#' - param_alloc: A vector with a resume of all the parameter given for the allocation.
#'
#' @examples
#'
#'# Load example data
#'data(beat.example)
#'
## Example 1
# Allocate the sample
#'allocation2st_1 <- beat.2st(file_strata=strata, errors=errors,
#'                            des_file=design, psu_file=PSU_strat,rho=rho)
# The total ammount of sample size is 191 PSU (36 SR + 155 NSR) and 15147 SSU. 
#'
## Example 2
# Assume 13000 SSUs is the maximum sample size to stick to our budget.
# Look at the sensitivity is in DOM1 for REG1 and REG2 due to V1.
#'allocation2st_1$sensitivity
# We can relax the constraints increasing the expected coefficients of variation for X1 by 10%
#'errors1 <- errors 
#'errors1[1,2] <- errors[1,2]+errors[1,2]*0.1
# Try the new allocation 
#'allocation2st_2 <- beat.2st(file_strata=strata, errors=errors1,
#'                            des_file=design, psu_file=PSU_strat,rho=rho)
#'
## Example 3
# On the contrary, if we tighten the constraints decreasing the expected coefficients of variation 
# for X1 by 10%
#'errors2 <- errors 
#'errors2[1,2] <- errors[1,2]-errors[1,2]*0.1
# The new allocation leads to a larger sample than the first example (around 18000)
#'allocation2st_3 <- beat.2st(file_strata=strata, errors=errors2,
#'                           des_file=design, psu_file=PSU_strat,rho=rho)

## Example 4
# Sometimes some budget constraints concern the number of PSU involved in the survey.
# Tuning the PSUs number is possible modyfing the MINIMUM in des_file. 
# Assume to increase the MINIMUM from 48 to 60
#'design1 <- design 
#'design1[,4] <- 60
#'allocation2st_4 <- beat.2st(file_strata=strata, errors=errors2,
#'                            des_file=design1, psu_file=PSU_strat, rho=rho)

# The PSUs number is decreased, while the SSUs number increased 
# due to cluster intra-correlation effect.
# Under the same expected errors, to offset a slight reduction of PSUs (from 221 to 207) 
# an increase of SSUs involved is observed.
#'allocation2st_3$expected
#'allocation2st_4$expected

## Example 5
# On the contrary, assume to decrease the MINIMUM from 48 to 24.
# The SSUs number strongly decrease in the face of an increase of PSUs,
# always under the same expected errors.
#'design2 <- design 
#'design2[,4] <- 24
#'allocation2st_5 <- beat.2st(file_strata=strata, errors=errors2,
#'                           des_file=design2, psu_file=PSU_strat, rho=rho)
#'allocation2st_4$expected
#'allocation2st_5$expected

## Example 6
# Assume that the SSUs are in turn clusters, for instance households composed by individuals.
# In the previous examples we always derived optimal allocations
# for sample of SSUs (i.e. households, because
# DELTA = 1).
#'design
#'design1
#'design2
# For obtaining a sample in terms of the elements composing SSUs
# (i.e., individuals) is just sufficient to 
# modify the DELTA in des_file.
#'design3 <- design 
#'design3$DELTA <- 2.31
# DELTA_IND=2.31, the average size of household in Italy.
#'allocation2st_6 <- beat.2st(file_strata=strata, errors=errors,
#'                            des_file=design3, psu_file=PSU_strat, rho=rho)

## Example 7
# Complete workflow
#'library(readr)
#'pop <- read_rds("https://github.com/barcaroli/R2BEAT_workflows/blob/master/pop.RDS?raw=true")
#'library(R2BEAT)
#'cv <- as.data.frame(list(DOM=c("DOM1","DOM2"),
#'                         CV1=c(0.02,0.03),
#'                         CV2=c(0.03,0.06),
#'                         CV3=c(0.03,0.06),
#'                         CV4=c(0.05,0.08)))
#'cv
#'samp_frame <- pop
#'samp_frame$one <- 1
#'id_PSU <- "municipality"  
#'id_SSU <- "id_ind"        
#'strata_var <- "stratum"   
#'target_vars <- c("income_hh","active","inactive","unemployed")   
#'deff_var <- "stratum"     
#'domain_var <- "region"  
#'delta =  1       # households = survey units
#'minimum <- 50    # minimum number of SSUs to be interviewed in each selected PSU
#'deff_sugg <- 1.5 # suggestion for the deff value
#'
#'inp <- prepareInputToAllocation1(samp_frame,
#'                                id_PSU,
#'                                 id_SSU,
#'                                 strata_var,
#'                                 target_vars,
#'                                 deff_var,
#'                                 domain_var,
#'                                 minimum,
#'                                 delta,
#'                                 deff_sugg)
#'inp$desfile$MINIMUM <- 50
#'alloc <- beat.2st(file_strata = inp$strata, 
#'                  errors = cv, 
#'                  des_file = inp$des_file, 
#'                  psu_file = inp$psu_file, 
#'                  rho = inp$rho, 
#'                  deft_start = NULL,
#'                  effst = inp$effst, 
#'                  minPSUstrat = 2,
#'                  minnumstrat = 50)
#'                  
#' @export
#' 
beat.2st <- function (file_strata, errors, des_file, psu_file, rho, deft_start = NULL, 
                      effst = NULL, epsilon1 = 5, mmdiff_deft = 1, maxi = 20, 
                      epsilon = 10^(-11), minPSUstrat = 2, minnumstrat = 2, maxiter = 200, 
                      maxiter1 = 25) 
{
  stratif=file_strata
  diffx = 999
  iterx = 0
  test_stages <- NULL
  test_stages <<- 2
  param_alloc <- as.data.frame(t(c(epsilon, minPSUstrat, minnumstrat, 
                                   mean(des_file$MINIMUM), maxiter, maxiter1, epsilon1, 
                                   diffx, iterx, mmdiff_deft, maxi, 2)))
  names(param_alloc) = c("p_epsilon", "p_minPSUstrat", "p_minnumstrat", 
                         "minimum", "p_maxiter", "p_maxiter1", "p_epsilon1", 
                         "p_diffx", "p_iterx", "p_mmdiff_deft", "p_maxi", "num_stages")
  nvar = ncol(errors) - 1
  if (!is.null(effst)) {
    colnames(stratif) <- toupper(colnames(stratif))
    sapply(1:nvar, function(i) {
      stratif[, paste("S_ORIG_no_effst", i, sep = "")] <<- stratif[, 
                                                                   paste("S", i, sep = "")]
    })
    sapply(1:nvar, function(i) {
      stratif[, paste("S_ORIG", i, sep = "")] <<- stratif[, 
                                                          paste("S", i, sep = "")] * effst[, paste("EFFST", 
                                                                                                   i, sep = "")]
    })
    sapply(1:nvar, function(i) {
      stratif[, paste("S", i, sep = "")] <<- (stratif[, 
                                                      paste("S", i, sep = "")] * effst[, paste("EFFST", 
                                                                                               i, sep = "")])
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
    deft <- merge(deft, deft_start, by = "STRATUM", all.x = T)
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
  colnames(iterations) <- c("iterations", "PSU_SR", "PSU NSR", 
                            "PSU Total", "SSU")
  des_file <- des_file[order(des_file$STRATUM), ]
  des_file$CAMP <- n
  infoloop <- merge(deft, stratif[, c("STRATUM", sapply(1:nvar, 
                                                        function(i) (paste("S", i, sep = ""))))], by = "STRATUM")
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
      des_file$THRESHOLD = ((des_file$MINIMUM/des_file$F) * 
                              des_file$DELTA)
      psu_file1 <- merge(des_file[, c("STRATUM", "THRESHOLD")], 
                         psu_file, by = "STRATUM")
      psu_file1$AR <- 0
      psu_file1$AR[(psu_file1$PSU_MOS >= psu_file1$THRESHOLD)] <- 1
      table(psu_file1$AR)
      psu_file1 <- psu_file1[order(psu_file1$STRATUM, 
                                   -psu_file1$PSU_MOS), ]
      psu_file1$SUB <- 0
      psu_file1$partial <- 0
      psu_file1$SUB[1] <- 1
      psu_file1$THRESHOLD_NAR <- psu_file1$THRESHOLD * 
        minPSUstrat
      for (i in 2:nrow(psu_file1)) {
        if (psu_file1$STRATUM[i] == psu_file1$STRATUM[i - 
                                                      1] & psu_file1$AR[i] == 1) {
          psu_file1$SUB[i] = psu_file1$SUB[i - 1] + 
            1
        }
        if (psu_file1$STRATUM[i] == psu_file1$STRATUM[i - 
                                                      1] & psu_file1$AR[i] == 0 & psu_file1$AR[i - 
                                                                                               1] == 1) {
          psu_file1$SUB[i] = psu_file1$SUB[i - 1] + 
            1
          psu_file1$partial[i] <- psu_file1$PSU_MOS[i]
        }
        if (psu_file1$STRATUM[i] != psu_file1$STRATUM[i - 
                                                      1] & psu_file1$AR[i] == 0) {
          psu_file1$SUB[i] = 1
          psu_file1$partial[i] <- psu_file1$PSU_MOS[i]
        }
        if (psu_file1$STRATUM[i] == psu_file1$STRATUM[i - 
                                                      1] & psu_file1$AR[i] == 0 & psu_file1$AR[i - 
                                                                                               1] == 0 & psu_file1$partial[i - 1] <= psu_file1$THRESHOLD_NAR[i] & 
            (psu_file1$partial[i - 1] + psu_file1$PSU_MOS[i]) <= 
            psu_file1$THRESHOLD_NAR[i]) {
          psu_file1$SUB[i] = psu_file1$SUB[i - 1]
          psu_file1$partial[i] <- (psu_file1$partial[i - 
                                                       1] + psu_file1$PSU_MOS[i])
        }
        if (psu_file1$STRATUM[i] == psu_file1$STRATUM[i - 
                                                      1] & psu_file1$AR[i] == 0 & psu_file1$AR[i - 
                                                                                               1] == 0 & psu_file1$partial[i - 1] <= psu_file1$THRESHOLD_NAR[i] & 
            (psu_file1$partial[i - 1] + psu_file1$PSU_MOS[i]) > 
            psu_file1$THRESHOLD_NAR[i]) {
          psu_file1$SUB[i] = psu_file1$SUB[i - 1] + 
            1
          psu_file1$partial[i] <- psu_file1$PSU_MOS[i]
        }
      }
      psu_file1$SUBSTRAT <- paste(psu_file1$STRATUM, psu_file1$SUB, 
                                  sep = "-")
      psu_strat <- aggregate(PSU_ID ~ SUBSTRAT, psu_file1, 
                             length)
      colnames(psu_strat)[2] <- "PSU_strat"
      psu_file1 <- merge(psu_file1, psu_strat, by = "SUBSTRAT")
      psu_file1$AR <- ifelse(psu_file1$PSU_strat == minPSUstrat, 
                             1, psu_file1$AR)
      psu_file1 <- psu_file1[order(psu_file1$STRATUM, 
                                   -psu_file1$PSU_MOS), ]
      psu_file1$SUB[1] <- 1
      psu_file1$POPAR = psu_file1$PSU_MOS * psu_file1$AR
      psu_file1$POPNAR = psu_file1$PSU_MOS - psu_file1$POPAR
      popolaz <- aggregate(psu_file1[, c("POPAR", "POPNAR")], 
                           by = list(STRATUM = psu_file1$STRATUM), sum)
      des_file <- merge(des_file, popolaz, by = "STRATUM")
      des_file <- merge(des_file, rho, by = "STRATUM")
      des_file$CAMPAR = ceiling((des_file$F) * (des_file$POPAR))
      des_file$CAMPNAR = des_file$CAMP - des_file$CAMPAR
      des_file$BDIS_AR = 1 * des_file$DELTA
      des_file$BDIS_NAR = des_file$MINIMUM * des_file$DELTA
      sapply(1:nvar, function(i) {
        des_file[, paste("DEFF", i, sep = "")] <<- ((des_file$CAMP/(des_file$STRAT_MOS^2)) * 
                                                      ((((des_file$POPAR^2)/(des_file$CAMPAR + epsilon)) * 
                                                          (1 + ((des_file[, paste("RHO_AR", i, sep = "")]) * 
                                                                  (des_file$BDIS_AR - 1)))) + (((des_file$POPNAR^2)/(des_file$CAMPNAR + 
                                                                                                                       epsilon)) * (1 + ((des_file[, paste("RHO_NAR", 
                                                                                                                                                           i, sep = "")]) * (des_file$BDIS_NAR - 1))))))
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
                                             sapply(1:nvar, function(i) (paste("DEFT", i, 
                                                                               sep = ""))))], by = "STRATUM")
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
                                                                   paste("DEFT", i, sep = "")] - loopdeft[, paste("oldDEFT", 
                                                                                                                  i, sep = "")])
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
                                                          paste("S_ORIG", i, sep = "")] * deft_start[, 
                                                                                                     paste("DEFT", i, sep = "")])
        })
      }
      else {
        sapply(1:nvar, function(i) {
          stratif[, paste("S", i, sep = "")] <<- (stratif[, 
                                                          paste("S_ORIG", i, sep = "")] * des_file[, 
                                                                                                   paste("DEFT", i, sep = "")])
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
                                                                          c * 2))] <- paste("S", 1:nvar, "_", iterx - 
                                                                                              1, sep = "")
      }
      infoloop <- merge(infoloop, des_file[, c("STRATUM", 
                                               sapply(1:nvar, function(i) (paste("DEFT", i, 
                                                                                 sep = ""))))], by = "STRATUM", suffixes = rep("", 
                                                                                                                               2))
      infoloop <- merge(infoloop, stratif[, c("STRATUM", 
                                              sapply(1:nvar, function(i) (paste("S", i, sep = ""))))], 
                        by = "STRATUM", suffixes = rep("", 2))
      ob <- beat.1st(stratif, errors, minnumstrat)
      n <- ob$n
      psu_file1$SUBSTRAT <- paste(psu_file1$STRATUM, psu_file1$SUB, 
                                  sep = "-")
      psu_strat <- aggregate(PSU_ID ~ SUBSTRAT, psu_file1, 
                             length)
      colnames(psu_strat)[2] <- "PSU_strat"
      psu_file1 <- merge(psu_file1[, which(colnames(psu_file1) != 
                                             "PSU_strat")], psu_strat, by = "SUBSTRAT")
      psu_file1 <- psu_file1[order(psu_file1$STRATUM, 
                                   -psu_file1$PSU_MOS), ]
      psu_file1$AR <- ifelse(psu_file1$PSU_strat <= minPSUstrat, 
                             1, psu_file1$AR)
      PSU_SR <- sum(aggregate(AR ~ SUBSTRAT, data = psu_file1, 
                              sum)[2])
      PSU_NSR <- length(unique(psu_file1$SUBSTRAT[psu_file1$AR == 
                                                    0])) * minPSUstrat
      PSU_Total <- PSU_SR + PSU_NSR
      iter_new <- c(iterx, PSU_SR, PSU_NSR, PSU_Total, 
                    sum(n))
      iterations <- rbind(iterations, iter_new)
      diffx <- max(abs(n - des_file$CAMP))
      loopdeft$diff_n <- (abs(n - des_file$CAMP))
      deft <- des_file[, c("STRATUM", sapply(1:nvar, function(i) (paste("DEFT", 
                                                                        i, sep = ""))))]
      output_beth12 <- list(i_diff_DEFT_fin = loopdeft, 
                            i_S_DEFT_loop = infoloop, i_arnar_fin = des_file, 
                            n_1st = n_1st, psu_trs = psu_file1)
      output_beth12 <<- output_beth12
      des_file <- des_file[, c("STRATUM", "STRAT_MOS", 
                               "MINIMUM", "DELTA")]
      des_file$CAMP <- n
      stratif <- stratif[, c("STRATUM", "N", sapply(1:nvar, 
                                                    function(i) paste("M", i, sep = "")), sapply(1:nvar, 
                                                                                                 function(i) paste("S", i, sep = "")), sapply(1:nvar, 
                                                                                                                                              function(i) paste("S_ORIG", i, sep = "")), sapply(1:ndom, 
                                                                                                                                                                                                function(i) paste("DOM", i, sep = "")), "COST", 
                             "CENS")]
    }
    obb <- list(ob = ob, iterations = iterations)
    return(obb)
  }
  ob_fin <- newDeft(des_file, stratif, errors, rho, psu_file, 
                    deft)
  iterations <- ob_fin$iteration
  ob_fin <- ob_fin$ob
  n_2st <- ob_fin$n
  Bethel_sample_n12 <- cbind(stratif[, c("STRATUM", "N", sapply(1:ndom, 
                                                                function(i) paste("DOM", i, sep = "")))], n_1st, n_2st)
  test_stages <<- 0
  output_beth12$n2_st <- n_2st
  output_beth12$Bethel_sample_n12 <- (cbind(stratif[, c("STRATUM", 
                                                        "N", sapply(1:ndom, function(i) paste("DOM", i, sep = "")))], 
                                            n_1st, n_2st))
  psu_trs <- output_beth12$psu_trs
  psu_trs$AR <- ifelse(psu_trs$PSU_strat <= minPSUstrat, 1, 
                       psu_trs$AR)
  iterations$PSU_SR[nrow(iterations)] <- sum(aggregate(AR ~ 
                                                         SUBSTRAT, data = psu_trs, sum)[2])
  iterations$`PSU NSR`[nrow(iterations)] <- length(unique(psu_trs$SUBSTRAT[psu_trs$AR == 
                                                                             0])) * minPSUstrat
  iterations$`PSU Total`[nrow(iterations)] <- sum(aggregate(AR ~ 
                                                              SUBSTRAT, data = psu_trs, sum)[2]) + length(unique(psu_trs$SUBSTRAT[psu_trs$AR == 
                                                                                                                                    0])) * minPSUstrat
  print(iterations)
  output_beth12 <<- output_beth12
  file_strata <- ob_fin$file_strata
  alloc <- ob_fin$alloc
  alloc[, -1] <- round(alloc[, -1], digits = 0)
  cv <- ob_fin$sensitivity
  planned <- cv[, c("Type","Dom","Dom_label","Var","PlannedCV")]
  varn <- unique(cv[, "Var"])
  step <- length(varn)
  namesP <- planned[seq(from = 1, to = dim(planned)[1], by = step), 
                    1:2]
  valueP <- matrix(round(planned[, "PlannedCV"], 5), ncol = step, byrow = TRUE)
  planned <- cbind(namesP, valueP)
  colnames(planned) <- c(colnames(namesP), paste(varn))
  expected <- cv[, c("Type","Dom","Dom_label","Var","ExpectedCV")]
  step <- length(unique(expected[, "Var"]))
  namesE <- expected[seq(from = 1, to = dim(expected)[1], 
                         by = step), 1:2]
  valueE <- matrix(round(expected[, "ExpectedCV"], 5), ncol = step, byrow = TRUE)
  expected <- cbind(namesE, valueE)
  colnames(expected) <- c(colnames(namesE), paste(varn))
  sensitivity <- cv[, c("Type","Dom","Dom_label","Var","Sensitivity10%")]
  step <- length(unique(sensitivity[, "Var"]))
  namesS <- sensitivity[seq(from = 1, to = dim(sensitivity)[1], 
                            by = step), 1:2]
  valueS <- matrix(sensitivity[, "Sensitivity10%"], ncol = step, byrow = TRUE)
  sensitivity <- cbind(namesS, valueS)
  colnames(sensitivity) <- c(colnames(namesE), paste(varn))
  deft <- output_beth12$i_S_DEFT_loop
  deft_c <- cbind(STRATUM = file_strata$STRATUM, deft[, startsWith(colnames(deft), 
                                                                   "DEFT")])
  out <- list(iterations = iterations, file_strata = file_strata, 
              alloc = alloc, psu_trs = psu_trs, planned = planned, 
              expected = expected, sensitivity = sensitivity, deft_c = deft_c, 
              param_alloc = param_alloc, minimum = des_file$MINIMUM)
  return(out)
}

