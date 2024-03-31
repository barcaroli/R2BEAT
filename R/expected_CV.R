#' expected_CV
#'  
#' @description
#' Function to calculate expected coefficients of variation for target variables Ys in a 'strata' dataset 
#' given an allocation 'alloc'
#' 
#' @details
#' The 'strata' dataframe can be obtained by executing the function 'prepareInputToAllocation_beat.1st'
#' The 'strata' dataframe must contain a variable with the indication of the sampling units allocated in each stratum.
#' 
#' @param strata name of the dataframe containing information in the sampling strata.
#' @param alloc_variable name of the variable in the strata dataframe containing the allocation of sampling units.
#' @param errors name of the dataframe 
#' 
#' @return a dataframe containing the maximum expected coefficients of variation in each domain for each target variable
#' 
#' @examples
#' load("./data/sample.RData")
#' target_vars <- c("active","inactive","unemployed","income_hh")
#' strata <- R2BEAT:::prepareInputToAllocation_beat.1st(samp_frame = samp,
#'                                                      ID = "id_hh",
#'                                                      stratum = "stratum_label",
#'                                                      dom = "region",
#'                                                      target = target_vars)
#' strata$CENS <- as.numeric(strata$CENS)
#' strata$COST <- as.numeric(strata$COST)
#' strata$CENS <- 0
#' cv <- as.data.frame(list(DOM = c("DOM1","DOM2"),
#'                          CV1 = c(0.05,0.10),
#'                          CV2 = c(0.05,0.10),
#'                          CV3 = c(0.05,0.10),
#'                          CV4 = c(0.05,0.10)))
#' allocation <- beat.1st(strata,cv)
#' 
#' strata$alloc <- allocation$alloc$ALLOC[-nrow(allocation$alloc)]
#' exp_cv <- expected_CV(strata,"alloc",cv)
#' exp_cv
expected_CV <- function (strata, alloc_variable, errors) 
{
  nvars <- ncol(errors) - 1
  colnames(errors)[2:(nvars + 1)] <- paste("cv(Y", c(1:nvars), 
                                           ")", sep = "")
  eval(parse(text = paste0("if (is.null(strata$", alloc_variable, 
                           ")) stop('There is no allocation of units in strata')")))
  ndom <- length(grep("DOM", colnames(strata)))
  exp_cv <- NULL
  for (d in c(1:ndom)) {
    eval(parse(text = paste0("vett <- strata$DOM", d)))
    vett <- as.integer(factor(as.character(vett)))
    strata2 <- aggrStrata(strata, nvars, vett, dominio = 1, 
                          censiti = 0)
    eval(parse(text = paste0("strata2$", alloc_variable, 
                             " <- tapply(strata$", alloc_variable, ",vett,FUN=sum)")))
    M_h <- S_h <- NULL
    for (i in (1:nrow(strata2))) {
      exp_cv2 <- errors[1, ]
      exp_cv2[, c(2:(nvars + 1))] <- NA
      strata3 <- strata2[i, ]
      for (j in 1:nvars) {
        eval(parse(text = paste0("n_h <- strata3$", 
                                 alloc_variable)))
        N_h <- strata3$N
        eval(parse(text = paste0("S_h <- strata3$S", 
                                 j)))
        eval(parse(text = paste0("M_h <- strata3$M", 
                                 j)))
        Y_h <- N_h * M_h
        Var_h <- (N_h^2) * (1 - n_h/N_h) * ((S_h^2)/n_h)
        CV <- sqrt(Var_h)/Y_h
        exp_cv2[1, j + 1] <- CV
        exp_cv2$DOM[1] <- paste0("DOM", d)
      }
      exp_cv <- rbind(exp_cv, exp_cv2)
    }
  }
  cv_columns <- grep("cv", names(exp_cv), value = TRUE)
  split_data <- split(exp_cv, exp_cv$DOM)
  calculate_max <- function(data_list) {
    sapply(cv_columns, function(column_name) {
      max(data_list[[column_name]], na.rm = TRUE)
    })
  }
  max_values <- lapply(split_data, calculate_max)
  max_values <- lapply(max_values,function(x) round(x,4))
  max_values_df <- do.call(rbind, max_values)
  row.names(max_values_df) <- unique(exp_cv$DOM)
  return(max_values_df)
}
