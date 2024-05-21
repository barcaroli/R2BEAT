#' expected_CV
#'  
#' @description
#' Function to report the expected coefficients of variation for target variables Ys in a 'strata' dataset 
#' given an allocation 'alloc' and the current set of precision constraints
#' 
#'  
#' @param strata name of the dataframe containing information in the sampling strata.
#' @param alloc vector containing the allocation of sampling units.
#' @param errors name of the dataframe 
#' 
#' @return a dataframe containing the maximum expected coefficients of variation in each domain level for each target variable
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
#' alloc <- allocation$alloc$ALLOC[-nrow(allocation$alloc)]
#' exp_cv <- expected_CV(strata,cv,alloc)
#' exp_cv
#' 


expected_CV <- function (strata, errors, alloc) 
{
  data <- beat.1cv(strata, errors, alloc)
  split_strings <- strsplit(data$`DOMAIN/VAR`, split = "/")
  data$Dom <- sapply(split_strings, `[`, 1)
  data$Var <- sapply(split_strings, `[`, 2)
  data <- data[, c("TYPE", "Dom", "Var", "ACTUAL_CV")]
  colnames(data) <- c("Type", "Dom", "Var", "Actual CV")
  ndom <- length(unique(data$Type))
  nvar <- length(unique(data$Var))
  unique_combinations <- unique(data[, c("Type", "Dom")])
  exp_cv <- as.data.frame(list(Type = unique_combinations[, 
                                                          1], DOM = unique_combinations[, 2], CV1 = rep(0, nrow(unique_combinations))))
  for (i in c(2:(nvar))) {
    eval(parse(text = paste0("exp_cv$CV", i, " <- 0")))
  }
  for (i in 1:nrow(unique_combinations)) {
    current_combination <- unique_combinations[i, ]
    exp_cv$Type[i] <- unique_combinations[i, 1]
    exp_cv$DOM[i] <- unique_combinations[i, 2]
    subset_data <- subset(data, Type == current_combination$Type & 
                            Dom == current_combination$Dom)
    for (j in c(1:nvar)) {
      eval(parse(text = paste0("exp_cv$CV", j, "[", i, 
                               "] <- subset_data$`Actual CV`[", j, "]")))
    }
  }
  exp_cv
  cv_columns <- grep("CV", names(exp_cv), value = TRUE)
  split_data <- split(exp_cv, exp_cv$Type)
  calculate_max <- function(data_list) {
    sapply(cv_columns, function(column_name) {
      max(data_list[[column_name]], na.rm = TRUE)
    })
  }
  max_values <- lapply(split_data, calculate_max)
  max_values <- lapply(max_values, function(x) round(x, 4))
  max_values_df <- do.call(rbind, max_values)
  row.names(max_values_df) <- unique(exp_cv$Type)
  max_values_df <- as.data.frame(max_values_df)
  max_values_df$DOM <- c(row.names(max_values_df))
  max_values_df <- max_values_df[,c(ncol(max_values_df),c(1:(ncol(max_values_df)-1))),]
  return(max_values_df)
}
