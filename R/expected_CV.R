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
