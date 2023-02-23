select_SSU <- function (df, PSU_code, SSU_code, PSU_sampled, verbose = TRUE) 
{
  select <- function(df, PSU_code, PSU_sampled) {
    st <- paste0("PSU_ID <- df$", PSU_code, "[1]")
    eval(parse(text = st))
    st <- paste0("PSU_size <- PSU_sampled$PSU_final_sample_unit[PSU_sampled$PSU_ID==df$", 
                 PSU_code, "[1]]")
    eval(parse(text = st))
    s <- sample(1:nrow(df), size = PSU_size)
    s <- df[s, ]
    st <- paste0("s$Pik <- PSU_sampled$Pik[PSU_sampled$PSU_ID==df$", 
                 PSU_code, "[1]]")
    eval(parse(text = st))
    s$Prob <- nrow(s)/nrow(df)
    if (verbose == TRUE) {
      cat("\nPSU = ", as.character(PSU_ID), " *** Selected SSU = ", 
          nrow(s))
    }
    return(s)
  }
  PSU_size <- NULL
  PSU_ID <- NULL
  df.split <- NULL
  a <- NULL
  st <- paste0("df <- df[df$", PSU_code, " %in% PSU_sampled$PSU_ID, ]")
  eval(parse(text = st))
  st <- paste0("df$municipality <- as.factor(df$", PSU_code, 
               ")")
  eval(parse(text = st))
  st <- paste0("df.split <- split(df,df$", PSU_code, ",drop=TRUE)")
  eval(parse(text = st))
  samp <- sapply(df.split, function(df) select(df, PSU_code, 
                                               PSU_sampled))
  samp <- data.frame(t(samp))
  st <- paste0(SSU_code," <- unlist(samp$", SSU_code, ")")
  eval(parse(text = st))
  Pik <- unlist(samp$Pik)
  Prob <- unlist(samp$Prob)
  st <- paste0("a <- as.data.frame(list(", SSU_code, ",Pik,Prob))")
  eval(parse(text = st))
  st <- paste0("colnames(a) <- c('", SSU_code, "','Prob_1st','Prob_2st')")
  eval(parse(text = st))
  samp <- merge(df, a)
  samp$Prob_tot <- samp$Prob_1st * samp$Prob_2st
  samp$weight <- 1/samp$Prob_tot
  if (verbose == TRUE) {
    cat("\n--------------------------------")
    cat("\nTotal PSU = ", nrow(PSU_sampled))
    cat("\nTotal SSU = ", nrow(samp))
    cat("\n--------------------------------")
  }
  # PSU <- PSU_sampled[, c("PSU_ID", "SR", "nSR", "stratum")]
  # colnames(PSU)[4] <- "stratum_2"
  # st <- paste0("samp <- merge(samp,PSU,by.x=c('", PSU_code, 
  #              "'),by.y='PSU_ID')")
  # eval(parse(text = st))
  return(samp)
}
