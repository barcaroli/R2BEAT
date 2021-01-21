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
      cat("\nPSU = ", PSU_ID, " *** Selected SSU = ", 
          nrow(s))
    }
    return(s)
  }
  st <- paste0("df <- df[df$", PSU_code, " %in% PSU_sampled$PSU_ID, ]")
  eval(parse(text = st))
  st <- paste0("df$municipality <- as.factor(df$", PSU_code, 
               ")")
  eval(parse(text = st))
  st <- paste0("df.split <- split(df,df$", PSU_code, 
               ",drop=TRUE)")
  eval(parse(text = st))
  samp <- sapply(df.split, function(df) select(df, PSU_code, 
                                                  PSU_sampled))
  samp <- data.frame(t(samp))
  st <- paste0("id_ind <- unlist(samp$", SSU_code, ")")
  eval(parse(text = st))
  Pik <- unlist(samp$Pik)
  Prob <- unlist(samp$Prob)
  st <- paste0("a <- as.data.frame(list(", SSU_code, ",Pik,Prob))")
  eval(parse(text = st))
  st <- paste0("colnames(a) <- c('", SSU_code, "','Prob_1st','Prob_2st')")
  eval(parse(text = st))
  samp <- merge(pop, a)
  samp$Prob_tot <- samp$Prob_1st * samp$Prob_2st
  samp$weight <- 1/samp$Prob_tot
  if (verbose == TRUE) {
    cat("\n--------------------------------")
    cat("\nTotal PSU = ", nrow(PSU_sampled))
    cat("\nTotal SSU = ", nrow(samp))
    cat("\n--------------------------------")
  }
  return(samp)
}
#-----------------------------------------
# Function to select 2nd stage units (SSU)
#-----------------------------------------

select_SSU <- function(df,PSU_code,SSU_code,PSU_sampled,verbose=TRUE) {
  select <- function(df,PSU_code,PSU_sampled) {
    st <- paste0("PSU_ID <- df$",PSU_code,"[1]")
    eval(parse(text=st))
    st <- paste0("PSU_size <- PSU_sampled$PSU_final_sample_unit[PSU_sampled$PSU_ID==df$",PSU_code,"[1]]")
    eval(parse(text=st))
    s <- sample(1:nrow(df),size=PSU_size)
    s <- df[s,]
    st <- paste0("s$Pik <- PSU_sampled$Pik[PSU_sampled$PSU_ID==df$",PSU_code,"[1]]")
    eval(parse(text=st))
    s$Prob <- nrow(s) / nrow(df)
    if (verbose==TRUE) {
      cat("\nPSU = ",PSU_ID," *** Selected SSU = ",nrow(s))
    }
    return(s)
  }
  st <- paste0("frame <- frame[frame$",PSU_code," %in% PSU_sampled$PSU_ID, ]")
  eval(parse(text=st))
  st <- paste0("frame$municipality <- as.factor(frame$",PSU_code,")")
  eval(parse(text=st))
  st <- paste0("frame.split <- split(frame,frame$",PSU_code,",drop=TRUE)")
  eval(parse(text=st))  
  samp <- sapply(frame.split,function(df) select(df,PSU_code,PSU_sampled))
  samp <- data.frame(t(samp))
  st <- paste0("id_ind <- unlist(samp$",SSU_code,")")
  eval(parse(text=st))   
  Pik <- unlist(samp$Pik)
  Prob <- unlist(samp$Prob)
  st <- paste0("a <- as.data.frame(list(",SSU_code,",Pik,Prob))")
  eval(parse(text=st)) 
  st <- paste0("colnames(a) <- c('",SSU_code,"','Prob_1st','Prob_2st')")
  eval(parse(text=st)) 
  samp <- merge(pop,a)
  samp$Prob_tot <- samp$Prob_1st * samp$Prob_2st
  samp$weight <- 1/samp$Prob_tot
  if (verbose==TRUE) {
    cat("\n--------------------------------")
    cat("\nTotal PSU = ",nrow(PSU_sampled))
    cat("\nTotal SSU = ",nrow(samp))
    cat("\n--------------------------------")
  }  
  return(samp)
}


