select_SSU <- function (df, PSU_code, SSU_code, PSU_sampled) 
{
  eval(parse(text = paste0("df$key <- paste(df$", PSU_code, 
                           ",df$", SSU_code, ",sep='*')")))
  df$keynew <- c(1:nrow(df))
  frame <- df
  eval(parse(text = paste0("df <- df[df$", PSU_code, " %in% PSU_sampled$PSU_ID,c('",PSU_code,"','",SSU_code,"')]")))
  df$keynew <- c(1:nrow(df))
  test <- NULL
  eval(parse(text = paste0("test <- length(unique(df$keynew)) > length(unique(df$", 
                           SSU_code, "))")))
  if (test) {
    eval(parse(text = paste0("df2 <- df[!duplicated(df$", 
                             SSU_code, "),]")))
  }
  if (!test) {
    df2 <- df
  }
  eval(parse(text = paste0("df2 <- df2[order(df2$", PSU_code, 
                           "),]")))
  PSU_sampled <- PSU_sampled[order(PSU_sampled$PSU_ID), ]
  s <- strata(df2, stratanames = PSU_code, size = PSU_sampled$PSU_final_sample_unit, 
              method = "srswor")
  samp <- getdata(df2, s)
  eval(parse(text = paste0("samp$key <- paste(samp$", PSU_code, 
                           ",samp$", SSU_code, ",sep='*')")))
  s <- frame[frame$key %in% samp$key, ]
  samp <- merge(s, samp, all.x = TRUE)
  eval(parse(text = paste0("samp <- merge(samp, PSU_sampled,by.x='", 
                           PSU_code, "',by.y='PSU_ID')")))
  if (test) {
    samp$weight <- samp$pik_2st <- NULL
    samp$ones <- 1
    a <- NULL
    eval(parse(text = paste0("a <- aggregate(ones~", PSU_code, 
                             ",data=samp,FUN=sum)")))
    df$ones <- 1
    b <- NULL
    eval(parse(text = paste0("b <- aggregate(ones~", PSU_code, 
                             ",data=df,FUN=sum)")))
    c <- cbind(a, b[, 2])
    c$pik_2st <- c[, 2]/c[, 3]
    samp <- merge(samp, c[, c(1, 4)])
    samp$Pik <- samp$ones <- NULL
  }
  samp$weight <- (1/samp$pik_1st) * (1/samp$pik_2st)
  samp$ones <- samp$Pik <- NULL
  samp$Prob <- samp$ID_unit <- samp$Stratum <- samp$key <- samp$keynew <- samp$stratum.x <- NULL
  colnames(samp)[colnames(samp) == "STRATUM"] <- "stratum"
  colnames(samp)[colnames(samp) == "stratum.y"] <- "stratum_2"
  cat("\n--------------------------------")
  cat("\nTotal PSUs = ", nrow(PSU_sampled))
  eval(parse(text = paste0("cat('\nTotal SSUs = ', length(unique(samp$", 
                           SSU_code, ")))")))
  if (test) 
    cat("\nTotal Units = ", nrow(samp))
  cat("\n--------------------------------")
  samp
}
