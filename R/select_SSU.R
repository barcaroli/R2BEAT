select_SSU <- function (df, PSU_code, SSU_code, PSU_sampled) {
  require(sampling)
  frame <- df
  eval(parse(text=paste0("frame$key <- paste(frame$",PSU_code,",frame$",SSU_code,",sep='*')")))
  eval(parse(text=paste0("df <- df[df$",PSU_code," %in% PSU_sampled$PSU_ID,c(PSU_code,SSU_code)]")))
  df$keynew <- c(1:nrow(df))
  eval(parse(text=paste0("test <- length(unique(df$keynew)) > length(unique(df$",SSU_code,"))")))
  if (test) {
    eval(parse(text=paste0("df2 <- df[!duplicated(df$",SSU_code,"),]")))
  }
  if (!test) {
    df2 <- df
  }
  eval(parse(text=paste0("df2 <- df2[order(df2$",PSU_code,"),]")))
  PSU_sampled <- PSU_sampled[order(PSU_sampled$PSU_ID),]
  s <- strata(df2,stratanames=PSU_code,size=PSU_sampled$PSU_final_sample_unit,method="srswor")
  samp <- getdata(df2,s)
  samp$keynew <- NULL
  if (test) {
    eval(parse(text=paste0("samp$key <- paste(samp$",PSU_code,",samp$",SSU_code,",sep='*')")))
    eval(parse(text=paste0("s <- frame[frame$key %in% samp$key, c('",PSU_code,"','",SSU_code,"')]")))
    samp <- merge(s,samp,all.x=TRUE)
    samp$key <- NULL
  }
  eval(parse(text=paste0("samp <- merge(samp, PSU_sampled,by.x='",PSU_code,"',by.y='PSU_ID')")))
  if (test) {
    samp$weight_2st <- samp$weight <- NULL
    samp$ones <- 1
    eval(parse(text=paste0("a <- aggregate(ones~",PSU_code,",data=samp,FUN=sum)")))
    df$ones <- 1
    eval(parse(text=paste0("b <- aggregate(ones~",PSU_code,",data=df,FUN=sum)")))
    c <- cbind(a,b[,2])
    c$weight_2st <- c[,3] / c[,2]
    samp <- merge(samp,c[,c(1,4)])
    samp$ones <- NULL
  }
  samp$Prob <- samp$ID_unit <- samp$Stratum <- NULL
  samp$weight <- samp$weight_1st * samp$weight_2st
  cat("\n--------------------------------")
  cat("\nTotal PSUs = ", nrow(PSU_sampled))
  eval(parse(text = paste0("cat('\nTotal SSUs = ', length(unique(samp$", SSU_code, ")))")))
  if (test) cat("\nTotal Units = ", nrow(samp))
  cat("\n--------------------------------")
  samp
}