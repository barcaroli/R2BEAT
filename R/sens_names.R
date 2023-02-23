sens_names <- function(s,target_vars,strata) {
  s$Var <- target_vars
  s$Dom[s$Type=="DOM1"] <- "Total"
  ndom <- sum(grepl("DOM",colnames(strata))) - 1
  for (j in c(1:ndom)) {
    eval(parse(text=paste0("t <- s$Dom[s$Type=='DOM",j+1,"']")))
    # t <- as.factor(t)
    eval(parse(text=paste0("strata$DOM",j+1," <- as.factor(strata$DOM",j+1,")")))
    eval(parse(text=paste0("for (i in c(1:length(levels(strata$DOM",j+1,")))) t[t==i] <- levels(strata$DOM",j+1,")[i]")))
    eval(parse(text=paste0("s$Dom[s$Type=='DOM",j+1,"'] <- t")))
  }
  return(s)
}
