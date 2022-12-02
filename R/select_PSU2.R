select_PSU2 <- function(alloc, 
                       type="ALLOC", 
                       var_ord=NULL) {
  require(sampling)
  allocation <- alloc$alloc[-nrow(alloc$alloc),c("STRATUM",type)]
  allocation$PSUs <- round(allocation$ALLOC/24)
  alloc$PSUs <- ifelse(alloc$PSUs == 0,1,alloc$PSUs)
  psu_file <- alloc$psu_trs[,c(4,2,5)]
  psu_file$ones <- 1
  PSUs <- aggregate(cbind(ones,PSU_MOS)~STRATUM,data=psu_file,FUN=sum)
  PSUs <- merge(PSUs,allocation)
  colnames(PSUs) <- c("STRATUM","PSUs_available","SSUs_available",
                      "SSUs_allocated","PSUs_allocated")
  PSUs <- PSUs[,c(1,2,5,3,4)]
  psu_frame <- merge(psu_file,PSUs,by=c("STRATUM"))
  if (!is.null(var_ord)) {
    psu_frame <- merge(psu_frame,var_ord,by="PSU_ID")
    vars <- colnames(var_ord)[c(2:ncol(var_ord))]
    stmt <- NULL  
    for (v in c(1:length(vars))) {
      if (v == 1) {
        stmt <- paste0("psu_frame <- psu_frame[order(psu_frame$",vars[v],",")
      }
      if (v > 1 & v < length(vars)) {
        stmt <- paste0(stmt,"psu_frame",vars[v],",")
      }
      if (v == length(vars)) {
        stmt <- paste0(stmt,"psu_frame$",vars[v],"),]")
      }
    }
    eval(parse(text=stmt)) 
  }
  for (i in c(unique(psu_frame$STRATUM))) {
    psu_to_be_selected <- PSUs$PSUs_allocated[PSUs$STRATUM == i]
    psu_frame$pik_1st[psu_frame$STRATUM == i] <- inclusionprobabilities(psu_frame$PSU_MOS[psu_frame$STRATUM == i],psu_to_be_selected)
    try(psu_frame$sampled[psu_frame$STRATUM == i] <- UPsystematic(pik=psu_frame$pik_1st[psu_frame$STRATUM == i]))
    cat("\n Stratum: ",i,"  PSUs to be selected: ",psu_to_be_selected,"  PSUs selected: ",sum(psu_frame$sampled[psu_frame$STRATUM == i & psu_frame$sampled == 1]))
  }
  sampled_PSUs <- psu_frame[psu_frame$sampled==1,]
  sampled_PSUs$pik_2st <- 24 / sampled_PSUs$PSU_MOS
  sampled_PSUs$weight <- (1/sampled_PSUs$pik_1st) * (1/sampled_PSUs$pik_2st)
  sampled_PSUs$ones <- NULL
  sampled_PSUs$PSUs_available <- NULL
  sampled_PSUs$PSUs_allocated <- NULL
  sampled_PSUs$SSUs_available <- NULL
  sampled_PSUs$SSUs_allocated <- NULL
  sampled_PSUs$sampled <- NULL
  return(sampled_PSUs)
}