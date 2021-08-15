eval_2stage <- function (df,
                         PSU_code,
                         SSU_code,
                         domain_var,
                         target_vars,
                         PSU_sampled,
                         nsampl = 100, 
                         writeFiles = TRUE,
                         progress = TRUE) 
{
  if (writeFiles == TRUE) {
    dire <- getwd()
    direnew <- paste(dire,"/simulation",sep="")
    if(!dir.exists(direnew)) dir.create(direnew)
    setwd(direnew)
  }
  numY <- length(target_vars)
  numdom <- eval(parse(text=paste0("length(levels(as.factor(df$",domain_var,")))")))
  param <- array(0, c(numY, numdom))
  for (i in c(1:numY)) {
    eval(parse(text = paste("param[i,] <- aggregate(df$",target_vars[i], ",by=list(df$",domain_var,"),FUN=mean)[,2]", 
                            sep = "")))
  }
  estim <- array(0, c(numdom, nsampl, numY))
  differ <- array(0, c(numdom, nsampl, numY))
  # create progress bar
  if (progress == TRUE) pb <- txtProgressBar(min = 0, max = nsampl, style = 3)
  for (j in (1:nsampl)) {
    if (progress == TRUE) Sys.sleep(0.01)
    # update progress bar
    if (progress == TRUE) setTxtProgressBar(pb, j)
    samp <- select_SSU(df,
                       PSU_code,
                       SSU_code,
                       PSU_sampled,
                       verbose=FALSE)
    
    for (k in 1:numY) {
      eval(parse(text = paste0("estim[,j,k] <- aggregate(samp$",target_vars[k],"*samp$weight,by=list(samp$",domain_var,"),FUN=sum)[,2] / aggregate(samp$weight,by=list(samp$",domain_var,"),FUN=sum)[,2]")))
    }
  }
  if (progress == TRUE) close(pb)
  for (j in (1:nsampl)) {
    for (i in (1:numdom)) {
      for (k in 1:numY) {
        differ[i, j, k] <- estim[i, j, k] - param[k, 
                                                  i]
      }
    }
  }
  cv <- array(0, c(numdom, numY))
#  bias <- array(0, c(numdom, numY))
  for (k in 1:numY) {
    for (i in (1:numdom)) {
      cv[i, k] <- sd(estim[i, , k])/mean(estim[i, , k])
#      bias[i, k] <- mean(differ[i, , k]/mean(estim[i, , k]))
    }
  }
  cv <- as.data.frame(cv,stringsAsFactors = TRUE)
  for (i in 1:numY) {
    stmt <- paste0("colnames(cv)[",i,"] <- c('CV",i,"')")
    eval(parse(text = stmt))
  }
  cv$dom <- paste("DOM", c(1:numdom), sep = "")
  if (writeFiles == TRUE) 
    write.table(cv, "expected_cv.csv", sep = ",", row.names = FALSE, 
                col.names = TRUE, quote = FALSE)
  cv1 <- NULL
  cv1$domainvalue <- rep((1:numdom), numY)
  cv1$cv <- rep(0, (numY * numdom))
  cv1$val <- rep(0, (numY * numdom))
  cv1 <- as.data.frame(cv1,stringsAsFactors = TRUE)
  for (k in (1:numY)) {
    for (j in (1:numdom)) {
      cv1$cv[(k - 1) * numdom + j] <- k
      stmt <- paste("cv1$val[(k-1)*numdom + j] <- cv$CV", 
                    k, "[j]", sep = "")
      eval(parse(text = stmt))
    }
  }
  diff <- NULL
  diff$dom <- rep(0, numdom * nsampl)
  diff$samp <- rep(0, numdom * nsampl)
  for (i in (1:numY)) {
    stmt <- paste("diff$diff", i, " <- rep(0,numdom*nsampl)", 
                  sep = "")
    eval(parse(text = stmt))
  }
  diff <- as.data.frame(diff,stringsAsFactors = TRUE)
  for (i in (1:numdom)) {
    for (j in (1:nsampl)) {
      diff$dom[(i - 1) * nsampl + j] <- i
      diff$samp[(i - 1) * nsampl + j] <- j
      for (k in (1:numY)) {
        stmt <- paste("diff$diff", k, "[(i-1)*nsampl+j] <- differ[i,j,k]", 
                      sep = "")
        eval(parse(text = stmt))
      }
    }
  }
  if (writeFiles == TRUE) 
    write.table(diff, "differences.csv", sep = ",", row.names = FALSE, 
                col.names = TRUE, quote = FALSE)
  ############################################   
# New code for bias  
  bias <- NULL
  for (i in (3:ncol(diff))) {
    stmt <- paste("bias$y",i-2," <- tapply(diff$diff",i-2,",diff$dom,mean)",sep="")
    eval(parse(text=stmt))
  } 
  bias$dom <- paste("DOM", c(1:numdom), sep = "")
  bias <- as.data.frame(bias,stringsAsFactors = TRUE)
  
  Y <- array(0, c(numdom, numY))
  
  for (i in c(1:numY)) {
    eval(parse(text=paste0("Y[,i] <- aggregate(df$",target_vars[i],",by=list(df$",domain_var,"),mean)$x")))
  }

  
  bias[,c(1:numY)] <- round(bias[,c(1:numY)]/Y,4)
  if (writeFiles == TRUE)
    write.table(bias, "expected_rel_bias.csv", sep = ",",
                row.names = FALSE, col.names = TRUE, quote = FALSE)
  if (numdom > 1) {
    if (writeFiles == TRUE) 
      # pdf("cv.pdf", width = 7, height = 5)
      png("cv.png")
    boxplot(val ~ cv, data = cv1, col = "orange", main = "Distribution of CV's in the domains", 
            xlab = "Variables Y", ylab = "Value of CV")
    if (writeFiles == TRUE) 
      dev.off()
    if (writeFiles == TRUE) 
      # pdf("rel_bias.pdf", width = 7, height = 5)
      png("rel_bias.png")
    boxplot(bias[,-ncol(bias)], col = "orange", main = "Distribution of relative bias in the domains",
            xlab = "Variables Y", ylab = "Relative bias")
    # boxplot(bias, col = "orange", main = "Distribution of relative bias in the domains", 
    #         xlab = "Variables Y", ylab = "Relative bias")
    if (writeFiles == TRUE) 
      dev.off()
    
  }
  # if (writeFiles == TRUE) 
  # # pdf("differences.pdf", width = 14, height = 10)
  # png("differences.png")
  # k <- ceiling(numY/4)
  # for (j in 1:k) {
  #   split.screen(c(2, 2))
  #   for (i in 1:4) {
  #     if (i + 4 * (j - 1) <= numY) {
  #       stmt <- paste("screen(", i, ")", sep = "")
  #       eval(parse(text = stmt))
  #       stmt <- paste("boxplot(diff", i, "~dom,data=diff,ylab='Differences',xlab='Domain',col = 'orange')", 
  #                     sep = "")
  #       eval(parse(text = stmt))
  #       stmt <- paste("mtext(expression(Y", i, "), side=3, adj=0, cex=1.0, line=1)", 
  #                     sep = "")
  #       eval(parse(text = stmt))
  #     }
  #   }
  #   close.screen(all.screens = TRUE)
  # }
  # if (writeFiles == TRUE) 
  #   dev.off()
  # results <- list(coeff_var = cv1, bias = bias1)
  est <- matrix(NA,nrow=numdom*nsampl,ncol=numY)
  est <- as.data.frame(est,stringsAsFactors = TRUE) 
  colnames(est) <- c(paste("Y",c(1:numY),sep=""))
  est$dom <- rep(c(1:numdom),each=nsampl)
  for (i in (1:numdom)) {
    est[est$dom == i,c(1:(numY))] <- estim[i,,]
  }
  if (writeFiles == TRUE) {
    write.table(est,"estimates.csv",sep=",",row.names=F,col.names=F)
  }

  cv[,c(1:numY)] <- round(cv[,c(1:numY)],4)
  results <- list(coeff_var = cv, rel_bias = bias, est = est)
  cv <- cbind(c(1:nrow(cv)),cv)
  colnames(cv) <- c("domain",paste("cv(Y",c(1:numY),")",sep=""))
  bias <- cbind(c(1:nrow(bias)),bias)
  colnames(bias) <- c("domain",paste("bias(Y",c(1:numY),")",sep=""))
  cv <- formattable(cv,list(area(col = 2:(numY+1)) ~ color_tile("#DeF7E9", "#71CA97")))
  bias <- formattable(bias,list(area(col = 2:(numY+1)) ~ color_tile("#DeF7E9", "#71CA97")))
  if (writeFiles == TRUE) {
    setwd(dire)
  }
  return(results)
}
