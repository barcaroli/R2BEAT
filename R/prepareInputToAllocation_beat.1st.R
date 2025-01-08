#' The function returns a dataframe, starting from the sampling frame (either universe or sample of a previous survey) with strata information.
#'
#' @param samp_frame: A dataframe with the sampling frame (either universe or sample of another survey).
#' @param ID: name of the variable which is the identifier of the units in the sampling frame.
#' @param stratum: either name of the variable in samp_frame which is taken as the stratum, or the name of the variables which have to be concatenated to obtain the stratum. In this last case, the variables used to build the stratum are retained..
#' @param dom: name of the variable(s) in samp_frame which are the domain(s) of the estimates.
#' @param target: variable(s) in samp_frame which is the target of the estimate.
#' @param samp_weight: name of the variable in the sample of another survey specifying the sampling weights. Default is NULL, meaning that population is taken as the reference.
#' @return A dataframe with strata information. The output dataframe can be used as input dataframe stratif for R2BEAT one-stage sample design (beat.1st)
prepareInputToAllocation_beat.1st <- function (samp_frame, ID, stratum, dom,  target, samp_weight=NULL){
  buildFrame <- function (df, id, X, Y, domainvalue) 
  {
    vars <- colnames(df)
    if (!id %in% vars) 
      stop("Id name not in the frame variables")
    for (i in (1:length(X))) {
      var <- X[i]
      if (!var %in% vars) {
        msg = paste("X variable ", var, " not in the frame variables", 
                    sep = "")
        stop(msg)
      }
    }
    for (i in (1:length(Y))) {
      var <- Y[i]
      if (!var %in% vars) {
        msg <- paste("Y variable ", var, " not in the frame variables", 
                    sep = "")
        stop(msg)
      }
    }
    for (i in (1:length(domainvalue))){
      var <- domainvalue[i]
      if (!var %in% vars){
        msg=paste("Domain name", var,"not in the frame variables")
        stop(msg)
      }
    }
    
    dframe <- NULL
    stmt <- paste("dframe$id <- df$", as.character(id), sep = "")
    eval(parse(text = stmt))
    for (i in (1:length(X))) {
      aux <- paste("dframe$X", i, " <- df$", as.character(X[i]), 
                   sep = "")
      eval(parse(text = aux))
    }
    for (i in (1:length(Y))) {
      target <- paste("dframe$Y", i, " <- df$", as.character(Y[i]), 
                      sep = "")
      eval(parse(text = target))
    }
    for (i in (1:length(domainvalue))){
      stmt <- paste("dframe$DOMAINVALUE", i, "<- df$", as.character(domainvalue[i]), 
                    sep = "")
      eval(parse(text = stmt))
      #stmt <- paste("dframe$DOMAINVALUE",i,"<-as.integer(dframe$DOMAINVALUE", i,")", sep="")
      stmt <- paste("dframe$DOMAINVALUE",i,"<-as.factor(dframe$DOMAINVALUE", i,")", sep="")
      eval(parse(text=stmt))
      
      stmt <- paste("dframe$domain_cat_labels", i, "<- df$", as.character(domainvalue[i]), 
                    sep = "")
      eval(parse(text = stmt))
    }
    
    #colnames(df)=toupper(colnames(df))
    #if (length(grep("WEIGHT", names(df))) == 1){
    if(!is.null(samp_weight)){
      txt <- paste("dframe$WEIGHT <- df$", as.character(samp_weight), sep="")
      eval(parse(text=txt))
      #dframe$WEIGHT=df$WEIGHT
    }
    dframe <- as.data.frame(dframe, stringsAsFactors = TRUE)
    
    return(dframe)
  }
  
  buildStrata <- function (dataset, progress = TRUE, verbose = TRUE) 
  {
    numdom <- NULL
    stratatot <- NULL
    stratadom <- NULL
    colnames(dataset) <- toupper(colnames(dataset))
    nvarX <- length(grep("X", names(dataset)))
    nvarY <- length(grep("Y", names(dataset)))
    
    stdev1 <- function(x, w) {
      mx <- sum(x * w)/sum(w)
      sqrt(sum(w * (x - mx)^2)/(sum(w) - 1))
    }
    stdev2 <- function(x, w) {
      mx <- sum(x * w)/sum(w)
      sqrt(sum(w * (x - mx)^2)/(sum(w)))
    }
    stdev_prop <- function(x, w) {
      mx <- sum(x * w)/sum(w)
      sqrt(mx*(1-mx))
    }
   
    avgf <- function(x, w){
      sum(x * w)/sum(w)
    }
    
    if (length(grep("WEIGHT", names(dataset))) == 1) {
      if (verbose == TRUE) {
        cat("\nComputations are being done on sampling data\n")
      }
      
      avg <- "avgf"
      stdev <- NULL
      for (i in 1:nvarY){
        if (eval(parse(text=paste("nlevels(as.factor(dataset$Y",i,"))==2", sep="")))){
        stdev[[i]] <- "stdev_prop"
        if(eval(parse(text=paste("all(levels(as.factor(dataset$Y", i, "))!=c('0','1'))", sep="")))){
          #check levels of the dummy variable and convert levels to compute the proportion
          
          cat(paste("\nConverting levels of factor variable Y", i," into 0/1 \n", sep=""))
          txt=paste("dataset$Y", i, " <- as.factor(dataset$Y",i,")", sep="")
          eval(parse(text=txt)) #to factor
          txt=paste("levels(dataset$Y", i,") <- c('0', '1')", sep="")
          eval(parse(text=txt)) #change levels
          txt=paste("dataset$Y", i, " <- as.numeric(levels(dataset$Y",i,"))[dataset$Y", i, "]", sep="")
          eval(parse(text=txt)) #to numeric
        }
        }
        else{
          stdev[[i]] <- "stdev1" 
        }

      }
    }
    if (length(grep("WEIGHT", names(dataset))) == 0) {
      if (verbose == TRUE) {
        cat("\nComputations are being done on population data\n")
      }
      dataset$WEIGHT <- rep(1, nrow(dataset))
      
      avg="avgf"
      stdev <- NULL
      for (i in 1:nvarY){
        if (eval(parse(text=paste("nlevels(as.factor(dataset$Y",i,"))==2", sep="")))){
          stdev[[i]] <- "stdev_prop"
          if(eval(parse(text=paste("all(levels(as.factor(dataset$Y", i, "))!=c('0','1'))", sep="")))){
            #check levels of the dummy variable and convert levels to compute the proportion
            
            cat(paste("\nConverting levels of factor variable Y", i," into 0/1 \n", sep=""))
            txt=paste("dataset$Y", i, " <- as.factor(dataset$Y",i,")", sep="")
            eval(parse(text=txt)) #to factor
            txt=paste("levels(dataset$Y", i,") <- c('0', '1')", sep="")
            eval(parse(text=txt)) #cahnge levels
            txt=paste("dataset$Y", i, " <- as.numeric(levels(dataset$Y",i,"))[dataset$Y", i, "]", sep="")
            eval(parse(text=txt)) #to numeric
          }
        }
        else{
         
         stdev[[i]] <- "stdev2" 
        }
       
      }
    }
    #define the number of domains:
    
    ndom <- sum(grepl("DOMAINVALUE", colnames(dataset)))
    
    
    for (idx_d in 1:ndom){
    #   if (verbose == TRUE) {
    #     cat(paste("\nDOM", idx_d, "\n", sep=""))
    #   }
   # dataset$DOMAINVALUE<- factor(dataset$DOMAINVALUE)
      txt <- paste("dataset$DOMAINVALUE", idx_d, "<- factor(dataset$DOMAINVALUE", idx_d, ")", sep="")
      eval(parse(text=txt)) #domains variables are converted to factors
    }
      # txt2 <- paste("numdom[idx_d] <- length(levels(dataset$DOMAINVALUE", idx_d, "))", sep="")
      # eval(parse(text=txt2)) #save the number of domain categories into the vector numdom
      # 
      # if (progress == TRUE){ 
      #   pb <- txtProgressBar(min = 0, max = numdom[idx_d], style = 3)
      # }
      # doms <- 0
      # #for each category d in domain idx_d:
      # for (d in (eval(parse(text=paste("levels(dataset$DOMAINVALUE", idx_d,")", sep=""))))){
      #   
        # if (verbose == TRUE) {
        #   cat(paste("\nDOM", idx_d, " category: ", d,"\n", sep=""))
        # }
        # 
        # doms <- doms + 1
        # if (progress == TRUE){
        #   Sys.sleep(0.1)
        #   setTxtProgressBar(pb, doms)
        # }
        
        # if (progress == TRUE) 
        #   setTxtProgressBar(pb, doms)
        # dom <- d
        # txt=paste("domain <- dataset[dataset$DOMAINVALUE", idx_d,"== dom,]", sep="")
        # eval(parse(text=txt)) #select only rows of dataset having category dom of variable idx_d
        # 
    
    listX <- NULL #save list of X variables 
    namesX <- NULL #save names of X variables 
    for (i in 1:nvarX) {
      name <- paste("X", i, sep = "")
      namesX <- cbind(namesX, name)
      if (i < nvarX) 
        listX <- paste(listX, "dataset$X", i, ",", sep = "")
      else listX <- paste(listX, "dataset$X", i, sep = "")
    }
    listM <- NULL
    listS <- NULL
    for (i in 1:nvarY) {
      if (i<nvarY){
        listM <- paste(listM, "M", i, ",", sep = "")
        listS <- paste(listS, "S", i, ",", sep = "")
      } 
      else{
        listM <- paste(listM, "M", i,  sep = "")
        listS <- paste(listS, "S", i,  sep = "")
        
      }
    }
        #to obtain stratum, paste the categories of the stratification variables, using as separator '*'
        stmt <- paste("dataset$STRATUM <- as.factor(paste(", listX, 
                      ",sep='*'))", sep = "")
        eval(parse(text = stmt))
        
        
        
        # if(idx_d==1){
          if (!is.null(dataset$COST)){
            cost <- tapply(dataset$WEIGHT * dataset$COST, dataset$STRATUM,sum)/tapply(dataset$WEIGHT, 
                                                                                   dataset$STRATUM, sum)
          } #cost function to be used later if the variable COST in dataset is specified.
          
          
          
          for (i in 1:nvarY) {
            WEIGHT <- NULL
            STRATUM <- NULL
            Y <- NULL
            stmt <- paste("Y <- dataset$Y", i, "[!is.na(dataset$Y", 
                          i, ")]", sep = "")
            eval(parse(text = stmt)) #save in Y the values of variable Yi without NA  
            if(length(Y)==0){
              cat(paste("\nVariable Y", i, 
                        #" in the category ", d," of domain DOM", idx_d,
                        " is completely missing. By default, M",i, 
                        " and S", i, 
                        " are set equal to 0 in all the strata \n", sep=""))
            }
            stmt <- paste("WEIGHT <- dataset$WEIGHT[!is.na(dataset$Y", 
                          i, ")]", sep = "")
            eval(parse(text = stmt)) #save in WEIGHT the values of variable WEIGHT (corresponding to non missing Y values) 
            stmt <- paste("STRATUM <- dataset$STRATUM[!is.na(dataset$Y", i, ")]", sep = "")
            eval(parse(text = stmt)) #save in STRATUM the values of variable STRATUM (corresponding to non missing Y values) 
            STRATUM <- factor(STRATUM)
            
            
            samp <- NULL
            stmt <- paste("samp <- dataset[!is.na(dataset$Y", 
                          i, "),]", sep = "")
            eval(parse(text = stmt)) #define samp as dataframe 'domain' (including only rows with non-missing values of the variable Yi)
            l.split <- split(samp, samp$STRATUM, drop = TRUE) #split divides the data samp into the groups defined by STRATUM.
            stmt <- paste("M", i, " <- sapply(l.split, function(df,x,w) ",  avg, "(df[,x],df[,w]), x='Y", i, "', w='WEIGHT')", 
                          sep = "")
            eval(parse(text = stmt))  # separately per STRATUM, save the mean value of the variable Yi
            
            stmt <- paste("S", i, " <- sapply(l.split, function(df,x,w) ",  stdev[[i]], "(df[,x],df[,w]), x='Y", i, "', w='WEIGHT')", 
                          sep = "")
            eval(parse(text = stmt)) # separately per STRATUM, save the standard deviation value of the variable Yi
            
            stmt <- paste("stratirid <- names(M", i,")", sep="")
            eval(parse(text = stmt)) #save the identifier of the stratum (only stratum corresponding to non-missing values of Yi)
            strati <- data.frame(X1 = levels(dataset$STRATUM), 
                                 stringsAsFactors = TRUE) #create a dataframe with variable X1: STRATUM
            stmt <- paste("m <- data.frame(cbind(X1=stratirid,X2=M", 
                          i, "), stringsAsFactors = TRUE)", sep = "")
            eval(parse(text = stmt)) #create a dataframe with variable x1: STRATUM, x2: mean of the target variable Yi
            m <- merge(strati, m, by = c("X1"), all = TRUE) #merge by strato
            
            m$X2 <- ifelse(is.na(m$X2), m$X2, as.character(m$X2)) #convert X2 as character
            m$X2 <- ifelse(is.na(m$X2), m$X2, as.numeric(m$X2))  #convert X2 as numeric
            m$X2 <- ifelse(is.na(m$X2), 0, m$X2) #replace NA in X2 with 0
            stmt <- paste("M", i, " <- m$X2", sep = "") # update Mi
            eval(parse(text = stmt))
            stmt <- paste("s <- data.frame(cbind(X1=stratirid,X2=S", 
                          i, "), stringsAsFactors = TRUE)", sep = "")
            eval(parse(text = stmt))  #create a dataframe with variable x1: STRATUM, x2: standard deviation of the target variable Yi
            s <- merge(strati, s, by = c("X1"), all = TRUE) #merge by stratum
            s$X2 <- ifelse(is.na(s$X2), s$X2, as.character(s$X2)) #convert X2 as character
            s$X2 <- ifelse(is.na(s$X2), s$X2, as.numeric(s$X2)) #convert X2 as numeric
            s$X2 <- ifelse(is.na(s$X2), 0, s$X2) #replace NA in X2 with 0
            stmt <- paste("S", i, " <- s$X2", sep = "")
            eval(parse(text = stmt)) # update Si
          }
          
          N <- round(tapply(dataset$WEIGHT, dataset$STRATUM, sum))
          #save number of units per stratum
          STRATUM <- dataset$STRATUM 
          #define COST variable
          if (is.null(dataset$COST)){
            COST <- rep(1, length(levels(dataset$STRATUM)))
          }
          if (!is.null(dataset$COST)){
            COST <- cost
          }
          
          
          #define CENS variable:
          if(!is.null(dataset$CENS)){
            CENS <- dataset$CENS[match(levels(STRATUM), STRATUM)]
          } else CENS <- rep(0, length(levels(dataset$STRATUM)))
          
          
          

          dataset_dom=data.frame(matrix(nrow=length(levels(STRATUM)), ncol=ndom))
          colnames(dataset_dom) = paste("DOM", 1:ndom, sep="")
          for (idx_d in 1:ndom){
            stmt=paste("DOM",idx_d, " <- dataset$DOMAINVALUE", idx_d, "[match(levels(STRATUM), STRATUM)]", sep="")
            eval(parse(text=stmt))  #define DOMidx_d variable
            stmt=paste("dataset_dom$DOM",idx_d, " <- DOM", idx_d, sep="")
            eval(parse(text=stmt))  #define DOMidx_d variable
            
          }
            
            # if (verbose == TRUE) {
            #   cat(paste("\nDOM", idx_d, "\n", sep=""))
            # }
            # txt <- paste("dataset$DOMAINVALUE", idx_d, "<- factor(dataset$DOMAINVALUE", idx_d, ")", sep="")
            # eval(parse(text=txt)) #domains variables are converted to factors
            #txt2 <- paste("numdom[idx_d] <- length(levels(dataset$DOMAINVALUE", idx_d, "))", sep="")
            #eval(parse(text=txt2)) #save the number of domain categories into the vector numdom

            # if (progress == TRUE){
            #   pb <- txtProgressBar(min = 0, max = numdom[idx_d], style = 3)
            # }
          #   doms <- 0
          #   #for each category d in domain idx_d:
          #   for (d in (eval(parse(text=paste("levels(dataset$DOMAINVALUE", idx_d,")", sep=""))))){
          # 
          #     if (verbose == TRUE) {
          #       cat(paste("\nDOM", idx_d, " category: ", d,"\n", sep=""))
          #     }
          # 
          #     doms <- doms + 1
          #     # if (progress == TRUE){
          #     #   Sys.sleep(0.1)
          #     #   setTxtProgressBar(pb, doms)
          #     # }
          # 
          #     # if (progress == TRUE)
          #     #   setTxtProgressBar(pb, doms)
          #     dom <- d
          # 
          #     txt=paste("domain <- dataset[dataset$DOMAINVALUE", idx_d,"== dom,]", sep="")
          #     eval(parse(text=txt)) #select only rows of dataset having category dom of variable idx_d
          #     domain$STRATUM<- factor(domain$STRATUM)
          #     if(doms==1){
          #       stmt=paste("DOM",idx_d, " <- rep(as.character(dom), length(levels(domain$STRATUM)))", sep="")
          #       eval(parse(text=stmt))  #define DOMidx_d variable
          #       #strat <- domain$STRATUM
          #     }
          #     else{
          #     stmt <- paste("DOM",idx_d, " <- c(DOM", idx_d, ",rep(as.character(dom), length(levels(domain$STRATUM))))", sep="")
          #     eval(parse(text=stmt))  #define DOMidx_d variable
          #     #strat=c(strat, domain$STRATUM)
          #     #dom_strat=cbind()
          #     }
          # 
          #   }
          #   if(idx_d==1){
          #     dom_df=data.frame(DOM1=DOM1, stringsAsFactors = TRUE)
          #   }
          #   else{
          #   stmt <- paste("dom_df <- as.data.frame(cbind(dom_df,DOM",idx_d,"),stringsAsFactors = TRUE)", sep="" )
          #   eval(parse(text = stmt)) #create a dataframe with variable STRATUM, N, M1:Mn, S1:Sn, COST, CENS, DOMidx_d
          #   }
          # }


          
          
          stmt <- paste("strata <- as.data.frame(cbind(STRATUM=levels(STRATUM),N,", listM,",", listS,",COST,CENS, dataset_dom),stringsAsFactors = TRUE)", sep="" )
          eval(parse(text = stmt)) #create a dataframe with variable STRATUM, N, M1:Mn, S1:Sn, COST, CENS, DOMidx_d
          #add domains
          
          for (i in 1:nvarX) {
            stmt <- paste("strata$X", i, " <- rep(0, length(levels(dataset$STRATUM)))", 
                          sep = "")
            eval(parse(text = stmt)) #initialize with null vector X variables
          }
          strata$STRATUM <- as.character(strata$STRATUM) #to character
          if(length(namesX)!=1){
            for (i in 1:nrow(strata)) {
              strata[i, c(namesX)] <- unlist(strsplit(strata$STRATUM[i], 
                                                      "\\*"))
              #fill X variables with their values (taken from STRATUM)
            }
          }
          else{
            strata[, c(namesX)] <-strata$STRATUM
            # if stratification variables are not provided, then returns stratum
          }
        #}
        #for all the other domain, as I have already save the values of CENS, COST, M1:Mn, S1:Sn, X1:Xk, I just need to save the stratum and the domain categories
      #   if(idx_d!=1){
      #     stmt=paste("DOM",idx_d, " <- rep(as.character(dom), length(levels(domain$STRATUM)))", sep="")
      #     eval(parse(text=stmt)) #define DOMidx_d variable
      #     STRATUM <- domain$STRATUM
      #     stmt <- paste("strata <- as.data.frame(cbind(STRATUM=levels(STRATUM), DOM",idx_d,"),stringsAsFactors = TRUE)", sep="" )
      #     eval(parse(text = stmt)) #create a dataframe with variable STRATUM and DOMidx_d
      #   }
      #   stratatot <- rbind(stratatot, strata) #rbind strata dataframe to obtain stratatot with the all categories of domain DOMidx_d
      #   
      # }
      # stratadom[[idx_d]] <- stratatot #save 
      # if (idx_d==1){
      #   strataok <- stratatot
      #   stratatot <- NULL
      # }
      # if (idx_d!=1){ 
      #   #
      #   strataok <- merge(strataok, stratatot, by="STRATUM") #merge the datasets created for each domain DOMidx_d (idx_d:1:ndom)
      #   stratatot <- NULL
      # }
      
    #}
    
    # #define the number of domains:
    # 
    # ndom <- sum(grepl("DOMAINVALUE", colnames(dataset)))
    # 
    # 
    # for (idx_d in 1:ndom){
    #   if (verbose == TRUE) {
    #     cat(paste("\nDOM", idx_d, "\n", sep=""))
    #   }
    #   txt <- paste("dataset$DOMAINVALUE", idx_d, "<- factor(dataset$DOMAINVALUE", idx_d, ")", sep="")
    #   eval(parse(text=txt)) #domains variables are converted to factors
    #   txt2 <- paste("numdom[idx_d] <- length(levels(dataset$DOMAINVALUE", idx_d, "))", sep="")
    #   eval(parse(text=txt2)) #save the number of domain categories into the vector numdom
    # 
    #   if (progress == TRUE){
    #     pb <- txtProgressBar(min = 0, max = numdom[idx_d], style = 3)
    #   }
    #   doms <- 0
    #   #for each category d in domain idx_d:
    #   for (d in (eval(parse(text=paste("levels(dataset$DOMAINVALUE", idx_d,")", sep=""))))){
    # 
    #     if (verbose == TRUE) {
    #       cat(paste("\nDOM", idx_d, " category: ", d,"\n", sep=""))
    #     }
    # 
    #     doms <- doms + 1
    #     if (progress == TRUE){
    #       Sys.sleep(0.1)
    #       setTxtProgressBar(pb, doms)
    #     }
    # 
    #     # if (progress == TRUE)
    #     #   setTxtProgressBar(pb, doms)
    #     dom <- d
    #     txt=paste("domain <- dataset[dataset$DOMAINVALUE", idx_d,"== dom,]", sep="")
    #     eval(parse(text=txt)) #select only rows of dataset having category dom of variable idx_d
    #     listX <- NULL #save list of X variables
    #     namesX <- NULL #save names of X variables
    #     for (i in 1:nvarX) {
    #       name <- paste("X", i, sep = "")
    #       namesX <- cbind(namesX, name)
    #       if (i < nvarX)
    #         listX <- paste(listX, "domain$X", i, ",", sep = "")
    #       else listX <- paste(listX, "domain$X", i, sep = "")
    #     }
    #     listM <- NULL
    #     listS <- NULL
    #     for (i in 1:nvarY) {
    #       if (i<nvarY){
    #         listM <- paste(listM, "M", i, ",", sep = "")
    #         listS <- paste(listS, "S", i, ",", sep = "")
    #       }
    #       else{
    #         listM <- paste(listM, "M", i,  sep = "")
    #         listS <- paste(listS, "S", i,  sep = "")
    # 
    #       }
    #     }
    #     #to obtain stratum, paste the categories of the stratification variables, using as separator '*'
    #     stmt <- paste("domain$STRATUM <- as.factor(paste(", listX,
    #                   ",sep='*'))", sep = "")
    #     eval(parse(text = stmt))
    # 
    # 
    # 
    #     if(idx_d==1){
    #       if (!is.null(dataset$COST)){
    #         cost <- tapply(domain$WEIGHT * domain$COST, domain$STRATUM,sum)/tapply(domain$WEIGHT,
    #                                                                                 domain$STRATUM, sum)
    #       } #cost function to be used later if the variable COST in dataset is specified.
    # 
    # 
    # 
    #       for (i in 1:nvarY) {
    #         WEIGHT <- NULL
    #         STRATUM <- NULL
    #         Y <- NULL
    #         stmt <- paste("Y <- domain$Y", i, "[!is.na(domain$Y",
    #                       i, ")]", sep = "")
    #         eval(parse(text = stmt)) #save in Y the values of variable Yi without NA
    #         if(length(Y)==0){
    #           cat(paste("\nVariable Y", i, " in the category ", d," of domain DOM", idx_d,
    #                     " is completely missing. By default, M",i,
    #                     " and S", i,
    #                     " are set equal to 0 in all the strata \n", sep=""))
    #         }
    #         stmt <- paste("WEIGHT <- domain$WEIGHT[!is.na(domain$Y",
    #                       i, ")]", sep = "")
    #         eval(parse(text = stmt)) #save in WEIGHT the values of variable WEIGHT (corresponding to non missing Y values)
    #         stmt <- paste("STRATUM <- domain$STRATUM[!is.na(domain$Y", i, ")]", sep = "")
    #         eval(parse(text = stmt)) #save in STRATUM the values of variable STRATUM (corresponding to non missing Y values)
    #         STRATUM <- factor(STRATUM)
    # 
    # 
    #         samp <- NULL
    #         stmt <- paste("samp <- domain[!is.na(domain$Y",
    #                       i, "),]", sep = "")
    #         eval(parse(text = stmt)) #define samp as dataframe 'domain' (including only rows with non-missing values of the variable Yi)
    #         l.split <- split(samp, samp$STRATUM, drop = TRUE) #split divides the data samp into the groups defined by STRATUM.
    #         stmt <- paste("M", i, " <- sapply(l.split, function(df,x,w) ",  avg, "(df[,x],df[,w]), x='Y", i, "', w='WEIGHT')",
    #         sep = "")
    #         eval(parse(text = stmt))  # separately per STRATUM, save the mean value of the variable Yi
    # 
    #         stmt <- paste("S", i, " <- sapply(l.split, function(df,x,w) ",  stdev[[i]], "(df[,x],df[,w]), x='Y", i, "', w='WEIGHT')",
    #                       sep = "")
    #         eval(parse(text = stmt)) # separately per STRATUM, save the standard deviation value of the variable Yi
    # 
    #         stmt <- paste("stratirid <- names(M", i,")", sep="")
    #         eval(parse(text = stmt)) #save the identifier of the stratum (only stratum corresponding to non-missing values of Yi)
    #         strati <- data.frame(X1 = levels(domain$STRATUM),
    #                              stringsAsFactors = TRUE) #create a dataframe with variable X1: STRATUM
    #         stmt <- paste("m <- data.frame(cbind(X1=stratirid,X2=M",
    #                       i, "), stringsAsFactors = TRUE)", sep = "")
    #         eval(parse(text = stmt)) #create a dataframe with variable x1: STRATUM, x2: mean of the target variable Yi
    #         m <- merge(strati, m, by = c("X1"), all = TRUE) #merge by strato
    # 
    #         m$X2 <- ifelse(is.na(m$X2), m$X2, as.character(m$X2)) #convert X2 as character
    #         m$X2 <- ifelse(is.na(m$X2), m$X2, as.numeric(m$X2))  #convert X2 as numeric
    #         m$X2 <- ifelse(is.na(m$X2), 0, m$X2) #replace NA in X2 with 0
    #         stmt <- paste("M", i, " <- m$X2", sep = "") # update Mi
    #         eval(parse(text = stmt))
    #         stmt <- paste("s <- data.frame(cbind(X1=stratirid,X2=S",
    #                       i, "), stringsAsFactors = TRUE)", sep = "")
    #         eval(parse(text = stmt))  #create a dataframe with variable x1: STRATUM, x2: standard deviation of the target variable Yi
    #         s <- merge(strati, s, by = c("X1"), all = TRUE) #merge by stratum
    #         s$X2 <- ifelse(is.na(s$X2), s$X2, as.character(s$X2)) #convert X2 as character
    #         s$X2 <- ifelse(is.na(s$X2), s$X2, as.numeric(s$X2)) #convert X2 as numeric
    #         s$X2 <- ifelse(is.na(s$X2), 0, s$X2) #replace NA in X2 with 0
    #         stmt <- paste("S", i, " <- s$X2", sep = "")
    #         eval(parse(text = stmt)) # update Si
    #       }
    # 
    #       N <- round(tapply(domain$WEIGHT, domain$STRATUM, sum))
    #       #save number of units per stratum
    #       STRATUM <- domain$STRATUM
    #       #define COST variable
    #       if (is.null(dataset$COST)){
    #         COST <- rep(1, length(levels(domain$STRATUM)))
    #       }
    #       if (!is.null(dataset$COST)){
    #         COST <- cost
    #       }
    # 
    # 
    #       #define CENS variable:
    #       if(!is.null(dataset$CENS)){
    #         CENS <- domain$CENS
    #       } else CENS <- rep(0, length(levels(domain$STRATUM)))
    # 
    # 
    #       stmt <- paste("DOM",idx_d, " <- rep(as.character(dom), length(levels(domain$STRATUM)))", sep="")
    #       eval(parse(text=stmt))  #define DOMidx_d variable
    # 
    # 
    #       stmt <- paste("strata <- as.data.frame(cbind(STRATUM=levels(STRATUM),N,", listM,",", listS,",COST,CENS,DOM",idx_d,"),stringsAsFactors = TRUE)", sep="" )
    #       eval(parse(text = stmt)) #create a dataframe with variable STRATUM, N, M1:Mn, S1:Sn, COST, CENS, DOMidx_d
    #       for (i in 1:nvarX) {
    #         stmt <- paste("strata$X", i, " <- rep(0, length(levels(domain$STRATUM)))",
    #                       sep = "")
    #         eval(parse(text = stmt)) #initialize with null vector X variables
    #       }
    #       strata$STRATUM <- as.character(strata$STRATUM) #to character
    #       if(length(namesX)!=1){
    #         for (i in 1:nrow(strata)) {
    #           strata[i, c(namesX)] <- unlist(strsplit(strata$STRATUM[i],
    #                                                   "\\*"))
    #           #fill X variables with their values (taken from STRATUM)
    #         }
    #       }
    #       else{
    #         strata[, c(namesX)] <-strata$STRATUM
    #         # if stratification variables are not provided, then returns stratum
    #       }
    #     }
    #     #for all the other domain, as I have already save the values of CENS, COST, M1:Mn, S1:Sn, X1:Xk, I just need to save the stratum and the domain categories
    #     if(idx_d!=1){
    #       stmt=paste("DOM",idx_d, " <- rep(as.character(dom), length(levels(domain$STRATUM)))", sep="")
    #       eval(parse(text=stmt)) #define DOMidx_d variable
    #       STRATUM <- domain$STRATUM
    #       stmt <- paste("strata <- as.data.frame(cbind(STRATUM=levels(STRATUM), DOM",idx_d,"),stringsAsFactors = TRUE)", sep="" )
    #       eval(parse(text = stmt)) #create a dataframe with variable STRATUM and DOMidx_d
    #     }
    #     stratatot <- rbind(stratatot, strata) #rbind strata dataframe to obtain stratatot with the all categories of domain DOMidx_d
    # 
    #   }
    #   stratadom[[idx_d]] <- stratatot #save
    #   if (idx_d==1){
    #     strataok <- stratatot
    #     stratatot <- NULL
    #   }
    #   if (idx_d!=1){
    #     #
    #     strataok <- merge(strataok, stratatot, by="STRATUM") #merge the datasets created for each domain DOMidx_d (idx_d:1:ndom)
    #     stratatot <- NULL
    #   }
    # 
    # }
    # strataok=strata
    # stratatot <- strataok
    stratatot=strata
    #strataok <- NULL
    strata <- NULL
    colnames(stratatot) <- toupper(colnames(stratatot)) #to upper
    
    for (i in 1:ndom){
      txt <- paste("stratatot$DOM", i, "<- as.factor(stratatot$DOM", i, ")", sep="")
      eval(parse(text=txt)) #convert DOMs variables to factors
    }
    options(scipen = 100)
    for (i in 1:nvarY){
      stmt <- paste("stratatot$M", i, " <- as.numeric(as.character(stratatot$M", i, 
                    "))", sep = "")
      eval(parse(text = stmt)) #convert Mi to numeric
      
      stmt <- paste("stratatot$S", i, " <- as.numeric(as.character(stratatot$S", i, 
                    "))", sep = "")
      eval(parse(text = stmt)) #convert Si to numeric
    }
    stratatot$N <- as.numeric(as.character(stratatot$N))
    for (j in (1:nrow(stratatot))) {
      stmt <- paste("stratatot$M", i, "[j] <- ifelse(stratatot$M", 
                    i, "[j] == 0,0.000000000000001,stratatot$M", i, 
                    "[j])", sep = "")
      eval(parse(text = stmt)) # replace 0 in Mi with 0.000000000000001
    }
    if (verbose == TRUE) {
      cat("\nNumber of strata: ", nrow(stratatot))
      cat("\n... of which with only one unit: ", sum(stratatot$N == 1))
    }
    stratatot <- stratatot[,c("STRATUM", paste("X", 1:nvarX, sep=""),
                           paste("DOM", 1:ndom, sep=""), 
                           "N", paste("M", 1:nvarY, sep=""), 
                           paste("S", 1:nvarY, sep=""), 
                           "CENS", "COST")] #reorder columns
    return(list(stratatot=stratatot, numdom=ndom, nvarX=nvarX, nvarY=nvarY))
    
  }
 
  #library(SamplingStrata)
  # if (is.null(samp_frame$one)) 
  #   samp_frame$one <- 1
  cat("\nCalculating strata...")
  frame <- buildFrame(df = samp_frame,
                        id = ID,
                        X = stratum, 
                        Y = target,
                        domainvalue = dom)
 # nvarY <- length(grep("Y", colnames(frame)))
  stratif <- buildStrata(dataset=frame, progress = TRUE, verbose=TRUE)
  nvarX <- stratif$nvarX
  nvarY <- stratif$nvarY
  numdom <- stratif$numdom
  stratif <- stratif$strata
  #append also national domain
  for (i in seq(numdom, 1, by=-1)){
    txt <- paste("stratif$DOM", i+1, "<- stratif$DOM", i, sep="")
    eval(parse(text=txt))
  }
  stratif$DOM1 <- 1 #national domain
  stratif <- stratif[,c("STRATUM", paste("X", 1:nvarX, sep=""),
                         paste("DOM", 1:(numdom+1), sep=""), 
                         "N", paste("M", 1:nvarY, sep=""), 
                         paste("S", 1:nvarY, sep=""), 
                         "CENS", "COST")] #reorder columns
  
  #adjust colnames
  if(length(grepl("X", names(stratif))) > 1){
    pos <- which(grepl("X", colnames(stratif)))
    colnames(stratif)[pos] <- stratum
  }
  
  stratif$CENS <- as.numeric(stratif$CENS)
  stratif$CENS <- 0
  stratif$COST <- as.numeric(stratif$COST)
  
  
  
 
  return(stratif)

}

    
 