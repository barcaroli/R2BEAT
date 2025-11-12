#' The function prepares the main input for the function *beat.1st* in R2BEAT R-package. It returns a list with the following objects: a dataframe, *file_strata*, with the strata information that is one of the inputs for performing the multivariate and multidomainain optimal allocation implemented in the function *beat.1st*; a vector, *var_label*, listing the target variable(s) which will be used for leading the optimal allocation.
#'
#' @param frame: A dataframe of population data (i.e. population frame or register) containing necessarily the identifier of the units, strata and domainain variables and optionally the target variable(s). 
#' @param ID: Name of the identifier of the units in the sampling frame.
#' @param stratum: Either name of the variable in data which is taken as the stratum, or the name of the variables which have to be concatenated to obtain the stratum. In this latter case, the variables used to build the stratum are retained.
#' @param domain: Name of the variable(s) in data which are the domain(s) of the estimates. Domain(s) are aggregations of strata for which the efficiency of the estimates of the target variable(s) will be kept under control when determining the optimal allocation. 
#' @param target: Name of variable(s) in data which will be the variables leading the optimal allocation, namely the efficiency of their estimates at domain(s) level will represent a constraint in the optimal allocation. Categorical variables must be encoded as factor previously.
#' @param sample: A dataframe of sample survey data, containing necessarily the strata and domains variables, the target variable(s) and the sampling weights. In this way, statistical summaries of the target variables will be estimated on the sample. Strata and domains variables must be consistent with those defined for frame dataframe. Default is NULL, meaning that population data are used.
#' @param weights: Name of the variable specifying the sampling weights whether sample survey data are used. Default is NULL, meaning that population data are used.
#' @return It returns a list with the following objects: a dataframe, *file_strata*, with the strata information that is one of the inputs for performing the multivariate and multidomain optimal allocation implemented in the function *beat.1st* in R2BEAT R-package; a vector, *var_label*, listing the target variable(s) which will be used for leading the optimal allocation.
prepareInputToAllocation_beat.1st <- function(frame, ID, stratum, domain, target, sample=NULL, weights=NULL){
###
# load dependences
###
suppressWarnings({library(dplyr)})
###
# checks
###
if (!is.null(sample)) {
  cat("Processing sample data...")
  data=sample
  var_list <- unique(c(ID,stratum))
  for (i in (1:length(var_list))) {
    var <- var_list[i]
    if (!var %in% colnames(frame)) {
      msg = paste(var, " is not in the frame. Please add it", sep="")
      stop(msg)
    }
  }
  var_list <- unique(c(stratum,domain,target,weights))
  for (i in (1:length(var_list))) {
    var <- var_list[i]
    if (!var %in% colnames(sample)) {
      msg = paste(var, " is not in the sample dataframe. Please add it", sep="")
      stop(msg)
    }
  }  
}else{
  cat("Processing population data...")
  data=frame
  var_list <- unique(c(ID,stratum,domain,target))
  for (i in (1:length(var_list))) {
    var <- var_list[i]
    if (!var %in% colnames(frame)) {
      msg = paste(var, " is not in the frame. Please add it", sep="")
      stop(msg)
    }
  }
  }
data <- data[,colnames(data) %in% c(stratum,domain,target,weights)]
###
# target variables classification
###
type <- c() 
for (i in 1:length(target)){
  if(is.null(weights) & is.numeric(data[,target[i]])){
    type[i] <- "numericPop"}
  if(!is.null(weights) & is.numeric(data[,target[i]])){
    type[i] <- "numericSamp"}
  if(!is.numeric(data[,target[i]]) & nlevels(data[,target[i]])==2){
    if(prod(levels(data[,target[i]])==c("0","1"))==1){
      type[i] <- "factor2"
      data[,target[i]] <- as.numeric(data[,target[i]])-1}
   if(prod(levels(data[,target[i]])==c("0","1"))==0){
     type[i] <- "factorP"}}
  if(!is.numeric(data[,target[i]]) & nlevels(data[,target[i]])>2){
    type[i] <- "factorP"}
}
# dichotomization of polytomous variables
if(sum(grepl("factorP", type))>0){
     Mfactor <- target[grepl("factorP", type)]
     target <- target[-which(grepl("factorP", type))]
     type <- type[-which(grepl("factorP", type))]
       for (i in 1:length(Mfactor)){
         adj <- as.matrix(fastDummies::dummy_cols(data.frame(x=data[,Mfactor[i]]), 
                          remove_selected_columns = TRUE, 
                          remove_first_dummy = FALSE))
         varname <- Mfactor[i]
         colnames(adj) <- paste(varname,levels(data[,Mfactor[i]]),sep="_")
         y <- colnames(adj)
         data <- cbind(data,adj)
         target <- c(target,y)
         type <- c(type,rep("factor2",length(y)))}
    }
###
# file_strata dataframe
###
data$WEIGHTS <- 1
if(!is.null(weights)){
  data$WEIGHTS <- data[[weights]]
}
# STRATUM variable
if(length(stratum)==1){
  data$STRATUM <- as.factor(data[[stratum]])
  } else {
  data$STRATUM <- as.factor(do.call(paste, c(data[stratum], sep = "_")))}
# population size (N) from sample frame
if(!is.null(sample)){
  frame$N=1
  # STRATUM variable
  if(length(stratum)==1){
    frame$STRATUM <- as.factor(frame[[stratum]])
  } else {
    frame$STRATUM <- as.factor(do.call(paste, c(frame[stratum], sep = "_")))}
}
# strata means
file_strata <- data.frame(
           data %>%
           group_by(across(all_of(c("STRATUM",stratum)))) %>%
           summarise(across(all_of(target),
               ~ weighted.mean(.x, !!sym("WEIGHTS"), na.rm = TRUE),
           .names = "{.col}"),
           .groups = "drop")
           )
colnames(file_strata)[colnames(file_strata) %in% target] <- paste("M",1:length(target),sep="")
# domains variables
data$DOM1 <- as.factor("Total")
DOM <-  eval(parse(text=paste("aggregate(WEIGHTS~STRATUM+DOM1+",paste(domain,collapse="+"),",data,sum)",sep="")))
colnames(DOM)[-1] <- c(paste("DOM",1:(length(domain)+1),sep=""),"N") 
# population size (N) from sample frame
if(!is.null(sample)){
  N_frame <- aggregate(N~STRATUM,frame,sum)
  DOM <- merge(DOM[,-ncol(DOM)], N_frame, by="STRATUM", all.x=TRUE)
}
if (nrow(file_strata)!=nrow(DOM)) {
  msg <- "Domains cannot split strata. Check the domain of interest (domain) and the strata variables."
  stop(msg)
}
ID_stratum <- frame[,which(colnames(frame) %in% c(ID,"STRATUM"))]
file_strata <- merge(file_strata,DOM,by="STRATUM",all.x=TRUE)
# population standard deviations by strata
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
S <- data.frame(STRATUM=file_strata$STRATUM)
for (i in 1:length(target)){
  if(type[i]=="numericSamp"){
    s <- data %>%
      group_by(across(all_of(c("STRATUM")))) %>%
      summarise(
        s = stdev1(.data[[target[i]]], !!sym("WEIGHTS")),
        .groups = "drop"
      )
  } 
  if(type[i]=="numericPop"){
    s <- data %>%
      group_by(across(all_of(c("STRATUM")))) %>%
      summarise(
        s = stdev2(.data[[target[i]]], !!sym("WEIGHTS")),
        .groups = "drop"
      )
  } 
  if(type[i]=="factor2"){
  s <- data %>%
    group_by(across(all_of(c("STRATUM")))) %>%
    summarise(
      s = stdev_prop(.data[[target[i]]], !!sym("WEIGHTS")),
      .groups = "drop"
    )
} 
colnames(s)[2] <- paste("S",i,sep="")
S <- merge(S,s,by="STRATUM",all.x=TRUE)
}
file_strata <- merge(file_strata,S,by="STRATUM",all.x=TRUE)  
file_strata <- file_strata[,unique(c(which(colnames(file_strata) %in% c("STRATUM",stratum)),
                            which(startsWith(colnames(file_strata),prefix="DOM")),
                            which(colnames(file_strata)=="N"),
                            which(startsWith(colnames(file_strata),prefix="M")),
                            which(startsWith(colnames(file_strata),prefix="S"))))]
# CENS variable by default
file_strata$CENS <- 0
# COST variable by default
file_strata$COST <- 1
###
# output preparation
###
out <- list(file_strata=file_strata,
            var_label=paste(1:length(target),target,sep="."),
            ID_stratum=ID_stratum)
return(out)
}



