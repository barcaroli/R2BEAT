#' adjustCVs
#'  
#' @description
#' Function to increase or decrease the precision constraints in order to obtain the desired sample size 
#'  
#' 
#' @param target_size desired sample size.
#' @param strata the 'strata' dataset.
#' @param errors the 'errors' dataset containing the current precision constraints
#' @param adj_rate the rate of adjustment (default=0.01): 
#' the smaller, the higher the precision in reaching the target size; the higher, the quicker is the adjustment
#' 
#' @return the new 'errors' dataset containing the modified precision constraints
#' 
#' @examples
#' data(beat.example)
#' errors
#' a <- beat.1st(strata,errors)
#' sum(a$alloc$ALLOC[-nrow(a$alloc)])
#' errors_new <- adjust_CVs(9000,strata,errors,adj_rate=0.005)
#' errors_new

adjust_CVs <- function(target_size,strata,errors,adj_rate=0.01) {
  a <- beat.1st(stratif=strata,errors=errors)
  cvnew <- errors
  size <- sum(a$alloc$ALLOC[-nrow(a$alloc)])
  if (size < target_size) {
    repeat {
      for (k in c(2:ncol(cvnew))) {
        cvnew[,k] <- cvnew[,k] - 0.01*cvnew[,k]
      }
      b <- beat.1st(stratif=strata,errors=cvnew)
      if (sum(b$alloc$ALLOC[-nrow(b$alloc)]) > target_size) break
    } 
  }
  if (size >= target_size) {
    repeat {
      for (k in c(2:ncol(cvnew))) {
        cvnew[,k] <- cvnew[,k] + 0.01*cvnew[,k]
      }
      b <- beat.1st(stratif=strata,errors=cvnew)
      if (sum(b$alloc$ALLOC[-nrow(b$alloc)]) < target_size) break
    } 
  }
  cat("\n Size: ",sum(b$alloc$ALLOC[-nrow(b$alloc)]))
  return(cvnew)
}