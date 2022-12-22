build_dummy_variables <- function(frame,domain_var,initial_target_vars,cv) {
  new_target_vars <- initial_target_vars
  eval(parse(text=paste0("mod_mat <- as.data.frame(model.matrix(~",domain_var,"-1,data=pop))")))
  for (k in c(1:length(initial_target_vars))) {
    tarvar <- initial_target_vars[k]
    eval(parse(text=paste0(tarvar," <- mod_mat * frame$",tarvar)))
    cat("\n Variable: ",tarvar)
    eval(parse(text=paste0(tarvar," <- mod_mat * frame$",tarvar)))
    st1 <- paste0("colnames(",tarvar,") <- paste('") 
    st2 <- paste0(tarvar,"',colnames(",tarvar,"),sep='_')")
    st <- paste0(st1,st2)
    eval(parse(text=st))
    eval(parse(text=paste0("new <- colnames(",tarvar,")")))
    new_target_vars <- c(new_target_vars,new)
    st <- paste0("frame <- cbind(frame,",tarvar,")")
    eval(parse(text=st))
  }
  cv2 <- as.data.frame(list(DOM = c("DOM1","DOM2")))
  for (i in c(1:length(new_target_vars))) {
    st <- paste0("cv2$CV",i," <- cv$CV1")
    eval(parse(text=st))
  }
  cv2[1,c((length(initial_target_vars)+1):(length(new_target_vars)+1))] <- 1.0
  cv2[2,c(2:(length(initial_target_vars)+1))] <- 1.0
  out <- list(frame = frame, new_target_vars = new_target_vars, cvnew = cv2)
}