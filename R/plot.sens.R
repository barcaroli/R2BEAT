plot.sens <- function(x,search,min,max) {
  grid <- seq(min, max, (max - min)/10)
  par(mfrow = c(1, 2))
  plot(grid, x[[1]], type = "l", ylab = "Number of PSUs", 
       xlab = paste0("Varying ",search))
  title("First Stage Size")
  plot(grid, x[[2]], type = "l", ylab = "Number of SSUs", 
       xlab = paste0("Varying ",search))
  title("Second Stage Size")
  par(mfrow = c(1, 1))
}

