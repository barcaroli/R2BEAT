plot_sens <- function(x,min,max) {
  grid <- seq(min, max, (max - min)/10)
  par(mfrow = c(1, 2))
  plot(grid, x[[1]], type = "l", ylab = "Number of PSUs", 
       xlab = "Varying minimum SSUs per PSU")
  title("First Stage Size")
  plot(grid, x[[2]], type = "l", ylab = "Number of SSUs", 
       xlab = "Varying minimum SSUs per PSU")
  title("Second Stage Size")
  par(mfrow = c(1, 1))
}

