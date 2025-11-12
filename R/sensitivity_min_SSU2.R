sensitivity_min_SSU2 <- function (strata, 
                         des_file,
                         psu_file,
                         rho,
                         effst,
                         errors, 
                         min, 
                         max, 
                         plot = TRUE) 
{
  grid <- seq(min, max, (max - min)/10)
  k <- 0
  PSU <- rep(NA, 10)
  SSU <- rep(NA, 10)
  minimum <- grid[1]

  k <- 0
  for (i in grid) {
    k <- k + 1
    cat("\n", k)
    des_file$MINIMUM = i
    alloc <- beat.2st(file_strata = strata, 
                      errors = errors, 
                      des_file = des_file, 
                      psu_file = psu_file, 
                      rho = rho, 
                      deft_start = NULL, 
                      effst = effst, 
                      # epsilon1 = 5, mmdiff_deft = 1, maxi = 15, epsilon = 10^(-11), 
                      minPSUstrat = 2, 
                      minnumstrat = 2)
    PSU[k] <- alloc$iterations$`PSU Total`[length(alloc$iterations$iterations)]
    SSU[k] <- alloc$iterations$SSU[length(alloc$iterations$iterations)]
  }
  two_stage_allocation <- list(PSU, SSU)
  if (plot == TRUE) 
    plot_sens(two_stage_allocation, min, max)
  return(two_stage_allocation)
}
