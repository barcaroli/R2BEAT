sensitivity <- function (samp_frame, errors, id_PSU, id_SSU, strata_var, target_vars, 
                         deff_var, domain_var, minimum = 50, delta = 1, f = 0.05, 
                         deff_sugg = 1.5, search = c("deff", "min_SSU", "sample_fraction"), 
                         min, max, plot=TRUE) 
{
  if (!(search %in% c("deff", "min_SSU", "sample_fraction"))) {
    stop("Value for search grid not acceptable: must be one of 'deff','min_SSU','sample_fraction'")
  }
  grid <- seq(min, max, (max - min)/10)
  k <- 0
  PSU <- rep(NA, 10)
  SSU <- rep(NA, 10)
  for (i in grid) {
    k <- k + 1
    cat("\n", k)
    inp <- prepareInputToAllocation(samp_frame, id_PSU, 
                                    id_SSU, strata_var, target_vars, deff_var, domain_var, 
                                    minimum, delta, f, deff_sugg = i)
    alloc <- beat.2st(stratif = inp$strata, errors = errors, 
                      des_file = inp$des_file, psu_file = inp$psu_file, 
                      rho = inp$rho, deft_start = NULL, effst = inp$effst, 
                      epsilon1 = 5, mmdiff_deft = 1, maxi = 15, epsilon = 10^(-11), 
                      minnumstrat = 2, maxiter = 200, maxiter1 = 25)
    PSU[k] <- alloc$iterations$`PSU Total`[length(alloc$iterations$iterations)]
    SSU[k] <- alloc$iterations$SSU[length(alloc$iterations$iterations)]
  }
  two_stage_allocation <- list(PSU, SSU)
  if (plot==TRUE) plot.sens(two_stage_allocation,search,min,max)
  return(two_stage_allocation)
}
