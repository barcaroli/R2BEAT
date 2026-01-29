\name{pareto}
\alias{pareto}
\title{Pareto function}
\usage{
pareto(strata,current_cvs,target_size,tolerance,cv_caps,t0,inactive_cv,
                   min_cv,eps_inact,bracket_fac,max_iter,max_same_iter,minnumstrat,
                   beat1cv_fun,verbose,plot,plot_dir,plot_prefix,save_pdf,show_targets,     
                   show_caps,plot_convergence)
}
\arguments{
\item{strata}{Data frame with strata information}
\item{current_cvs}{Data frame with initial CV targets per domain.}
\item{target_size}{Desired total sample size.}
\item{tolerance}{Allowed shortfall from target_size (scalar, >= 0).}    
\item{cv_caps}{Optional data frame with columns DOM, VAR, MAX_CV to cap CVs by domain/variable. VAR can be V#, CV#, or integer index.}     
\item{t0}{Initial scaling factor for CV targets.}    
\item{inactive_cv}{CV values >= inactive_cv are treated as inactive targets.} 
\item{min_cv}{Lower bound for active CV targets after scaling.} 
\item{bracket_fac}{Multiplicative step used to expand the bracketing interval.} 
\item{max_iter}{Maximum bisection iterations.} 
\item{max_same_iter}{Stop bracketing if n(t) repeats this many times (caps plateau).} 
\item{minnumstrat}{Minimum sample size per stratum passed to beat1cv_fun.}
\item{beat1cv_fun}{Optional function to compute expected CVs (e.g., beat.1cv_new).}
\item{verbose}{If TRUE, print progress messages.}
\item{plot}{If TRUE, plot CVs by domain and allocations.}
\item{plot_dir}{Output directory for plots.}
\item{plot_prefix}{Filename prefix for plot outputs.}
\item{save_pdf}{If TRUE, save plots to a PDF file.}
\item{show_targets}{If TRUE, overlay initial target CVs on plots.}
\item{show_caps}{If TRUE, overlay CV caps on plots (when provided).}
\item{plot_convergence}{If TRUE, create convergence plot showing bracketing and bisection}
}
\examples{
out <- pareto(
  strata = strata,
  current_cvs = cv_equal,
  target_size = target_value,
  tolerance = 5,
  cv_caps = caps,
  beat1cv_fun = beat.1cv_new,
  max_iter = 50,
  max_same_iter = 20,
  plot = TRUE,
  plot_dir = "outputs",
  plot_prefix = "italy_domains",
  show_targets = TRUE,     
  show_caps = TRUE      
)

}
