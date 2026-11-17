#' Launch the R2BEAT Shiny application
#'
#' Starts the interactive R2BEAT application, which gives access to three
#' workflows: CV optimisation for one-stage allocation (Pareto + genetic
#' algorithm), two-stage design starting from a sampling frame, and two-stage
#' design starting from a previous survey (ReGenesees objects).
#'
#' If any of the packages required by the app (shiny, shinyjs, shinydashboard,
#' DT, ggplot2, genalg) is not installed, it is installed from CRAN first.
#'
#' @param launch.browser Passed to \code{shiny::runApp}.
#' @param ... Further arguments passed to \code{shiny::runApp}.
#' @return Called for its side effect (runs the app).
#' @export
runR2BEATapp <- function(launch.browser = TRUE, ...) {
  pkgs <- c("shiny", "shinyjs", "shinydashboard", "DT", "ggplot2", "genalg")
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss) > 0) {
    message("Installing the packages required by the app: ",
            paste(miss, collapse = ", "))
    repos <- getOption("repos")
    if (is.null(repos) || identical(unname(repos["CRAN"]), "@CRAN@"))
      repos <- c(CRAN = "https://cloud.r-project.org")
    utils::install.packages(miss, repos = repos)
    still <- miss[!vapply(miss, requireNamespace, logical(1), quietly = TRUE)]
    if (length(still) > 0)
      stop("Could not install the following packages: ",
           paste(still, collapse = ", "), call. = FALSE)
  }
  appDir <- system.file("shiny", "R2BEAT_app", package = "R2BEAT")
  if (appDir == "") stop("Could not find the app directory. Re-install R2BEAT.", call. = FALSE)
  shiny::runApp(appDir, launch.browser = launch.browser, ...)
}
