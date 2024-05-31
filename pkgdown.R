require(devtools)
# use_readme_rmd()
# use_news_md()
#use_vignette("R2BEAT")
#use_github_links()
# use_travis()
# use_cran_badge()

# install.packages("tvthemes") # v1.1.0
# library(tvthemes)


# devtools::install_github("hadley/pkgdown", version=1.6.1)
setwd("D:/Google Drive/Sampling/R2BEAT/R2BEAT_1.0.6/")
# usethis::use_github_action("pkgdown")
# setwd("C:\\Users\\UTENTE\\Google Drive\\Sampling\\R2BEAT\\R2BEAT_1.0.5")
library(pkgdown)
# usethis::use_pkgdown()
# Sys.setenv(R_WIN_USE_SAFE_CALLBACKS = "0")
init_site(pkg = ".")
devtools::build_site()




