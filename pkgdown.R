# require(devtools)
# use_readme_rmd()
# use_news_md()
# use_vignette("R2BEAT")
# use_github_links()
# use_travis()
# use_cran_badge()

# install.packages("tvthemes") # v1.1.0
# library(tvthemes)


# devtools::install_github("hadley/pkgdown")
# setwd("C:/Users/UTENTE/Google Drive/Sampling/R2BEAT/R2BEAT_1.0.34/R2BEAT")
setwd("D:\\Google Drive\\Sampling\\R2BEAT\\R2BEAT_1.0.4")

library(pkgdown)
# usethis::use_pkgdown()
# build_favicon(pkg = ".")
# init_site(pkg = ".")
# build_site()
init_site()
build_home()
build_reference()
build_articles()
# build_tutorials()
build_news()
build_reference_index(pkg = "D:\\Documenti\\R\\win-library\\4.0\\R2BEAT")
#build_reference_index(pkg = "C:\\Users\\UTENTE\\Documents\\R\\win-library")

