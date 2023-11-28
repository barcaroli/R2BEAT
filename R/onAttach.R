.onAttach<-function(libname, pkgname){
  if (!"ReGenesees" %in% rownames(installed.packages())) {
    packageStartupMessage("\n", appendLF = FALSE)
    packageStartupMessage("Installing ReGenesees as it is required by R2BEAT", appendLF = FALSE)
    packageStartupMessage("\n", appendLF = FALSE)
    if (!"devtools" %in% rownames(installed.packages())) {
      install.packages("devtools")
    }
    devtools::install_github('DiegoZardetto/ReGenesees')
  }
  library(ReGenesees)
  packageStartupMessage("\n", appendLF = FALSE)
  packageStartupMessage("Report issues at https://github.com/barcaroli/R2BEAT/issues", appendLF = FALSE)
  packageStartupMessage("\n", appendLF = FALSE)
  packageStartupMessage("Get a complete documentation on https://barcaroli.github.io/R2BEAT/", appendLF = FALSE)
  packageStartupMessage("\n", appendLF = FALSE)
}


      