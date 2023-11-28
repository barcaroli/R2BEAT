.onAttach<-function(libname, pkgname){
  if (!"ReGenesees" %in% rownames(installed.packages())) {
    packageStartupMessage("\n", appendLF = FALSE)
    packageStartupMessage("Please install package ReGenesees, as it is required by R2BEAT, by executing: ", appendLF = FALSE)
    packageStartupMessage("\n", appendLF = FALSE)
    packageStartupMessage("devtools::install_github('DiegoZardetto/ReGenesees')", appendLF = FALSE)
  }
  packageStartupMessage("\n", appendLF = FALSE)
  packageStartupMessage("Report issues at https://github.com/barcaroli/R2BEAT/issues", appendLF = FALSE)
  packageStartupMessage("\n", appendLF = FALSE)
  packageStartupMessage("Get a complete documentation on https://barcaroli.github.io/R2BEAT/", appendLF = FALSE)
  packageStartupMessage("\n", appendLF = FALSE)
}


      