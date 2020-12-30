
<!-- README.md is generated from README.Rmd. Please edit README.Rmd file -->

[![License](http://img.shields.io/badge/license-GPL%20%28%3E=%202%29-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-3.0.html)
[![CRAN](http://www.r-pkg.org/badges/version/R2BEAT)](https://cran.r-project.org/package=R2BEAT)
[![Downloads](http://cranlogs.r-pkg.org/badges/R2BEAT?color=brightgreen)](http://www.r-pkg.org/pkg/R2BEAT)

# R2BEAT: Multistage Sampling Design and PSUs selection

Multivariate optimal allocation for different domains in one and two
stages stratified sample design. R2BEAT extends the Neyman (1934) –
Tschuprow (1923) allocation method to the case of several variables,
adopting a generalization of the Bethel’s proposal (1989). R2BEAT
develops this methodology but, moreover, it allows to determine the
sample allocation in the multivariate and multi-domains case of
estimates for two-stage stratified samples. It also allows to perform
Primary Stage Units selection.

## Installation

You can install the released version of R2BEAT from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("R2BEAT")
```

or the last version of R2BEAT from github with:

``` r
install.packages("devtools")
devtools::install_github("barcaroli/R2BEAT")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(R2BEAT)
# Load example data
data(beat.example)
head(strata)
#   STRATUM       N DOM1  DOM2         M1         M2 COST CENS        S1        S2
# 1       1 2302353 REG1 PROV1 0.02097409 0.01732549    1    0 0.1432975 0.1304811
# 2       2  179562 REG1 PROV2 0.02304681 0.01997698    1    0 0.1500522 0.1399210
# 3       3  371802 REG2 PROV3 0.02284502 0.01894097    1    0 0.1494093 0.1363166
# 4       4  592303 REG2 PROV3 0.02375591 0.01077033    1    0 0.1522878 0.1032198
# 5       5  221687 REG2 PROV3 0.03215470 0.01088398    1    0 0.1764108 0.1037570
# 6       6  440613 REG2 PROV3 0.04496313 0.01967690    1    0 0.2072232 0.1388874
head(errors)
#    DOM  CV1  CV2
# 1 DOM1 0.10 0.99
# 2 DOM2 0.99 0.99
head(design)
#   STRATUM STRAT_MOS DELTA MINIMUM
# 1       1   2302353     1      48
# 2       2    179562     1      48
# 3       3    371802     1      48
# 4       4    592303     1      48
# 5       5    221687     1      48
# 6       6    440613     1      48
head(PSU_strat)
#   STRATUM PSU_MOS PSU_ID
# 1       1    2591      1
# 2       1    3808      2
# 3       1     465      3
# 4       1    1778      4
# 5       1     713      5
# 6       1    6378      6
head(rho)
#   STRATUM RHO_AR1   RHO_NAR1 RHO_AR2     RHO_NAR2
# 1       1       1 0.07563272       1 0.0005484989
# 2       2       1 0.08327200       1 0.0157297282
# 3       3       1 0.03444641       1 0.0001944330
# 4       4       1 0.05433433       1 0.0021201253
# 5       5       1 0.03926152       1 0.0020177360
# 6       6       1 0.02452771       1 0.0293768074
# Allocate the sample
allocation <- beat.2st(stratif=strata, 
                            errors=errors,
                            des_file=design, 
                            psu_file=PSU_strat,
                            rho=rho)
#    iterations PSU_SR PSU NSR PSU Total   SSU
# 1           0      0       0         0  6858
# 2           1     16      79        95 17170
# 3           2     42     161       203 14736
# 4           3     33     153       186 15273
# 5           4     36     155       191 15147
# 6           5     33     156       189 15274
# 7           6     36     155       191 15146
# 8           7     33     156       189 15273
# 9           8     36     155       191 15147
# 10          9     33     156       189 15274
# 11         10     36     155       191 15146
# 12         11     33     156       189 15273
# 13         12     36     155       191 15147
# 14         13     33     156       189 15274
# 15         14     36     155       191 15146
# 16         15     33     156       189 15273
# 17         16     36     155       191 15147
# 18         17     33     156       189 15274
# 19         18     36     155       191 15146
# 20         19     33     156       189 15273
# 21         20     36     155       191 15147
str(allocation)
# Selection of PSUs
sample_PSU <- StratSel(dataPop= PSU_strat,
                       idpsu= ~ PSU_ID, 
                       dom= ~ STRATUM, 
                       final_pop= ~ PSU_MOS, 
                       size= ~ PSU_MOS, 
                       PSUsamplestratum= 1, 
                       min_sample= 2, 
                       min_sample_index= FALSE, 
                       dataAll=allocation$alloc,
                       domAll= ~ factor(STRATUM), 
                       f_sample= ~ ALLOC, 
                       planned_min_sample= NULL, 
                       launch= F)
```
