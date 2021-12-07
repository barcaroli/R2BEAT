
<!-- README.md is generated from README.Rmd. Please edit README.Rmd file -->

[![CRAN](http://www.r-pkg.org/badges/version/R2BEAT)](https://cran.r-project.org/package=R2BEAT)
[![Downloads](http://cranlogs.r-pkg.org/badges/R2BEAT?color=brightgreen)](http://www.r-pkg.org/pkg/R2BEAT)
[![Mentioned in Awesome Official
Statistics](https://awesome.re/mentioned-badge.svg)](http://www.awesomeofficialstatistics.org)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/barcaroli/R2BEAT_workflows/HEAD?filepath=R2BEAT_workflows.ipynb)

# R2BEAT: Multistage Sampling Design and Sample Selection

Multivariate optimal allocation for different domains in one and two
stages stratified sample design. R2BEAT extends the Neyman (1934) –
Tschuprow (1923) allocation method to the case of several variables,
adopting a generalization of the Bethel’s proposal (1989). R2BEAT
develops this methodology but, moreover, it allows to determine the
sample allocation in the multivariate and multi-domains case of
estimates for two-stage stratified samples. It also allows to perform
both Primary Stage Units and Secondary Stage Units selection.

*R2BEAT* easily manages all the complexity due to the optimal sample
allocation in two-stage sampling design and provides several outputs for
evaluating the allocation. Its name stands for “R ‘to’ Bethel Extended
Allocation for Two-stage”. It is an extension of another open-source
software called *Mauss-R* (Multivariate Allocation of Units in Sampling
Surveys), implemented by ISTAT researchers
(<https://www.istat.it/en/methods-and-tools/methods-and-it-tools/design/design-tools/mauss-r>).
*Mauss-R* determines the optimal sample allocation in multivariate and
multi-domains estimation, for one-stage stratified samples.

To complete the suite of tools developed by Istat in order to cover the
stratified sample design, we cite *SamplingStrata*
(<https://CRAN.R-project.org/package=SamplingStrata>), that allows to
jointly optimize both the stratification of the sampling frame and the
allocation, still in the multivariate multidomain case (only for
one-stage designs), and *MultiWay.Sample.Allocation*, that allows to
determine the optimal sample allocation for multi-way stratified
sampling designs and incomplete stratified sampling designs
(<https://www.istat.it/en/methods-and-tools/methods-and-it-tools/design/design-tools/multiwaysampleallocation>).

For a complete illustration of the methodology see the vignette “R2BEAT
methodology and use”
(<https://barcaroli.github.io/R2BEAT/articles/R2BEAT_methodology.html>).

A complete example is illustrated in the vignette “Two-stage sampling
design workflow”
(<https://barcaroli.github.io/R2BEAT/articles/R2BEAT_workflow.html>).

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

``` r
library(R2BEAT)

#-------------------------------------------------------------------------------
# Read sampling frame
library(readr)
pop <- read_rds("https://github.com/barcaroli/R2BEAT_workflows/blob/master/pop.RDS?raw=true")
str(pop)
# 'data.frame': 2258507 obs. of  13 variables:
#   $ region       : Factor w/ 3 levels "north","center",..: 1 1 1 1 1 1 1 1 1 1 ...
# $ province     : Factor w/ 6 levels "north_1","north_2",..: 1 1 1 1 1 1 1 1 1 1 ...
# $ municipality : num  1 1 1 1 1 1 1 1 1 1 ...
# $ id_hh        : Factor w/ 963018 levels "H1","H10","H100",..: 1 1 1 2 3 3 3 3 1114 1114 ...
# $ id_ind       : int  1 2 3 4 5 6 7 8 9 10 ...
# $ stratum      : Factor w/ 24 levels "1000","2000",..: 12 12 12 12 12 12 12 12 12 12 ...
# $ stratum_label: chr  "north_1_6" "north_1_6" "north_1_6" "north_1_6" ...
# $ sex          : int  1 2 1 2 1 1 2 2 1 1 ...
# $ cl_age       : Factor w/ 8 levels "(0,14]","(14,24]",..: 3 7 8 5 4 6 6 4 4 1 ...
# $ active       : num  1 1 0 1 1 1 1 1 1 0 ...
# $ income_hh    : num  30488 30488 30488 21756 29871 ...
# $ unemployed   : num  0 0 0 0 0 0 0 0 0 0 ...
# $ inactive     : num  0 0 1 0 0 0 0 0 0 1 ...
#-------------------------------------------------------------------------------
# Precision constraints
cv <- as.data.frame(list(DOM=c("DOM1","DOM2"),
                         CV1=c(0.02,0.03),
                         CV2=c(0.03,0.06),
                         CV3=c(0.03,0.06),
                         CV4=c(0.05,0.08)))
cv
#    DOM  CV1  CV2  CV3  CV4
# 1 DOM1 0.02 0.03 0.03 0.05
# 2 DOM2 0.03 0.06 0.06 0.08
#-------------------------------------------------------------------------------
# Preparation
samp_frame <- pop
samp_frame$one <- 1
id_PSU <- "municipality"  
id_SSU <- "id_ind"        
strata_var <- "stratum"   
target_vars <- c("income_hh","active","inactive","unemployed")   
deff_var <- "stratum"     
domain_var <- "region"  
delta =  1       # households = survey units
minimum <- 50    # minimum number of SSUs to be interviewed in each selected PSU
deff_sugg <- 1.5 # suggestion for the deff value
inp <- prepareInputToAllocation1(samp_frame,
                                 id_PSU,
                                 id_SSU,
                                 strata_var,
                                 target_vars,
                                 deff_var,
                                 domain_var,
                                 minimum,
                                 delta,
                                 deff_sugg)
#-------------------------------------------------------------------------------
# Optimal allocation
alloc <- beat.2st(stratif = inp$strata, 
                  errors = cv, 
                  des_file = inp$des_file, 
                  psu_file = inp$psu_file, 
                  rho = inp$rho, 
                  deft_start = NULL,
                  effst = inp$effst, 
                  minPSUstrat = 2,
                  minnumstrat = 50)
#   iterations PSU_SR PSU NSR PSU Total  SSU
# 1          0      0       0         0 7887
# 2          1     31     104       135 8328
# 3          2     39     104       143 8317
# 4          3     38     104       142 8320
#-------------------------------------------------------------------------------
# First stage selection
sample_1st <- select_PSU(alloc, type="ALLOC", pps=TRUE, plot=TRUE)
sample_1st$PSU_stats
#    STRATUM PSU PSU_SR PSU_NSR  SSU SSU_SR SSU_NSR
# 1     1000   2      2       0  286    286       0
# 2     2000   9      3       6  452    152     300
# 3     3000   4      0       4  200      0     200
# 4     4000   2      0       2  100      0     100
# 5     5000   2      2       0  219    219       0
# 6     6000   2      0       2  100      0     100
# 7     7000   2      0       2  100      0     100
# 8     8000   2      0       2  100      0     100
# 9     9000   1      1       0  557    557       0
# 10   10000   6      6       0  587    587       0
# 11   11000  26      2      24 1300    100    1200
# 12   12000   8      0       8  400      0     400
# 13   13000   1      1       0  703    703       0
# 14   14000   4      4       0  577    577       0
# 15   15000  27      9      18 1361    461     900
# 16   16000  18      0      18  900      0     900
# 17   17000   1      1       0  154    154       0
# 18   18000   4      2       2  200    100     100
# 19   19000   7      1       6  350     50     300
# 20   20000   4      0       4  200      0     200
# 21   21000   1      1       0  125    125       0
# 22   22000   3      3       0  150    150       0
# 23   23000   4      0       4  200      0     200
# 24   24000   2      0       2  100      0     100
# 25   Total 142     38     104 9421   4221    5200
#-------------------------------------------------------------------------------
# Second stage selection
sample_2st <- select_SSU(df=pop,
                         PSU_code="municipality",
                         SSU_code="id_ind",
                         PSU_sampled=sample_1st$sample_PSU)
head(sample_2st)
#   municipality id_ind region province  id_hh stratum stratum_label sex  cl_age active income_hh
# 1           10   1570  north  north_1 H49617   12000     north_1_6   2 (64,74]      1  25393.18
# 2           10   1632  north  north_1 H49638   12000     north_1_6   1 (64,74]      1  24261.16
# 3           10   1638  north  north_1 H49639   12000     north_1_6   1 (24,34]      1  60761.11
# 4           10   1690  north  north_1 H49656   12000     north_1_6   2 (34,44]      1  40065.25
# 5           10   1709  north  north_1 H49665   12000     north_1_6   1 (44,54]      1  13320.58
# 6           10   1767  north  north_1 H49682   12000     north_1_6   1 (34,44]      1  27106.55
#   unemployed inactive  Prob_1st   Prob_2st    Prob_tot weight SR nSR stratum_2
# 1          0        0 0.2363018 0.02665245 0.006298022 158.78  0   1   12000-1
# 2          0        0 0.2363018 0.02665245 0.006298022 158.78  0   1   12000-1
# 3          0        0 0.2363018 0.02665245 0.006298022 158.78  0   1   12000-1
# 4          0        0 0.2363018 0.02665245 0.006298022 158.78  0   1   12000-1
# 5          0        0 0.2363018 0.02665245 0.006298022 158.78  0   1   12000-1
# 6          0        0 0.2363018 0.02665245 0.006298022 158.78  0   1   12000-1
```
