
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

# Methodology

## Multivariate allocation in two-stage sampling design

Two-stage sampling design with stratification of PSUs are very common
for households surveys in Official Statistics. However, their sample
allocations are usually a non-easy task, because the households surveys
have multi-domains and multi-purpose objectives, so they have to provide
accurate estimates for different variables and different domains
(e.g. geographical areas such as national, regional, etc.).

The whole sample size of a survey, *n*, is often an exogenous
information, because it is usually defined by budget and, sometimes,
also by logistic constraints.

Instead, the allocation of the sample *n* among the strata
(*h* = 1, …, *L*), generally the lower level at which estimates are
required, can be defined in different ways (mainly with uniform,
proportional or optimal allocation). In Official Statistics, optimal
allocation is the most widely used method, especially because
information on the size of the strata and the variance of target
variables in the strata can be obtained. In particular, the variances of
the target variables or at least proxy variables for each stratum can be
computed, for instance, from register data or previous occasions of the
same survey.

*R2BEAT* easily manages all the complexity due to the optimal sample
allocation in two-stage sampling design and provides several outputs for
evaluating the allocation. Its name stands for “R ‘to’ Bethel Extended
Allocation for Two-stage”. It is an extension of another open-source
software called *Mauss-R* (Multivariate Allocation of Units in Sampling
Surveys), implemented by ISTAT researchers
(<a href="https://www.istat.it/en/methods-and-tools/methods-and-it-tools/design/design-tools/mauss-r" class="uri">https://www.istat.it/en/methods-and-tools/methods-and-it-tools/design/design-tools/mauss-r</a>).
*Mauss-R* determines the optimal sample allocation in multivariate and
multi-domains estimation, for one-stage stratified samples. It extends
the Neyman (1934) Tschuprow (1923) allocation method to the case of
several variables, adopting a generalization of the Bethel’s proposal
(1989).

To complete the suite of tools developed by Istat in order to cover the
stratified sample design, we cite *SamplingStrata*
(<a href="https://CRAN.R-project.org/package=SamplingStrata" class="uri">https://CRAN.R-project.org/package=SamplingStrata</a>),
that allows to jointly optimize both the stratification of the sampling
frame and the allocation, still in the multivariate multidomain case
(only for one-stage designs), and *MultiWay.Sample.Allocation*, that
allows to determine the optimal sample allocation for multi-way
stratified sampling designs and incomplete stratified sampling designs
(<a href="https://www.istat.it/en/methods-and-tools/methods-and-it-tools/design/design-tools/multiwaysampleallocation" class="uri">https://www.istat.it/en/methods-and-tools/methods-and-it-tools/design/design-tools/multiwaysampleallocation</a>).

The main idea for optimal allocation is that strata with larger size and
larger variability with respect to the target variables need more sample
size to provide better estimates. Since, in our case, data from previous
surveys occasions were available, an optimal allocation could be
deployed.

*R2BEAT* develops Bethel’s proposal but, besides, enables the
determination of optimal sample allocation for two-stage stratified
samples, by using the function *beat.2st.*

Before using this function, a preliminary step is needed. In fact, it is
necessary to define some parameters for performing the function (see the
package documentation for the input formats). In particular:

1.  population size of each stratum in terms of households or
    individuals depending on the aim of the survey *N*);
2.  strata to be censused (*CENS*);
3.  cost of interviews per stratum (*COST*);
4.  number of PSUs to be selected in each stratum (*MINCS*);
5.  average dimension of SSUs in each stratum (*DELTA*, *Δ*);
6.  minimum number of SSUs to be interviewed in each selected PSU
    (*MINIMUM*);
7.  mean of each target variables in each strata (*M*, *p̂* or *x̂*);
8.  estimated population standard deviation of each target variable in
    each stratum (*S*, *Ŝ*);
9.  intra-class correlation coefficients for each target variable in
    each stratum (*RHO*, *ρ̂*);
10. estimator effect for each target variable in each stratum (*EFFST*,
    *e**f**f**s**t*).

To better understand the optimal allocation, a few words on points 7-10
are necessary. In particular, the mean can be estimated using the
function *svystatTM* of the package *ReGenesees*. For the estimated
population standard deviation, a distinction between binary and
quantitative variable is needed. In the case of binary variables, the
estimated population standard deviation of variable *i* is equal to

$$Var(\\hat{Y\_{g}})=\\sum\_{h=1}^{H}N\_{h}^{2} (1- \\frac{ 
  n\_{h}}
  {N\_{h}}) \\frac{ 
  S\_{h,g}^{2}}
  {n\_{h}} \\;\\;\\;  g=1,...,G$$

$$ 
\\hat{S} = \\sqrt{\\hat{p}\_i \\times (1 - \\hat{p}\_i)} 
$$

where *p̂*<sub>*i*</sub> is an estimate of the proportion of individuals
with a given characteristic in the population strata, that is the mean
of a binary variable. While, in the case of a quantitative variable, it
is equal to where
$\\hat{M}\_i^2 = \\frac{\\sum\_{k \\in s} y\_{ik}^2 \\ w\_k}{N}$ and
$\\hat{\\bar{Y}}\_i = \\frac{ \\sum\_{k \\in s} y\_{ik} \\ w\_k}{N}$ are
the quadratic mean and the arithmetic mean estimated on sample data of
the previous occasion(s) of the survey, respectively, with
*y*<sub>*i**k*</sub> the value of the target variable *i* observed on
the unit *k*, *w*<sub>*k*</sub> the corresponding sampling weight and
*N* the population size.

A crucial statistic for the optimal sample allocation in two-stage
sampling design is the intra-class correlation coefficient, *ρ̂*. Due to
logistics, in household surveys the sample is not wide-spread in the
whole country, but often \`\`clusterised" in PSUs and, then, in SSUs
(clusters of family members). If the units of clusters are too similar
to each another, it is not efficient to collect too many units from the
same cluster. Anyway, sometimes a small lost in efficiency of the
estimates is accepted, because it is paid off by a gain in the
organization of the survey and, moreover, by cost reduction.

Therefore, to optimize the sample size, it is important to compute the
intra-class correlation coefficient, because it provides a measure of
data clustering in PSUs and SSUs. In general, if the value of *ρ̂* is
close to 1 the clustering is high and it is convenient to collect only a
few units in the cluster. On the contrary, if *ρ̂* is close to 0,
collection of units from the same cluster does not affect the efficiency
of the estimates.

To compute *ρ̂*, another important statistic must be introduced: the
design effect (*DEFF*, *d**e**f**f*). The design effect measures how
much the sampling variance under the adopted sampling design is inflated
with respect to simple random sample (*s**r**s*), with the same sample
size. In formula: where *b* is the average cluster (i.e. PSU) size in
terms of the final sampling units and *ρ̂*<sub>*i*</sub> the intra-class
correlation within the cluster (PSU) for the variable *i*.

The package *R2BEAT* takes into account a more general expression of
*d**e**f**f*. This expression refers to a typical situation for
households surveys in which PSUs are assigned to Self-Representing (SR)
strata, that is they are included for sure in the sample, or to
Not-Self-Representing (NSR) strata, where they are selected by chance.
In practice, this assignment is usually performed by comparing their
measure of size (*MOS*) with respect to the threshold: where *m̄* is the
minimum number of SSUs to be interviewed in each selected PSU
(*MINIMUM*), *f* = *n*/*N* is the sampling fraction and *Δ* (*DELTA*) is
the average dimension of the SSU in terms of elementary survey units.
Then *DELTA* must be set equal to 1 if, for the survey, the selection
units are the same as the elementary units (that is, household-household
or individuals-individuals), whereas it must be set equal to the average
dimension of the households if the elementary units are individuals,
while the selection units are the households.

PSUs with *MOS* exceeding the threshold are identified as SR, while the
remaining PSUs are NSR.

Then the extended expression of *d**e**f**f* is

where, for *S**R* and *N**S**R* strata,
Of course, if there are no *S**R* strata the expression recurs. The
design effect is equal to 1 under the *s**r**s* design and increases for
each additional stage of selection, due to intra-class correlation
coefficient which is, usually, positive. It can be computed using the
function *svystatTM* from *ReGenesees*, setting up the parameter
*deff=TRUE*.

The intra-class correlation coefficient for *N**S**R* can be derived
from the expression of *d**e**f**f*, that is While it is not necessary
to compute the intra-class correlation coefficient for *S**R* strata
because, just one PSU is selected and the intra-class correlation is 1
by definition.

The last statistic to be defined is the estimator effect (*EFFST*,
*e**f**f**s**t*), that is, how much the sampling variance of the applied
estimator under the adopted design is inflated or deflated with respect
to the HT estimator, on the same sample. The *e**f**f**s**t* is equal to
It is an optional parameter for *R2BEAT*, but it useful to take into
account, from the allocation phase, of the impact on the estimates of a
different estimator other than the Horvitz-Thompson estimator (*H**T*).
Indeed, the most applied estimator for the households survey is the
calibrated estimator that, through the use of auxiliary variables,
provides better estimates than HT. Then, a reduction of the final sample
size can be expected when applying the calibration estimator in place of
the HT.

The optimal allocation is defined by *R2BEAT* solving the minimum
optimization problem
$$
\\left\\lbrace
\\begin{array}{l} 
    C = \\text{min} \\\\
    \\sigma\\left(\\hat{\\bar{Y}}\_{i,h}\\right) \\leq \\delta \\left(\\hat{\\bar{Y}}\_{i,h} \\right)
    \\hspace{0.5cm} \\begin{array}{c}
        i = 1, \\dots, J \\\\
        h = 1, \\dots, L 
    \\end{array}
\\end{array}
\\right.,
$$
where *C* is the global cost of the survey (if *COST* is equal to 1 in
all the strata, *C* = *n*) and
$\\sigma\\left(\\hat{\\bar{Y}}\_{i,h}\\right)$ is the error we expect to
observe estimating the *Y*<sub>*i*</sub> variable in the stratum *h*
(expected errors), that must be lower than the precision constraints
defined by the user.

The solution is obtained with an iterative algorithm. Then, *Ŝ* (see,
e.g, expression and ) is multiplied by the new design effect and the
estimator effect and a new allocation is computed. The algorithm stops
when the difference between two consecutive iterations is lower than a
predefined threshold. At each step, an allocation is provided and the
design effect is updated following the expression .

The package *R2BEAT* provides several outputs that help the evaluation
of the allocation.

It is important to point out that, as stated at the beginning of this
section, the *n* is usually given. Then, to find the optimal allocation
it is necessary to tune the precision constraints until the desired
sample size is matched.

## First stage selection

Once performed the allocation step, the function *StratSel* allows to
perform the selection of the Primary Stage Units (PSUs).

The function carries out the stratification of the PSUs, for each
estimation domain, according to a size threshold . PSUs with measure of
size exceeding the threshold are identified as SR (see also description
in the previous section). The remaining NSR PSUs are ordered by their
measure of size and divided into strata whose sizes are approximately
equal to the threshold multiplied by the number of PSUs to be selected
in each stratum (*MINCS*). In this way, strata are composed by PSUs
having as equal as possible size. For each PSU, *StratSel* determines
the size of final sample units on the basis of the size planned in
input. Then the selection of a fixed number of PSUs per stratum is
carried out using the Sampford’s method (unequal probabilities, without
replacement, fixed sample size), implemented by the *UPsampford*
function of the R package *sampling*.

The function produces an output list, reporting the overall information
on selected PSUs. In particular, in the fourth output, it provides the
first order inclusion probability for each PSU, that is, for the first
selection stage. The inclusion probability for the second stage can be
easily obtained dividing the number of SSUs to be selected in the PSU by
the measure of size of the PSU. Then, the design weights for each unit
in the sample are equal to the inverse of the product between the first
order and the second stage inclusion probabilities.

The inclusion probability of the first and second stage of the sampling
scheme are such that the overall inclusion probabilities are almost
constant for each unit in the sample, thus obtaining self-weighting
samples.

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
allocat <- allocation$alloc[-nrow(allocation$alloc),]
str(allocat)
allocat$STRATUM <- as.character(as.numeric(allocat$STRATUM)*100)
unique(allocat$STRATUM)
PSU_strat$STRATUM <- PSU_strat$STRATUM*100
unique(PSU_strat$STRATUM)
sampled_PSU <- StratSel(dataPop= PSU_strat,
                       idpsu= ~ PSU_ID,
                       dom= ~ STRATUM,
                       final_pop= ~ PSU_MOS,
                       size= ~ PSU_MOS,
                       PSUsamplestratum= 1,
                       min_sample= 12,
                       min_sample_index= FALSE,
                       dataAll=allocat,
                       domAll= ~ factor(STRATUM),
                       f_sample= ~ ALLOC,
                       planned_min_sample= NULL,
                       launch= F)
sampled_PSU[[2]]
#    Domain SRdom nSRdom SRdom+nSRdom SR_PSU_final_sample_unit NSR_PSU_final_sample_unit
# 1     100   126     53          179                     8240                       660
# 2     200    20     15           35                      305                        99
# 3     300     2      8           10                       28                        37
# 4     400     6     16           22                      199                       288
# 5     500     1      6            7                       68                       203
# 6     600     6      8           14                       11                        81
# 7     700     1      6            7                      939                       381
# 8     800     4      5            9                      140                       473
# 9     900     2      8           10                      423                       498
# 10   1000     3      8           11                      751                       186
# 11   1100     1      3            4                       51                        97
# 12   1200    10     23           33                       88                       192
# 13   1300     3     17           20                       37                        70
# 14   1400     1      7            8                      113                       104
# 15   1500    26     30           56                       26                        70
# 16   1600     7     39           46                       98                        59
# 17   1700    25     39           64                       42                        95
# 18  Total   244    291          535                    11559                      3593
# 19   Mean                                                680                       211
```
