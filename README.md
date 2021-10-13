
<!-- README.md is generated from README.Rmd. Please edit README.Rmd file -->

[![CRAN](http://www.r-pkg.org/badges/version/R2BEAT)](https://cran.r-project.org/package=R2BEAT)
[![Downloads](http://cranlogs.r-pkg.org/badges/R2BEAT?color=brightgreen)](http://www.r-pkg.org/pkg/R2BEAT)
[![Mentioned in Awesome Official
Statistics](https://awesome.re/mentioned-badge.svg)](http://www.awesomeofficialstatistics.org)

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
