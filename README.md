# RMSTSens

<!-- badges: start -->
[![Project Status](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active/)
[![Package version](https://img.shields.io/badge/GitHub-1.0.0-orange.svg)](https://github.com/seungjae2525/RMSTSens/)
[![minimal R version](https://img.shields.io/badge/R-v4.0.0+-blue.svg)](https://cran.r-project.org/)
[![last_update](https://img.shields.io/badge/last%20update--blue.svg)](https://github.com/seungjae2525/RMSTSens/)  
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![R-CMD-check](https://github.com/seungjae2525/RMSTSens/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/seungjae2525/RMSTSens/actions/workflows/R-CMD-check.yaml)
[![Travis](https://travis-ci.org/seungjae2525/RMSTSens.svg)](https://travis-ci.org/seungjae2525/RMSTSens/)
<!-- badges: end -->


## Description
This is the source code for the `RMSTSens` package in R. 
`RMSTSens` is a package aimed at providing a novel sensitivity analysis method for the RMST difference when unmeasured confounding is suspected.
Given a user-specified sensitivity parameter, the sensitivity range of the difference in adjusted RMST is calculated with the percentile bootstrap confidence interval for the population sensitivity range.
 
### References
Not yet...


## Installation
### Current GitHub release:
Installation using R package remotes:

```r
install.packages("remotes") # if devtools not already installed
remotes::install_github("seungjae2525/RMSTSens")
```


### Bug Reports:
You can also report bugs on GitHub under [Issues](https://github.com/seungjae2525/RMSTSens/issues/).
