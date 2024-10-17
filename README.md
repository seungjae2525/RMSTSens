# RMSTSens

<!-- badges: start -->
[![Project Status](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active/)
[![Package version](https://img.shields.io/badge/GitHub-1.0.0-orange.svg)](https://github.com/seungjae2525/RMSTSens/)
[![minimal R version](https://img.shields.io/badge/R-v4.0.0+-blue.svg)](https://cran.r-project.org/)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![R-CMD-check](https://github.com/seungjae2525/RMSTSens/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/seungjae2525/RMSTSens/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

# Sensitivity analysis for unmeasured confounding in estimating the difference in restricted mean survival time

## Description
This is the source code for the `RMSTSens` package in R. 
`RMSTSens` is a package aimed at providing a novel sensitivity analysis method for the RMST difference when unmeasured confounding is suspected.
Given a user-specified sensitivity parameter, the sensitivity range of the difference in adjusted RMST is calculated with the percentile bootstrap confidence interval for the population sensitivity range. See reference for details.
 
### Reference
Lee, S., Park, J. H., and Lee, W. Sensitivity analysis for unmeasured confounding in estimating the difference in restricted mean survival time. *Statistical Methods in Medical Research*. 2024. <https://doi.org/10.1177/09622802241280782>


## Installation
### Current GitHub release:
Installation using R package `remotes`:

```r
if (!require("devtools")) { install.packages("devtools") } # if devtools not already installed
remotes::install_github("seungjae2525/RMSTSens")
library(RMSTSens)
```

### Bug Reports:
You can also report bugs on GitHub under [Issues](https://github.com/seungjae2525/RMSTSens/issues/).
