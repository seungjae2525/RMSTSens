#' RMSTSens: Sensitivity analysis for unmeasured confounding in estimating the difference in restricted mean survival time
#'
#' @description
#' R package \bold{RMSTSens} is sensitivity analysis package for the difference in adjusted restricted mean survival time (RMST)
#' to possible unmeasured confounding.
#'
#' @name RMSTSens-package
#'
#' @author Seungjae Lee \email{seungjae2525@@gmail.com} and Woojoo Lee \email{lwj221@@gmail.com}
#'
#' @references
#' Lee, S., Park, J. H., and Lee, W.
#' Sensitivity analysis for unmeasured confounding in estimating the difference in restricted mean survival time.
#' \emph{Statistical Methods in Medical Research}. 2024.
#' \doi{10.1177/09622802241280782}
#'
#' @keywords package
#'
#' @import parallel
#' @import optimParallel
#' @import ggplot2
#' @import doParallel
#' @import foreach
#' @import doRNG
#'
#' @importFrom stats binomial glm optim predict qlogis quantile spline
#' @importFrom splines ns
"_PACKAGE"
