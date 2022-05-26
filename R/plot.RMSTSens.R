#' @title Sensitivity analysis for RMST
#'
#' @description Function for sensitivity analysis of unmeasured confounding for restricted mean survival time using propensity score
#'
#' @importFrom graphics plot
#' @rdname autoplot.RMSTSens
plot.RMSTSen <- function(x, ...) {
  print(ggplot2::autoplot(x, ...))
}
