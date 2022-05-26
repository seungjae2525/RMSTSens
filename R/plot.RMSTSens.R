#' @title Sensitivity analysis for RMST
#'
#' @description Function for sensitivity analysis of unmeasured confounding for restricted mean survival time using propensity score
#'
#' @rdname autoplot.RMSTSens
#' @param ... further arguments passed to or from other methods.
#'
#' @importFrom graphics plot
#' @export
plot.RMSTSen <- function(x, ...) {
  print(autoplot(x, ...))
}
