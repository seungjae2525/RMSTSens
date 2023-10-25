#' @title Print for \code{RMSTSens} objects
#'
#' @description Print the sensitivity analysis results for object of class \code{RMSTSens}.
#'
#' @param x An object for class \code{RMSTSens}. Either object \code{RMSTSens} or \code{RMSTSens.ci} is allowed.
#' If you want to see the results of several \code{RMSTSens} objects together, use the \code{RMSTSens.merge} function. See \code{RMSTSens.merge}.
#' @param digits Print digits. Default: max(1L, getOption("digits") - 3L).
#' @param ... Further arguments (currently not used).
#'
#' @details Print the sensitivity analysis results for object (\code{RMSTSens} or \code{RMSTSens.ci}) of class \code{RMSTSens}.
#'
#' @examples
#' dat <- gbsg
#'
#' ## Estimation of propensity score
#' denom.fit <- glm(hormon~age+meno+size+factor(grade)+nodes+pgr+er,
#'                  data=dat, family=binomial(link="logit"))
#' dat$Ps <- predict(denom.fit, type="response")
#'
#' ## Performing the sensitivity analysis - sensitivity range
#' # Using direct optimization method
#' results.optim <- RMSTSens(time="rfstime", status="status", exposure="hormon",
#'                           level.exposed="1", ps="Ps", data=dat, methods="Optim",
#'                           use.multicore=TRUE, n.core=2,
#'                           lambda=1.2, tau=365.25*5, ini.par=1, verbose=FALSE)
#' print(results.optim)
#'
#' # Using approximate optimization method
#' results.approx <- RMSTSens(time="rfstime", status="status", exposure="hormon",
#'                            level.exposed="1", ps="Ps", data=dat, methods="Approx",
#'                            use.multicore=TRUE, n.core=2,
#'                            lambda=1.2, tau=365.25*5, ini.par=1, verbose=FALSE)
#' print(results.approx)
#'
#' ## Performing the sensitivity analysis - sensitivity range with multiple lambda
#' # Using approximate optimization method
#' results.approx2 <- RMSTSens(time="rfstime", status="status", exposure="hormon",
#'                             level.exposed="1", ps="Ps", data=dat, methods="Approx",
#'                             use.multicore=TRUE, n.core=2,
#'                             lambda=c(1,1.2), tau=365.25*5, ini.par=1, verbose=FALSE)
#' print(results.approx2)
#'
#' # Percentile bootstrap CI for population sensitivity range
#' re.ap.boot <- RMSTSens.ci(x=results.approx2,
#'               B=50, # Set B=50 to reduce computation time for R checks
#'               level=0.95, seed=220524,
#'               formula=hormon~age+meno+size+factor(grade)+nodes+pgr+er,
#'               model="logistic", use.multicore=TRUE, n.core=2)
#' print(re.ap.boot)
#'
#' @keywords print
#'
#' @seealso
#'  \code{\link[RMSTSens]{RMSTSens}}, \code{\link[RMSTSens]{RMSTSens.ci}}, \code{\link[RMSTSens]{RMSTSens.merge}}
#'
#' @export

print.RMSTSens <- function (x, digits = max(1L, getOption("digits") - 3L), ...){

  if (!inherits(x, "RMSTSens")){
    stop("Argument 'x' must be an object of class \"RMSTSens\".")
  }

  xx <- x$result.df

  cat("\n")
  cat("N =", xx$N[1])
  cat(", N in exposed group =", xx$N.exposed[1])
  cat(", N in unexposed group =", xx$N.unexposed[1], "\n")

  cat("\n")
  cat("Number of events =", xx$N.event.exposed[1] + xx$N.event.unexposed[1], "\n")
  cat("Nunbmer of event in exposed group =", xx$N.event.exposed[1], "\n")
  cat("Nunbmer of event in unexposed group =", xx$N.event.unexposed[1], "\n")

  cat("\n")
  cat("Censoring rate =", round(xx$cen.rate[1], digits), "\n")
  cat("Censoring rate in exposed group =", round(xx$cen.rate.exposed[1], digits), "\n")
  cat("Censoring rate in unexposed group =", round(xx$cen.rate.unexposed[1], digits), "\n")

  lambda.digit <- sapply(xx$Lambda, function(x) match(TRUE, round(x, 1:20) == x))
  spr <- sprintf(paste0("%.",max(lambda.digit),"f"), xx$Lambda)

  cat("\n")
  cat("Parameters: \n")
  cat("Optimization method =", xx$Method[1])
  if(length(xx$Lambda) == 1){
    cat(", Lambda =", spr)
  } else {
    cat(",", paste0("Lambda = (",paste0(spr, collapse=", "),")"))
  }
  cat(", Tau =", xx$Tau[1], "\n")

  # savedig <- options(digits = digits)
  # on.exit(options(savedig))

  cat("\n")
  cat("Results: \n")
  cat("Adjusted RMST for exposed group: \n")
  tmp1 <- round(cbind(xx$min.exposed, xx$max.exposed), digits)
  rownames(tmp1) <- paste0(rep(expression(Lambda), times=nrow(xx)), '=', spr)
  colnames(tmp1) <- c("    Minimum", "    Maximum")
  print(tmp1)

  cat("Adjusted RMST for unexposed group: \n")
  tmp2 <- round(cbind(xx$min.unexposed, xx$max.unexposed), digits)
  rownames(tmp2) <- paste0(rep(expression(Lambda), times=nrow(xx)), '=', spr)
  colnames(tmp2) <- c("    Minimum", "    Maximum")
  print(tmp2)
  # printCoefmat(tmp2, digits = digits)

  cat("Sensitivity range of difference in adjusted RMST: \n")
  tmp3 <- round(cbind(xx$RMST.diff.min, xx$RMST.diff.max), digits)
  rownames(tmp3) <- paste0(rep(expression(Lambda), times=nrow(xx)), '=', spr)
  colnames(tmp3) <- c("Lower bound", "Upper bound")
  print(tmp3)

  if("RMST.diff.max.upper" %in% colnames(xx)){
    cat("Confidence interval for population sensitivity range: \n")
    tmp4 <- round(cbind(xx$RMST.diff.min.lower, xx$RMST.diff.max.upper), digits)
    rownames(tmp4) <- paste0(rep(expression(Lambda), times=nrow(xx)), '=', spr)
    colnames(tmp4) <- c("Lower bound", "Upper bound")
    print(tmp4)
  }

  # invisible(x)
  invisible(NULL)
}
