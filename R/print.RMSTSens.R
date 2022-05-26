#' @title Print for "RMSTSens" objects
#' @description Print for objects of class "RMSTSens".
#' @param x an object of class "RMSTSens"
#' @param digits print digits, Default: max(1L, getOption("digits") - 3L)
#' @param ... further arguments passed to or from other methods.
#' @details print.RMSTSens tries to be smart about results for sensitivity analysis
#' @examples
#' if(interactive()){
#'  dat <- gbsg
#'  dat$size2 <- ifelse(dat$size <= 20, 0,
#'                      ifelse(dat$size > 20 & dat$size <= 50, 1, 2))
#'  dat$age2 <- dat$age/100
#'  dat$er2 <- dat$er/1000
#'
#'  ## Estimation of propensity score
#'  denom.fit <- glm(hormon~(age2)^3+(age2)^3*log(age2)+meno+factor(size2)+sqrt(nodes)+er2,
#'                   data=dat, family=binomial(link='logit'))
#'  dat$Ps <- predict(denom.fit, type='response')
#'
#'  ## Between-group difference in adjusted RMST based on shifted propensity score
#'  ## Adjusted RMST with tau equal to 5-year
#'  # Using direct optimization method
#'  results.optim <- RMSTsensitivity(time='rfstime', status='status', exposure='hormon',
#'                                   exposed.ref.level=1, ps='Ps', data=dat, methods='Optim',
#'                                   use.multicore=TRUE, n.core=parallel::detectCores()/2,
#'                                   lambda=1.5, tau=365.25*5, ini.par=1, verbose=FALSE)
#'  print(results.optim)
#'
#'  # Using approximate optimization method
#'  results.approx <- RMSTsensitivity(time='rfstime', status='status', exposure='hormon',
#'                                    exposed.ref.level=1, ps='Ps' ,data=dat, methods='Approx',
#'                                    use.multicore=TRUE, n.core=2,
#'                                    lambda=1.5, tau=365.25*5, ini.par=1, verbose=FALSE)
#'  print(results.approx)
#'
#'  ## Adjusted RMST with not specified tau and with multiple lambda
#'  # Using approximate optimization method
#'  results.approx2 <- RMSTsensitivity(time='rfstime', status='status', exposure='hormon',
#'                                     exposed.ref.level=1, ps='Ps' ,data=dat, methods='Approx',
#'                                     use.multicore=TRUE, n.core=2,
#'                                     lambda=c(1,1.5), tau=365.25*5, ini.par=1, verbose=FALSE)
#'  print(results.approx2)
#'
#'  re.ap.boot <- boot.ci.RMST(x=results.approx2, B=20, level=0.95, seed=220524,
#'                formula=hormon~(age2)^3+(age2)^3*log(age2)+meno+factor(size2)+sqrt(nodes)+er2,
#'                model="logistic", use.multicore=TRUE, n.core=2, verbose=TRUE)
#'  print(re.ap.boot)
#' }
#'
#' @export

print.RMSTSens <- function (x, digits = max(1L, getOption("digits") - 3L), ...){

  if (!inherits(x, "RMSTSens")){
    stop("Argument 'x' must be an object of class \"RMSTSens\".")
  }

  xx <- x$result.df

  cat("\n")
  cat("N=", xx$N[1])
  cat(", N in exposed group=", xx$N.exposed[1])
  cat(", N in unexposed group=", xx$N.unexposed[1], "\n")

  cat("\n")
  cat("Number of events=", xx$N.event.exposed[1] + xx$N.event.unexposed[1], "\n")
  cat("Nunbmer of event in exposed group=", xx$N.event.exposed[1], "\n")
  cat("Nunbmer of event in unexposed group=", xx$N.event.unexposed[1], "\n")

  cat("\n")
  cat("Censoring rate=", round(xx$cen.rate[1], digits), "\n")
  cat("Censoring rate in exposed group=", round(xx$cen.rate.exposed[1], digits), "\n")
  cat("Censoring rate in unexposed group=", round(xx$cen.rate.unexposed[1], digits), "\n")

  lambda.digit <- sapply(xx$Lambda, function(x) match(TRUE, round(x, 1:20) == x))
  spr <- sprintf(paste0("%.",max(lambda.digit),"f"), xx$Lambda)

  cat("\n")
  cat("Parameters: \n")
  cat("Optimization method=", xx$Method[1])
  if(length(xx$Lambda) == 1){
    cat(", Lambda=", spr)
  } else {
    cat(",", paste0("Lambda=(",paste0(spr, collapse=", "),")"))
  }
  cat(", Tau=", xx$Tau[1], "\n")

  # savedig <- options(digits = digits)
  # on.exit(options(savedig))

  cat("\n")
  cat("Results: \n")
  cat("RMST for exposed group: \n")
  tmp1 <- round(cbind(xx$min.exposed, xx$max.exposed), digits)
  rownames(tmp1) <- paste0(rep(expression(Lambda),nrow(xx)), '=', spr)
  colnames(tmp1) <- c("    Minimum", "    Maximum")
  print(tmp1)

  cat("RMST for unexposed group: \n")
  tmp2 <- round(cbind(xx$min.unexposed, xx$max.unexposed), digits)
  rownames(tmp2) <- paste0(rep(expression(Lambda),nrow(xx)), '=', spr)
  colnames(tmp2) <- c("    Minimum", "    Maximum")
  print(tmp2)
  # printCoefmat(tmp2, digits = digits)

  cat("Range of RMST difference: \n")
  tmp3 <- round(cbind(xx$RMST.diff.min, xx$RMST.diff.max), digits)
  rownames(tmp3) <- paste0(rep(expression(Lambda),nrow(xx)), '=', spr)
  colnames(tmp3) <- c("Lower bound", "Upper bound")
  print(tmp3)

  if("RMST.diff.max.upper" %in% colnames(xx)){
    cat("Confidence interval of RMST difference: \n")
    tmp4 <- round(cbind(xx$RMST.diff.min.lower, xx$RMST.diff.max.upper), digits)
    rownames(tmp4) <- paste0(rep(expression(Lambda),nrow(xx)), '=', spr)
    colnames(tmp4) <- c("Lower bound", "Upper bound")
    print(tmp4)
  }

  # invisible(x)
  invisible(NULL)
}
