#' @title Print for "senRMST" objects
#' @description Print for objects of class "senRMST".
#' @param x an object of class "senRMST"
#' @param digits print digits, Default: max(1L, getOption("digits") - 3L)
#' @param ... PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname print.RMSTsensitivity
#' @export

print.RMSTsensitivity <- function (x, digits = max(1L, getOption("digits") - 3L), ...){

  if (!inherits(x, "senRMST")){
    stop("Argument 'x' must be an object of class \"senRMST\".")
  }

  cat("\n")
  cat("N=", x$N)
  cat(", N in exposed group=", x$N.exposed)
  cat(", N in unexposed group=", x$N.unexposed, "\n")

  cat("\n")
  cat("Number of events=", x$N.event.exposed + x$N.event.unexposed, "\n")
  cat("Nunbmer of event in exposed group=", x$N.event.exposed, "\n")
  cat("Nunbmer of event in unexposed group=", x$N.event.unexposed, "\n")

  cat("\n")
  cat("Censoring rate=", round(x$cen.rate, digits), "\n")
  cat("Censoring rate in exposed group=", round(x$cen.rate.exposed, digits), "\n")
  cat("Censoring rate in unexposed group=", round(x$cen.rate.unexposed, digits), "\n")

  cat("\n")
  cat("Parameters: \n")
  cat("Optimization method=", x$Method)
  cat(", Lambda=", x$Lambda)
  cat(", Tau=", x$Tau, "\n")

  savedig <- options(digits = digits)
  on.exit(options(savedig))

  cat("\n")
  cat("Results: \n")
  tmp <- rbind(x$min.exposed, x$max.exposed, x$min.unexposed, x$max.unexposed,
               x$RMST.diff.min, x$RMST.diff.max)
  rownames(tmp) <- c("Minimum of Shifted RMST in exposed group",
                     "Maximum of Shifted RMST in exposed group",
                     "Minimum of Shifted RMST in unexposed group",
                     "Maximum of Shifted RMST in unexposed group",
                     "Minimum of Between-group difference in Shifted RMST",
                     "Maximum of Between-group difference in Shifted RMST")
  colnames(tmp) <- c("Estimate")
  printCoefmat(tmp, digits = digits, ...)

  invisible(x)
}
