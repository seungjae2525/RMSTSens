#' @title Function to merge \code{RMSTSens} objects about results of sensitivity analysis
#'
#' @description Merge the \code{RMSTSens} object, which is the result of the previously performed analysis, and the additional \code{RMSTSens} objects, which are the results of the additional analyses.
#'
#' @param x A list of objects for class \code{RMSTSens}. For example, list(result1, result2, result3, ...). See examples.
#'
#' @rdname merge_object
#'
#' @return New merged \code{RMSTSens} object
#'
#' @examples
#' dat <- gbsg
#' dat$size2 <- ifelse(dat$size <= 20, 0,
#'                     ifelse(dat$size > 20 & dat$size <= 50, 1, 2))
#' dat$age2 <- dat$age/100
#' dat$er2 <- dat$er/1000
#'
#' ## Estimation of propensity score
#' denom.fit <- glm(hormon~(age2)^3+(age2)^3*log(age2)+meno+factor(size2)+sqrt(nodes)+er2,
#'                  data=dat, family=binomial(link='logit'))
#' dat$Ps <- predict(denom.fit, type='response')
#'
#' ## Between-group difference in adjusted RMST based on shifted propensity score
#' ## Adjusted RMST with not specified tau and with multiple lambda
#' # Using approximate optimization method
#' results.approx2 <- RMSTSens(time='rfstime', status='status', exposure='hormon',
#'                             exposed.ref.level=1, ps='Ps', data=dat, methods='Approx',
#'                             use.multicore=TRUE, n.core=2,
#'                             lambda=c(1,1.5,2.0), tau=365.25*5, ini.par=1, verbose=FALSE)
#' merge_object(x=list(results.approx2))
#'
#' results.approx3 <- RMSTSens(time='rfstime', status='status', exposure='hormon',
#'                             exposed.ref.level=1, ps='Ps', data=dat, methods='Approx',
#'                             use.multicore=TRUE, n.core=2,
#'                             lambda=c(1.7), tau=365.25*5, ini.par=1, verbose=FALSE)
#' merge_object(x=list(results.approx2, results.approx3))
#'
#' @export
merge_object <- function(x = list()){
  if(length(x) == 1){ # an object
    if (!inherits(x[[1]], "RMSTSens")){
      stop("Argument 'object' must be an object of class \"RMSTSens\".")
    }
    xx.data.frame <- x[[1]]$result.df
    x.out <- x[[1]]
  } else { # a list of objects
    xx.data.frame <- x[[1]]$result.df

    for(i in 1:length(x)){
      if (!inherits(x[[i]], "RMSTSens")){
        stop("Argument 'object' must be an object of class \"RMSTSens\".")
      }
      xx.data.frame <- rbind(xx.data.frame, x[[i]]$result.df)
    }

    ## Remove duplicates
    xx.data.frame <- unique(xx.data.frame)
    xx.data.frame <- xx.data.frame[order(xx.data.frame$Lambda),]
    x.out <- x[[1]]
  }

  x.out$result.df <- xx.data.frame

  class(x.out) <- "RMSTSens"

  return(x.out)
}
