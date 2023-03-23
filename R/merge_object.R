#' @title Merging \code{RMSTSens} objects for the results of sensitivity analysis
#'
#' @description \code{RMSTSens} is the function to merge the \code{RMSTSens} object, which includes the result of the previously performed sensitivity analysis,
#' and the additional \code{RMSTSens} object, which includes the result of the another performed sensitivity analysis.
#'
#' @param x A list of objects for class \code{RMSTSens}. For example, list(result1, result2, ...). See Examples.
#'
#' @rdname merge_object
#'
#' @return New merged object for class \code{RMSTSens}
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
#'                  data=dat, family=binomial(link="logit"))
#' dat$Ps <- predict(denom.fit, type="response")
#'
#' ## Between-group difference in adjusted RMST based on shifted propensity score
#' ## Adjusted RMST with not specified tau and with multiple lambda
#' # Using approximate optimization method
#' results.approx2 <- RMSTSens(time="rfstime", status="status", exposure="hormon",
#'                             level.exposed="1", ps="Ps", data=dat, methods="Approx",
#'                             use.multicore=TRUE, n.core=2,
#'                             lambda=c(1,1.5,2.0), tau=365.25*5, ini.par=1, verbose=FALSE)
#' merge_object(x=list(results.approx2))
#'
#' results.approx3 <- RMSTSens(time="rfstime", status="status", exposure="hormon",
#'                             level.exposed="1", ps="Ps", data=dat, methods="Approx",
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

      if(i != 1){
        if(!identical(as.numeric(x[[i-1]]$data[,x$argument[5]]), as.numeric(x[[i]]$data[,x$argument[5]]))){
          stop("Each object of class \"RMSTSens\" has a different propensity score.")
        }
        if(!identical(x[[i-1]]$argument[1:5], x[[i]]$argument[1:5])){
          stop("Arguments ('time', 'status', 'exposure', 'level.exposed', or 'ps') in each object of class \"RMSTSens\" are different. \n These must be the same for all objects.")
        }
        if(!identical(x[[i-1]]$result.df$Tau[1], x[[i]]$result.df$Tau[1])){
          stop("'Tau' in one or more objects of class \"RMSTSens\" is not the same.")
        }
      }
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
