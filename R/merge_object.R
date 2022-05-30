#' @title Function to make object for plotting results of sensitivity analysis
#'
#' @param xxx An object or list of objects for class \code{RMSTSens}.
#'
#' @rdname make_object
#'
#' @return Dataframe of objects
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
#'                             lambda=c(1,1.5), tau=365.25*5, ini.par=1, verbose=FALSE)
#' merge_object(results.approx2)
#'
#' results.approx3 <- RMSTSens(time='rfstime', status='status', exposure='hormon',
#'                             exposed.ref.level=1, ps='Ps', data=dat, methods='Approx',
#'                             use.multicore=TRUE, n.core=2,
#'                             lambda=c(1.7), tau=365.25*5, ini.par=1, verbose=FALSE)
#' merge_object(list(results.approx2, results.approx3))
#'
#' @export
merge_object <- function(xxx){
  if(is.data.frame(xxx$data)){ # an object
    if (!inherits(xxx, "RMSTSens")){
      stop("Argument 'object' must be an object of class \"RMSTSens\".")
    }
    xx.data.frame <- xxx$result.df
    xxx.out <- xxx
  } else { # a list of objects
    xx.data.frame <- xxx[[1]]$result.df

    for(i in 1:length(xxx)){
      if (!inherits(xxx[[i]], "RMSTSens")){
        stop("Argument 'object' must be an object of class \"RMSTSens\".")
      }
      xx.data.frame <- rbind(xx.data.frame, xxx[[i]]$result.df)
    }

    ## Remove duplicates
    xx.data.frame <- unique(xx.data.frame)
    xx.data.frame <- xx.data.frame[order(xx.data.frame$Lambda),]
    xxx.out <- xxx[[1]]
  }

  xxx.out$result.df <- xx.data.frame

  class(xxx.out) <- "RMSTSens"

  return(xxx.out)
}
