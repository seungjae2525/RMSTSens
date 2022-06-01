#' @title Confidence interval for bias-adjusted restricted mean survival time
#'
#' @description Function for constructing the confidence interval(s) for restricted mean survival time using propensity score.
#'
#' @param x An object for class \code{RMSTSens}. If you want to input several \code{RMSTSens} objects, use the \code{merge_object} function. See \code{merge_object}.
#' @param B The number of bootstrap replicates, Default: 1000.
#' @param level The confidence level required (i.e., \eqn{1-\alpha}), Default: 0.95.
#' @param seed The seed number. If the propensity score was estimated using methods in the \code{randomForest} or \code{gbm} package, then enter the seed number used at that time.
#' @param formula The formula for estimating propensity score. See Examples.
#' @param model The method for estimating propensity score, Default: "logistic". "logistic", "rf", or "gbm" can be available. If model is not "logistic", \code{randomForest} or \code{gbm} package are available.
#' @param trunc.prop Optional truncation percentile (0-0.5), Default: 0. E.g. when trunc.prop=0.01, the left tail is truncated to the 1st percentile, and the right tail is truncated to the 99th percentile. When specified, truncated propensity scores are returned.
#' @param use.multicore Logical scalar indicating whether to parallelize our optimization problem, Default: TRUE.
#' @param n.core The number of cores to use, Default: parallel::detectCores()/2.
#' @param verbose According to the verbose level, whether or not print the completion message for each \code{B}/100\emph{th} bootstrap, Default: TRUE.
#' @param \dots Additional arguments passed on to \code{randomForest} or \code{gbm} function in \code{randomForest} or \code{gbm} package.
#'
#' @return The object is a data.frame with class \code{RMSTSens}. The function returns following components:
#' \item{N}{Total number of subjects}
#' \item{N.exposed}{The number of subjects in exposed group}
#' \item{N.unexposed}{The number of subjects in unexposed group}
#' \item{N.event.exposed}{The number of events in exposed group}
#' \item{N.event.unexposed}{The number of events in unexposed group}
#' \item{cen.rate}{Total censoring rate}
#' \item{cen.rate.exposed}{Censoring rate in exposed group}
#' \item{cen.rate.unexposed}{Censoring rate in unexposed group}
#' \item{Lambda}{A scalar or vector of sensitivity parameter \eqn{\Lambda} used}
#' \item{Tau}{User-specific time point \eqn{\tau}, If tau not specified (NULL), use the minimum between the largest observed event time in each groups}
#' \item{Method}{A optimization method used}
#' \item{min.exposed}{The minimum of adjusted RMST based on the shifted propensity score for exposed group}
#' \item{max.exposed}{The maximum of adjusted RMST based on the shifted propensity score for exposed group}
#' \item{min.unexposed}{The minimum of adjusted RMST based on the shifted propensity score for unexposed group}
#' \item{max.unexposed}{The maximum of adjusted RMST based on the shifted propensity score for unexposed group}
#' \item{RMST.diff.min}{Lower bound of point estimate for between-group difference in adjusted RMST based on shifted propensity score}
#' \item{RMST.diff.max}{Upper bound of point estimate for between-group difference in adjusted RMST based on shifted propensity score}
#' \item{RMST.diff.min.lower}{Lower bound of (\eqn{1-\alpha})% confidence interval for between-group difference in adjusted RMST based on shifted propensity score}
#' \item{RMST.diff.max.upper}{Upper bound of (\eqn{1-\alpha})% confidence interval for between-group difference in adjusted RMST based on shifted propensity score}
#' \item{min.exposed.lower}{Lower (\eqn{\alpha/2})-quantile of adjusted RMST based on the shifted propensity score for exposed group}
#' \item{max.exposed.upper}{Upper (\eqn{\alpha/2})-quantile of adjusted RMST based on the shifted propensity score for exposed group}
#' \item{min.unexposd.lower}{Lower (\eqn{\alpha/2})-quantile of adjusted RMST based on the shifted propensity score for unexposed group}
#' \item{max.unexposed.upper}{Upper (\eqn{\alpha/2})-quantile of adjusted RMST based on the shifted propensity score for unexposed group}
#' The results for the \code{\link{RMSTSens.ci}} are printed with the \code{\link{print.RMSTSens}} functions.
#' To generate result plot comparing sensitivity parameters \eqn{\Lambda} with confidence interval and range of adjusted RMST based on shifted propensity score, use the \code{\link{autoplot.RMSTSens}} function.
#'
#' @details To assess details of method for sensitivity analysis, see Lee et al. (2022).
#' When estimating the range of the adjusted RMST based on the shifted propensity score using \code{RMSTSens} function, there is no need for the formula and model used to estimate the propensity score.
#' However, when estimating the confidence interval, the formula and model are absolutely necessary.
#' Note that to use \code{RMSTSens.ci} function in our package, propensity score should be estimated by using 1) \code{glm} that has a binomial distribution with logit link function or 2) \code{randomForest} or \code{gbm} function in \code{randomForest} or \code{gbm} package.
#' Also, if propensity score was estimated using methods in the \code{randomForest} or \code{gbm} package, then you should enter the seed number used at that time in "seed" argument.
#'
#' @examples
#' \dontrun{
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
#'
#' ## Between-group difference in adjusted RMST based on shifted propensity score
#' ## Adjusted RMST with not specified tau and with multiple lambda
#' # Using approximate optimization method
#' results.approx2 <- RMSTSens(time="rfstime", status="status", exposure="hormon",
#'                             level.exposed="1", ps="Ps", data=dat, methods="Approx",
#'                             use.multicore=TRUE, n.core=2,
#'                             lambda=c(1,1.5), tau=365.25*5, ini.par=1, verbose=FALSE)
#' re.ap.boot <- RMSTSens.ci(x=results.approx2, B=40, level=0.95, seed=220524,
#'               formula=hormon~(age2)^3+(age2)^3*log(age2)+meno+factor(size2)+sqrt(nodes)+er2,
#'           model="logistic", trunc.prop=0, use.multicore=TRUE, n.core=2, verbose=TRUE)
#' re.ap.boot
#'
#'
#' ## Estimate of propensity score using random forest
#' library(randomForest)
#' set.seed(220528)
#' model.rf <- randomForest(formula=factor(hormon)~age2+meno+size2+nodes+er2, data=dat)
#' dat$Ps.rf <- as.numeric(predict(model.rf, type= "prob")[,2])
#' sum(dat$Ps.rf== 0) + sum(dat$Ps.rf== 1) # There are some 0 or 1 values
#' # Trimming the propensity score
#' trunc.prop <- 0.01
#' dat$Ps.rf <- ifelse(dat$Ps.rf > quantile(dat$Ps.rf, 1-trunc.prop),
#'                     quantile(dat$Ps.rf, 1-trunc.prop),
#'                     ifelse(dat$Ps.rf < quantile(dat$Ps.rf, trunc.prop),
#'                            quantile(dat$Ps.rf, trunc.prop), dat$Ps.rf))
#' sum(dat$Ps.rf== 0) + sum(dat$Ps.rf== 1)
#' # Range of RMST
#' results.approx.rf <- RMSTSens(time="rfstime", status="status", exposure="hormon",
#'                               level.exposed="1", ps="Ps.rf", data=dat, methods="Approx",
#'                               use.multicore=TRUE, n.core=2,
#'                               lambda=c(1,1.5,2), tau=365.25*5, ini.par=1, verbose=FALSE)
#' # CI of RMST
#' re.rf <- RMSTSens.ci(x=results.approx.rf, B=40, level=0.95, seed=220528,
#'               formula=factor(hormon)~age2+meno+size2+nodes+er2,
#'                      model="rf", trunc.prop=trunc.prop, use.multicore=TRUE, n.core=2, verbose=TRUE,
#'               trContol=ctrl)
#' re.rf
#'
#'
#' ## Estimate of propensity score using generalized boosted regression modeling
#' library(gbm)
#' set.seed(220528)
#' model.gbm <- gbm(formula=hormon~age2+meno+size2+nodes+er2,
#'                  data=dat, distribution= "bernoulli", verbose= FALSE)
#' dat$Ps.gbm <- as.numeric(predict(model.gbm, type= "response"))
#' sum(dat$Ps.gbm== 0) + sum(dat$Ps.gbm== 1)
#' # Range of RMST
#' results.approx.gbm <- RMSTSens(time="rfstime", status="status", exposure="hormon",
#'                                level.exposed="1", ps="Ps.gbm", data=dat, methods="Approx",
#'                              use.multicore=TRUE, n.core=2,
#'                              lambda=c(1,1.5,2), tau=365.25*5, ini.par=1, verbose=FALSE)
#' # CI of RMST
#' re.gbm <- RMSTSens.ci(x=results.approx.gbm, B=40, level=0.95, seed=220528,
#'                   formula=hormon~age2+meno+size2+nodes+er2,
#'                       model="gbm", trunc.prop=0, use.multicore=TRUE, n.core=2, verbose=TRUE,
#'                   trContol=ctrl)
#' re.gbm
#' }
#'
#' @author Seungjae Lee \email{seungjae2525@@gmail.com}
#'
#' @seealso
#'  \code{\link[RMSTSens]{RMSTSens}}, \code{\link[RMSTSens]{print.RMSTSens}}, \code{\link[RMSTSens]{autoplot.RMSTSens}}
#'
#' @references
#' Bakbergenuly I, Hoaglin DC, Kulinskaya E (2020):
#' Methods for estimating between-study variance and overall
#' effect in meta-analysis of odds-ratios.
#' \emph{Research Synthesis Methods},
#' DOI: 10.1002/jrsm.1404
#'
#' @keywords methods
#'
#' @export
RMSTSens.ci <- function(x, B=1000, level=0.95, seed=NULL, formula, model="logistic", trunc.prop=0,
                        use.multicore=TRUE, n.core=parallel::detectCores()/2, verbose=TRUE, ...){

  if (!inherits(x, "RMSTSens")){
    stop("Argument 'x' must be an object of class \"RMSTSens\".")
  }

  if (use.multicore==TRUE){
    if(is.null(n.core)){
      n.core <- parallel::detectCores()/2
    }
  } else {
    n.core <- 0
  }

  if (use.multicore==TRUE & n.core > parallel::detectCores()) {
    stop("\n Error: \"n.core\" is more than the number of available cores.")
  }

  if (trunc.prop < 0 | trunc.prop > 0.5){
    stop("\n Error: \"trunc.prop\" should be between 0 and 0.5.")
  }

  ## arguments and parameters used in previous range calculation
  time <- x$argument[1]
  status <- x$argument[2]
  exposure <- x$argument[3]
  level.exposed <- x$argument[4]
  data <- x$data
  methods <- x$result.df$Method[1]
  lambda <- x$result.df$Lambda
  tau <- x$result.df$Tau[1]
  ini.par <- as.numeric(x$argument[6])
  propensity <- as.numeric(data[,x$argument[5]])

  ## Check propensity score
  if(model=="logistic"){
    ps.present <- as.numeric(glm(formula, data=data, family=binomial(link="logit"))$fitted.values)

  } else if (model == "rf") {
    if(!requireNamespace("randomForest", quietly=TRUE)) {
      stop("\n Error: If model is \"rf\", then \"randomForest\" package needed for this function to work. Please install it.")
    }
    set.seed(seed)
    model.ps <- randomForest::randomForest(formula=formula, data=data, do.trace=FALSE, ...)
    ps.present <- as.numeric(predict(model.ps, type = "prob")[,2])

  } else if (model == "gbm"){
    if(!requireNamespace("gbm", quietly=TRUE)) {
      stop("\n Error: If model is \"gbm\", then \"gbm\" package needed for this function to work. Please install it.")
    }
    set.seed(seed)
    model.ps <- gbm::gbm(formula=formula, data=data, distribution = "bernoulli", verbose=FALSE, ...)
    ps.present <- suppressMessages(as.numeric(predict(model.ps, type = "response")))

  } else {
    stop("\n Error: model must be \"logistic\", \"rf\", or \"gbm\".")
  }

  if (trunc.prop != 0) {
    ps.present <- ifelse(ps.present > quantile(ps.present, 1-trunc.prop), quantile(ps.present, 1-trunc.prop),
                         ifelse(ps.present < quantile(ps.present, trunc.prop), quantile(ps.present, trunc.prop), ps.present))
  }

  stopcond <- identical(ps.present, propensity)

  if (stopcond==FALSE) {
    stop("\n Error: Formula used here may be different from the formula used to calculate range of RMST.")
  }

  ## Empty vector
  data$Ps <- NULL
  min.exposd <- max.exposed <- c()
  min.unexposd <- max.unexposed <- c()
  RMST.diff.min <- RMST.diff.max <- c()
  avaiable.mintau <- c()

  ## Percentile bootstrap
  set.seed(seed)
  for(i in 1:B){
    dat.temp <- data[sample.int(nrow(data), nrow(data), replace=TRUE), ]
    if(model == "logistic"){
      dat.temp$PS <- glm(formula, data=dat.temp, family=binomial(link="logit"))$fitted.values
    } else if (model == "rf") {
      model.ps.temp <- randomForest::randomForest(formula=formula, data=dat.temp, do.trace=FALSE, ...)
      dat.temp$PS <- as.numeric(predict(model.ps.temp, type = "prob")[,2])
    } else if (model == "gbm"){
      model.ps.temp <- gbm::gbm(formula=formula, data=dat.temp, distribution = "bernoulli", verbose=FALSE, ...)
      dat.temp$PS <- suppressMessages(as.numeric(predict(model.ps.temp, type = "response")))
    }

    if(sum(dat.temp$PS == 0 | dat.temp$PS == 1) != 0){
      warning(paste0("There is a propensity score of 0 or 1 in ",B,"-th replicate. 0 or 1 will be truncated by \"trunc.prop\". \n"))

      if (trunc.prop==0) {
        stop("\n Error: Because there is a propensity score of 0 or 1, \"trunc.prop\" should not be 0.")
      }

      dat.temp$PS <- ifelse(dat.temp$PS > quantile(dat.temp$PS, 1-trunc.prop), quantile(dat.temp$PS, 1-trunc.prop),
                            ifelse(dat.temp$PS < quantile(dat.temp$PS, trunc.prop), quantile(dat.temp$PS, trunc.prop), dat.temp$PS))
    }

    for(j in 1:length(lambda)){
      ##
      max.exposed.event.time <- max(dat.temp[,time][(dat.temp[,status] == 1) & (dat.temp[,exposure] == level.exposed)])
      max.unexposed.event.time <- max(dat.temp[,time][(dat.temp[,status] == 1) & (dat.temp[,exposure] != level.exposed)])
      avaiable.mintau[(i-1)*length(lambda)+j] <- min(c(max.exposed.event.time, max.unexposed.event.time))

      ## Optimization
      temp.re <- RMSTSens(time=time, status=status, exposure=exposure, level.exposed=level.exposed,
                          ps="PS", data=dat.temp,
                          methods=methods, use.multicore=use.multicore, n.core=n.core,
                          lambda=lambda[j], tau=tau, ini.par=ini.par)

      min.exposd[(i-1)*length(lambda)+j] <- temp.re$result.df$min.exposed
      max.exposed[(i-1)*length(lambda)+j] <- temp.re$result.df$max.exposed
      min.unexposd[(i-1)*length(lambda)+j] <- temp.re$result.df$min.exposed
      max.unexposed[(i-1)*length(lambda)+j] <- temp.re$result.df$max.unexposed
      RMST.diff.min[(i-1)*length(lambda)+j] <- temp.re$result.df$RMST.diff.min
      RMST.diff.max[(i-1)*length(lambda)+j] <- temp.re$result.df$RMST.diff.max

      if (verbose & ((i-1)*length(lambda)+j) %% ((B/100)*length(lambda)) == 0) {
        cat(paste0("[", Sys.time(), "]"), ((i-1)*length(lambda)+j)/length(lambda),"th end! \n")
      }
    }
  }

  # which(avaiable.mintau < tau)

  ## Dataframe for summary results
  results.summary.dat <- data.frame(min.exposd = min.exposd,
                                    max.exposed = max.exposed,
                                    min.unexposd = min.unexposd,
                                    max.unexposed = max.unexposed,
                                    RMST.diff.min = RMST.diff.min,
                                    RMST.diff.max = RMST.diff.max,
                                    lambda = rep(lambda, times=B))

  ## Confidence interval
  aa <- (1 - level)/2 # alpha/2
  # (alpha/2) quantile of infimum
  x$result.df$RMST.diff.min.lower <- sapply(x$result.df$Lambda,
                                            function(x) quantile(results.summary.dat$RMST.diff.min[results.summary.dat$lambda == x], aa))
  # (1-alpha/2) quantile of supremum
  x$result.df$RMST.diff.max.upper <- sapply(x$result.df$Lambda,
                                            function(x) quantile(results.summary.dat$RMST.diff.max[results.summary.dat$lambda == x], 1-aa))

  x$result.df$min.exposed.lower <- sapply(x$result.df$Lambda,
                                          function(x) quantile(results.summary.dat$min.exposd[results.summary.dat$lambda == x], aa))
  x$result.df$max.exposed.upper <- sapply(x$result.df$Lambda,
                                          function(x) quantile(results.summary.dat$max.exposed[results.summary.dat$lambda == x], 1-aa))
  x$result.df$min.unexposd.lower <- sapply(x$result.df$Lambda,
                                           function(x) quantile(results.summary.dat$min.unexposd[results.summary.dat$lambda == x], aa))
  x$result.df$max.unexposed.upper <- sapply(x$result.df$Lambda,
                                            function(x) quantile(results.summary.dat$max.unexposed[results.summary.dat$lambda == x], 1-aa))

  return(x)
}
