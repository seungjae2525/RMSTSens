#' @title Confidence interval for population sensitivity range
#'
#' @description \code{RMSTSens.ci()} is the main function of \code{RMSTSens} package and
#' constructs the percentile bootstrap confidence interval(s) for population sensitivity range.
#'
#' @param x An object for class \code{RMSTSens}. If you want to input several \code{RMSTSens} objects,
#' use the \code{RMSTSens.merge} function. See \code{\link[RMSTSens]{RMSTSens.merge}}.
#' @param B The number of bootstrap replicates. Default: 1000.
#' @param level The confidence level required (i.e., \eqn{1-\alpha}). Default: 0.95.
#' @param seed The seed number. If the propensity score was estimated using methods in the \code{randomForest} or \code{gbm} package,
#' then enter the seed number used at that time.
#' @param formula The formula for estimating propensity score. See Examples.
#' @param model The method for estimating propensity score. "logistic", "rf", or "gbm" can be available.
#' If model is not "logistic", \code{randomForest} or \code{gbm} package is needed. Default: "logistic"
#' @param trunc.prop Optional truncation percentile (Values between 0 and 0.5 are allowed).
#' For example, when trunc.prop=0.01, the left tail is truncated to the 1st percentile,
#' and the right tail is truncated to the 99th percentile. When specified, truncated propensity scores are returned. Default: 0.
#' @param use.multicore Logical scalar indicating whether to parallelize our optimization problem. Default: TRUE.
#' @param n.core The number of CPU cores to use. Default: parallel::detectCores()/2.
#' @param \dots Additional arguments passed on to \code{\link[randomForest]{randomForest}} or \code{\link[gbm]{gbm}} function in \code{randomForest} or \code{gbm} package.
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
#' \item{Lambda}{A scalar or vector of sensitivity parameter(s) \eqn{\Lambda} used}
#' \item{Tau}{User-specific time point \eqn{\tau}; If tau not specified (NULL), the minimum value of last event times in each group}
#' \item{Method}{A optimization method used}
#' \item{min.exposed}{The minimum value of adjusted RMST for exposed group}
#' \item{max.exposed}{The maximum value of adjusted RMST for exposed group}
#' \item{min.unexposed}{The minimum value of adjusted RMST for unexposed group}
#' \item{max.unexposed}{The maximum value of adjusted RMST for unexposed group}
#' \item{RMST.diff.min}{Lower bound of the sensitivity range for the difference in adjusted RMST}
#' \item{RMST.diff.max}{Upper bound of the sensitivity range for the difference in adjusted RMST}
#' \item{min.exposed.lower}{Lower (\eqn{\alpha/2})-quantile of adjusted RMST for exposed group}
#' \item{max.exposed.upper}{Upper (\eqn{\alpha/2})-quantile of adjusted RMST for exposed group}
#' \item{min.unexposed.lower}{Lower (\eqn{\alpha/2})-quantile of adjusted RMST for unexposed group}
#' \item{max.unexposed.upper}{Upper (\eqn{\alpha/2})-quantile of adjusted RMST for unexposed group}
#' \item{RMST.diff.min.lower}{Lower bound of (\eqn{1-\alpha})% percentile bootstrap confidence interval for the population sensitivity range}
#' \item{RMST.diff.max.upper}{Upper bound of (\eqn{1-\alpha})% percentile bootstrap confidence interval for the population sensitivity range}
#' The results for the \code{\link{RMSTSens.ci}} are printed with the \code{\link{print.RMSTSens}} function.
#' To generate the plot of results for the \code{\link{RMSTSens.ci}}, use the \code{\link{autoplot.RMSTSens}} function.
#'
#' @details To assess details of method for sensitivity analysis, see Lee et al. (2024).
#' When estimating the sensitivity range of the difference in adjusted RMST using \code{\link{RMSTSens}} function,
#' there is no need for the argument of formula and model to estimate the propensity score.
#' However, when estimating the confidence interval for the population sensitivity range, the argument of formula and model are absolutely necessary.
#' Note that to use \code{\link{RMSTSens.ci}} function in this package, propensity score should be re-estimated by using
#' 1) \code{\link[stats]{glm}} that has a binomial distribution with logit link function or
#' 2) \code{\link[randomForest]{randomForest}} or \code{\link[gbm]{gbm}} function in \code{randomForest} or \code{gbm} package.
#'
#' Also, if propensity score was re-estimated using method either \code{randomForest} or \code{gbm} package,
#' then you should enter the seed number used at that time in "seed" argument.
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
#' # Using approximate optimization method
#' results.approx2 <- RMSTSens(time="rfstime", status="status", exposure="hormon",
#'                             level.exposed="1", ps="Ps", data=dat, methods="Approx",
#'                             use.multicore=TRUE, n.core=2,
#'                             lambda=c(1,1.1), tau=365.25*5, ini.par=1, verbose=FALSE)
#'
#' # Additional sensitivity analysis when lambda=1.2
#' results.approx3 <- RMSTSens(time="rfstime", status="status", exposure="hormon",
#'                             level.exposed="1", ps="Ps", data=dat, methods="Approx",
#'                             use.multicore=TRUE, n.core=2,
#'                             lambda=c(1.2), tau=365.25*5, ini.par=1, verbose=FALSE)
#'
#' # Percentile bootstrap CI for population sensitivity range
#' re.ap.boot <- RMSTSens.ci(x=RMSTSens.merge(x=list(results.approx2, results.approx3)),
#'               B=50, # Set B=50 to reduce computation time for R checks
#'               level=0.95, seed=220524,
#'               formula=hormon~age+meno+size+factor(grade)+nodes+pgr+er,
#'               model="logistic", trunc.prop=0, use.multicore=TRUE, n.core=2)
#' re.ap.boot
#'
#' @seealso
#'  \code{\link[RMSTSens]{RMSTSens}}, \code{\link[RMSTSens]{print.RMSTSens}}, \code{\link[RMSTSens]{autoplot.RMSTSens}}
#'
#' @references{
#' Lee, S., Park, J. H., and Lee, W.
#' Sensitivity analysis for unmeasured confounding in estimating the difference in restricted mean survival time.
#' \emph{Statistical Methods in Medical Research}. 2024.
#' \doi{10.1177/09622802241280782}
#'
#' Zhao, Q., Small, D. S., and Bhattacharya, B. B.
#' Sensitivity analysis for inverse probability weighting estimators via the percentile bootstrap.
#' \emph{Journal of the Royal Statistical Society Series B: Statistical Methodology}. 2019.
#' \doi{10.1111/rssb.12327}
#' }
#'
#'
#' @keywords methods
#'
#' @export
RMSTSens.ci <- function(x, B=1000, level=0.95, seed=NULL, formula, model="logistic", trunc.prop=0,
                        use.multicore=TRUE, n.core=parallel::detectCores()/2, ...){

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
    ps.present <- as.numeric(predict(model.ps, type="prob")[,2])

  } else if (model == "gbm"){
    if(!requireNamespace("gbm", quietly=TRUE)) {
      stop("\n Error: If model is \"gbm\", then \"gbm\" package needed for this function to work. Please install it.")
    }
    set.seed(seed)
    model.ps <- gbm::gbm(formula=formula, data=data, distribution="bernoulli", verbose=FALSE, ...)
    ps.present <- suppressMessages(as.numeric(predict(model.ps, type="response")))

  } else {
    stop("\n Error: model must be one of \"logistic\", \"rf\", or \"gbm\".")
  }

  if (trunc.prop == 0 & sum(propensity == 0) + sum(propensity == 1) != 0) {
    stop("\n Error: Because there is a propensity score of 0 or 1, \"trunc.prop\" should not be 0.")
  }

  if (trunc.prop != 0) {
    ps.present <- ifelse(ps.present > quantile(ps.present, 1-trunc.prop), quantile(ps.present, 1-trunc.prop),
                         ifelse(ps.present < quantile(ps.present, trunc.prop), quantile(ps.present, trunc.prop), ps.present))
  }

  stopcond <- identical(ps.present, propensity)

  if (stopcond == FALSE) {
    stop("\n Error: Formula used here may be different from the formula used to calculate range of RMST.")
  }

  ## Empty vector
  data$Ps <- NULL
  min.exposd <- max.exposed <- c()
  min.unexposd <- max.unexposed <- c()
  RMST.diff.min <- RMST.diff.max <- c()

  # Setting for parallel processing
  cores <- n.core
  chunk.size <- B/(cores)

  # split data by ourselves
  cl <- makeCluster(cores)
  registerDoParallel(cl, cores=cores)

  ## Percentile bootstrap
  registerDoRNG(seed)
  mat.total <- foreach(iii=1:cores,
                       .combine='rbind',
                       .export=c('RMSTSens','optim_data','optim_f'),
                       .packages = c("parallel"),
                       .multicombine=TRUE) %dopar%
    {
      # local data for results
      mat <- matrix(0, nrow=chunk.size, ncol=(7)*length(lambda))

      for(xx in ((iii-1)*chunk.size+1):(iii*chunk.size)){
        dat.temp <- data[sample.int(nrow(data), nrow(data), replace=TRUE), ]
        if(model == "logistic"){
          dat.temp$PS.boots <- glm(formula, data=dat.temp, family=binomial(link="logit"))$fitted.values
        } else if (model == "rf") {
          model.ps.temp <- randomForest::randomForest(formula=formula, data=dat.temp, do.trace=FALSE, ...)
          dat.temp$PS.boots <- as.numeric(predict(model.ps.temp, type="prob")[,2])
        } else if (model == "gbm"){
          model.ps.temp <- gbm::gbm(formula=formula, data=dat.temp, distribution="bernoulli", verbose=FALSE, ...)
          dat.temp$PS.boots <- suppressMessages(as.numeric(predict(model.ps.temp, type="response")))
        }

        if(sum(dat.temp$PS.boots == 0) + sum(dat.temp$PS.boots == 1) != 0){
          warning(paste0("There is a propensity score of 0 or 1 in ",iii,"-th replicate. 0 or 1 will be truncated by \"trunc.prop\". \n"))
          dat.temp$PS.boots <- ifelse(dat.temp$PS.boots > quantile(dat.temp$PS.boots, 1-trunc.prop),
                                      quantile(dat.temp$PS.boots, 1-trunc.prop),
                                ifelse(dat.temp$PS.boots < quantile(dat.temp$PS.boots, trunc.prop),
                                       quantile(dat.temp$PS.boots, trunc.prop), dat.temp$PS.boots))
        }

        ## Optimization
        temp.re <- RMSTSens(time=time, status=status, exposure=exposure, level.exposed=level.exposed,
                            ps="PS.boots", data=dat.temp,
                            methods=methods, use.multicore=FALSE, n.core=NULL,
                            lambda=lambda, tau=tau, ini.par=ini.par)

        ## save
        mat[xx - (iii-1)*chunk.size, (length(lambda)*0+1):(length(lambda)*1)] <- temp.re$result.df$min.exposed
        mat[xx - (iii-1)*chunk.size, (length(lambda)*1+1):(length(lambda)*2)] <- temp.re$result.df$max.exposed
        mat[xx - (iii-1)*chunk.size, (length(lambda)*2+1):(length(lambda)*3)] <- temp.re$result.df$min.unexposed
        mat[xx - (iii-1)*chunk.size, (length(lambda)*3+1):(length(lambda)*4)] <- temp.re$result.df$max.unexposed
        mat[xx - (iii-1)*chunk.size, (length(lambda)*4+1):(length(lambda)*5)] <- temp.re$result.df$RMST.diff.min
        mat[xx - (iii-1)*chunk.size, (length(lambda)*5+1):(length(lambda)*6)] <- temp.re$result.df$RMST.diff.max
        mat[xx - (iii-1)*chunk.size, (length(lambda)*6+1):(length(lambda)*7)] <- xx
      }
      # return local results
      return(mat)
    }
  ## Shut down the workers
  stopImplicitCluster()
  stopCluster(cl)

  ##
  colnames(mat.total) <- c(paste0(rep(c("min.exposd", "max.exposed", "min.unexposd", "max.unexposed",
                                        "RMST.diff.min", "RMST.diff.max"), each=length(lambda)), "." ,lambda),
                           paste0("xx", ".", lambda))
  results.boots <- data.frame(mat.total)
  rm(mat.total)

  ###
  min.exposd <- as.vector(unlist(results.boots[, (length(lambda)*0+1):(length(lambda)*1)]))
  max.exposed <- as.vector(unlist(results.boots[, (length(lambda)*1+1):(length(lambda)*2)]))
  min.unexposd <- as.vector(unlist(results.boots[, (length(lambda)*2+1):(length(lambda)*3)]))
  max.unexposed <- as.vector(unlist(results.boots[, (length(lambda)*3+1):(length(lambda)*4)]))
  RMST.diff.min <- as.vector(unlist(results.boots[, (length(lambda)*4+1):(length(lambda)*5)]))
  RMST.diff.max <- as.vector(unlist(results.boots[, (length(lambda)*5+1):(length(lambda)*6)]))

  ## Dataframe for summary results
  results.summary.dat <- data.frame(min.exposd=min.exposd,
                                    max.exposed=max.exposed,
                                    min.unexposd=min.unexposd,
                                    max.unexposed=max.unexposed,
                                    RMST.diff.min=RMST.diff.min,
                                    RMST.diff.max=RMST.diff.max,
                                    lambda=rep(lambda, each=B))

  ##
  aa <- (1 - level)/2 # alpha/2
  x$result.df$min.exposed.lower <- sapply(lambda,
                                          function(x) quantile(results.summary.dat$min.exposd[results.summary.dat$lambda == x], aa))
  x$result.df$max.exposed.upper <- sapply(lambda,
                                          function(x) quantile(results.summary.dat$max.exposed[results.summary.dat$lambda == x], 1-aa))
  x$result.df$min.unexposed.lower <- sapply(lambda,
                                            function(x) quantile(results.summary.dat$min.unexposd[results.summary.dat$lambda == x], aa))
  x$result.df$max.unexposed.upper <- sapply(lambda,
                                            function(x) quantile(results.summary.dat$max.unexposed[results.summary.dat$lambda == x], 1-aa))

  ## Confidence interval
  # (alpha/2) quantile of infimum
  x$result.df$RMST.diff.min.lower <- sapply(lambda,
                                            function(x) quantile(results.summary.dat$RMST.diff.min[results.summary.dat$lambda == x], aa))
  # (1-alpha/2) quantile of supremum
  x$result.df$RMST.diff.max.upper <- sapply(lambda,
                                            function(x) quantile(results.summary.dat$RMST.diff.max[results.summary.dat$lambda == x], 1-aa))

  return(x)
}
