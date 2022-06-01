RMSTSens <- function(...) UseMethod("RMSTSens")

#' @title Sensitivity analysis for RMST
#'
#' @description Function for sensitivity analysis of unmeasured confounding for restricted mean survival time using propensity score.
#'
#' @param time The name of variable for time to event.
#' @param status The name of variable for censoring indicator (0 if censored, 1 if event).
#' @param exposure The name of variable for exposure (0 if unexposed, 1 if exposed).
#' @param level.exposed Level for exposed group in exposure variable, Default: "1".
#' @param ps The name of variable for propensity score variable or the vector for propensity score i.e., P(A=1|L).
#' @param data A data frame in which contains the follow-up time (time), the event (status), the exposure (exposure), and the propensity score (ps).
#' @param methods A character with the methods how to calculate the adjusted RMST ("Optim", "Approx", "LP1", "LP2"), Default: "Approx". See Details.
#' @param use.multicore Logical scalar indicating whether to parallelize our optimization problem, Default: TRUE.
#' @param n.core The number of cores to use, Default: parallel::detectCores()/2.
#' @param lambda A scalar or vector for sensitivity parameter \eqn{\Lambda}, Default: 2.
#' @param tau User-specific time point, If tau not specified (NULL), use the minimum of the largest observed event time in both groups, Default: NULL.
#' @param ini.par Initial parameter for direct optimization method, Default: 1.
#' @param verbose Conditional on the verbose level, print the message that each optimization (minimization or maximization) for each group was ended, Default: FALSE.
#'
#' @return An object of class \code{RMSTSens}. The object is a data.frame with the following components:
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
#' The results for the \code{RMSTSens} are printed with the \code{\link{print.RMSTSens}} functions.
#' To generate result plot comparing sensitivity parameters \eqn{\Lambda} with range of adjusted RMST based on shifted propensity score, use the \code{\link{autoplot.RMSTSens}} functions.
#'
#' @details There are four possible methods for our sensitivity analysis.
#' In general survival analysis setting, if the censoring rates in both exposure group are less than 0.7, the approximate optimization method (methods="Approx") is not inferior to the direct optimization method in terms of both bias and computational time.
#' Otherwise, one may use the direct optimization method (methods="Optim") to perform the sensitivity analysis, although it takes slightly more computational time than the approximate optimization method.
#' In special settings, to obtain the solution of optimization problem, only need to compute the objective function of linear fractional programming (methods="LP1" or methods="LP2) for at most m (the number of event) candidates of optimization parameters.
#' To assess detail of methods and special settings, see Lee et al. (2022).
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
#' ## Adjusted RMST with tau equal to 5-year
#' # Using direct optimization method
#' results.optim <- RMSTSens(time="rfstime", status="status", exposure="hormon",
#'                           level.exposed="1", ps="Ps", data=dat, methods="Optim",
#'                           use.multicore=TRUE, n.core=2,
#'                           lambda=1.5, tau=365.25*5, ini.par=1, verbose=FALSE)
#' results.optim
#'
#' # Using approximate optimization method
#' results.approx <- RMSTSens(time="rfstime", status="status", exposure="hormon",
#'                            level.exposed="1", ps="Ps", data=dat, methods="Approx",
#'                            use.multicore=TRUE, n.core=2,
#'                            lambda=1.5, tau=365.25*5, ini.par=1, verbose=FALSE)
#' results.approx
#'
#' ## Adjusted RMST with not specified tau and with multiple lambda
#' # Using approximate optimization method
#' results.approx2 <- RMSTSens(time="rfstime", status="status", exposure="hormon",
#'                             level.exposed="1", ps="Ps", data=dat, methods="Approx",
#'                             use.multicore=TRUE, n.core=2,
#'                             lambda=c(1,1.5), tau=365.25*5, ini.par=1, verbose=FALSE)
#' results.approx2
#'
#' @author Seungjae Lee \email{seungjae2525@@gmail.com}
#'
#' @seealso
#'  \code{\link[RMSTSens]{print.RMSTSens}}, \code{\link[RMSTSens]{autoplot.RMSTSens}}, \code{\link[RMSTSens]{RMSTSens.ci}}
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
RMSTSens <- function(time, status,
                     exposure, level.exposed="1", ps, data,
                     methods="Approx", use.multicore=TRUE, n.core=parallel::detectCores()/2,
                     lambda=2, tau=NULL, ini.par=1, verbose=FALSE) {
  if (!(methods %in% c("Optim", "Approx", "LP1", "LP2"))) {
    stop("\n Error: Method must be one of \"Optim\", \"Approx\", \"LP1\", or \"LP2\".")
  }

  if (sum(data[, time] < 0) > 0) {
    stop("\n Error: Time must be positive.")
  }

  if (sum(data[, status] != 0 & data[, status] != 1) > 0) {
    stop("\n Error: Status must be 0 or 1.")
  }

  if (length(unique(data[, exposure])) != 2) {
    stop("\n Error: Exposure must be binary.")
  }

  if (sum(unique(data[, exposure]) %in% c(0, 1)) != 2) {
    stop("\n Error: Exposure must be 0 or 1.")
  }

  if (is.character(ps)) {
    if (sum(data[,ps] < 0) + sum(data[,ps] > 1) > 0) {
      stop("\n Error: Propensity score must be between 0 and 1.")
    }
    if (sum(data[,ps] == 0) + sum(data[,ps] == 1) > 0) {
      stop("\n Error: There is a propensity score of 0 or 1.")
    }
  } else if (is.vector(ps)) {
    if (sum(ps < 0) + sum(ps > 1) > 0) {
      stop("\n Error: Propensity score must be between 0 and 1.")
    }
    if (sum(ps == 0) + sum(ps == 1) > 0) {
      stop("\n Error: There is a propensity score of 0 or 1.")
    }
    data$ps <- ps
    ps <- "ps"
  } else {
    stop("\n Error: Propensity scores are not specified.")
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

  if (any(lambda < 1)) {
    stop("\n Error: \"lambda\" must be larger than or equal to 1.")
  } else if (length(lambda) > 1 & sum(lambda == 1) == 0) {
    stop("\n Error: If \"lambda\" is a vector, then it must include 1.")
  }

  if(length(lambda) != length(unique(lambda))){
    warning("There is a duplicate in \"lambda\". We will remove the duplicate \"lambda\".")
    lambda <- unique(lambda)
  }

  ## Shut down an R parallel cluster
  if(!is.null(getDefaultCluster())) try(parallel::setDefaultCluster(NULL), silent = TRUE)
  if(requireNamespace("doParallel", quietly=TRUE)) doParallel::stopImplicitCluster()

  ## Sorting Lambda
  lambda <- sort(lambda)

  ## Set tau
  group <- unique(data[, exposure])
  .tau <- min(c(max(data[, time][data[, status] == 1 &
                                   data[, exposure] == group[group != level.exposed]]),
                max(data[, time][data[, status] == 1 &
                                   data[, exposure] == level.exposed])))
  if (is.null(tau)) {
    tau <- .tau
  } else if (tau > .tau) {
    ## If tau not specified, use the minimum of the largest observed event time in both groups.
    warning("\"tau\" may have to smaller than the minimum of the largest observed event time in both groups.")
  }

  ## Make dataset for optimization
  dat.exposed <- optim_data(data, time=time, status=status, exposure=exposure,
                            ps=ps, group=level.exposed)

  dat.unexposed <- optim_data(data, time=time, status=status, exposure=exposure,
                              ps=ps, group=group[group != level.exposed])

  dat.exposed.event <- dat.exposed[dat.exposed$status == 1, ]
  dat.unexposed.event <- dat.unexposed[dat.unexposed$status == 1, ]

  len.exposed <- nrow(dat.exposed.event)
  len.unexposed <- nrow(dat.unexposed.event)

  ## Fit model
  if(!is.vector(lambda)){
    if(lambda == 1){
      opt1_min_result1 <- opt1_max_result1 <- optim_f(1, dat.exposed, minmax="min", lambda=1, tau=tau)
      opt1_min_result0 <- opt1_max_result0 <- optim_f(1, dat.unexposed, minmax="min", lambda=1, tau=tau)

      ## Return results
      result.df <- data.frame(N=nrow(data),
                              N.exposed=nrow(dat.exposed), N.unexposed=nrow(dat.unexposed),
                              N.event.exposed=sum(dat.exposed$status == 1),
                              N.event.unexposed=sum(dat.unexposed$status == 1),
                              cen.rate=sum(data[, status] == 0)/nrow(data),
                              cen.rate.exposed=sum(dat.exposed$status == 0)/nrow(dat.exposed),
                              cen.rate.unexposed=sum(dat.unexposed$status == 0)/nrow(dat.unexposed),
                              Lambda=lambda, Tau=tau, Method=methods,
                              min.exposed=opt1_min_result1,
                              max.exposed=opt1_max_result1,
                              min.unexposed=opt1_min_result0,
                              max.unexposed=opt1_max_result0,
                              RMST.diff.min=opt1_min_result1 - opt1_max_result0,
                              RMST.diff.max=opt1_max_result1 - opt1_min_result0)

    } else {
      if (methods == "Optim") {
        if (ini.par > lambda | ini.par < (1/lambda)) {
          stop("\n Error: \"ini.par\" must be between 1/Lambda and Lambda.")
        } else {
          ## If the optimParallel package is not available
          if (use.multicore == FALSE) {
            ## Not using multi-core

            ## Set initial optimization parameter zi
            par.exposed <- rep(ini.par, times=len.exposed)
            par.unexposed <- rep(ini.par, times=len.unexposed)

            # Minimum of shifted RMST for exposed group
            opt1_min_result1 <- stats::optim(par=par.exposed, data=dat.exposed, fn=optim_f,
                                             lower=1/lambda, upper=lambda, method="L-BFGS-B",
                                             minmax="min", lambda=lambda, tau=tau,
                                             control=list(fnscale=1, maxit=1000), hessian=FALSE)
            if(verbose==TRUE) cat('End: Minimum of shifted RMST for exposed group \n')

            # Maximum of shifted RMST for exposed group
            opt1_max_result1 <- stats::optim(par=par.exposed, data=dat.exposed, fn=optim_f,
                                             lower=1/lambda, upper=lambda, method="L-BFGS-B",
                                             minmax="max", lambda=lambda, tau=tau,
                                             control=list(fnscale=-1, maxit=1000), hessian=FALSE)
            if(verbose==TRUE) cat('End: Maximum of shifted RMST for exposed group \n')

            # Minimum of shifted RMST for unexposed group
            opt1_min_result0 <- stats::optim(par=par.unexposed, data=dat.unexposed, fn=optim_f,
                                             lower=1/lambda, upper=lambda, method="L-BFGS-B",
                                             minmax="min", lambda=lambda, tau=tau,
                                             control=list(fnscale=1, maxit=1000), hessian=FALSE)
            if(verbose==TRUE) cat('End: Minimum of shifted RMST for unexposed group \n')

            # Maximum of shifted RMST for unexposed group
            opt1_max_result0 <- stats::optim(par=par.unexposed, data=dat.unexposed, fn=optim_f,
                                             lower=1/lambda, upper=lambda, method="L-BFGS-B",
                                             minmax="max", lambda=lambda, tau=tau,
                                             control=list(fnscale=-1, maxit=1000), hessian=FALSE)
            if(verbose==TRUE) cat('End: Maximum of shifted RMST for unexposed group \n')

          } else {
            ## Using multi-core
            if (tolower(.Platform$OS.type) != "windows") {
              cl <- makeCluster(spec=n.core, type="FORK", outfile="")
            } else {
              cl <- makeCluster(spec=n.core, outfile="")
            }
            setDefaultCluster(cl=cl)

            ## Set initial parameters
            par.exposed <- rep(ini.par, times=len.exposed)
            par.unexposed <- rep(ini.par, times=len.unexposed)

            # Minimum of shifted RMST for exposed group
            opt1_min_result1 <- optimParallel(par=par.exposed, data=dat.exposed, fn=optim_f,
                                              lower=1/lambda, upper=lambda, method="L-BFGS-B",
                                              minmax="min", lambda=lambda, tau=tau,
                                              control=list(fnscale=1, maxit=1000), hessian=FALSE,
                                              parallel=list(loginfo=FALSE, forward=TRUE))
            if(verbose==TRUE) cat('End: Minimum of shifted RMST for exposed group \n')

            # Maximum of shifted RMST for exposed group
            opt1_max_result1 <- optimParallel(par=par.exposed, data=dat.exposed, fn=optim_f,
                                              lower=1/lambda, upper=lambda, method="L-BFGS-B",
                                              minmax="max", lambda=lambda, tau=tau,
                                              control=list(fnscale=-1, maxit=1000), hessian=FALSE,
                                              parallel=list(loginfo=FALSE, forward=TRUE))
            if(verbose==TRUE) cat('End: Maximum of shifted RMST for exposed group \n')

            # Minimum of shifted RMST for unexposed group
            opt1_min_result0 <- optimParallel(par=par.unexposed, data=dat.unexposed, fn=optim_f,
                                              lower=1/lambda, upper=lambda, method="L-BFGS-B",
                                              minmax="min", lambda=lambda, tau=tau,
                                              control=list(fnscale=1, maxit=1000), hessian=FALSE,
                                              parallel=list(loginfo=FALSE, forward=TRUE))
            if(verbose==TRUE) cat('End: Minimum of shifted RMST for unexposed group \n')

            # Maximum of shifted RMST for unexposed group
            opt1_max_result0 <- optimParallel(par=par.unexposed, data=dat.unexposed, fn=optim_f,
                                              lower=1/lambda, upper=lambda, method="L-BFGS-B",
                                              minmax="max", lambda=lambda, tau=tau,
                                              control=list(fnscale=-1, maxit=1000), hessian=FALSE,
                                              parallel=list(loginfo=FALSE, forward=TRUE))
            if(verbose==TRUE) cat('End: Maximum of shifted RMST for unexposed group \n')

            setDefaultCluster(cl=NULL)
            stopCluster(cl)
          }
        }

        ## Return results
        result.df <- data.frame(N=nrow(data),
                                N.exposed=nrow(dat.exposed), N.unexposed=nrow(dat.unexposed),
                                N.event.exposed=sum(dat.exposed$status == 1),
                                N.event.unexposed=sum(dat.unexposed$status == 1),
                                cen.rate=sum(data[, status] == 0)/nrow(data),
                                cen.rate.exposed=sum(dat.exposed$status == 0)/nrow(dat.exposed),
                                cen.rate.unexposed=sum(dat.unexposed$status == 0)/nrow(dat.unexposed),
                                Lambda=lambda, Tau=tau, Method=methods,
                                min.exposed=opt1_min_result1$value,
                                max.exposed=opt1_max_result1$value,
                                min.unexposed=opt1_min_result0$value,
                                max.unexposed=opt1_max_result0$value,
                                RMST.diff.min=opt1_min_result1$value - opt1_max_result0$value,
                                RMST.diff.max=opt1_max_result1$value - opt1_min_result0$value)

      } else if (methods == "Approx") {
        ## Set initial optimization parameter zi
        par.exposed.min <- rep(1/lambda, times=len.exposed)
        par.exposed.max <- rep(lambda, times=len.exposed)

        par.unexposed.min <- rep(1/lambda, times=len.unexposed)
        par.unexposed.max <- rep(lambda, times=len.unexposed)

        ## Make the list of parameters
        aa1 <- aa2 <- aa3 <- aa4 <- list()
        for (j in 1:len.exposed) {
          aa1[[j]] <- par.exposed.min
          par.exposed.min[j] <- lambda

          aa2[[j]] <- par.exposed.max
          par.exposed.max[j] <- 1/lambda
        }
        for (j in 1:len.unexposed) {
          aa3[[j]] <- par.unexposed.min
          par.unexposed.min[j] <- lambda

          aa4[[j]] <- par.unexposed.max
          par.unexposed.max[j] <- 1/lambda
        }

        if (use.multicore == FALSE){
          ## Not using multi-core
          # Minimum of shifted RMST for exposed group
          re1 <- sapply(aa1, function(x) optim_f(x, dat.exposed, minmax="min", lambda=lambda, tau=tau))
          if(verbose==TRUE) cat('End: Minimum of shifted RMST for exposed group \n')

          # Maximum of shifted RMST for exposed group
          re2 <- sapply(aa2, function(x) optim_f(x, dat.exposed, minmax="max", lambda=lambda, tau=tau))
          if(verbose==TRUE) cat('End: Maximum of shifted RMST for exposed group \n')

          # Minimum of shifted RMST for unexposed group
          re3 <- sapply(aa3, function(x) optim_f(x, dat.unexposed, minmax="min", lambda=lambda, tau=tau))
          if(verbose==TRUE) cat('End: Minimum of shifted RMST for unexposed group \n')

          # Maximum of shifted RMST for unexposed group
          re4 <- sapply(aa4, function(x) optim_f(x, dat.unexposed, minmax="max", lambda=lambda, tau=tau))
          if(verbose==TRUE) cat('End: Maximum of shifted RMST for unexposed group \n')

        } else {
          ## Using multi-core
          if (tolower(.Platform$OS.type) != "windows") {
            cl <- makeCluster(spec=getOption("cl.cores", n.core), type="FORK", outfile="")
          } else {
            cl <- makeCluster(spec=getOption("cl.cores", n.core), outfile="")
          }
          setDefaultCluster(cl=cl)
          clusterExport(cl=cl, envir=environment(),
                        varlist=c("optim_f", "dat.exposed", "dat.unexposed", "lambda", "tau"))

          ## Optimization using "parSapply" function
          # Minimum of shifted RMST for exposed group
          re1 <- parSapply(cl, aa1, function(x) optim_f(x, dat.exposed, minmax="min", lambda=lambda, tau=tau))
          if(verbose==TRUE) cat('End: Minimum of shifted RMST for exposed group \n')

          # Maximum of shifted RMST for exposed group
          re2 <- parSapply(cl, aa2, function(x) optim_f(x, dat.exposed, minmax="max", lambda=lambda, tau=tau))
          if(verbose==TRUE) cat('End: Maximum of shifted RMST for exposed group \n')

          # Minimum of shifted RMST for unexposed group
          re3 <- parSapply(cl, aa3, function(x) optim_f(x, dat.unexposed, minmax="min", lambda=lambda, tau=tau))
          if(verbose==TRUE) cat('End: Minimum of shifted RMST for unexposed group \n')

          # Maximum of shifted RMST for unexposed group
          re4 <- parSapply(cl, aa4, function(x) optim_f(x, dat.unexposed, minmax="max", lambda=lambda, tau=tau))
          if(verbose==TRUE) cat('End: Maximum of shifted RMST for unexposed group \n')

          setDefaultCluster(cl=NULL)
          stopCluster(cl)
        }

        ## Return results
        result.df <- data.frame(N=nrow(data),
                                N.exposed=nrow(dat.exposed), N.unexposed=nrow(dat.unexposed),
                                N.event.exposed=sum(dat.exposed$status == 1),
                                N.event.unexposed=sum(dat.unexposed$status == 1),
                                cen.rate=sum(data[, status] == 0)/nrow(data),
                                cen.rate.exposed=sum(dat.exposed$status == 0)/nrow(dat.exposed),
                                cen.rate.unexposed=sum(dat.unexposed$status == 0)/nrow(dat.unexposed),
                                Lambda=lambda, Tau=tau, Method=methods,
                                min.exposed=min(re1),
                                max.exposed=max(re2),
                                min.unexposed=min(re3),
                                max.unexposed=max(re4),
                                RMST.diff.min=min(re1) - max(re4),
                                RMST.diff.max=max(re2) - min(re3))

      } else if (methods == "LP1"){
        ## Minimum censoring time
        min.cen.time <- min(data[, time][data[, status] == 0])
        ## Maximum event time
        max.event.time <- max(data[, time][data[, status] == 1])

        ## When a closed cohort that the study entry times are the same for all subjects and
        #  there is no censoring apart from administrative censoring at the end of follow-up,
        #  the optimization problem can be reduced to linear programme.
        if (min.cen.time < max.event.time) {
          stop("\n Error: There must not be loss to follow-up nor other early censoring in the data in both groups (only administrative censoring at the end of follow-up).")
        } else {
          ## Minimum of shifted RMST for exposed group
          re1 <- optim_LP(dat.exposed, minmax="min", lambda=lambda, tau=tau)
          if(verbose==TRUE) cat('End: Minimum of shifted RMST for exposed group \n')

          ## Maximum of shifted RMST for exposed group
          re2 <- optim_LP(dat.exposed, minmax="max", lambda=lambda, tau=tau)
          if(verbose==TRUE) cat('End: Maximum of shifted RMST for exposed group \n')

          ## Minimum of shifted RMST for unexposed group
          re3 <- optim_LP(dat.unexposed, minmax="min", lambda=lambda, tau=tau)
          if(verbose==TRUE) cat('End: Minimum of shifted RMST for unexposed group \n')

          ## Maximum of shifted RMST for unexposed group
          re4 <- optim_LP(dat.unexposed, minmax="max", lambda=lambda, tau=tau)
          if(verbose==TRUE) cat('End: Maximum of shifted RMST for unexposed group \n')
        }

        ## Return results
        result.df <- data.frame(N=nrow(data),
                                N.exposed=nrow(dat.exposed), N.unexposed=nrow(dat.unexposed),
                                N.event.exposed=sum(dat.exposed$status == 1),
                                N.event.unexposed=sum(dat.unexposed$status == 1),
                                cen.rate=sum(data[, status] == 0)/nrow(data),
                                cen.rate.exposed=sum(dat.exposed$status == 0)/nrow(dat.exposed),
                                cen.rate.unexposed=sum(dat.unexposed$status == 0)/nrow(dat.unexposed),
                                Lambda=lambda, Tau=tau, Method=methods,
                                min.exposed=re1,
                                max.exposed=re2,
                                min.unexposed=re3,
                                max.unexposed=re4,
                                RMST.diff.min=re1 - re4,
                                RMST.diff.max=re2 - re3)

      } else if (methods == "LP2"){
        ## Minimum censoring time
        min.cen.time <- min(data[, time][data[, status] == 0])
        ## Maximum event time
        max.event.time <- max(data[, time][data[, status] == 1])

        ## When a closed cohort that the study entry times are the same for all subjects,
        #  if the minimum censoring time is longer than or equal to the pre-specified time point,
        #  then, the optimization problem can be reduced to linear programme.
        if(min.cen.time < tau){
          stop("\n Error: The minimum censoring time must be larger or equal to the pre-specified time point.")
        } else {
          ## Minimum of shifted RMST for exposed group
          re1 <- optim_LP(dat.exposed, minmax="min", lambda=lambda, tau=tau)
          if(verbose==TRUE) cat('End: Minimum of shifted RMST for exposed group \n')

          ## Maximum of shifted RMST for exposed group
          re2 <- optim_LP(dat.exposed, minmax="max", lambda=lambda, tau=tau)
          if(verbose==TRUE) cat('End: Maximum of shifted RMST for exposed group \n')

          ## Minimum of shifted RMST for unexposed group
          re3 <- optim_LP(dat.unexposed, minmax="min", lambda=lambda, tau=tau)
          if(verbose==TRUE) cat('End: Minimum of shifted RMST for unexposed group \n')

          ## Maximum of shifted RMST for unexposed group
          re4 <- optim_LP(dat.unexposed, minmax="max", lambda=lambda, tau=tau)
          if(verbose==TRUE) cat('End: Maximum of shifted RMST for unexposed group \n')
        }

        ## Return results
        result.df <- data.frame(N=nrow(data),
                                N.exposed=nrow(dat.exposed), N.unexposed=nrow(dat.unexposed),
                                N.event.exposed=sum(dat.exposed$status == 1),
                                N.event.unexposed=sum(dat.unexposed$status == 1),
                                cen.rate=sum(data[, status] == 0)/nrow(data),
                                cen.rate.exposed=sum(dat.exposed$status == 0)/nrow(dat.exposed),
                                cen.rate.unexposed=sum(dat.unexposed$status == 0)/nrow(dat.unexposed),
                                Lambda=lambda, Tau=tau, Method=methods,
                                min.exposed=re1,
                                max.exposed=re2,
                                min.unexposed=re3,
                                max.unexposed=re4,
                                RMST.diff.min=re1 - re4,
                                RMST.diff.max=re2 - re3)

      }
    }

  } else {
    ## If Lmabda is a vector
    result.df <- data.frame(matrix(NA, ncol = 17))
    colnames(result.df) <- c("N","N.exposed","N.unexposed","N.event.exposed","N.event.unexposed",
                             "cen.rate","cen.rate.exposed","cen.rate.unexposed",
                             "Lambda","Tau","Method",
                             "min.exposed","max.exposed","min.unexposed","max.unexposed",
                             "RMST.diff.min","RMST.diff.max")

    for(i in 1:length(lambda)){

      if (lambda[i] == 1) {
        opt1_min_result1 <- opt1_max_result1 <- optim_f(1, dat.exposed, minmax="min", lambda=1, tau=tau)
        opt1_min_result0 <- opt1_max_result0 <- optim_f(1, dat.unexposed, minmax="min", lambda=1, tau=tau)

        ## Return results
        result.df[i,] <- data.frame(N=nrow(data),
                                    N.exposed=nrow(dat.exposed), N.unexposed=nrow(dat.unexposed),
                                    N.event.exposed=sum(dat.exposed$status == 1),
                                    N.event.unexposed=sum(dat.unexposed$status == 1),
                                    cen.rate=sum(data[, status] == 0)/nrow(data),
                                    cen.rate.exposed=sum(dat.exposed$status == 0)/nrow(dat.exposed),
                                    cen.rate.unexposed=sum(dat.unexposed$status == 0)/nrow(dat.unexposed),
                                    Lambda=lambda[i], Tau=tau, Method=methods,
                                    min.exposed=opt1_min_result1,
                                    max.exposed=opt1_max_result1,
                                    min.unexposed=opt1_min_result0,
                                    max.unexposed=opt1_max_result0,
                                    RMST.diff.min=opt1_min_result1 - opt1_max_result0,
                                    RMST.diff.max=opt1_max_result1 - opt1_min_result0)

      } else {
        if (methods == "Optim") {
          if (ini.par > lambda[i] | ini.par < (1/lambda[i])) {
            stop("\n Error: \"ini.par\" must be between 1/Lambda and Lambda.")
          } else {
            ## If the optimParallel package is not available
            if (use.multicore == FALSE) {
              ## Not using multi-core

              ## Set initial optimization parameter zi
              par.exposed <- rep(ini.par, times=len.exposed)
              par.unexposed <- rep(ini.par, times=len.unexposed)

              # Minimum of shifted RMST for exposed group
              opt1_min_result1 <- optim(par=par.exposed, data=dat.exposed, fn=optim_f,
                                        lower=1/lambda[i], upper=lambda[i], method="L-BFGS-B",
                                        minmax="min", lambda=lambda[i], tau=tau,
                                        control=list(fnscale=1, maxit=1000), hessian=FALSE)
              if(verbose==TRUE) cat('End: Minimum of shifted RMST for exposed group \n')

              # Maximum of shifted RMST for exposed group
              opt1_max_result1 <- optim(par=par.exposed, data=dat.exposed, fn=optim_f,
                                        lower=1/lambda[i], upper=lambda[i], method="L-BFGS-B",
                                        minmax="max", lambda=lambda[i], tau=tau,
                                        control=list(fnscale=-1, maxit=1000), hessian=FALSE)
              if(verbose==TRUE) cat('End: Maximum of shifted RMST for exposed group \n')

              # Minimum of shifted RMST for unexposed group
              opt1_min_result0 <- optim(par=par.unexposed, data=dat.unexposed, fn=optim_f,
                                        lower=1/lambda[i], upper=lambda[i], method="L-BFGS-B",
                                        minmax="min", lambda=lambda[i], tau=tau,
                                        control=list(fnscale=1, maxit=1000), hessian=FALSE)
              if(verbose==TRUE) cat('End: Minimum of shifted RMST for unexposed group \n')

              # Maximum of shifted RMST for unexposed group
              opt1_max_result0 <- optim(par=par.unexposed, data=dat.unexposed, fn=optim_f,
                                        lower=1/lambda[i], upper=lambda[i], method="L-BFGS-B",
                                        minmax="max", lambda=lambda[i], tau=tau,
                                        control=list(fnscale=-1, maxit=1000), hessian=FALSE)
              if(verbose==TRUE) cat('End: Maximum of shifted RMST for unexposed group \n')

            } else {
              ## Using multi-core
              if (tolower(.Platform$OS.type) != "windows") {
                cl <- makeCluster(spec=n.core, type="FORK", outfile="")
              } else {
                cl <- makeCluster(spec=n.core, outfile="")
              }
              setDefaultCluster(cl=cl)

              ## Set initial parameters
              par.exposed <- rep(ini.par, times=len.exposed)
              par.unexposed <- rep(ini.par, times=len.unexposed)

              # Minimum of shifted RMST for exposed group
              opt1_min_result1 <- optimParallel(par=par.exposed, data=dat.exposed, fn=optim_f,
                                                lower=1/lambda[i], upper=lambda[i], method="L-BFGS-B",
                                                minmax="min", lambda=lambda[i], tau=tau,
                                                control=list(fnscale=1, maxit=1000), hessian=FALSE,
                                                parallel=list(loginfo=FALSE, forward=TRUE))
              if(verbose==TRUE) cat('End: Minimum of shifted RMST for exposed group \n')

              # Maximum of shifted RMST for exposed group
              opt1_max_result1 <- optimParallel(par=par.exposed, data=dat.exposed, fn=optim_f,
                                                lower=1/lambda[i], upper=lambda[i], method="L-BFGS-B",
                                                minmax="max", lambda=lambda[i], tau=tau,
                                                control=list(fnscale=-1, maxit=1000), hessian=FALSE,
                                                parallel=list(loginfo=FALSE, forward=TRUE))
              if(verbose==TRUE) cat('End: Maximum of shifted RMST for exposed group \n')

              # Minimum of shifted RMST for unexposed group
              opt1_min_result0 <- optimParallel(par=par.unexposed, data=dat.unexposed, fn=optim_f,
                                                lower=1/lambda[i], upper=lambda[i], method="L-BFGS-B",
                                                minmax="min", lambda=lambda[i], tau=tau,
                                                control=list(fnscale=1, maxit=1000), hessian=FALSE,
                                                parallel=list(loginfo=FALSE, forward=TRUE))
              if(verbose==TRUE) cat('End: Minimum of shifted RMST for unexposed group \n')

              # Maximum of shifted RMST for unexposed group
              opt1_max_result0 <- optimParallel(par=par.unexposed, data=dat.unexposed, fn=optim_f,
                                                lower=1/lambda[i], upper=lambda[i], method="L-BFGS-B",
                                                minmax="max", lambda=lambda[i], tau=tau,
                                                control=list(fnscale=-1, maxit=1000), hessian=FALSE,
                                                parallel=list(loginfo=FALSE, forward=TRUE))
              if(verbose==TRUE) cat('End: Maximum of shifted RMST for unexposed group \n')

              setDefaultCluster(cl=NULL)
              stopCluster(cl)
            }
          }

          ## Return results
          result.df[i,] <- data.frame(N=nrow(data),
                                      N.exposed=nrow(dat.exposed), N.unexposed=nrow(dat.unexposed),
                                      N.event.exposed=sum(dat.exposed$status == 1),
                                      N.event.unexposed=sum(dat.unexposed$status == 1),
                                      cen.rate=sum(data[, status] == 0)/nrow(data),
                                      cen.rate.exposed=sum(dat.exposed$status == 0)/nrow(dat.exposed),
                                      cen.rate.unexposed=sum(dat.unexposed$status == 0)/nrow(dat.unexposed),
                                      Lambda=lambda[i], Tau=tau, Method=methods,
                                      min.exposed=opt1_min_result1$value,
                                      max.exposed=opt1_max_result1$value,
                                      min.unexposed=opt1_min_result0$value,
                                      max.unexposed=opt1_max_result0$value,
                                      RMST.diff.min=opt1_min_result1$value - opt1_max_result0$value,
                                      RMST.diff.max=opt1_max_result1$value - opt1_min_result0$value)

        } else if (methods == "Approx") {
          ## Set initial optimization parameter zi
          par.exposed.min <- rep(1/lambda[i], times=len.exposed)
          par.exposed.max <- rep(lambda[i], times=len.exposed)

          par.unexposed.min <- rep(1/lambda[i], times=len.unexposed)
          par.unexposed.max <- rep(lambda[i], times=len.unexposed)

          ## Make the list of parameters
          aa1 <- aa2 <- aa3 <- aa4 <- list()
          for (j in 1:len.exposed) {
            aa1[[j]] <- par.exposed.min
            par.exposed.min[j] <- lambda[i]

            aa2[[j]] <- par.exposed.max
            par.exposed.max[j] <- 1/lambda[i]
          }
          for (j in 1:len.unexposed) {
            aa3[[j]] <- par.unexposed.min
            par.unexposed.min[j] <- lambda[i]

            aa4[[j]] <- par.unexposed.max
            par.unexposed.max[j] <- 1/lambda[i]
          }

          if (use.multicore == FALSE){
            ## Not using multi-core
            # Minimum of shifted RMST for exposed group
            re1 <- sapply(aa1, function(x) optim_f(x, dat.exposed, minmax="min", lambda=lambda[i], tau=tau))
            if(verbose==TRUE) cat('End: Minimum of shifted RMST for exposed group \n')

            # Maximum of shifted RMST for exposed group
            re2 <- sapply(aa2, function(x) optim_f(x, dat.exposed, minmax="max", lambda=lambda[i], tau=tau))
            if(verbose==TRUE) cat('End: Maximum of shifted RMST for exposed group \n')

            # Minimum of shifted RMST for unexposed group
            re3 <- sapply(aa3, function(x) optim_f(x, dat.unexposed, minmax="min", lambda=lambda[i], tau=tau))
            if(verbose==TRUE) cat('End: Minimum of shifted RMST for unexposed group \n')

            # Maximum of shifted RMST for unexposed group
            re4 <- sapply(aa4, function(x) optim_f(x, dat.unexposed, minmax="max", lambda=lambda[i], tau=tau))
            if(verbose==TRUE) cat('End: Maximum of shifted RMST for unexposed group \n')

          } else {
            ## Using multi-core
            if (tolower(.Platform$OS.type) != "windows") {
              cl <- makeCluster(spec=getOption("cl.cores", n.core), type="FORK", outfile="")
            } else {
              cl <- makeCluster(spec=getOption("cl.cores", n.core), outfile="")
            }
            setDefaultCluster(cl=cl)
            clusterExport(cl=cl, envir=environment(),
                          varlist=c("optim_f", "dat.exposed", "dat.unexposed", "lambda", "tau"))

            ## Optimization using "parSapply" function
            # Minimum of shifted RMST for exposed group
            re1 <- parSapply(cl, aa1, function(x) optim_f(x, dat.exposed, minmax="min", lambda=lambda[i], tau=tau))
            if(verbose==TRUE) cat('End: Minimum of shifted RMST for exposed group \n')

            # Maximum of shifted RMST for exposed group
            re2 <- parSapply(cl, aa2, function(x) optim_f(x, dat.exposed, minmax="max", lambda=lambda[i], tau=tau))
            if(verbose==TRUE) cat('End: Maximum of shifted RMST for exposed group \n')

            # Minimum of shifted RMST for unexposed group
            re3 <- parSapply(cl, aa3, function(x) optim_f(x, dat.unexposed, minmax="min", lambda=lambda[i], tau=tau))
            if(verbose==TRUE) cat('End: Minimum of shifted RMST for unexposed group \n')

            # Maximum of shifted RMST for unexposed group
            re4 <- parSapply(cl, aa4, function(x) optim_f(x, dat.unexposed, minmax="max", lambda=lambda[i], tau=tau))
            if(verbose==TRUE) cat('End: Maximum of shifted RMST for unexposed group \n')

            setDefaultCluster(cl=NULL)
            stopCluster(cl)
          }

          ## Return results
          result.df[i,] <- data.frame(N=nrow(data),
                                      N.exposed=nrow(dat.exposed), N.unexposed=nrow(dat.unexposed),
                                      N.event.exposed=sum(dat.exposed$status == 1),
                                      N.event.unexposed=sum(dat.unexposed$status == 1),
                                      cen.rate=sum(data[, status] == 0)/nrow(data),
                                      cen.rate.exposed=sum(dat.exposed$status == 0)/nrow(dat.exposed),
                                      cen.rate.unexposed=sum(dat.unexposed$status == 0)/nrow(dat.unexposed),
                                      Lambda=lambda[i], Tau=tau, Method=methods,
                                      min.exposed=min(re1),
                                      max.exposed=max(re2),
                                      min.unexposed=min(re3),
                                      max.unexposed=max(re4),
                                      RMST.diff.min=min(re1) - max(re4),
                                      RMST.diff.max=max(re2) - min(re3))

        } else if (methods == "LP1"){
          ## Minimum censoring time
          min.cen.time <- min(data[, time][data[, status] == 0])
          ## Maximum event time
          max.event.time <- max(data[, time][data[, status] == 1])

          ## When a closed cohort that the study entry times are the same for all subjects and
          #  there is no censoring apart from administrative censoring at the end of follow-up,
          #  the optimization problem can be reduced to linear programme.
          if (min.cen.time < max.event.time) {
            stop("\n Error: There must not be loss to follow-up nor other early censoring in the data in both groups (only administrative censoring at the end of follow-up).")
          } else {
            ## Minimum of shifted RMST for exposed group
            re1 <- optim_LP(dat.exposed, minmax="min", lambda=lambda[i], tau=tau)
            if(verbose==TRUE) cat('End: Minimum of shifted RMST for exposed group \n')

            ## Maximum of shifted RMST for exposed group
            re2 <- optim_LP(dat.exposed, minmax="max", lambda=lambda[i], tau=tau)
            if(verbose==TRUE) cat('End: Maximum of shifted RMST for exposed group \n')

            ## Minimum of shifted RMST for unexposed group
            re3 <- optim_LP(dat.unexposed, minmax="min", lambda=lambda[i], tau=tau)
            if(verbose==TRUE) cat('End: Minimum of shifted RMST for unexposed group \n')

            ## Maximum of shifted RMST for unexposed group
            re4 <- optim_LP(dat.unexposed, minmax="max", lambda=lambda[i], tau=tau)
            if(verbose==TRUE) cat('End: Maximum of shifted RMST for unexposed group \n')
          }

          ## Return results
          result.df[i,] <- data.frame(N=nrow(data),
                                      N.exposed=nrow(dat.exposed), N.unexposed=nrow(dat.unexposed),
                                      N.event.exposed=sum(dat.exposed$status == 1),
                                      N.event.unexposed=sum(dat.unexposed$status == 1),
                                      cen.rate=sum(data[, status] == 0)/nrow(data),
                                      cen.rate.exposed=sum(dat.exposed$status == 0)/nrow(dat.exposed),
                                      cen.rate.unexposed=sum(dat.unexposed$status == 0)/nrow(dat.unexposed),
                                      Lambda=lambda[i], Tau=tau, Method=methods,
                                      min.exposed=re1,
                                      max.exposed=re2,
                                      min.unexposed=re3,
                                      max.unexposed=re4,
                                      RMST.diff.min=re1 - re4,
                                      RMST.diff.max=re2 - re3)

        } else if (methods == "LP2"){
          ## Minimum censoring time
          min.cen.time <- min(data[, time][data[, status] == 0])
          ## Maximum event time
          max.event.time <- max(data[, time][data[, status] == 1])

          ## When a closed cohort that the study entry times are the same for all subjects,
          #  if the minimum censoring time is longer than or equal to the pre-specified time point,
          #  then, the optimization problem can be reduced to linear programme.
          if(min.cen.time < tau){
            stop("\n Error: The minimum censoring time must be larger or equal to the pre-specified time point.")
          } else {
            ## Minimum of shifted RMST for exposed group
            re1 <- optim_LP(dat.exposed, minmax="min", lambda=lambda[i], tau=tau)
            if(verbose==TRUE) cat('End: Minimum of shifted RMST for exposed group \n')

            ## Maximum of shifted RMST for exposed group
            re2 <- optim_LP(dat.exposed, minmax="max", lambda=lambda[i], tau=tau)
            if(verbose==TRUE) cat('End: Maximum of shifted RMST for exposed group \n')

            ## Minimum of shifted RMST for unexposed group
            re3 <- optim_LP(dat.unexposed, minmax="min", lambda=lambda[i], tau=tau)
            if(verbose==TRUE) cat('End: Minimum of shifted RMST for unexposed group \n')

            ## Maximum of shifted RMST for unexposed group
            re4 <- optim_LP(dat.unexposed, minmax="max", lambda=lambda[i], tau=tau)
            if(verbose==TRUE) cat('End: Maximum of shifted RMST for unexposed group \n')
          }

          ## Return results
          result.df[i,] <- data.frame(N=nrow(data),
                                      N.exposed=nrow(dat.exposed), N.unexposed=nrow(dat.unexposed),
                                      N.event.exposed=sum(dat.exposed$status == 1),
                                      N.event.unexposed=sum(dat.unexposed$status == 1),
                                      cen.rate=sum(data[, status] == 0)/nrow(data),
                                      cen.rate.exposed=sum(dat.exposed$status == 0)/nrow(dat.exposed),
                                      cen.rate.unexposed=sum(dat.unexposed$status == 0)/nrow(dat.unexposed),
                                      Lambda=lambda[i], Tau=tau, Method=methods,
                                      min.exposed=re1,
                                      max.exposed=re2,
                                      min.unexposed=re3,
                                      max.unexposed=re4,
                                      RMST.diff.min=re1 - re4,
                                      RMST.diff.max=re2 - re3)

        }
      }
    }
  }

  data <- data
  argument <- c(time, status, exposure, level.exposed, ps, ini.par)

  final.result <- list(data=data, argument=argument, result.df=result.df)

  class(final.result) <- c("RMSTSens")
  return(final.result)
}

