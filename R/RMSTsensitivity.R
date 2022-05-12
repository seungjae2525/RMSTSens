#' @title Sensitivity analysis for RMST
#'
#' @description Function for sensitivity analysis of unmeasured confounding for restricted mean survival time using propensity score
#'
#' @param time The name of the variable for time to event
#' @param status The name of the variable for status (0 if censored, 1 if event)
#' @param exposure The name of the variable for exposure (0 if unexposed, 1 if exposed)
#' @param exposed.ref.level Reference level in exposure variable, Default: 1
#' @param ps The name of the variable for propensity score variable i.e., P(A=1|L)>0
#' @param data A data frame in which contains the follow-up time (time), the event (status), the exposure (exposure), and the propensity score (ps)
#' @param methods A character with the methods how to calculate the adjusted RMST ("Optim", "Approx", "LP1", "LP2"), Default: 'Approx'
#' @param use.multicore Logical scalar indicating whether to parallelize our optimization problem, Default: TRUE
#' @param n.core The number of usable cores, Default: parallel::detectCores()/2
#' @param lambda The sensitivity parameter, Default: 1.5
#' @param tau User-specific time point, If tau not specified (NULL), use the minimum of the largest observed event time in both groups.
#' @param ini.par Initial parameter for direct optimization method, Default: 1
#' @param verbose Conditional on the verbose level, print the message that each optimization (minimization or maximization) for each group was ended, Default: FALSE
#'
#' @return An object of class \code{senRMST}. The object is a data.frame with the following components:
#' \item{N}{Total number of subjects}
#' \item{N.exposed}{The number of subjects in exposed group}
#' \item{N.unexposed}{The number of subjects in unexposed group}
#' \item{N.event.exposed}{The number of events in exposed group}
#' \item{N.event.unexposed}{The number of events in unexposed group}
#' \item{cen.rate}{Total censoring rate}
#' \item{cen.rate.exposed}{Censoring rate in exposed group}
#' \item{cen.rate.unexposed}{Censoring rate in unexposed group}
#' \item{Lambda}{A used sensitivity parameter }
#' \item{Tau}{User-specific time point, If tau not specified (NULL), use the minimum of the largest observed event time in both groups}
#' \item{Method}{A used method}
#' \item{min.exposed}{The minimum of adjusted RMST based on the shifted propensity score for exposed group}
#' \item{max.exposed}{The maximum of adjusted RMST based on the shifted propensity score for exposed group}
#' \item{min.unexposed}{The minimum of adjusted RMST based on the shifted propensity score for unexposed group}
#' \item{max.unexposed}{The maximum of adjusted RMST based on the shifted propensity score for unexposed group}
#' \item{RMST.diff.min}{The minimum of between-group difference in adjusted RMST based on shifted propensity score}
#' \item{RMST.diff.max}{The maximum of between-group difference in adjusted RMST based on shifted propensity score}
#' The results for the RMST_sensitivity are printed with the \code{\link{print.RMSTsensitivity}} functions.
#' To generate graphs comparing Lambda with interval of adjusted RMST based on shifted propensity score use the aaa functions.
#'
#' @details To assess details of method of sensitivity analysis, see Lee et al. (2022) for details.
#'
#' @examples
#' if(interactive()){
#'  ## EXAMPLE
#'  library(survival)
#'
#'  dat <- gbsg
#'  dat$size2 <- ifelse(dat$size <= 20, 0,
#'                      ifelse(dat$size > 20 & dat$size <= 50, 1, 2))
#'  dat$age2 <- dat$age/100
#'  dat$er2 <- dat$er/1000
#'
#'  ## Estimation of propensity score
#'  denom.fit <- glm(hormon ~ (age2)^3 + (age2)^3*log(age2) + meno + factor(size2) + sqrt(nodes) + er2,
#'                   data=dat, family=binomial(link='logit'))
#'  dat$Ps <- predict(denom.fit, type='response')
#'
#'  ## Between-group difference in adjusted RMST based on shifted propensity score
#'  # Using approximate optimization method
#'  results.approx <- RMST_sensitivity(time='rfstime', status='status', exposure='hormon',
#'                                     exposed.ref.level=1, ps='Ps' ,data=dat, methods="Approx",
#'                                     use.multicore=TRUE, n.core=parallel::detectCores()/2,
#'                                     lambda=1.5, tau=365.25*5)
#'  results.approx
#'  ## Adjusted RMST with tau equal to 5-year
#'  # Using direct optimization method
#'  results.optim <- RMST_sensitivity(time='rfstime', status='status', exposure='hormon',
#'                                    exposed.ref.level=1, data=dat, ps='Ps', methods="Optim",
#'                                    use.multicore=TRUE, n.core=parallel::detectCores()/2,
#'                                    lambda=1.5, tau=365.25*5)
#'  results.optim
#'  }
#'
#' @author Seungjae Lee \email{seungjae2525@@gmail.com}
#'
#' @seealso
#'  \code{\link[parallel]{detectCores}}, \code{\link[parallel]{makeCluster}}
#'
#' @rdname RMSTsensitivity
#'
#' @references
#' Bakbergenuly I, Hoaglin DC, Kulinskaya E (2020):
#' Methods for estimating between-study variance and overall
#' effect in meta-analysis of odds-ratios.
#' \emph{Research Synthesis Methods},
#' DOI: 10.1002/jrsm.1404
#'
#' @import parallel
#' @import optimParallel
#'
#' @export
RMSTsensitivity <- function(time, status,
                            exposure, exposed.ref.level=1, ps, data,
                            methods="Approx", use.multicore=TRUE, n.core=parallel::detectCores()/2,
                            lambda=2, tau=NULL, ini.par=1, verbose=FALSE) {
  if (!(methods %in% c("Optim", "Approx", "LP1", "LP2"))) {
    stop("\n Error: Method must be \"Optim\", \"Approx\", \"LP1\", or \"LP2\".")
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

  if (sum(data[,ps] < 0) + sum(data[,ps] > 1) > 0) {
    stop("\n Error: Propensity score must be between 0 and 1.")
  }

  if (use.multicore==TRUE){
    if(is.null(n.core)){
      n.core <- parallel::detectCores()/2
    }
  } else {
    if(is.null(n.core)){
      n.core <- 0
    }
  }

  if (use.multicore==TRUE & n.core > parallel::detectCores()) {
    stop("\n Error: \"n.core\" is more than the number of available cores.")
  }

  if (lambda < 1) {
    stop("\n Error: \"lambda\" must be larger than or equal to 1.")
  }

  ## Shut down an R parallel cluster
  if(!is.null(getDefaultCluster())) try(parallel::setDefaultCluster(NULL), silent = TRUE)
  if(!requireNamespace("doParallel", quietly=TRUE)) doParallel::stopImplicitCluster()

  f.dat <- data

  ## Set tau
  group <- unique(data[, exposure])
  .tau <- min(c(max(f.dat[, time][f.dat[, status] == 1 &
                                    f.dat[, exposure] == group[group != exposed.ref.level]]),
                max(f.dat[, time][f.dat[, status] == 1 &
                                    f.dat[, exposure] == exposed.ref.level])))
  if (is.null(tau)) {
    tau <- .tau
  } else if (tau > .tau) {
    ## If tau not specified, use the minimum of the largest observed event time in both groups.
    warning("\"tau\" may have to smaller than the minimum of the largest observed event time in both groups.\n")
  }

  ## Make dataset for optimization
  dat.exposed <- optim_data(data, time=time, status=status, exposure=exposure,
                            ps=ps, group=exposed.ref.level)

  dat.unexposed <- optim_data(data, time=time, status=status, exposure=exposure,
                              ps=ps, group=group[group != exposed.ref.level])

  dat.exposed.event <- dat.exposed[dat.exposed$status == 1, ]
  dat.unexposed.event <- dat.unexposed[dat.unexposed$status == 1, ]

  len.exposed <- nrow(dat.exposed.event)
  len.unexposed <- nrow(dat.unexposed.event)

  ## Fit model
  if (methods == "Optim") {
    if (ini.par > lambda | ini.par < (1/lambda)) {
      stop("\n Error: \"ini.par\" must be between 1/Lambda and Lambda.")
    } else {
      ## If the optimParallel package is not available
      if (use.multicore == FALSE) {
        ## Not using multi-core

        ## Set initial optimization parameter zi
        par.exposed <- rep(ini.par, len.exposed)
        par.unexposed <- rep(ini.par, len.unexposed)

        # Minimum of shifted RMST for exposed group
        opt1_min_result1 <- optim(par=par.exposed, data=dat.exposed, fn=optim_f,
                                  lower=1/lambda, upper=lambda, method="L-BFGS-B",
                                  minmax="min", lambda=lambda, tau=tau,
                                  control=list(fnscale=1, maxit=1000), hessian=FALSE)
        if(verbose==TRUE) cat('End: Minimum of shifted RMST for exposed group \n')

        # Maximum of shifted RMST for exposed group
        opt1_max_result1 <- optim(par=par.exposed, data=dat.exposed, fn=optim_f,
                                  lower=1/lambda, upper=lambda, method="L-BFGS-B",
                                  minmax="max", lambda=lambda, tau=tau,
                                  control=list(fnscale=-1, maxit=1000), hessian=FALSE)
        if(verbose==TRUE) cat('End: Maximum of shifted RMST for exposed group \n')

        # Minimum of shifted RMST for unexposed group
        opt1_min_result0 <- optim(par=par.unexposed, data=dat.unexposed, fn=optim_f,
                                  lower=1/lambda, upper=lambda, method="L-BFGS-B",
                                  minmax="min", lambda=lambda, tau=tau,
                                  control=list(fnscale=1, maxit=1000), hessian=FALSE)
        if(verbose==TRUE) cat('End: Minimum of shifted RMST for unexposed group \n')

        # Maximum of shifted RMST for unexposed group
        opt1_max_result0 <- optim(par=par.unexposed, data=dat.unexposed, fn=optim_f,
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
        par.exposed <- rep(ini.par, len.exposed)
        par.unexposed <- rep(ini.par, len.unexposed)

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
    result.df <- data.frame(N=nrow(f.dat),
                            N.exposed=nrow(dat.exposed), N.unexposed=nrow(dat.unexposed),
                            N.event.exposed=sum(dat.exposed$status == 1),
                            N.event.unexposed=sum(dat.unexposed$status == 1),
                            cen.rate=sum(f.dat[, status] == 0)/nrow(f.dat),
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
    par.exposed.min <- rep(1/lambda, len.exposed)
    par.exposed.max <- rep(lambda, len.exposed)

    par.unexposed.min <- rep(1/lambda, len.unexposed)
    par.unexposed.max <- rep(lambda, len.unexposed)

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
    result.df <- data.frame(N=nrow(f.dat),
                            N.exposed=nrow(dat.exposed), N.unexposed=nrow(dat.unexposed),
                            N.event.exposed=sum(dat.exposed$status == 1),
                            N.event.unexposed=sum(dat.unexposed$status == 1),
                            cen.rate=sum(f.dat[, status] == 0)/nrow(f.dat),
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
    min.cen.time <- min(f.dat[, time][f.dat[, status] == 0])
    ## Maximum event time
    max.event.time <- max(f.dat[, time][f.dat[, status] == 1])

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
    result.df <- data.frame(N=nrow(f.dat),
                            N.exposed=nrow(dat.exposed), N.unexposed=nrow(dat.unexposed),
                            N.event.exposed=sum(dat.exposed$status == 1),
                            N.event.unexposed=sum(dat.unexposed$status == 1),
                            cen.rate=sum(f.dat[, status] == 0)/nrow(f.dat),
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
    min.cen.time <- min(f.dat[, time][f.dat[, status] == 0])
    ## Maximum event time
    max.event.time <- max(f.dat[, time][f.dat[, status] == 1])

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
    result.df <- data.frame(N=nrow(f.dat),
                            N.exposed=nrow(dat.exposed), N.unexposed=nrow(dat.unexposed),
                            N.event.exposed=sum(dat.exposed$status == 1),
                            N.event.unexposed=sum(dat.unexposed$status == 1),
                            cen.rate=sum(f.dat[, status] == 0)/nrow(f.dat),
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

  class(result.df) <- c("senRMST")
  return(result.df)
}


