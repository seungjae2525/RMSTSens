## Optimization function
optim_f <- function(par, data, minmax="min", lambda=2, tau=NULL) {
  time <- data$time
  status <- data$status
  exp_logit_Pr <- data$exp_logit_Pr
  lambda <- ifelse(minmax == "min", 1/lambda, lambda)

  ## If tau not specified, use the minimum of the largest observed event time in both groups
  if (is.null(tau)) {
    tau <- max(data$time[data$status == 1])
  }

  ## Calculate adjusted Kaplan-Meier curve
  tj <- c(0, unique(time[status == 1]))
  djyj <- sapply(tj, function(x) {
    dj <- sum((1 + par * exp_logit_Pr[status == 1])[time[status == 1] == x])
    yj <- sum((1 + par * exp_logit_Pr[status == 1])[time[status == 1] >= x]) +
      sum((1 + lambda * exp_logit_Pr[status == 0])[time[status == 0] >= x])

    c(dj=dj, yj=yj)
  })

  dj <- djyj[1, ]
  yj <- djyj[2, ]
  st <- cumprod(1 - (dj/yj))

  ## RMST
  rtime <- tj <= tau
  tj_r <- sort(c(tj[rtime], tau))
  st_r <- st[rtime]
  time_diff <- diff(c(0, tj_r))
  areas <- time_diff * c(1, st_r)
  rmst <- sum(areas)

  return(rmst)
}


## Optimization function via liear programming
optim_LP <- function(data, minmax="min", lambda=3, tau=NULL) {
  time <- data$time
  status <- data$status
  exp_logit_Pr <- data$exp_logit_Pr

  ## If tau not specified, use the minimum of the largest observed event time in both groups
  if (is.null(tau)) {
    tau <- max(data$time[data$status == 1])
  }

  ## Redefine time using pre-specified tau
  time <- ifelse(status == 1 & time < tau, time, tau)

  ## Optimization using Proposition 2 in Zhao et al. (2019)
  if (!(minmax %in% c("min", "max"))) {
    cat("Error: minmax must be \"min\" or \"max\"")
  } else if (minmax == "min") {
    zi1 <- c(rep(lambda, sum(status == 1)), rep(1/lambda, sum(status == 0)))
    zi2 <- c(rep(1/lambda, nrow(data)))

    ## Numerator
    num.up <- (1 + zi1 * exp_logit_Pr) * time
    num.low <- (1 + zi2 * exp_logit_Pr) * time

    ## Denominator
    den.up <- (1 + zi1 * exp_logit_Pr)
    den.low <- (1 + zi2 * exp_logit_Pr)

    num.min <- c(0, cumsum(num.up)) + c(rev(cumsum(rev(num.low))), 0)
    den.min <- c(0, cumsum(den.up)) + c(rev(cumsum(rev(den.low))), 0)
    optimum <- min(num.min/den.min)
  } else if (minmax == "max") {
    zi1 <- c(rep(lambda, nrow(data)))
    zi2 <- c(rep(1/lambda, sum(status == 1)), rep(lambda, sum(status == 0)))

    ## Numerator
    num.up <- (1 + zi1 * exp_logit_Pr) * time
    num.low <- (1 + zi2 * exp_logit_Pr) * time

    ## Denominator
    den.up <- (1 + zi1 * exp_logit_Pr)
    den.low <- (1 + zi2 * exp_logit_Pr)

    num.max <- c(0, cumsum(num.low)) + c(rev(cumsum(rev(num.up))), 0)
    den.max <- c(0, cumsum(den.low)) + c(rev(cumsum(rev(den.up))), 0)
    optimum <- max(num.max/den.max)
  }

  return(optimum)
}

