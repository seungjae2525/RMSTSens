#' @title Adjusted RMST
#' @description Function for RMST estimate from adjusted Kaplan-Meier curve via SIPW (or IPW)
#'
#' @param data A data frame in which contains the follow-up time (time), the event (status), the exposure (exposure), and the optional propensity score (ps)
#' @param time The name of the variable for time to event
#' @param status The name of the variable for status (0 if censored, 1 if event)
#' @param exposure The name of the variable for exposure (0 if unexposed, 1 if exposed)
#' @param ps The name of the variable for propensity score variable P(A=1|L)>0, Default: NULL
#' @param stabilize A logical value. If true, stabilized ipw., Default: TRUE
#' @param tau pre-specified time point to estimate of restricted mean survival time (RMST), Default: NULL
#' @param var.est A logical value. variance of RMST via xie (2005) is calculated., Default: FALSE
#' @param alpha A significance level, Default: 0.05
#'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  # EXAMPLE
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
#'  ## Adjusted RMST with tau qual to 5-year
#'  rr <- adj_rmst(data=dat, time='rfstime', status='status', exposure='hormon',
#'                 ps='Ps', stabilize=FALSE, tau=365.25*5, alpha=.05, var.est = TRUE)
#'  round(c(rr$rmst_diff, rr$rmst_diff_low, rr$rmst_diff_upp, rr$rmst_diff_pval), 3)
#'  }
#' }
#' @rdname adj_rmst
#' @import stats
#' @export

adj_rmst <- function(time, status, exposure, data,
                     ps=NULL, stabilize=TRUE, tau=NULL, var.est=FALSE, alpha=0.05) {
  # Based on "akm_RMST.R" from https://github.com/s-conner/akm-rmst
  # Author: Conner, S. C., et al.
  # Ref: Adjusted restricted mean survival times in observational studies (2019).

  if (sum(data[, time] < 0) > 0) {
    stop("\n Error: Time must be positive.")
  } else if (sum(data[,ps] < 0) + sum(data[,ps] > 1) > 0) {
    stop("\n Error: Propensity score must be between 0 and 1.")
  } else if (sum(data[, status] != 0 & data[, status] != 1) > 0) {
    stop("\n Error: Status must be 0 or 1.")
  } else {
    if (is.null(ps)) {
      ## Crude Kaplan-Meier curve
      weight <- rep(1, length(data[, time]))
    } else {
      ps <- data[, ps] # P(A=1|L)
      pn <- as.numeric(glm(formula(paste0(eval(exposure),"~1")),
                           family=binomial(), data=data)$fitted.value[1]) # P(A=1)
      if(stabilize == TRUE){
        weight <- ifelse(data[, exposure] == 0, (1-pn)/(1-ps), pn/ps) # sipw
      } else {
        weight <- ifelse(data[, exposure] == 0, 1/(1-ps), 1/ps) # ipw
      }
    }

    data <- data.frame(cbind(data[, time], data[, status], data[, exposure], weight))
    colnames(data) <- c("time", "status", "exposure", "weight")
    data <- data[!is.na(data$exposure) & !is.na(data$time), ]
    data <- data[order(data$exposure), ]
    data$exposure <- factor(data$exposure)

    ## If tau not specified, use the minimum of the largest observed event time in both groups
    j <- length(unique(data$exposure))
    if (is.null(tau)) {
      taui <- c()
      for (i in (1:j)) {
        exposureval <- (levels(data$exposure)[i])
        dat_exposure <- data[which(data$exposure == (exposureval)), ]
        taui[i] <- max(dat_exposure$time[dat_exposure$status == 1])
      }
      tau <- min(taui)
    }

    ## Calculate restricted mean survival time via adjusted Kaplan-Meier curve
    rmst <- c()
    exposureval <- c()
    rmst_var <- c()
    rmst_se <- c()

    for (i in 1:j) {
      exposureval[i] <- (levels(data$exposure)[i])
      dat_exposure <- data[which(data$exposure == (exposureval[i])), ]

      ## Adjusted Kaplan-Meier curve Based on 'ipw.survival' function in "RISCA"
      ## Author: Le Borgne, F., et al.
      tj <- c(0, sort(unique(dat_exposure$time[dat_exposure$status == 1])))
      dj <- sapply(tj, function(x) {
        sum(dat_exposure$weight[dat_exposure$time == x & dat_exposure$status == 1])
      })
      yj <- sapply(tj, function(x) {
        sum(dat_exposure$weight[dat_exposure$time >= x])
      })
      st <- cumprod(1 - (dj / yj))
      ft <- data.frame(tj, yj, dj, st, i)

      ## Restricted mean survival time
      rtime <- ft$tj <= tau
      tj_r <- sort(c(ft$tj[rtime], tau))
      st_r <- ft$st[rtime]
      yj_r <- ft$yj[rtime]
      dj_r <- ft$dj[rtime]
      time_diff <- diff(c(0, tj_r))
      areas <- time_diff * c(1, st_r)
      rmst[i] <- sum(areas)

      ## Variance of RMST
      if (var.est == TRUE) {
        m <- sapply(tj, function(x) {
          sum((dat_exposure$weight[dat_exposure$time >= x])^2)
        })
        mj <- ((yj^2)/m)
        mj_r <- mj[rtime]
        var_r <- ifelse((yj_r - dj_r) == 0, 0, dj_r / (mj_r * (yj_r - dj_r)))
        var_r <- c(var_r, 0)
        rmst_var[i] <- sum(cumsum(rev(areas[-1]))^2 * rev(var_r)[-1])
        rmst_se[i] <- sqrt(rmst_var[i])
      }
    }

    if (var.est == TRUE) {
      ## Variance of adjusted Kaplan-Meier curve Based on 'rmst2 function' in "survRM2"
      ## Author: Uno, H., et al.
      results <- data.frame(exposureval, rmst, rmst_se, rmst_var)
      label_diff <- paste0("exposures ", results[2, ]$exposureval,
                           " vs. ", results[1, ]$exposureval)
      rmst_diff <- (results[2, ]$rmst - results[1, ]$rmst)
      rmst_diff_se <- sqrt(results[2, ]$rmst_var + results[1, ]$rmst_var)
      rmst_diff_low <- rmst_diff - qnorm(1 - alpha/2) * rmst_diff_se
      rmst_diff_upp <- rmst_diff + qnorm(1 - alpha/2) * rmst_diff_se
      rmst_diff_pval <- 2 * (1 - pnorm(abs(rmst_diff)/rmst_diff_se))

      final <- c(tau, results[2, ], results[1, ])
      names(final) <- c("tau", "exposed", "rmst1", "rmst_se1", "rmst_var1",
                        "unexposed", "rmst0", "rmst_se0", "rmst_var0")
      final$label_diff <- label_diff
      final$rmst_diff <- rmst_diff
      final$rmst_diff_se <- rmst_diff_se
      final$rmst_diff_low <- rmst_diff_low
      final$rmst_diff_upp <- rmst_diff_upp
      final$rmst_diff_pval <- rmst_diff_pval
    } else {
      results <- data.frame(exposureval, rmst)
      final <- c(tau, results[2, ], results[1, ])
      names(final) <- c("tau", "exposed", "rmst1", "unexposed", "rmst0")
    }
  }
  return(final)
}
