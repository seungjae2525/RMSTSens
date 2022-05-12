## Make dataset for optimization
optim_data <- function(time, status, exposure, ps, data, group=1) {
  time <- data[, time]                    # survival time
  status <- data[, status]                # status
  exposure <- data[, exposure]            # exposure
  ps <- data[, ps]                        # propensity score: P(A=1|L)
  .Pr <- ifelse(exposure == 0, 1-ps, ps)  # P(A=a|L)
  exp_logit_Pr <- exp(-qlogis(.Pr))       # exp(-logit(P(A=a|L)))

  time <- time[exposure == group]
  status <- status[exposure == group]
  exp_logit_Pr <- exp_logit_Pr[exposure == group]

  data <- data.frame(time=time, status=status, exp_logit_Pr=exp_logit_Pr)
  data <- data[order(-data$status, data$time), ]
  rownames(data) <- NULL

  return(data)
}
