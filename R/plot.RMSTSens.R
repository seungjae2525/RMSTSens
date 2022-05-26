#' @title Plot for sensitivity analysis
#'
#' @description Plot for sensitivity analysis either or both of range and confidence interval for bias-adjusted RMST
#'
#' @param x an object of "autoplot.RMSTSens"
#' @param ... further arguments passed to or from other methods.
#'
#' @importFrom graphics plot
#'
#' @rdname plot.RMSTSens
#'
#' @export
plot.RMSTSens <- function(x, ...) {
  print(autoplot(x, ...))
}



#' @title Plot for sensitivity analysis
#'
#' @param x an object of class "RMSTSens"
#' @param smooth.degree Degree of smooth for geom_ribbon, Default: 11
#' @param alpha.ci It refers to the opacity of confidence interval. Values of alpha range from 0 to 1, with lower values corresponding to more transparent colors, Default: 0.9
#' @param alpha.range It refers to the opacity of range.Values of alpha range from 0 to 1, with lower values corresponding to more transparent colors, Default: 0.4
#' @param yscale 1, 10, 100, 1000, Default: 100
#' @param ytickdiff Distance between y-axis tick, Default: 100
#' @param point.size Size of estimate of lower and upper for bias-adjusted RMST
#' @param h.width Horizon lines width. By default, set to 1
#' @param axis.title.size Size of x and y axis title
#' @param axis.text.size Size of x and y axis text
#' @param save.plot When TRUE, it will save image, Default: FALSE
#' @param save.plot.name File name to create on disk, Default: 'Plot'
#' @param save.plot.device Device to use. Can either be a device function (e.g. png), or one of "eps", "ps", "tex" (pictex), "pdf", "jpeg", "tiff", "png", "bmp", "svg" or "wmf" (windows only).
#' @param save.plot.width Plot width size in units
#' @param save.plot.height Plot height size in units
#' @param save.plot.dpi Resolution of plot. By default, set to 300. Also accepts a string input: "retina" (320), "print" (300), or "screen" (72). Applies only to raster output types.
#'
#' @return Results for sensitivity analysis plot.
#'
#' @details If the object contains results of bootstrap confidence interval, then it will plot of range and confidence interval for bias-adjusted RMST, otherwise it will only plot of range for bias-adjusted RMST.
#'
#' @examples
#' if(interactive()){
#'  dat <- gbsg
#'  dat$size2 <- ifelse(dat$size <= 20, 0,
#'                      ifelse(dat$size > 20 & dat$size <= 50, 1, 2))
#'  dat$age2 <- dat$age/100
#'  dat$er2 <- dat$er/1000
#'
#'  ## Estimation of propensity score
#'  denom.fit <- glm(hormon~(age2)^3+(age2)^3*log(age2)+meno+factor(size2)+sqrt(nodes)+er2,
#'                   data=dat, family=binomial(link='logit'))
#'  dat$Ps <- predict(denom.fit, type='response')
#'
#'  ## Between-group difference in adjusted RMST based on shifted propensity score
#'  ## Adjusted RMST with not specified tau and with multiple lambda
#'  # Using approximate optimization method
#'  results.approx2 <- RMSTsensitivity(time='rfstime', status='status', exposure='hormon',
#'                                     exposed.ref.level=1, ps='Ps' ,data=dat, methods='Approx',
#'                                     use.multicore=TRUE, n.core=2,
#'                                     lambda=c(1,1.5), tau=365.25*5, ini.par=1, verbose=FALSE)
#'  plot(x=results.approx2, smooth.degree=11, alpha.ci=0.9, alpha.range=0.4,
#'       yscale=100, ytickdiff=100, point.size=1.4, h.width=1,
#'       axis.title.size=15, axis.text.size=12,
#'       save.plot=FALSE, save.plot.name="Plot", save.plot.device="png",
#'       save.plot.width=10, save.plot.height=6, save.plot.dpi=300)
#'
#'  re.ap.boot <- boot.ci.RMST(x=results.approx2, B=20, level=0.95, seed=220524,
#'                formula=hormon~(age2)^3+(age2)^3*log(age2)+meno+factor(size2)+sqrt(nodes)+er2,
#'                model="logistic", use.multicore=TRUE, n.core=2, verbose=TRUE)
#'  plot(x=re.ap.boot, smooth.degree=11, alpha.ci=0.9, alpha.range=0.4,
#'       yscale=100, ytickdiff=100, point.size=1.4, h.width=1,
#'       axis.title.size=15, axis.text.size=12,
#'       save.plot=FALSE, save.plot.name="Plot", save.plot.device="png",
#'       save.plot.width=10, save.plot.height=6, save.plot.dpi=300)
#' }
#'
#' @seealso
#'  \code{\link[splines]{ns}} \code{\link[ggplot2]{theme}}
#'
#' @importFrom ggplot2 autoplot
#'
#' @rdname plot.RMSTSens
#'
#' @export
autoplot.RMSTSens <- function(x, smooth.degree=11, alpha.ci=0.9, alpha.range=0.4,
                         yscale=100, ytickdiff=100, point.size=1.4, h.width=1,
                         axis.title.size=15, axis.text.size=12,
                         save.plot=FALSE, save.plot.name="Plot", save.plot.device="png",
                         save.plot.width=10, save.plot.height=6, save.plot.dpi=300) {
  if (!inherits(x, "RMSTSens")){
    stop("Argument 'x' must be an object of class \"RMSTSens\".")
  }

  xx <- x$result.df

  if (length(xx$Lambda) == 1) {
    stop("\n Error: To plot the results, \"lambda\" must be a vector.")
  }

  Lambda <- RMST.diff.min.lower <- RMST.diff.max.upper <- RMST.diff.min <- RMST.diff.max <- NULL
  if("RMST.diff.max.upper" %in% colnames(xx)){
    ylabel <- seq(floor(min(xx$RMST.diff.min.lower)/yscale)*yscale,
                  ceiling(max(xx$RMST.diff.max.upper)/yscale)*yscale, ytickdiff)

    g1 <- ggplot(xx) +
      stat_smooth(aes(x=Lambda, y=RMST.diff.min.lower, colour = "min"),
                  method = "glm", formula = y ~ splines::ns(x,2), n=smooth.degree) +
      stat_smooth(aes(x=Lambda, y=RMST.diff.max.upper, colour = "max"),
                  method = "glm", formula = y ~ splines::ns(x,2), n=smooth.degree)
    gg1 <- ggplot_build(g1)

    g2 <- ggplot(xx) +
      stat_smooth(aes(x=Lambda, y=RMST.diff.min, colour = "min"),
                  method = "glm", formula = y ~ splines::ns(x,2), n=smooth.degree) +
      stat_smooth(aes(x=Lambda, y=RMST.diff.max, colour = "max"),
                  method = "glm", formula = y ~ splines::ns(x,2), n=smooth.degree)
    gg2 <- ggplot_build(g2)

    df <- data.frame(Lambda = gg1$data[[1]]$x,
                     RMST.diff.min.lower = gg1$data[[1]]$y,
                     RMST.diff.max.upper = gg1$data[[2]]$y,
                     RMST.diff.min = gg2$data[[1]]$y,
                     RMST.diff.max = gg2$data[[2]]$y)
    df <- rbind(df, xx[,c("Lambda","RMST.diff.min.lower","RMST.diff.max.upper",
                          "RMST.diff.min","RMST.diff.max")])
    df <- df[order(df$Lambda),]
    df <- df[!duplicated(round(df, 10), fromLast = TRUE),]

    g <- ggplot(df) +
      geom_ribbon(aes(x=Lambda, ymin=RMST.diff.min.lower, ymax=RMST.diff.max.upper),
                  alpha=alpha.ci, inherit.aes=F, fill="pink") +
      geom_ribbon(aes(x=Lambda, ymin=RMST.diff.min, ymax=RMST.diff.max),
                  alpha=alpha.range, inherit.aes=F, fill="red") +
      geom_line(aes(x=Lambda, y=xx$RMST.diff.min[1]), colour="red", size=h.width) +
      geom_point(data=xx, aes(x=Lambda, y=RMST.diff.min), size=point.size, colour="blue") +
      geom_point(data=xx, aes(x=Lambda, y=RMST.diff.max), size=point.size, colour="blue") +
      geom_line(aes(x=Lambda, y=0), colour="black", size=h.width, linetype = "dashed", alpha=0.5) +
      xlab(expression(Lambda)) + ylab("Between-group difference in bias-adjusted RMST")  +
      theme_bw() +
      scale_x_continuous(breaks = xx$Lambda,
                         labels = sprintf("%.1f", xx$Lambda), expand = c(0.005,0.005)) +
      scale_y_continuous(breaks = ylabel, labels = ylabel) +
      theme(axis.title.x=element_text(face="bold", size=axis.title.size),
            axis.title.y=element_text(face="bold", size=axis.title.size, angle=90),
            axis.text=element_text(face="bold", size=axis.text.size),
            panel.grid.minor.x = element_blank(),
            panel.grid.minor.y = element_blank())

  } else {
    ylabel <- seq(floor(min(xx$RMST.diff.min)/yscale)*yscale,
                  ceiling(max(xx$RMST.diff.max)/yscale)*yscale, ytickdiff)

    g2 <- ggplot(xx) +
      stat_smooth(aes(x=Lambda, y=RMST.diff.min, colour = "min"),
                  method = "glm", formula = y ~ splines::ns(x,2), n=smooth.degree) +
      stat_smooth(aes(x=Lambda, y=RMST.diff.max, colour = "max"),
                  method = "glm", formula = y ~ splines::ns(x,2), n=smooth.degree)
    gg2 <- ggplot_build(g2)

    df <- data.frame(Lambda = gg2$data[[1]]$x,
                     RMST.diff.min = gg2$data[[1]]$y,
                     RMST.diff.max = gg2$data[[2]]$y)
    df <- rbind(df, xx[,c("Lambda","RMST.diff.min","RMST.diff.max")])
    df <- df[order(df$Lambda),]
    df <- df[!duplicated(round(df, 10), fromLast = TRUE),]

    g <- ggplot(df) +
      geom_ribbon(aes(x=Lambda, ymin=RMST.diff.min, ymax=RMST.diff.max),
                  alpha=alpha.range, inherit.aes=F, fill="red") +
      geom_line(aes(x=Lambda, y=xx$RMST.diff.min[1]), colour="red", size=h.width) +
      geom_point(data=xx, aes(x=Lambda, y=RMST.diff.min), size=point.size, colour="blue") +
      geom_point(data=xx, aes(x=Lambda, y=RMST.diff.max), size=point.size, colour="blue") +
      geom_line(aes(x=Lambda, y=0), colour="black", size=h.width, linetype = "dashed", alpha=0.5) +
      xlab(expression(Lambda)) + ylab("Between-group difference in bias-adjusted RMST")  +
      theme_bw() +
      scale_x_continuous(breaks = xx$Lambda,
                         labels = sprintf("%.1f", xx$Lambda), expand = c(0.005,0.005)) +
      scale_y_continuous(breaks = ylabel, labels = ylabel) +
      theme(axis.title.x=element_text(face="bold", size=axis.title.size),
            axis.title.y=element_text(face="bold", size=axis.title.size, angle=90),
            axis.text=element_text(face="bold", size=axis.text.size),
            panel.grid.minor.x = element_blank(),
            panel.grid.minor.y = element_blank())
  }

  ##
  if(save.plot == TRUE){
    file.name <- paste0(save.plot.name, ".", save.plot.device)
    ggsave(filename=save.plot.name, plot=g,
           width=save.plot.width, height=save.plot.height,
           device=save.plot.device, dpi=save.plot.dpi)
    # dev.off()
  }

  g
}

