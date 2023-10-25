#' @title Plot for sensitivity analysis results
#'
#' @param object An object for class \code{RMSTSens}. Either object \code{RMSTSens} or \code{RMSTSens.ci} is allowed.
#' If you want to see the results of several \code{RMSTSens} objects together, use the \code{RMSTSens.merge} function. See \code{RMSTSens.merge}.
#' @param alpha.range A opacity of sensitivity range. Values of \code{alpha.range} range from 0 to 1, with lower values corresponding to more transparent colors.
#' Default: 0.4.
#' @param alpha.ci A opacity of confidence interval. Values of \code{alpha.ci} range from 0 to 1, with lower values corresponding to more transparent colors.
#' Default: 0.9.
#' @param smooth.par Smooth parameter. It corresponds to argument \code{n} in \code{spline} function in \code{stats} package.
#' Default: 100.
#' @param ytickdiff Distance between y-axis tick. Default: 100.
#' @param point.size Point size of the the lower and upper bound of the sensitivity range. Default: 1.4.
#' @param h.width Width of horizon line which represents the estimated difference in adjusted RMST. Default: 1.
#' @param axis.title.size Size of x and y axis title. Default: 15.
#' @param axis.text.size Size of x and y axis text. Default: 12.
#' @param save.plot When TRUE, it will save the resulting plot. To save the resulting plot to a specific directory,
#' the working directory will be probably need to set using \code{setwd()}. Default: FALSE.
#' @param save.plot.name File name to create on disk. Default: "Plot".
#' @param save.plot.device Device to use. Can either of "png", "eps", "ps", "tex" (pictex), "pdf", "jpeg", "tiff", bmp", "svg", or "wmf" (windows only). Default: "png".
#' @param save.plot.width Plot width size in units. Default: 10.
#' @param save.plot.height Plot height size in units. Default: 6.
#' @param save.plot.dpi Resolution of plot. Also accepts a string input: "retina" (320), "print" (300), or "screen" (72). Applies only to raster output types.
#' Default: 300.
#' @param ... Further arguments (currently not used).
#'
#' @details If the object contains the results of confidence interval for the population sensitivity range as well as the sensitivity range,
#' then it will return a plot with sensitivity range with confidence interval, otherwise it will only return a plot of sensitivity range.
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
#'                             lambda=c(1,1.1,1.2), tau=365.25*5, ini.par=1, verbose=FALSE)
#' autoplot(object=results.approx2, alpha.range=0.4, alpha.ci=0.9, smooth.par=100,
#'          ytickdiff=100, point.size=1.4, h.width=1,
#'          axis.title.size=15, axis.text.size=12,
#'          save.plot=FALSE, save.plot.name="Plot", save.plot.device="png",
#'          save.plot.width=10, save.plot.height=6, save.plot.dpi=300)
#' # Additional sensitivity analysis when lambda=1.15
#' results.approx3 <- RMSTSens(time="rfstime", status="status", exposure="hormon",
#'                             level.exposed="1", ps="Ps", data=dat, methods="Approx",
#'                             use.multicore=TRUE, n.core=2,
#'                             lambda=c(1.15), tau=365.25*5, ini.par=1, verbose=FALSE)
#' # After merging two results, plot the analysis results.
#' autoplot(object=RMSTSens.merge(list(results.approx2, results.approx3)), smooth.par=100,
#'          alpha.range=0.4, alpha.ci=0.9,
#'          ytickdiff=100, point.size=1.4, h.width=1,
#'          axis.title.size=15, axis.text.size=12,
#'          save.plot=FALSE, save.plot.name="Plot", save.plot.device="png",
#'          save.plot.width=10, save.plot.height=6, save.plot.dpi=300)
#'
#' # Percentile bootstrap CI for population sensitivity range
#' re.ap.boot <- RMSTSens.ci(x=RMSTSens.merge(list(results.approx2, results.approx3)),
#'               B=50, # Set B=50 to reduce computation time for R checks
#'               level=0.95, seed=220524,
#'               formula=hormon~age+meno+size+factor(grade)+nodes+pgr+er,
#'               model="logistic", use.multicore=TRUE, n.core=2)
#' autoplot(object=re.ap.boot, alpha.range=0.4, alpha.ci=0.9, smooth.par=100,
#'          ytickdiff=100, point.size=1.4, h.width=1,
#'          axis.title.size=15, axis.text.size=12,
#'          save.plot=FALSE, save.plot.name="Plot", save.plot.device="png",
#'          save.plot.width=10, save.plot.height=6, save.plot.dpi=300)
#'
#' @seealso
#'  \code{\link[RMSTSens]{RMSTSens}}, \code{\link[RMSTSens]{RMSTSens.ci}}, \code{\link[RMSTSens]{RMSTSens.merge}}
#'
#' @import ggplot2
#'
#' @rdname autoplot.RMSTSens
#'
#' @export
autoplot.RMSTSens <- function(object=object,
                              alpha.range=0.4, alpha.ci=0.9, smooth.par=100,
                              ytickdiff=100, point.size=1.4, h.width=1,
                              axis.title.size=15, axis.text.size=12,
                              save.plot=FALSE, save.plot.name="Plot", save.plot.device="png",
                              save.plot.width=10, save.plot.height=6, save.plot.dpi=300, ...){
  autoplot_RMSTSens(xxx=object,
                    alpha.range=alpha.range, alpha.ci=alpha.ci, smooth.par=smooth.par,
                    ytickdiff=ytickdiff, point.size=point.size, h.width=h.width,
                    axis.title.size=axis.title.size, axis.text.size=axis.text.size,
                    save.plot=save.plot, save.plot.name=save.plot.name,
                    save.plot.device=save.plot.device, save.plot.width=save.plot.width,
                    save.plot.height=save.plot.height, save.plot.dpi=save.plot.dpi, ...=...)
}




#' @title Plot for sensitivity analysis results
#'
#' @description Plot for sensitivity analysis results of either or both the sensitivity range of the difference in adjusted RMST and
#' the confidence interval for the population sensitivity range.
#'
#' @param x An object for class \code{RMSTSens}.
#'
#' @examples
#' \dontrun{
#' plot(results.approx2, alpha.range=0.4, alpha.ci=0.9, smooth.par=100,
#'      ytickdiff=100, point.size=1.4, h.width=1,
#'      axis.title.size=15, axis.text.size=12,
#'      save.plot=FALSE, save.plot.name="Plot", save.plot.device="png",
#'      save.plot.width=10, save.plot.height=6, save.plot.dpi=300)
#' }
#'
#' @keywords plot
#'
#' @importFrom graphics plot
#'
#' @rdname autoplot.RMSTSens
#'
#' @export
plot.RMSTSens <- function(x, ...) {
  print(autoplot(x, ...))
}



#########---------------------------------------------------------------------------------------------
autoplot_RMSTSens <- function(xxx=NULL,
                              alpha.range=NULL, alpha.ci=NULL, smooth.par=NULL,
                              ytickdiff=NULL, point.size=NULL, h.width=NULL,
                              axis.title.size=NULL, axis.text.size=NULL,
                              save.plot=FALSE, save.plot.name=NULL, save.plot.device=NULL,
                              save.plot.width=NULL, save.plot.height=NULL, save.plot.dpi=NULL, ...) {
  ## Function for counting the number of decimal places
  decimalplaces <- function(x) {
    rr <- c()
    for(i in 1:length(x)){
      if(i == 1){
        rr[i] <- "1.0"
      } else {
        if (abs(x[i] - round(x[i])) > .Machine$double.eps^0.5) {
          rrr <- nchar(strsplit(sub('0+$', '', as.character(x[i])), ".", fixed = TRUE)[[1]][[2]])
        } else {
          rrr <- 0
        }
        rr[i] <- sprintf(paste0("%.",rrr,"f"), x[i])
      }
    }
    return(rr)
  }

  ##
  xx <- xxx$result.df

  ## Make the line smooth
  # smooth.par <- 2*length(xx$Lambda)-1

  ##
  if (length(xx$Lambda) == 1) {
    stop("\n Error: To plot the sensitivity analysis results, \"lambda\" must be a vector.")
  }

  ## Plot
  Lambda <- RMST.diff.min.lower <- RMST.diff.max.upper <- RMST.diff.min <- RMST.diff.max <- NULL
  yscale <- 10^(floor(log10(ytickdiff)))
  if("RMST.diff.max.upper" %in% colnames(xx)){
    ylabel <- seq(floor(min(xx$RMST.diff.min.lower)/yscale)*yscale,
                  ceiling(max(xx$RMST.diff.max.upper)/yscale)*yscale, ytickdiff)

    if(length(xx$Lambda) == 2){
      x0 <- seq(xx$Lambda[1], xx$Lambda[2], by = 0.01)

      aa1 <- (xx$RMST.diff.min.lower[2] - xx$RMST.diff.min.lower[1])/(xx$Lambda[2] - xx$Lambda[1])^2
      f1 <- function(x) -aa1*(x-xx$Lambda[2])^2 + xx$RMST.diff.min.lower[2]

      aa2 <- (xx$RMST.diff.min[2] - xx$RMST.diff.min[1])/(xx$Lambda[2] - xx$Lambda[1])^2
      f2 <- function(x) -aa2*(x-xx$Lambda[2])^2 + xx$RMST.diff.min[2]

      aa3 <- (xx$RMST.diff.max.upper[2] - xx$RMST.diff.max.upper[1])/(xx$Lambda[2] - xx$Lambda[1])^2
      f3 <- function(x) -aa3*(x-xx$Lambda[2])^2 + xx$RMST.diff.max.upper[2]

      aa4 <- (xx$RMST.diff.max[2] - xx$RMST.diff.max[1])/(xx$Lambda[2] - xx$Lambda[1])^2
      f4 <- function(x) -aa4*(x-xx$Lambda[2])^2 + xx$RMST.diff.max[2]

      df <- data.frame(Lambda = x0,
                       RMST.diff.min.lower = f1(x = x0),
                       RMST.diff.max.upper = f3(x = x0),
                       RMST.diff.min = f2(x = x0),
                       RMST.diff.max = f4(x = x0))
      df <- rbind(df, xx[,c("Lambda","RMST.diff.min.lower","RMST.diff.max.upper",
                            "RMST.diff.min","RMST.diff.max")])

      g <- ggplot(df) +
        geom_ribbon(aes(x=Lambda, ymin=RMST.diff.min.lower, ymax=RMST.diff.max.upper),
                    alpha=alpha.ci, inherit.aes=F, fill="pink") +
        geom_ribbon(aes(x=Lambda, ymin=RMST.diff.min, ymax=RMST.diff.max),
                    alpha=alpha.range, inherit.aes=F, fill="red") +
        geom_line(aes(x=Lambda, y=xx$RMST.diff.min[1]), colour="red", linewidth=h.width) +
        geom_point(data=xx, aes(x=Lambda, y=RMST.diff.min), size=point.size, colour="blue") +
        geom_point(data=xx, aes(x=Lambda, y=RMST.diff.max), size=point.size, colour="blue") +
        geom_line(aes(x=Lambda, y=0), colour="black", linewidth=h.width, linetype = "dashed", alpha=0.5) +
        xlab(expression(bold(Lambda))) + ylab("Difference in RMST")  +
        theme_bw() +
        scale_x_continuous(breaks = xx$Lambda,
                           labels = decimalplaces(xx$Lambda), expand = c(0.005,0.005)) +
        scale_y_continuous(breaks = ylabel, labels = ylabel) +
        theme(axis.title.x=element_text(face="bold", size=axis.title.size),
              axis.title.y=element_text(face="bold", size=axis.title.size, angle=90),
              axis.text=element_text(face="bold", size=axis.text.size),
              panel.grid.minor.x = element_blank(),
              panel.grid.minor.y = element_blank())

    } else {
      s1 <- spline(xx$Lambda, xx$RMST.diff.min.lower, n=smooth.par)
      s2 <- spline(xx$Lambda, xx$RMST.diff.max.upper, n=smooth.par)
      s3 <- spline(xx$Lambda, xx$RMST.diff.min, n=smooth.par)
      s4 <- spline(xx$Lambda, xx$RMST.diff.max, n=smooth.par)

      df <- data.frame(Lambda = s1$x,
                       RMST.diff.min.lower = s1$y,
                       RMST.diff.max.upper = s2$y,
                       RMST.diff.min = s3$y,
                       RMST.diff.max = s4$y)

      df <- rbind(df, xx[,c("Lambda","RMST.diff.min.lower","RMST.diff.max.upper",
                            "RMST.diff.min","RMST.diff.max")])
      df <- df[order(df$Lambda),]
      df <- df[!duplicated(round(df, 10), fromLast = TRUE),]

      g <- ggplot(df) +
        geom_ribbon(aes(x=Lambda, ymin=RMST.diff.min.lower, ymax=RMST.diff.max.upper),
                    alpha=alpha.ci, inherit.aes=F, fill="pink") +
        geom_ribbon(aes(x=Lambda, ymin=RMST.diff.min, ymax=RMST.diff.max),
                    alpha=alpha.range, inherit.aes=F, fill="red") +
        geom_line(aes(x=Lambda, y=xx$RMST.diff.min[1]), colour="red", linewidth=h.width) +
        geom_point(data=xx, aes(x=Lambda, y=RMST.diff.min), size=point.size, colour="blue") +
        geom_point(data=xx, aes(x=Lambda, y=RMST.diff.max), size=point.size, colour="blue") +
        geom_line(aes(x=Lambda, y=0), colour="black", linewidth=h.width, linetype = "dashed", alpha=0.5) +
        xlab(expression(bold(Lambda))) + ylab("Difference in RMST")  +
        theme_bw() +
        scale_x_continuous(breaks = xx$Lambda,
                           labels = decimalplaces(xx$Lambda), expand = c(0.005,0.005)) +
        scale_y_continuous(breaks = ylabel, labels = ylabel) +
        theme(axis.title.x=element_text(face="bold", size=axis.title.size),
              axis.title.y=element_text(face="bold", size=axis.title.size, angle=90),
              axis.text=element_text(face="bold", size=axis.text.size),
              panel.grid.minor.x = element_blank(),
              panel.grid.minor.y = element_blank())

    }
  } else {
    ylabel <- seq(floor(min(xx$RMST.diff.min)/yscale)*yscale,
                  ceiling(max(xx$RMST.diff.max)/yscale)*yscale, ytickdiff)

    if(length(xx$Lambda) == 2){
      x0 <- seq(xx$Lambda[1], xx$Lambda[2], by = 0.01)

      aa2 <- (xx$RMST.diff.min[2] - xx$RMST.diff.min[1])/(xx$Lambda[2] - xx$Lambda[1])^2
      f2 <- function(x) -aa2*(x-xx$Lambda[2])^2 + xx$RMST.diff.min[2]

      aa4 <- (xx$RMST.diff.max[2] - xx$RMST.diff.max[1])/(xx$Lambda[2] - xx$Lambda[1])^2
      f4 <- function(x) -aa4*(x-xx$Lambda[2])^2 + xx$RMST.diff.max[2]

      df <- data.frame(Lambda = x0,
                       RMST.diff.min = f2(x = x0),
                       RMST.diff.max = f4(x = x0))
      df <- rbind(df, xx[,c("Lambda","RMST.diff.min","RMST.diff.max")])

      g <- ggplot(df) +
        geom_ribbon(aes(x=Lambda, ymin=RMST.diff.min, ymax=RMST.diff.max),
                    alpha=alpha.range, inherit.aes=F, fill="red") +
        geom_line(aes(x=Lambda, y=xx$RMST.diff.min[1]), colour="red", linewidth=h.width) +
        geom_point(data=xx, aes(x=Lambda, y=RMST.diff.min), size=point.size, colour="blue") +
        geom_point(data=xx, aes(x=Lambda, y=RMST.diff.max), size=point.size, colour="blue") +
        geom_line(aes(x=Lambda, y=0), colour="black", linewidth=h.width, linetype = "dashed", alpha=0.5) +
        xlab(expression(bold(Lambda))) + ylab("Difference in RMST")  +
        theme_bw() +
        scale_x_continuous(breaks = xx$Lambda,
                           labels = decimalplaces(xx$Lambda), expand = c(0.005,0.005)) +
        scale_y_continuous(breaks = ylabel, labels = ylabel) +
        theme(axis.title.x=element_text(face="bold", size=axis.title.size),
              axis.title.y=element_text(face="bold", size=axis.title.size, angle=90),
              axis.text=element_text(face="bold", size=axis.text.size),
              panel.grid.minor.x = element_blank(),
              panel.grid.minor.y = element_blank())
    } else {
      s3 <- spline(xx$Lambda, xx$RMST.diff.min, n=smooth.par)
      s4 <- spline(xx$Lambda, xx$RMST.diff.max, n=smooth.par)

      df <- data.frame(Lambda = s3$x,
                       RMST.diff.min = s3$y,
                       RMST.diff.max = s4$y)

      df <- rbind(df, xx[,c("Lambda","RMST.diff.min","RMST.diff.max")])
      df <- df[order(df$Lambda),]
      df <- df[!duplicated(round(df, 10), fromLast = TRUE),]

      g <- ggplot(df) +
        geom_ribbon(aes(x=Lambda, ymin=RMST.diff.min, ymax=RMST.diff.max),
                    alpha=alpha.range, inherit.aes=F, fill="red") +
        geom_line(aes(x=Lambda, y=xx$RMST.diff.min[1]), colour="red", linewidth=h.width) +
        geom_point(data=xx, aes(x=Lambda, y=RMST.diff.min), size=point.size, colour="blue") +
        geom_point(data=xx, aes(x=Lambda, y=RMST.diff.max), size=point.size, colour="blue") +
        geom_line(aes(x=Lambda, y=0), colour="black", linewidth=h.width, linetype = "dashed", alpha=0.5) +
        xlab(expression(bold(Lambda))) + ylab("Difference in RMST")  +
        theme_bw() +
        scale_x_continuous(breaks = xx$Lambda,
                           labels = decimalplaces(xx$Lambda), expand = c(0.005,0.005)) +
        scale_y_continuous(breaks = ylabel, labels = ylabel) +
        theme(axis.title.x=element_text(face="bold", size=axis.title.size),
              axis.title.y=element_text(face="bold", size=axis.title.size, angle=90),
              axis.text=element_text(face="bold", size=axis.text.size),
              panel.grid.minor.x = element_blank(),
              panel.grid.minor.y = element_blank())
    }
  }

  ##
  if(save.plot == TRUE){
    file.name <- paste0(save.plot.name, ".", save.plot.device)
    ggsave(filename=file.name, plot=g,
           width=save.plot.width, height=save.plot.height,
           device=save.plot.device, dpi=save.plot.dpi)
    # dev.off()
  }

  return(g)
}



