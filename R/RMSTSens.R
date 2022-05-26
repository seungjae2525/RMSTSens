RMSTSens <- function(x, ...) UseMethod("RMSTSens")

#' @import parallel
#' @import optimParallel
RMSTSens.default <- function(time, status,
                             exposure, exposed.ref.level=1, ps, data,
                             methods="Approx", use.multicore=TRUE, n.core=parallel::detectCores()/2,
                             lambda=2, tau=NULL, ini.par=1, verbose=FALSE){
  est <- RMSTsensitivity(time=time, status=status,
                         exposure=exposure, exposed.ref.level=exposed.ref.level, ps=ps, data=data,
                         methods=methods, use.multicore=use.multicore, n.core=n.core,
                         lambda=lambda, tau=tau, ini.par=ini.par, verbose=FALSE)
  class(est) <- "RMSTSens"
  est
}
