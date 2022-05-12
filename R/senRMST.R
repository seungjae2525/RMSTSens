RMST_sen <- function(x, ...) UseMethod("senRMST")

#' @import parallel
#' @import optimParallel
RMST_sen.default <- function(time, status, exposure, exposed.ref.level=1, ps, data,
                             methods="Approx", use.multicore=TRUE, n.core=parallel::detectCores()/2,
                             lambda=2, tau=NULL, ini.par=1, verbose=FALSE){
  est <- RMST_sensitivity(time=time, status=status, ps=ps,
                          exposure=exposure, exposed.ref.level=exposed.ref.level,
                          data, methods=methods,
                          use.multicore=use.multicore, n.core=n.core,
                          lambda=lambda, tau=tau, ini.par=ini.par, verbose=FALSE)
  class(est) <- "senRMST"
  est
}
