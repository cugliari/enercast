#' Title
#'
#' @param formulaXa 
#' @param formulaXf 
#' @param data 
#' @param ini 
#' @param for.each 
#' @param cores 
#'
#' @return
#' @export
#'
#' @examples
ssm <- function(formulaXa, formulaXf = NULL, data, ini = NULL,
                for.each = 'hour', cores = parallel::detectCores() - 1){
  cls <- parallel::makeCluster(cores)
  parallel::registerDoParallel(cls)
  on.exit(parallel::stopCluster(cls)) ##
  h <- unique(data[,for.each])
  k <- length(h)
  res <- foreach::foreach(j = 1:k, .final = function(x) setNames(x,h),
                          .packages = c('FKF', 'enercast')) %dopar% {
    data.j <- data[data[, for.each] == h[j], names(data) != for.each]
    ini.j <- NULL
    if (!is.null(ini)) ini.j <- ini[[j]]
    ssm.mod(formulaXa, formulaXf, data = data.j, ini = ini.j)
    }

    res$formulas <- list(Xa = formulaXa, Xf = formulaXf)
  return(res)
}
