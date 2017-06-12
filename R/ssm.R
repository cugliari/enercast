#' ssm
#'
#' Fit an univariate State Space Model
#'
#' @param formulaXa Formula with random regressors (\code{y~x1+x2+x3+x2:x3})
#' @param formulaXf Formula with fixed regressors *without* response variable (\code{~x4+x5})
#' @param data data.frame
#' @param ini List of intial value parameters of length of 'for.each', 24 in hourly data
#' @param for.each Variable of data.frame that indicate frequency
#' @param cores For parallel
#' 
#' @return   fkf.optim      fk list(dt=dt,ct=ct,Tt=Tt,Zt=Zt,a0=a0,P0=P0)
#' @export
#' 
#' @references #' Dordonnat, V. and Koopman, S. J. and Ooms, M. and Dessertaine, 
#' A. and Collet, J. (2008)  An hourly periodic state space model for modelling 
#' French national electricity load. \emph{International Journal of  Forecasting}, 
#' \bold{24(4)}:566--587.
#' 
#' Dordonnat, V. and Koopman, S. J. and Ooms, M. (2012) 
#' Dynamic factors in periodic time-varying regressions with an application to 
#' hourly electricity load modelling. 
#' \emph{Computational Statistics & Data Analysis}, \bold{56(11)}:3134--3152.}
#' 
#' @author Jairo Cugliari, Andres Castrillejo, Fernando Massa, Ignacio Ramirez
#' 
#' @examples
#' @export
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
