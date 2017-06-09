#' rf
#'
#' Sets up random forest model for prediction
#'
#' @param y Character. Variable name.
#' @param data data.frame
#' @param for.each Character. Variable name.
#' @param cores Parallel execution
#'
#' @return
#' @export
#'
#' @author Jairo Cugliari, Andres Castrillejo, Fernando Massa, Ignacio Ramirez
#' @examples
# \usage{rf(y = "load", data, for.each = "hour", cores = detectCores() - 1)}

rf <- function(y = 'load', data, for.each = 'hour', 
               cores = parallel::detectCores() - 1){
  cls <- parallel::makeCluster(cores)
  parallel::registerDoParallel(cls)
  on.exit(parallel::stopCluster(cls))
  K <- unique(data[,for.each])
  RES <- foreach::foreach(j = K ,.packages = 'randomForest') %dopar% {
      data.j <- data[data[,for.each] == j, names(data) != for.each]
      XX <- data.j[ , names(data.j) != eval(y)]
      YY <- data.j[, eval(y)] 
      rf <- randomForest::randomForest(y = YY, x = XX, na.action = na.omit)}
  names(RES) <- K
  return(RES)
}

