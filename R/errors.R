#' performance
#'
#' Performance measure for a predicted set of values.  
#'
#' @param y 
#' @param yhat 
#' @param measure 
#'
#' @return single value
#' @export
#'
#' @examples performance(rep(0, 10), rnorm(10))
#' @export
performance <- function(y, yhat, measure = "rmse") {
  switch(measure,
         mape = mean(abs(y - yhat) / abs(y)),
         mae  = mean(abs(y - yhat)),
         rmse = sqrt(mean((y - yhat)^2)),
         pdad = mean(abs(y - yhat)/yhat) * 100)
}
#mape <- function(y, yhat){mean(abs(y - yhat)/abs(y))}
#mae <- function(y, yhat){mean(abs(y - yhat))}
#rmse <- function(y, yhat){sqrt(mean((y - yhat)^2))}
#pdad <- function(y, yhat){mean(abs(y - yhat)/yhat) * 100
