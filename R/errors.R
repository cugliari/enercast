#' performance
#'
#' Performance measure for a predicted set of values.  
#'
#' @param y vector of observed load
#' @param yhat vector of predicted load
#' @param measure 
#' \describe{
#'   \item{rmse} {Root Mean Square Error: \deqn{rmse(y,\hat{y})=\sqrt{\frac{1}{N}\sum\limits_{t=1}^N (y_t - \hat{y}_t)^2}}}
#'   \item{mape} {Mean Absolute Percentage Error: \deqn{mape(y,\hat{y})=\frac{1}{N}\sum\limits_{t=1}^N \frac{|y_t - \hat{y}_t|}{|y_t|}}}
#'   \item{mae} {Mean Absolute Error:\deqn{mae(y,\hat{y}) =\frac{1}{N}\sum\limits_{t=1}^N |y_t - \hat{y}_t|}}
#'   \item{pdad} {Promedio de Desviacion Absoluta de la Demanda: \deqn{pdad(y,\hat{y})=\frac{1}{N}\sum\limits_{t=1}^N \frac{|y_t - \hat{y}_t|}{\hat{y}}*100}}
#' }
#' 
#' @return single value
#' @export
#'
#' @examples performance(rep(0, 10), rnorm(10))
#' @export
#' @author Jairo Cugliari, Andres Castrillejo, Fernando Massa, Ignacio Ramirez
#' @references{  Hyndman, R.J. and Athanasopoulos, G. 2014 \emph{Forecasting: principles and practice}. OTexts.}
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
