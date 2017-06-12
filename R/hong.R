#' hong
#'
#' Fit Hong Vanilla Benchmark multiple regression model for electricity load
#' forecasting
#' 
#' @param data Data frame with columns:
#' \describe{
#'   \item{year}{year (...2012, 2013, 2014...)}
#'   \item{month}{month (01, 02, ..., 12)}
#'   \item{wday}{day of the week} (Mon Thu Wed Tue Fri Sat Sun)}  
#'   \item{instant}{(01, 02, ..., 24) for hourly data or 
#'                  (01, 02, ..., 48) or half hourly data}  
#'   \item{temp}{temperature} 
#'   \item{load}{electricity loads}  
#'}
#'
#' @return An object of class \code{lm}
#' @export
#' @details 
#' The Hong Vanilla Model:
#'  \deqn{E(Load) = \beta_0 + \beta_1 Trend +\beta_2 Day \times Hour  + 
#'                  \beta_3 Month + \beta_4 Month \times Temp + 
#'                  \beta_5 Month  \times Temp^2  +\beta_6 Month  \times Temp^3 + 
#'                  \beta_7 Hour \times Temp + \beta_8 Hour \times Temp^2 +
#'                  \beta_9  Hour \times Temp^3}
#'
#' @author Jairo Cugliari, Andres Castrillejo, Fernando Massa, Ignacio Ramirez
#' 
#' @references Hong, T. and Wang, P. and Lee Willis, H. (2011) A naive multiple 
#' linear regression benchmark for short term load forecasting. In \emph{Power 
#' and Energy Society General Meeting, 2011 IEEE}, 1--6.
#' 
#' Hong, T. and Pinson, P.and Fan, S. (2014) Global energy forecasting 
#' competition 2012. \emph{International Journal of Forecasting}, 
#' \bold{30(2)}:357--363
#'
#' @examples  1
#' @export
hong <- function(data){
  
  trend <- 1:nrow(data)
    
  if (!is.character(data$hour))   data$hour  <- as.character(data$hour)
  if (!is.character(data$month))  data$month <- as.character(data$month)
 
  hv <- lm(load ~ trend + wday * hour + month + month * temp + 
                  month * (temp^2) + month * (temp^3) +  hour * temp +  
                  hour  * (temp^2) +  hour * (temp^3),
           data = data)
    
  class(hv) <- c('lm', 'hv_class')
  return(hv)
}
