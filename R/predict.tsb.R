#' predict.tsb
#'
#' Predicted values for Time Series Benchmark Model
#' 
#' @param mod An ARIMA model object fitting from \code{tsb}
#' @param steps Number of steps to forecast
#' @param n_run Length of the univariate time series to parse to KalmanRun 
#'        (must be smaller than the total length of data in mod)
#' @return vector of predicted values
#' @author Jairo Cugliari, Andres Castrillejo, Fernando Massa, Ignacio Ramirez
#' @seealso \link{tsb} \link{makeARIMA}
#' @export
predict.tsb <- function(mod, steps = 24, n_run = 1000){

	dem <- mod$data
	n <- length(dem)

	if (n < n_run) {
	  n_run <- n - 1 
	  warning("Not enough data or n_run too large: n_run set as n - 1")
	}
	
	tita <- comb.coef(comb.coef(mod$ma$reg, mod$ma$est1, 24), mod$ma$est2, 168)
	fi <- comb.coef(mod$ar$reg, mod$ar$est1, 24)

	delta <- 1
	if (mod$d[1] > 0) 
	  for (i in 1:mod$d[1]) delta <- convolve(delta, c(-1, 1), type = 'open')
	
	if (mod$d[2] > 0) 
	  for (i in 1:mod$d[2]) delta <- convolve(delta, c(-1, rep(0, mod$frec[1] - 1), 1),
	                                         type = 'open')
	if (mod$d[3] > 0) 
	  for (i in 1:mod$d[3]) delta <- convolve(delta, c(-1, rep(0, mod$frec[2] - 1), 1),
	                                         type = 'open')
	delta[abs(delta) < 1e-10] <- 0
	delta <- delta[-1]
	
	if (is.null(fi))   fi   <- vector('numeric')
	if (is.null(tita)) tita <- vector('numeric')
        
	modelo <- makeARIMA(phi = fi, Delta = -delta, theta = tita, SSinit = "Rossignol2011")
#	res <- KalmanRun(dem[(n - n_run):n], mod = modelo)
#	modelo$a <- res$states[nrow(res$states), ]
	pred <- KalmanForecast(steps, modelo)
	return(pred$pred)
}

