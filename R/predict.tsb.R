#' Title
#'
#' @param mod 
#' @param steps 
#'
#' @return
#' @export 
#'
#' @examples 3
predict.tsb <- function(mod, steps = 24){

	dem <- mod$data
	n <- length(dem)

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
	res <- KalmanRun(dem[(n - 1000):n], mod = modelo)
	modelo$a <- res$states[nrow(res$states), ]
	pred <- KalmanForecast(steps, modelo)
	return(pred$pred)
}

