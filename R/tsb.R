#' tsb
#'
#' Time series benchmark
#'
#' @param yt load data sequence. type numeric.
#' @param p AR order. type numeric.
#' @param d Degree of differencing. type numeric.
#' @param q MA order. type numeric.
#' @param P1 AR order. type numeric.
#' @param D1 Degree of differencing. type numeric.
#' @param Q1 MA order. type numeric.
#' @param P2 AR order. type numeric.
#' @param D2 Degree of differencing. type numeric.
#' @param Q2 MA order. type numeric.
#' @param s1 Seasonal frequencies. type numeric.
#' @param s2 Seasonal frequencies. type numeric.
#' @param k number of lags
#'
#' @return A list with ARIMA estimated model
#' \describe{
#'   \item{ar}{(reg=ar,est1=AR1,est2=AR2)}
#'   \item{ma}{(reg=ma,est1=MA1,est2=MA2)}
#'   \item{d}{(d,D1,D2)}
#'   \item{frec}{}
#'   \item{dif}{}
#'   \item{iter}{}
#'   \item{data}{load data}
#' }
#' @author Jairo Cugliari, Andres Castrillejo, Fernando Massa, Ignacio Ramirez
#' @references{Brockwell, P.J.and Davis, R.A. 1991 \emph{Time Series: theory and methods} 2nd ed. Springer.}
#'
#' @seealso \code{\link{predict.tsb}}
#' @examples
tsb <- function(yt, p = 1, d = 1, q = 1, P1 = 1, D1 = 0, Q1 = 1, P2 = 0, D2 = 1,
                Q2 = 1, s1 = 24, s2 = 168, k = 200){
	if (p  == 0)  ar  = NULL else ar  = rep(0.1, p)
	if (P1 == 0)  AR1 = NULL else AR1 = rep(0.1, P1)
	if (P2 == 0)  AR2 = NULL else AR2 = rep(0.1, P2)
	if (q  == 0)  ma  = NULL else ma  = rep(0.1, q)
	if (Q1 == 0)  MA1 = NULL else MA1 = rep(0.1, Q1)
	if (Q2 == 0)  MA2 = NULL else MA2 = rep(0.1, Q2)

	param <- c(ar, AR1, AR2, ma, MA1, MA2)

	par <- nlminb(start = param, objective = loss, yt = yt, 
	              p = p, d = d, q = q, 
	              P1 = P1, D1 = D1, Q1 = Q1, P2 = P2, D2 = D2, Q2 = Q2,
	              s1 = s1, s2 = s2, k = k)
	dif  <- par$objective
	iter <- par$iterations
	par  <- par$par

	if (p  > 0)  {ar  <- par[1:p] ; par <- par[-seq(p)]}  else ar  <- NULL
	if (P1 > 0)  {AR1 <- par[1:P1]; par <- par[-seq(P1)]} else AR1 <- NULL
	if (P2 > 0)  {AR2 <- par[1:P2]; par <- par[-seq(P2)]} else AR2 <- NULL
	if (q  > 0)  {ma  <- par[1:q] ; par <- par[-seq(q)]}  else ma  <- NULL
	if (Q1 > 0)  {MA1 <- par[1:Q1]; par <- par[-seq(Q1)]} else MA1 <- NULL
	if (Q2 > 0)  {MA2 <- par[1:Q2]; par <- par[-seq(Q2)]} else MA2 <- NULL

	ma  <- maInvert(ma)
	MA1 <- maInvert(MA1)
	MA2 <- maInvert(MA2)
	return(list(ar = list(reg = ar, est1 = AR1, est2 = AR2),
	            ma = list(reg = ma, est1 = MA1, est2 = MA2),
	            d  = c(d, D1, D2), frec = c(s1, s2), 
	            dif = dif, iter = iter, data = yt))
}
