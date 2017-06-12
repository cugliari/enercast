comb.coef <- function(rv, Sv, s){
	rn <- length(rv); Sn <- length(Sv)
	if (rn == 1 & 0 %in% rv) { 
		v <- rep(0, s * Sn);
		for (i in 1:Sn) v[s * i] <- Sv[i]
		return(v)
	}
	if ( rn >= 1 & Sn >= 1 & sum(rv != 0) != 0 & sum(Sv != 0) != 0) {
		v <- rep(0, (s * Sn + rn + 1))
		for (i in 1:(Sn + 1)) {
			for (j in 1:(rn + 1)) {
				v[j + s * (i - 1)] <- c(1, rv)[j] * c(1, Sv)[i]
			}
		}
		v <- v[-1]
	return(v)
	}
}


###
maInvert <- function(ma){
    q <- length(ma)
    q0 <- max(which(c(1, ma) != 0)) - 1
    if (!q0) return(ma)
    roots <- polyroot(c(1, ma[1:q0]))
    ind <- Mod(roots) < 1
    if (all(!ind)) return(ma)
    if (q0 == 1) return(c(1/ma[1], rep(0, q - q0)))
    roots[ind] <- 1/roots[ind]
    x <- 1
    for (r in roots) x <- c(x, 0) - c(0, x)/r
    c(Re(x[-1]), rep(0, q - q0))
}


###
acft <- function(ar = NULL, AR1 = NULL, AR2 = NULL, ma = NULL,
                 MA1 = NULL, MA2 = NULL, s1 = 24, s2 = 168, k = 500){

	if (is.null(ar))  ar = 0
	if (is.null(AR1)) AR1 = 0
	if (is.null(AR2)) AR2 = 0
	if (is.null(ma))  ma = 0
	if (is.null(MA1)) MA1 = 0
	if (is.null(MA2)) MA2 = 0
	ma <- maInvert(ma)
	MA1 <- maInvert(MA1)
	MA2 <- maInvert(MA2)
	
	maT = ma
	arT = ar
	if (any(MA1 != 0)) maT <- comb.coef(rv = ma,  Sv = MA1, s = s1) else maT = ma
	if (any(AR1 != 0)) arT <- comb.coef(rv = ar,  Sv = AR1, s = s1) else arT = ar
	if (any(MA2 != 0)) maT <- comb.coef(rv = maT, Sv = MA2, s = s2) else maT = maT
	if (any(AR2 != 0)) arT <- comb.coef(rv = arT, Sv = AR2, s = s2) else arT = arT
	
	ARMAacf(ar = arT, ma = maT, lag.max = k)[-1]
}



###
loss <- function(par, yt, p = 0 , d = 0, q = 0,
                     P1 = 0, D1 = 0, Q1 = 0, P2 = 0, D2 = 0, Q2 = 0,
                     s1 = 24, s2 = 168, k = 200){
  
  ytc <- yt
    
	if (p > 0) {
		ar <- par[1:p]
		par <- par[-seq(p)]
	} else ar <- NULL
	if (P1 > 0) {
		AR1 <- par[1:P1]
		par <- par[-seq(P1)]
	} else AR1 < -NULL
	if (P2 > 0) {
		AR2 <- par[1:P2]
		par <- par[-seq(P2)]
	} else AR2 <- NULL
	if (q > 0) {
		ma  <- par[1:q]
		par <- par[-seq(q)]
	} else ma <- NULL
	if (Q1 > 0) {
		MA1 <- par[1:Q1]
		par <- par[-seq(Q1)]
	} else MA1 <- NULL
	if (Q2 > 0) {
		MA2 <- par[1:Q2]
		par <- par[-seq(Q2)]
	} else MA2 <- NULL
	if (d > 0)  ytc <- diff(ytc, differences = d)
	if (D1 > 0) ytc <- diff(ytc, differences = D1, lag = s1)
	if (D2 > 0) ytc <- diff(ytc, differences = D2, lag = s2)
	acf1 <- acf(ytc, lag.max = k, plot = FALSE)$acf[-1]
	acfteo <- acft(ar = ar, AR1 = AR1, AR2 = AR2, ma = ma, 
	               MA1 = MA1, MA2 = MA2, s1 = s1, s2 = s2, k = k)
	
	return(sum((acfteo - acf1)^2))
}
