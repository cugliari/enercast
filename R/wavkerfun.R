#' Calibrate the KWF model 
#' 
#' Using   
#'
#' @param obj 
#' @param dist 
#' @param r 
#' @param H 
#' @param EPS 
#'
#' @return fitted wavkerfun model
#' 
#' @export
#'
#' @examples
wavkerfun <-  function(obj, dist = 1, r, H, EPS = 1e-6){

  orig <- obj$X
  n    <- ncol(obj$X)
  p    <- obj$p
  x    <- obj$D0
  y    <- obj$D0
#  c0   <- obj$S0
  gr   <- obj$gr
  
 
 # Checking ...
 if (is.null(gr))     gr <- as.factor(rep(1, n))
 if (n != length(gr)) stop("Check dimensions: gr")

 grs   <-  which(gr[n] == gr[-n])
 ng    <-  length(grs)

 if (missing(r)) r <- floor(ng / 10) 
 if (r > ng)     stop("Check dimensions: r is too big")

 if (missing(H))  {
   sigma <- mean(sd(orig))
   H     <- numeric(2)
   H[1]  <- sigma / 8
   H[2]  <- 3 * sigma
 }

 f <- function(h, x, y, n, p_x, p_y, r) CVkerfon(x, y, p, r, h, EPS)

 if (H[1] == H[2]) {
   hopt <- H[1]
   } else {
   hopt <- optimize(f, H, x[, grs], y[, grs + 1], ng,
                    p_past, p_fut, r)$minimum
   }
   
 res <- list(X = obj$X, S0 = obj$S0, D0 = obj$D0, p = p,
             J = obj$J, gr = gr, h = hopt) 

 class(res) <- "wavkerfun"
 return(res)
}

