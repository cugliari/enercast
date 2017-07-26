#' Functional wavelet kernal
#'
#' @param wavkerfun object 
#' @param ... 
#' @param newdata 
#' @param sd 
#'
#' @return predicted values
#' 
#' @export
#'
#' @examples
predict.wavkerfun <- function(object, ..., newdata = NULL, sd){ 
  if ((!is.null(class(object))) && (class((object)) != "wavkerfun"))
      stop("object is not of class wavkerfun")

  if ( !is.null(newdata) )  stop("newdata is not used for prediction")

  h     <- object$h
  p     <- object$p
  J     <- object$J
  N     <- 2^J
  
  x0     <- object$D0
  y0     <- object$D0
  gr0    <- object$gr
  n0     <- ncol(x0)  

  w   <- wgt(x0, gr = gr0, h = h )
  Mw  <- matrix(rep(w, nrow(y0)), ncol = n0 - 1, byrow = TRUE)
  Mres   <- y0[,-1] * Mw
  predD0 <- apply( Mres, 1, sum )
  predS0 <- 0

  c0  <- object$S0[-1]
  lev <- c0[length(c0)] + sum( w[-1]*diff(c0) )
  predS0 <- lev 
  
   
  # Reconstruction by IDWT 
  empty   <- wavethresh::wd(1:N, filter.number = 6)
  empty$D <- predD0
  empty   <- wavethresh::putC(empty, level = 0, predS0)
  predX   <- spline(wavethresh::wr(empty), n = p)$y 

  res <- list(X = predX, S0 = predS0, D0 = predD0, p = p, J = J)
  class(res) <- "wkdata"
  return( res )
}

