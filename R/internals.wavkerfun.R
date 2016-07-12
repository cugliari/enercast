wgt <- function (X, gr = NULL, h){ 
   
   nr   <- nrow(X)
   nc   <- ncol(X)

  # Check groups
  if( is.null(gr) ) gr <- rep(1,nc)
  if( length(gr) != nc )  stop("Check dimensions: wrong length of gr")

  grE  <- which( gr[-nc] == gr[nc] )

  w    <- apply( X[,grE ], 2, "-", X[,nc])
  w    <- apply( w, 2, DistWav)

  # Check h for multi-dimensional optiminization
  if( length(h) > 1 ) h <- h[gr[nc]]

  # Normalization
  w    <- exp ( - 0.5*( w/h )^2 )
  sumw <- sum(w)

  wc   <- numeric(nc-1)
  wc[grE] <- w

  if( sumw == 0) { # If no similarity is founded
    w <- rep( 1 / (nc - 1), nc-1 ) 
  } else { 
    w <- wc/sumw 
  }  
  w
}
