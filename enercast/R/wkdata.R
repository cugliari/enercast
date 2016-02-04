wkdata <- function(X, p, gr, J){
  nX   <- length(X)
  ncol <- nX / p
  J2   <- floor(log(p, 2)) + 1 
  if( missing(J) )  J <- J2
  if(J != round(J))  stop("p and J aren't consistents")

  if(missing(gr)) gr <- NULL
    
  if(ncol != floor(ncol)) stop("Data is incomplete")

  X1  <- matrix(X, nrow = p, ncol = ncol, byrow = FALSE)
  X2  <- spline(X, method = "natural", n = ncol * 2^J)$y
  X2  <- matrix(X2, nrow = 2^J, ncol = ncol, byrow = FALSE)

  X2  <- apply(X2, 2, function(x){
      wdx <- wd(x, filter.number = 6)
      return(c(wdx$D, wdx$C[length(wdx$C)]))
    })

  res <- list(X = X1, S0 = X2[2^J, ], D0 = X2[- (2^J), ],
              p = p, J = J, gr = gr)
  class(res) <- "wkdata"
  return(res)
}

