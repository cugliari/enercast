## File : functional-clustering.r
## Description : Clustering of a functional dataset using wavelets.
## Reference : Clustering functional data using wavelets
##             Antoniadis, A., Brossat, X., Cugliari, J. and Poggi, J.-M. (2013)
##             Int. J. of Wavelets, Multiresolution and Information Processing,


## Function : toDWT 
## Description : Transforms a matrix of functional time series using the 
##               Descrete Wavelet Transform (DWT).
## Args : x  Vector with a discretized time series
##        filter.number, familiy  Options parsed to wc (see wavethresh)
## Returns : Vector with the details coefficients of the DWT  
toDWT <- function(x, filter.number = 6, family = "DaubLeAsymm"){ 
  require(wavethresh)
  x2   <- spline(x, n = 2^ceiling( log(length(x), 2) ),
                 method = 'natural')$y
  Dx2 <- wd(x2, family = family, filter.number = filter.number)$D
  Dx2
}


## Function : contrib
## Description : Computes the absolute or relative contributions of each 
##               scale of the DWT to the global energy of the vector.
## Args : x  Vector containig the details' coefficients of the DWT.
##        rel  Logical. Should the relative contributions be calculated.
## Returns : A vector containing either the absolute or the relative 
##          contributions.
contrib <- function(x, rel = FALSE, logit = rel) { 
  J   <- log( length(x)+1, 2)
  res <- numeric(J)
  t0  <- 1
  t1  <- 0
  for( j in 1:J ) {
    t1     <- t1 + 2^(J-j)
    res[j] <- sqrt( sum( x[t0:t1]^2 ) )
    t0     <- t1 + 1
  }

  if(rel) res <- res / sum(res)
  if(logit) res <- log( 1 / (1 - res) )
  
  return(res)
}



## TEST

# Create artificial data (100 trajectories of a unscaled  wiener process)
#mat <- matrix(rnorm(20 * 100), nrow = 20)
#mat <- t(apply(mat, 2, cumsum))                       # observations are on rows

# Construct the features for the clustering
#mat_dwt <- t(apply(mat, 1, toDWT))                    # DWT over the rows
#mat_contr_abs <- t(apply(mat_dwt, 1, contrib))         
#mat_contr_rel <- t(apply(mat_dwt, 1, contrib, rel = TRUE))

# Clustering using Partitioning around mediods from library cluster.
# Any other clustering technique can be used since we have a data matrix
# with observations on the rows and variables on the columns.
#library(cluster)  # we need the pam function

#mat_pam_abs <- pam(mat_contr_abs, k = 3)
#mat_pam_rel <- pam(mat_contr_rel, k = 3)

