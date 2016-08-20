# DEPENDENCY: Julien Mairal SPAMS for R
library("spams")

predict.sparse <- function(sp,new.data) {
  if(inherits(sp,'sparse_class',TRUE) != 1)
    stop('Object is not of sparse_class')
  
  m = dim(new.data)[2]
  n = dim(new.data)[1]
  #
  # construct "single day samples" as  the concatenation of the m daily variables
  T = 24 # period, could be a parameter in the future
  N <- n / T 
  M = T*m
  X <- matrix(nrow=M,ncol=N)
  for (i in 1:m) {
    Xi <- new.data[,i]
    #    for (j in 1:N) {
    #      Xi[,j] <- solve(sp$scaling[[i]],Xi[,j] - sp$centering[[i]])
    #    }
    dim(Xi) <- c(T,N)
    X[((i-1)*T+1):(i*T),] <- Xi
  }
  #
  # trick: normalize but by twice the variance
  # we don't know the variance of the part to be estimated, but that can be assumed to be the same
  # as that of the part we do have.
  cm <- apply(X,2,mean)
  cv <- apply(X,2,sd)
  for (j in 1:N) {
    X[,j] <- ( X[,j] - cm[j]) / (2.*cv[j])
  }
  A <- spams.lasso(X,sp$Dtoday,lambda1=sp$lambda,mode='PENALTY')
  X1 <- sp$Dtomorrow %*% A
  for (j in 1:N) {
    X1[,j] <- 2.*cv[j]*X1[,j] + cm[j]
  }
  
  for (j in 1:N) {
    Xj <- X1[,j]
    for (i in 1:m) {
      Xji <- Xj[((i-1)*T+1):(i*T)]
      X1[((i-1)*T+1):(i*T),j] <- Xji
    }
  }
  return(X1)
}