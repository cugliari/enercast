# DEPENDENCY: Julien Mairal SPAMS for R
library("spams")

predict.sparse <- function(sp,new.data){
  if(inherits(sp,'sparse_class',TRUE) != 1)
    stop('Object is not of sparse_class')
  
  m = dim(new.data)[2]
  n = dim(new.data)[1]
  #
  # construct "single day samples" as  the concatenation of the m daily variables
  T = 24 # period, could be a parameter in the future
  N <- n / T 
  M = T*m
  #
  # normalize this variable to be N(0,Identity)
  #
  Xi <- new.data[,1]
  dim(Xi) <- c(T,N)
  print(sp)
  for (j in 1:N) {
    Xi[,j] <- solve(sp$whitening[[1]],Xi[,j] - sp$centering[[1]])
  }
  X <- Xi  
  for (i in 2:m) {
    print(i)
    Xi <- new.data[,i]
    dim(Xi) <- c(T,N)
    for (j in 1:N) {
      Xi[,j] <- solve(sp$whitening[[i]],Xi[,j] - sp$centering[[i]])
    }
    dim(Xi) <- c(T,N)
    X <- rbind(X,Xi)
  }
  M0 <- sp$M
  D0 <- sp$dictionary[1:M0,]
  D1 <- sp$dictionary[(M0+1):(2*M0),]
  A <- spams.lasso(X,D0,lambda1=sp$lambda,mode='PENALTY')
  X1 <- D1 %*% A
  for (j in 1:N) {
    Xj <- X1[,j]
    for (i in 1:m) {
      Xji <- Xj[((i-1)*T+1):(i*T)]
      X1[((i-1)*T+1):(i*T),j] <- sp$whitening[[i]] %*% Xji + sp$centering[[i]]
    }
  }
  return(X1)
}
