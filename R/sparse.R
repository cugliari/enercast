library("spams")
#
# The predictive model is based on a "colinearity" argument between one day and the other
# Each traning sample consists of two consecutive days (Ltoday,Ltomorrow). A dictionary is then learnt to represent
# a set of training samples with that form. The dictionary will then be an ensemble of two-day atoms, the upper part 
# corresponding to "todays", and the lower part to "tomorrows", [Dtod^T Dtom^T]^T
# To predict tomorrow's load, w represent today's load Ltod in terms of Dtod and obtain a set of coefficients a
# Whe then estimate Ltom as Dtom*a. 
# The underlying hypothesis is that "a" would be essentially the same if we had asked to represent the whole vector [Ltod,Ltom] 
#
# INPUT
# 
# train.data ... an NxM matrix where each column is an hour by hour variable (load, temperature, etc.)
# nmode ........ normalization mode
# natoms ....... number of atoms in target dictionary
# penalty ...... penalty term to be used
#

sparse <- function(train.data, nmode,natoms, lambda) {
  
  #
  # we rearrange the data as 48xM vectors, each representing 24 hourly samples of the M input variables
  #
  m = dim(train.data)[2]
  n = dim(train.data)[1]
  # first, a "single day sample" is the concatenation of the m daily variables
  T = 24 # period, could be a parameter in the future
  N <- n / T 
  M = T*m
  #
  # normalize
  #
  X <- matrix(nrow=M,ncol=N)
  # scaling <- list()
  # centering <- list()
  for (i in 1:m) {
     Xi <- train.data[,i]
     dim(Xi) <- c(T,N)
  #   if (nmode == 'whitening') {
  #     ui <- rowMeans(Xi)
  #     Si <- var(t(Xi))
  #     Wi <- chol(Si)
  #     scaling[[i]] <- Wi
  #     centering[[i]] <- ui
  #     for (j in 1:N) {
  #       Xi[,j] <- solve(Wi,Xi[,j] - ui)
  #     }
  #   } else if (nmode == 'diag') {
  #     ui <- rowMeans(Xi)
  #     Wi <- diag(apply(Xi,1,var))
  #     scaling[[i]] <- Wi
  #     centering[[i]] <- ui
  #     for (j in 1:N) {
  #       Xi[,j] <- solve(Wi,Xi[,j] - ui)
  #     }
  #   } else if (nmode == 'max') {
  #     ui <- rowMeans(Xi)
  #     Wi <- diag(rep(abs(max(Xi)),T))
  #     scaling[[i]] <- Wi
  #     centering[[i]] <- ui
  #     for (j in 1:N) {
  #       Xi[,j] <- solve(Wi,Xi[,j] - ui)
  #     }
  #   } else if (nmode == 'none') {
  #     ui <- rep(0,T)
  #     Wi <- diag(rep(1,T))
  #     scaling[[i]] <- Wi
  #     centering[[i]] <- ui
  #     for (j in 1:N) {
  #       Xi[,j] <- solve(Wi,Xi[,j] - ui)
  #     }
  #   }
     X[((i-1)*T+1):(i*T),] <- Xi
  }
  # then two-day samples are created by concatenating X with itself, shifted one day back
  X <- rbind(X[,1:(N-1)],X[,2:N])
  # Train a sparse model of natoms elements using train.data 
  # as training data, lambda as the regularization parameter and
  # mode as the regularization mode (PENALTY, LAGRANGIAN, ERROR)
  # mode=PENALTY means that D will minimize ||X-D*A|| + lambda||A||_1
  # together with A (which is discarded for now)
  # possible values
  #  'L1COEFFS' = 0,
  #  'L2ERROR' = 1,
  #  'PENALTY' = 2,
  #  'SPARSITY' = 3,
  #  'L2ERROR2' = 4,
  #  'PENALTY2' = 5,
  #  'FISTAMODE' = 6
  
  D <- spams.trainDL(X, lambda1= lambda, K=natoms, mode='PENALTY',return_model= FALSE, verbose= FALSE)
  D0 <- D[1:M,]
  D1 <- D[(M+1):(2*M),]
  
  # store the resulting dictionary and penalty as object attributes
  #sparse <- list(Dtoday=D0,Dtomorrow=D1,M=M,lambda=lambda,nmode=nmode,centering=centering,scaling=scaling)
  sparse <- list(Dtoday=D0,Dtomorrow=D1,M=M,lambda=lambda)
  class(sparse) <- c('sparse_class')
  return(sparse)
}
