#' sparse
#'
#' The predictive model is based on a "colinearity" argument between one day and
#' the other
#' Each traning sample consists of two consecutive days (Ltoday,Ltomorrow). A 
#' dictionary is then learnt to represent a set of training samples with that 
#' form. The dictionary will then be an ensemble of two-day atoms, the upper 
#' part corresponding to "todays", and the lower part to "tomorrows", 
#' [Dtod^T Dtom^T]^T
#' To predict tomorrow's load, w represent today's load Ltod in terms of Dtod 
#' and obtain a set of coefficients a
#' Whe then estimate Ltom as Dtom*a. 
#' The underlying hypothesis is that "a" would be essentially the same if we had asked to represent the whole vector [Ltod,Ltom] 
#'
#' @param train.data an NxM matrix where each column is an hour by hour variable (load, temperature, etc.)
#' @param natoms normalization mode
#' @param lambda number of atoms in target dictionary
#' @param delta penalty term to be used
#'
#' @return a list with dictionaries
#' @examples 3
#' @export
#library("spams")
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

sparse <- function(train.data, natoms, lambda, delta = 24) {

  if (!requireNamespace("spams", quietly = TRUE)) {
    stop("Package spams is needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  #
  # Rearrange the data as 2 * delta x M vectors (e.g. each representing 
  # 24 hourly samples of the M input variables)
  m = dim(train.data)[2]
  n = dim(train.data)[1]
  # first, a "single day sample" is the concatenation of the m daily variables
#  T = 24 # period, could be a parameter in the future
  N <- n / delta 
  M <- delta * m
  #
  # normalize
  #
  X <- matrix(nrow = M, ncol = N)
  for (i in 1:m) {
    Xi <- train.data[,i]
    dim(Xi) <- c(delta, N)
    X[((i - 1) * delta + 1):(i * delta),] <- Xi
  }
  #
  # normalize each column
  #
  #
  # then two-day samples are created by concatenating X with itself, shifted one day back
  #
  X <- rbind(X[,1:(N - 1)], X[, 2:N])
  X <- scale(X)
  D <- spams.trainDL(X, lambda1 = lambda, K = natoms, mode = 'PENALTY',
                     return_model = FALSE, verbose = FALSE)
  D0 <- D[1:M, ]
  D1 <- D[(M + 1):(2 * M),]
  
  # store the resulting dictionary and penalty as object attributes
  #sparse <- list(Dtoday=D0,Dtomorrow=D1,M=M,lambda=lambda,nmode=nmode,centering=centering,scaling=scaling)
  sparse <- list(Dtoday = D0, Dtomorrow = D1, M = M, lambda = lambda)
  class(sparse) <- c('sparse_class')
  return(sparse)
}
