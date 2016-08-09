# DEPENDENCY: Julien Mairal SPAMS for R
library("spams")

predict.sparse <- function(sp,new.data){
  if(inherits(sp,'sparse_class',TRUE) != 1)
    stop('Object is not of sparse_class')

    A <- spams.lasso(Xvn,sp.D,lambda1=sp.lambda,mode='PENALTY')
    Xhat <- sp.D %*% A


}

