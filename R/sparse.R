library("spams")
#library("lattice") what was this one for??

sparse <- function(train.data,natoms, lambda) {

    # Train a sparse model of natoms elements using train.data 
    # as training data, lambda as the regularization parameter and
    # mode as the regularization mode (PENALTY, LAGRANGIAN, ERROR)
    # mode=PENALTY means that D will minimize ||X-D*A|| + lambda||A||_1
    # together with A (which is discarded for now)

    D <- spams.trainDL(Xtn, lambda1= lambda, K=natoms, mode='PENALTY',return_model= FALSE, verbose= FALSE)

    # store the resulting dictionary and penalty as object attributes
    sparse <- list(dictionary=D,lambda=lambda)
    class(sparse) <- c('sparse_class')
    return(sparse)
}
