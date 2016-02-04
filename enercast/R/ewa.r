#' Calibrate the EWA algorithm
#' 
#' \code{ewa} agregates expert predictions using exponential weights and a
#' learning rate parameter \eqn{\eta} . 
#' 
#' @param eta : numeric value to start the initial grid for \eqn{\eta}
#' @param y : numeric vector contianing the time series
#' @param experts : matrix of time points X individual predictions (experts) 
#' @param EPS     : Epsilon tolerance
#' @return   A list containg: \code{w} the matrix of fitted weights and \code{eta} the obtain
#'           value for \eqn{\eta}.
ewa <- function(eta, y, experts, EPS = 1e-5){

  experts <- as.matrix(experts)
  p <- ncol(experts)     # Number of experts
  n <- nrow(experts)
  
  if(length(y) != n)    stop("Length of y must be equal to nrow(experts).")


  regmat <- rep(0, p)
  w <- matrix(0, ncol = p, nrow = n)

  w[1, ] <- rep(1/p, p)
  
  # Recursively predict and update weights from past losses
  for(time in 1:(n - 1)) {
    # prediction
    pred_n         <- experts[time, ] %*% w[time, ]
    loss_pred_n    <- as.numeric(loss(pred_n      , y[time])) 
    loss_experts_n <- as.numeric(loss(experts[time, ], y[time]))
    
    # update weights
    regmat <- regmat + (loss_pred_n - loss_experts_n)
  
    w[time + 1, ] <- exp(eta * regmat) 
    w[time + 1, ] <- checkweights(w[time + 1, ], EPS) # Normalize weights and eps-correction
  }
  
  return(list(w = w, eta = eta))
}


checkweights <- function(w, EPS) {
    w[w < EPS] <- 0
    w / sum(w)
}  

loss <- function(x, y) (x - y)^2  

