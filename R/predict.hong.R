#' predict.hong
#'
#' Predict method for Hong Vanilla Model
#' 
#' @param obj.hv Object of class inheriting from use of \code{hong} function. 
#' @param new.data Data frame with predictors. If omitted, the fitted values are used. 
#'
#' @return predicted values
#' @export
#' @details
#' \code{new.data} must contain the same predictors that data in \link{hong}.
#' Aditionally \code{new.data} must contain \code{trend} variable. This variable was computed like a sequence from 1 to number of rows of input data.frame in \link{hong} model function.
#' In forecasting new consecutive periods, \code{trend} variable in new.data must be computed like a sequence from number of rows of original data.
#' @author Jairo Cugliari, Andres Castrillejo, Fernando Massa, Ignacio Ramirez
#' @seealso \code{\link[stats]{predict.lm}}\code{\link{hong}
#'
predict.hong <- function(obj.hv,new.data=NULL){
  if(inherits(obj.hv,'hv_class',TRUE)!=2)
    stop('Object is not of hv_class')
  pred <- predict(object=obj.hv,newdata=new.data)

}


