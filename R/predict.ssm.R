#' predict.ssm
#'
#' @param ssm.obj 
#' @param newdata 
#' @param for.each 
#'
#' @return predicted values
#' 
#' @export
#'
#' @examples
predict.ssm <- function(ssm.obj, newdata, for.each = 'hour') {
  h <- length(ssm.obj) - 1
  auxH <- newdata[, for.each]
  newdata <- newdata[, -which(colnames(newdata) %in% for.each)]
  list.dat <- by(newdata, auxH, subset)
  
  MF <- lapply(list.dat, FUN = model.frame, formula = ssm.obj$formulas$Xa)
  XA <- lapply(MF, FUN = model.matrix, object = ssm.obj$formulas$Xa)
  
  ct <- matrix(0, 1, 1)
  if (!is.null(ssm.obj$formulas$Xf)) {
    XF <- lapply(MF, FUN = model.matrix, object = ssm.obj$formulas$Xf)
    ct <- lapply(XF, as.matrix)
  }
  
  k <- length(unique(auxH))
  for (i in 1:k) {
    cual <- names(ssm.obj) %in% names(XA)[i]
    cual <- which(cual)
    ssm.obj[[cual]]$xa.pred <- XA[[i]]
    if (!is.null(ssm.obj$formulas$Xf)) {
      ssm.obj[[cual]]$xf.pred <- ct[[i]]
    }
  }
  cuales <- names(ssm.obj) %in% names(XA)
  PRED <- lapply(ssm.obj[cuales], FUN = pred.ssm)
  
  pred <- vector('numeric', length(auxH))
  se.pred <-  vector('numeric', length(auxH))
  
  aux2 <- unique(auxH)
  for (i in 1:length(aux2)) {
    aux3 <- PRED[[aux2[i]]]
    pred[auxH == aux2[i]] <- aux3$predt
    se.pred[auxH == aux2[i]] <- sqrt(diag(aux3$covt))
  }
  
  return(list(pred = pred, se.pred = se.pred))
}
