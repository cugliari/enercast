predict.rf <- function(listRF, newdata=NULL, for.each='hour'){
  if(is.null(newdata)){
    PRED.RF <- lapply(listRF,predict)}
  if(!is.null(newdata)){
    k <- which(names(newdata)%in%for.each)
    PRED.RF <- NULL
    for(i in 1:nrow(newdata)){
      wrf <- which(names(listRF)%in%newdata[i,k])
      pred.rf <- predict(listRF[[wrf]],newdata[i,])
      PRED.RF <- c(PRED.RF,pred.rf)
    }
  }
  return(PRED.RF)
}
