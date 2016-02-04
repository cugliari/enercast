
predict.hong <- function(obj.hv,new.data=NULL){
  if(inherits(obj.hv,'hv_class',TRUE)!=2)
    stop('Object is not of hv_class')
  pred <- predict(object=obj.hv,newdata=new.data)

}


