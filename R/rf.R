rf <- function(y='load', data, for.each='hour', cores= detectCores()-1){
  cls <- makeCluster(cores)
  registerDoParallel(cls)
  on.exit(stopCluster(cls))
  K <- unique(data[,for.each])
  RES <- foreach(j = K ,.packages='randomForest') %dopar% {
      data.j <- data[data[,for.each] == j, names(data)!=for.each]
      XX <- data.j[ , names(data.j)!=eval(y)]
      YY <- data.j[,eval(y)] 
      rf <- randomForest(y=YY, x=XX, na.action=na.omit)}
  names(RES) <- K
  return(RES)
}




