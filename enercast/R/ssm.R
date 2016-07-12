ssm <- function(formulaXa,formulaXf=NULL, data, ini=NULL,for.each='hour', cores= detectCores()-1){
  cls <- makeCluster(cores)
  registerDoParallel(cls)
  on.exit(stopCluster(cls)) ##
  h <- unique(data[,for.each])
  k <- length(h)
  RES <- foreach(j = 1:k,.final = function(x) setNames(x,h),.packages=c('FKF','enercast')) %dopar% {
    data.j <- data[data[,for.each] == h[j], names(data)!=for.each]
    ini.j <- NULL
    if(!is.null(ini))
        ini.j <- ini[[j]]
    ssm.mod(formulaXa,formulaXf,data=data.j,ini=ini.j)}
  #stopImplicitCluster()#stopCluster(cls)
  RES$formulas <- list(Xa=formulaXa,Xf=formulaXf)
  return(RES)
}
