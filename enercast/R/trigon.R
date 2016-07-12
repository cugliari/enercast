trigon <- function(data,frecH=1,n=nrow(data),w=24,q=3){

    is.leapyear <- function(year){
  return(((year %% 4 == 0) & (year %% 100 != 0)) | (year %% 400 == 0))}


    if(!is.numeric(data$hour))
        data$hour <- as.numeric(data$hour)
    if(!is.numeric(data$year))
        data$year <- as.numeric(data$year)
    
    if (frecH==1)
        xtime <- data$hour
    if (frecH==7)
        xtime <- rep(1:(w*frecH),length=nrow(data))
    if (frecH==365){
        anio <- 365 + is.leapyear(unique(data$year))
        xtime <- unlist(sapply(anio,function(x)seq(1,w*x)))}

    RES <- NULL
    for (i in 1:q){
        s1 <- sin(xtime*2*pi*i/(frecH*w))
        c1 <- cos(xtime*2*pi*i/(frecH*w))
        RES <- cbind(RES,s1,c1)
    }
    colnames(RES) <- paste(rep(c('sin','cos'),q),rep(1:q,each=2),sep='_')
    return(RES)
}
