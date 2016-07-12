rcs <- function(x,m=5,nod=NULL){
    if(is.null(nod))
        nod <- quantile(x,probs=seq(0,1,1/(m+1)))[-c(1,m+2)]
    if(is.null(m))
        m <- length(nod)
    Cs <- matrix(NA,nrow=length(x),ncol=m-2)
    
    for (i in 1:(m-2))
        Cs[,i] <- apply(cbind(0,(x-nod[i])^3),1,max) - apply(cbind(0,(x-nod[m-1])^3),1,max)*(nod[m]-nod[i])/(nod[m]-nod[m-1]) + apply(cbind(0,(x-nod[m])^3),1,max)*(nod[m-1]-nod[i])/(nod[m]-nod[m-1])
    
    Cs <- cbind(x,Cs)
    colnames(Cs) <- paste('C',0:(m-2),sep='')
    Cs <- sweep(Cs,MARGIN=2,FUN="/",STATS=apply(Cs,2,function(x)sqrt(sum(x^2))))
    return(Cs)
}

