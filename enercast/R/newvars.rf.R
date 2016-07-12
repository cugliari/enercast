newvars.rf <- function(data,type,t,...){
    xrezaga <- function(vec,t,nom='xtemp'){
        n <- length(vec) ;     mtz <- matrix(0,n,t) ;    diag(mtz) <- 1
        AUX <- apply(apply(mtz,2,cumsum),2,cumsum)
        M <- matrix(NA,nrow(AUX),ncol(AUX));  colnames(M) <- c(nom,paste(nom,seq(t-1),sep='.'))
        for (j in 1:ncol(mtz))
            M[,j] <- c(rep(NA,(j-1)),vec[AUX[,j]])
        return(M[,-1])
    }
    xsqcb <- function(vec){
        return(cbind(temp2=vec^2,temp3=vec^3))
    }
  
    xdesv <- function(vec,t,nom='dvtemp'){
        n <- length(vec) ;     mtz <- matrix(0,n,t) ;    diag(mtz) <- 1
        AUX <- apply(apply(mtz,2,cumsum),2,cumsum)
        AUX <- AUX[-seq(t-1), -1]
        M <- matrix(NA,nrow(AUX),ncol(AUX))
        M.upper <- matrix(NA,(t-1),ncol(AUX))
        for (j in 1:ncol(AUX)){
            z1 <- AUX[ ,1] ;      zj <- AUX[ ,j]
            M[ ,j] <- abs(((vec[z1]+vec[zj])/2)-z1)  }
        M <- rbind(M.upper,M);  colnames(M) <- paste(nom,seq(t-1),sep='.')
        return(M)
    }

    switch(type,
           load.lag  = xrezaga(data$load,t=t,nom='xload'),
           temp.lag  = xrezaga(data$temp,t=t,nom='xtemp'),
           temp.sqcb = xsqcb(data$temp),
           temp.desv = xdesv(data$temp,t=t),
           load.desv = xdesv(data$load,t=t,nom='dvload'))
}

