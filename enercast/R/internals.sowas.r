###########################################################
# CONTINUOUS WAVELET TRANSFORMATION OF A TIME SERIES OBJECT
###########################################################

cwt.ts <- function(ts,s0,noctave=5,nvoice=10,w0=2*pi){
  
  if (class(ts)!="ts"){
    
    cat("# This function needs a time series object as input. You may construct this by using the function ts(data,start,deltat). Try '?ts' for help.\n")
    
  }
  else{
    
    t=time(ts)
    dt=t[2]-t[1]
    
    s0unit=s0/dt*w0/(2*pi)   
    s0log=as.integer((log2(s0unit)-1)*nvoice+1.5)
    
    if (s0log<1){
      cat(paste("# s0unit = ",s0unit,"\n",sep=""))
      cat(paste("# s0log  = ",s0log,"\n",sep=""))
      cat("# s0 too small for w0! \n")
    }
    totnoct=noctave+as.integer(s0log/nvoice)+1
    
    totts.cwt=cwt(ts,totnoct,nvoice,w0,plot=0)
    
    ts.cwt=totts.cwt[,s0log:(s0log+noctave*nvoice)]
    
    #Normalization
    sqs <- sqrt(2^(0:(noctave*nvoice)/nvoice)*s0)
    smat <- matrix(rep(sqs,length(t)),nrow=length(t),byrow=TRUE)
    
    ts.cwt*smat
    
  }
  
}

#####################################
# WSP
#####################################

wsp <- function(ts,s0=1,noctave=5,nvoice=10,w0=2*pi,sw=0,tw=0,swabs=0,siglevel=0.95,critval=NULL,nreal=1000,arealsiglevel=0.9,kernel=0,markt=-999,marks=-999,logscale=FALSE,plot=TRUE,units="",device="screen",file="wsp",color=TRUE,pwidth=10,pheight=7,labsc=1,labtext="",sigplot=3){
  
  if (class(ts)!="ts"){
    
    cat("# This function needs a time series object as input. You may construct this by using the function ts(data,start,deltat). Try '?ts' for help.\n")
    
  }
  else{
    
    if (sw!=0 & swabs==0)
      swabs <- as.integer(sw*nvoice)
    if (swabs!=0 & sw==0)
      sw <- swabs/nvoice
    
    sllist <- checkarealsiglevel(sw,tw,w0,arealsiglevel,siglevel,0)
    arealsiglevel <- sllist$arealsiglevel
    siglevel <- sllist$siglevel
    
    at <- NULL
    
    t <- time(ts)
    dt <- deltat(ts)
    s0rem <- s0
    s0 <- adjust.s0(s0,dt)
    dnoctave <- as.integer(log(s0/s0rem-0.000000000001)/log(2))+1
    
    noctave <- adjust.noctave(length(ts),dt,s0,tw,noctave)
    scalevector <- 2^(0:(noctave*nvoice)/nvoice)*s0
    tsanom <- ts-mean(ts)
    
    #WAVELET TRANSFORMATION
    ts.cwt <- cwt.ts(tsanom,s0,noctave,nvoice,w0)
    
    #SMOOTHING
    wsp <- smooth.matrix(Mod(ts.cwt)^2,swabs)
    smwsp <- smooth.time(wsp,tw,dt,scalevector)
    
    #POINTWISE SIGNIFICANCE TEST
    if (is.null(critval)==FALSE){ # is critval empty?
      if (dim(critval)[2]!=dim(smwsp)[2]){ # is critval of the wrong dimension?
        if (siglevel[1]!=0 & nreal!=0) critval <-
          criticalvaluesWSP(tsanom,s0,noctave,nvoice,w0,swabs,tw,siglevel,nreal)
        #critval is of wrong dimension and siglevel and nreal are given
        else {
          critval <- NULL # no test possible
          arealsiglevel <- 0
          cat("# dimension of critval is wrong \n")
          cat("# areawise test only possible with pointwise test \n")
        }
      }
    }
    else{ # critval is empty, but nreal or siglevel is given
      if (siglevel[1]!=0 & nreal!=0) critval <-
        criticalvaluesWSP(tsanom,s0,noctave,nvoice,w0,swabs,tw,siglevel,nreal)
      else {
        critval <- NULL
        arealsiglevel <- 0
        cat("# areawise test only possible with pointwise test \n")
      }
    }
    
    #AREAL SIGNIFICANCE TEST
    if (arealsiglevel!=0){
      v <- critval[1,]
      v[v==0] <- 9999
      cvmat <- matrix(rep(v,length(t)),nrow=length(t),byrow=TRUE)
      atest <- arealtest(smwsp/cvmat,dt,s0,noctave,nvoice,w0,swabs,tw,siglevel,arealsiglevel,kernel,0)
      at <- atest$at
      kernel <- atest$kernel
    }
    
    if (s0rem<s0){
      smwsp <- addvalues(nvoice,dnoctave,smwsp,NA)
      critval <- addvalues(nvoice,dnoctave,critval,1)
      #at <- addvalues(nvoice,dnoctave,at,NA)
      noctave <- noctave+dnoctave
      s0 <- s0/2^dnoctave
      scalevector <- 2^(0:(noctave*nvoice)/nvoice)*s0
    }
    
    #PARAMETERS
    wclist <-
      list(modulus=smwsp,phase=NULL,time=t,s0=s0,noctave=noctave,nvoice=nvoice,w0=w0,scales=scalevector,critval=critval,at=at,kernel=kernel)
    
    class(wclist) <- "wt"
    
    #PLOTTING
    if (plot)
      plot(wclist,markt,marks,NULL,NULL,logscale,FALSE,units,"Wavelet Power Spectrum",device,file,FALSE,color,pwidth,pheight,labsc,labtext,sigplot)
    
    wclist
    
  }
  
}

#####################################
# WCO
#####################################

wco <-
  function(ts1,ts2,s0=1,noctave=5,nvoice=10,w0=2*pi,sw=0,tw=0,swabs=0,siglevel=0.95,arealsiglevel=0.9,kernel=0,markt=-999,marks=-999,sqr=FALSE,phase=TRUE,plot=TRUE,units="",device="screen",file="wcoh",split=FALSE,color=TRUE,pwidth=10,pheight=7,labsc=1,labtext="",sigplot=3){
    
    if (class(ts1)!="ts"){
      
      cat("# This function needs two time series objects as input. You may construct them by using the function ts(data,start,deltat). Try '?ts' for help.\n")
      
    }
    else{
      
      if (sw!=0 & swabs==0)
        swabs <- as.integer(sw*nvoice)
      if (swabs!=0 & sw==0)
        sw <- swabs/nvoice
      
      sllist <- checkarealsiglevel(sw,tw,w0,arealsiglevel,siglevel,1)
      arealsiglevel <- sllist$arealsiglevel
      siglevel <- sllist$siglevel
      
      if (sw==0 & tw==0 & swabs==0) {
        cat("# coherence without averaging makes no sense! \n")
        siglevel <- 0
        arealsiglevel <- 0
      }
      
      if (phase==FALSE) phs <- NULL
      
      at <- NULL
      
      tsadjust <- adjust.timeseries(ts1,ts2)
      ts1 <- tsadjust$ts1
      ts2 <- tsadjust$ts2
      
      t <- time(ts1)
      dt <- deltat(ts1)
      
      s0rem <- s0
      s0 <- adjust.s0(s0,dt)
      dnoctave <- as.integer(log(s0/s0rem-0.000000000001)/log(2))+1
      
      noctave <- adjust.noctave(length(ts1),dt,s0,tw,noctave)
      
      scalevector <- 2^(0:(noctave*nvoice)/nvoice)*s0
      
      ts1anom <- ts1-mean(ts1)
      ts2anom <- ts2-mean(ts2)
      
      ts1.cwt <- cwt.ts(ts1anom,s0,noctave,nvoice,w0)
      ts2.cwt <- cwt.ts(ts2anom,s0,noctave,nvoice,w0)
      
      cosp <- Re(ts1.cwt)*Re(ts2.cwt) + Im(ts1.cwt)*Im(ts2.cwt)
      quad <- Im(ts1.cwt)*Re(ts2.cwt) - Re(ts1.cwt)*Im(ts2.cwt)
      wsp1 <- Mod(ts1.cwt)^2
      wsp2 <- Mod(ts2.cwt)^2
      
      smcosp <- smooth.matrix(cosp,swabs)
      smquad <- smooth.matrix(quad,swabs)
      smwsp1 <- smooth.matrix(wsp1,swabs)
      smwsp2 <- smooth.matrix(wsp2,swabs)
      
      smsmcosp <- smooth.time(smcosp,tw,dt,scalevector)
      smsmquad <- smooth.time(smquad,tw,dt,scalevector)
      smsmwsp1 <- smooth.time(smwsp1,tw,dt,scalevector)
      smsmwsp2 <- smooth.time(smwsp2,tw,dt,scalevector)
      
      if (sqr==FALSE)
        wcoh <- sqrt((smsmcosp^2+smsmquad^2)/(smsmwsp1*smsmwsp2))
      else
        wcoh <- (smsmcosp^2+smsmquad^2)/(smsmwsp1*smsmwsp2)
      
      if (phase)
        phs <- atan2(smsmquad,smsmcosp)
      else phs <- NULL
      
      #POINTWISE SIGNIFICANCE TEST
      if (siglevel[1]!=0) critval <- criticalvaluesWCO(s0,noctave,nvoice,w0,swabs,tw,siglevel)
      else critval <- NULL
      if (sqr==TRUE & is.null(critval)==FALSE)
        critval <- critval^2
      
      #AREAWISE SIGNIFICANCE TEST
      if (arealsiglevel!=0){
        atest <- arealtest(wcoh/critval,dt,s0,noctave,nvoice,w0,swabs,tw,siglevel,arealsiglevel,kernel,1)
        at <- atest$at
        kernel <- atest$kernel
      }
      
      wcoh[1,1] <- 0
      wcoh[1,2] <- 1
      
      if (phase){
        phs[1,1] <- -pi
        phs[1,2] <- pi 
      }
      
      if (s0rem<s0){
        wcoh <- addvalues(nvoice,dnoctave,wcoh,NA)
        phs <- addvalues(nvoice,dnoctave,phs,NA)
        noctave <- noctave+dnoctave
        s0 <- s0/2^dnoctave
        scalevector <- 2^(0:(noctave*nvoice)/nvoice)*s0
      }
      
      wclist <- list(modulus=wcoh,phase=phs,s0=s0,noctave=noctave,nvoice=nvoice,w0=w0,time=t,scales=scalevector,critval=critval,at=at,kernel=kernel)
      
      class(wclist) <- "wt"
      
      if (plot) plot(wclist,markt,marks,NULL,NULL,FALSE,phase,units,"Wavelet Coherence",device,file,split,color,pwidth,pheight,labsc,labtext,sigplot)
      
      wclist
      
    }
    
  }

#####################################
# WCS
#####################################

wcs <- function(ts1,ts2,s0=1,noctave=5,nvoice=10,w0=2*pi,sw=0,tw=0,swabs=0,markt=-999,marks=-999,logscale=FALSE,phase=TRUE,plot=TRUE,units="",device="screen",file="wcsp",split=FALSE,color=TRUE,pwidth=10,pheight=7,labsc=1,labtext=""){
  
  if (class(ts1)!="ts"){
    
    cat("# This function needs two time series objects as input. You may construct them by using the function ts(data,start,deltat). Try '?ts' for help. \n")
    
  }
  else{
    
    if (sw!=0 & swabs==0)
      swabs <- as.integer(sw*nvoice)
    if (swabs!=0 & sw==0)
      sw <- swabs/nvoice
    
    tsadjust <- adjust.timeseries(ts1,ts2)
    ts1 <- tsadjust$ts1
    ts2 <- tsadjust$ts2
    
    t <- time(ts1)
    dt <- deltat(ts1)
    
    s0rem <- s0
    s0 <- adjust.s0(s0,dt)
    dnoctave <- as.integer(log(s0/s0rem-0.000000000001)/log(2))+1
    
    noctave <- adjust.noctave(length(ts1),dt,s0,tw,noctave)
    
    scalevector <- 2^(0:(noctave*nvoice)/nvoice)*s0
    
    ts1anom <- ts1-mean(ts1)
    ts2anom <- ts2-mean(ts2)
    
    ts1.cwt <- cwt.ts(ts1anom,s0,noctave,nvoice,w0)
    ts2.cwt <- cwt.ts(ts2anom,s0,noctave,nvoice,w0)
    
    cosp <- Re(ts1.cwt)*Re(ts2.cwt) + Im(ts1.cwt)*Im(ts2.cwt)
    quad <- Im(ts1.cwt)*Re(ts2.cwt) - Re(ts1.cwt)*Im(ts2.cwt)
    
    smcosp <- smooth.matrix(cosp,swabs)
    smquad <- smooth.matrix(quad,swabs)
    
    smsmcosp <- smooth.time(smcosp,tw,dt,scalevector)
    smsmquad <- smooth.time(smquad,tw,dt,scalevector)
    
    wcsp <- smsmcosp^2+smsmquad^2
    
    if (phase)
      phs <- atan2(smsmquad,smsmcosp)
    else phs <- NULL
    
    if (phase){
      phs[1,1] <- -pi
      phs[1,2] <- pi 
    }
    
    if (s0rem<s0){
      wcsp <- addvalues(nvoice,dnoctave,wcoh,NA)
      phs <- addvalues(nvoice,dnoctave,phs,NA)
      noctave <- noctave+dnoctave
      s0 <- s0/2^dnoctave
      scalevector <- 2^(0:(noctave*nvoice)/nvoice)*s0
    }
    
    wclist <- list(modulus=wcsp,phase=phs,s0=s0,noctave=noctave,nvoice=nvoice,w0=w0,time=t,scales=scalevector,critval=NULL,at=NULL,kernel=NULL)
    
    class(wclist) <- "wt"
    
    if (plot) plot(wclist,markt,marks,NULL,NULL,logscale,phase,units,"Wavelet Cross Spectrum",device,file,split,color,pwidth,pheight,labsc,labtext)
    
    wclist
    
  }
  
}

##########################################
# POINTWISE SIGNIFICANCE TEST
##########################################

rawWSP <- function(ts,s0=1,noctave=5,nvoice=20,w0=2*pi,swabs=0,tw=0,scalevector){
  
  tsanom <- ts-mean(ts)
  
  ts.cwt <- cwt.ts(tsanom,s0,noctave,nvoice,w0)
  
  wsp <- Mod(ts.cwt)^2
  
  smwsp <- smooth.matrix(wsp,swabs)
  smsmwsp <- smooth.time(smwsp,tw,deltat(ts),scalevector)
  
  smsmwsp
  
}

criticalvaluesWSP <- function(ts,s0=1,noctave=5,nvoice=10,w0=2*pi,swabs=0,tw=0,siglevel=0.95,nreal=1000){
  
  t=time(ts)
  dt=deltat(ts)
  
  s0 <- adjust.s0(s0,dt)
  noctave <- adjust.noctave(length(ts),dt,s0,tw,noctave)
  
  scalevector <- 2^(0:(noctave*nvoice)/nvoice)*s0
  
  cat("# calculating critical values by means of MC-simulations \n")
  
  s0unit=s0/dt*w0/(2*pi)                
  s0log=as.integer((log2(s0unit)-1)*nvoice+1.5)
  
  wsp0 <- rawWSP(ts,s0,noctave,nvoice,w0,swabs,tw,scalevector)
  
  S <- dim(wsp0)[2]
  
  n1 <- 1 + 2*as.integer( sqrt(2) * 2^((S-swabs+s0log-1)/nvoice+1) )
  n2 <- max(scalevector)*tw*2/dt*1.1
  n3 <- 2*tw*s0*2^noctave/dt+100
  n <- max(n1,n2,n3)
  
  center <- (n-1)/2+1
  cv <- matrix(rep(0,S*length(siglevel)),ncol=S)
  rmatrix <- matrix(0,nrow=nreal,ncol=S)
  
  # Fitting AR1-process (red noise) to data
  arts0 <- ar(ts,order.max=1)
  sdts0 <- sd(ts[1:length(ts)]) 
  
  if (arts0$order==0){
    se <- sqrt(arts0$var)
    arts0$ar <- 0.000001
  }
  else
    se <- sqrt(sdts0*sdts0*(1-arts0$ar*arts0$ar))
  
  tsMC <- ts(data=rep(0,n),start=t[1],frequency=1/dt)
  
  # MC Simulations
  for (i in 1:nreal){
    
    tsMC[1:n] <- arima.sim(list(ar = arts0$ar), n+1, sd = se)[2:(n+1)]
    
    rmatrix[i,] <- rawWSP(tsMC,s0,noctave,nvoice,w0,swabs,tw,scalevector)[center,]
    
  }
  
  for (s in (1+swabs):(S-swabs)) rmatrix[,s] <- sort(rmatrix[,s])
  
  for (i in 1:length(siglevel)){ 
    sigindex <- as.integer(nreal*siglevel[i])
    cvv <- rmatrix[sigindex,]
    cvv[is.na(cvv)] <- 0
    cv[i,] <- cvv
  }
  
  cv
  
}

###########

criticalvaluesWCO <- function(s0,noctave,nvoice,w0,swabs,tw,siglevel=0.95){
  
  cv=rep(0,length(siglevel))
  
  for (i in 1:length(siglevel)){
    
    if (siglevel[i]!=0.9 && siglevel[i]!=0.95 && siglevel[i]!=0.99) siglevel[i] <- 0.95
    
    if (siglevel[i]==0.9){
      cat("# significance level set to 0.90 \n")
      sw <- 1.0*swabs/nvoice
      cv[i] <- 0.00246872*w0^2*sw + 0.0302419*w0*sw^2 + 0.342056*sw^3 -
        0.000425853*w0^2 - 0.101749*w0*sw - 0.703537*sw^2 +
        0.00816219*w0 + 0.442342*sw + 0.970315
    }
    
    if (siglevel[i]==0.95){
      cat("# significance level set to 0.95 \n")
      sw <- swabs*100.0/3.0/nvoice
      cv[i] <- 0.0000823*w0^3 + 0.0000424*w0^2*sw + 0.0000113*w0*sw^2 +
        0.0000154*sw^3 - 0.0023*w0^2 - 0.00219*w0*sw - 0.000751*sw^2 +
        0.0205*w0 + 0.0127*sw + 0.95
    }
    
    if (siglevel[i]==0.99){
      cat("# significance level set to 0.99 \n")
      sw <- 1.0*swabs/nvoice
      cv[i] <- 0.00102304*w0^2*sw +  0.00745772*w0*sw^2 + 0.230035*sw^3 -
        0.000361565*w0^2 - 0.0502861*w0*sw - 0.440777*sw^2 +
        0.00694998*w0 + 0.29643*sw + 0.972964
    }
    
    if (cv[i]>1) cv[i] <- 1
    
    cat(paste("# significance testing, cv=",cv[i]," \n",sep=""))
    
  }
  
  cv
  
}

#############################
# AREAWISE SIGNIFICANCE TEST
#############################

slide <- function(data,kernellist,s0,noctave,nvoice,cv){
  
  # slides kernel over data set
  #---------------------------- 
  # data:       matrix containing data
  # kernellist: matrix containing kernel
  # s0:         lowest scale
  # noctave:    number of octaves
  # nvoice:     number of voices per octave
  # cv:         critical value, all values higher are set to one
  
  #Time:  (rows) n,i the kernel is to be scaled in this direction 
  #Scale: (cols) m,j
  
  data[data<cv] <- 0
  
  kernel <- kernellist$bitmap
  
  js <- kernellist$is
  
  sm <- s0*2^noctave
  
  dbin <- tobin(data)
  kbin <- tobin(kernel)
  
  dn <- nrow(dbin)
  dm <- ncol(dbin)
  
  kn <- nrow(kbin)
  km <- ncol(kbin)
  
  mark <- matrix(rep(0,dn*dm),nrow=dn)
  
  for (j in 1:(dm-km+1)){
    
    s <- s0*2^((j+js-1)/nvoice)
    kscn <- as.integer(kn*s/sm);
    if (kscn==0) kscn <- 1
    
    ksc <- scaleKernel(kbin,kscn)
    kscm <- km
    
    for (i in 1:(dn-kscn+1)){
      
      subbin <- dbin[i:(i+kscn-1),j:(j+kscm-1)]
      
      if (sum(ksc*subbin)==sum(ksc))
        mark[i:(i+kscn-1),j:(j+kscm-1)] <- mark[i:(i+kscn-1),j:(j+kscm-1)]+ksc
      
    }
    
  }
  
  mark <- tobin(mark)
  
  mark
  
}

arealtest <- function(wt,dt=1,s0=1,noctave=5,nvoice=20,w0=2*pi,swabs=0,tw=0,siglevel,arealsiglevel=0.9,kernel=0,type=0){
  
  slp <- slope(w0,swabs,tw,nvoice,siglevel,arealsiglevel,type)
  
  if (length(kernel)<2){ 
    maxarea <- s0*2^noctave*slp/10*nvoice/dt
    cat(paste("# calculating kernel (maxarea=",maxarea,")\n",sep=""))
    cvkernel <-
      kernelRoot(s0,w0,a=maxarea,noctave,nvoice,swabs,tw,dt)
    
    cat("# building kernel bitmap \n")
    kernel <-
      kernelBitmap(cvkernel,s0,w0,noctave,nvoice,swabs,tw,dt)
    
  }
  
  cat("# performing arealtest \n")
  sl <- slide(wt,kernel,s0,noctave,nvoice,1)  
  
  list(at=sl,kernel=kernel)
  
}

#############################
# PLOTTING
#############################

plotat <- function(t,wt,at,sigplot){
  
  if (length(at)>1){
    linewidth <- 1
    if (sigplot==3) 
      linewidth <- 5
    
    contour(x=t,z=at,levels=0.5,add=TRUE,drawlabels=FALSE,axes=FALSE,lwd=linewidth,col="black")
  }
  
}


plotcv <- function(t,wt,cv){
  
  if (length(dim(cv))==0)
    
    contour(x=t,z=wt,levels=c(cv),drawlabels=FALSE,axes=FALSE,add=TRUE,col="black",lwd=1)
  
  else{
    
    for(i in 1:nrow(cv)){
      
      v <- cv[i,]
      v[v==0] <- 9999
      m <- matrix(rep(v,length(t)),nrow=length(t),byrow=TRUE)
      
      contour(x=t,z=wt/m,levels=1,drawlabels=FALSE,axes=FALSE,add=TRUE,col="black",lwd=1)
      
    }
    
  }
  
}

plotcoi <- function(t,s0,noctave,w0){
  
  tv <- as.vector(t)
  tvl <- tv[tv-tv[1]<(tv[length(tv)]-tv[1])/2]
  tvr <- tv[tv-tv[1]>=(tv[length(tv)]-tv[1])/2]
  
  lines(tvl,log2(((tvl-tv[1])*4*pi/((w0+sqrt(2+w0*w0))*sqrt(2)))/s0)/noctave,col="black")
  lines(tvr,log2(((tv[length(tv)]-tvr)*4*pi/((w0+sqrt(2+w0*w0))*sqrt(2)))/s0)/noctave,col="black")
  
}

plotmarks <- function(t,s0,noctave,markt,marks){
  
  if (markt[1]!=-999){
    
    for (i in 1:length(markt)){
      lines(c(markt[i],markt[i]),c(0,1),lty="dotted")
    }
    
  }
  
  if (marks[1]!=-999){
    
    for (i in 1:length(marks)){
      lines(c(t[1],t[length(t)]),c(log2(marks[i]/s0)/noctave,log2(marks[i]/s0)/noctave),lty="dotted")
    }
    
  }
  
}


#####################
# PLOT.WT
#####################

plot.wt <- function(wt,markt=-999,marks=-999,t1=NULL,t2=NULL,logscale=FALSE,phase=FALSE,units="",plottitle="",device="screen",file="wt",split=FALSE,color=TRUE,pwidth=10,pheight=5,labsc=1.5,labtext="",sigplot=3,xax=NULL,xlab=NULL,yax=NULL,ylab=NULL){
  
  plotwt(wt$modulus,wt$phase,wt$time,wt$s0,wt$noctave,wt$w0,wt$critval,wt$at,markt,marks,t1,t2,logscale,phase,units,plottitle,device,file,split,color,pwidth,pheight,labsc,labtext,sigplot,xax,xlab,yax,ylab)
  
}

plotwt <-
  function(wt,phs,t,s0,noctave,w0,cv=NULL,at=NULL,markt=-999,marks=-999,t1=NULL,t2=NULL,logscale=FALSE,phase=FALSE,units="",plottitle="Wavelet Plot",device="screen",file="wavelet",split=FALSE,color=TRUE,pwidth=10,pheight=7,labsc=1,labtext="",sigplot=1,xax=NULL,xlab=NULL,yax=NULL,ylab=NULL){
    
    if (is.null(phs)) phase <- FALSE
    
    mgpv <- c(3,1,0)
    if (labsc>1) mgpv[1] <- 3-(labsc-1.5) 
    
    ncolors <- 64
    colors <- wtcolors(ncolors)
    if (color==FALSE) colors <- gray((0:ncolors)/ncolors*0.7+0.3)
    
    rangev <- (0:(ncolors-1))/(ncolors-1)
    rangebar <- matrix(rangev,nrow=2,ncol=64,byrow=TRUE) 
    
    s <- 2^(0:(noctave))*s0
    sn <- (0:(noctave))/noctave
    
    deltat <- (t[length(t)]-t[1])/(length(t)-1)
    cut <- FALSE
    if ((!is.null(t1)) | (!is.null(t2))){
      if (t1<t2 & t1>=t[1] & t2<=t[length(t)]){
        
        cut <- TRUE
        
        i1 <- (t1-t[1])/deltat+1
        i2 <- (t2-t[1])/deltat+1
        
        t <- t[t>=t1 & t<=t2]
        
        wt <- wt[i1:i2,]
        if (phase) phs <- phs[i1:i2,]
        if (length(at)>1) at <- at[i1:i2,]
        
      }
    }
    
    #----------------
    # Setting layout
    #----------------
    
    if (device=="ps" && split==FALSE){
      
      file <- paste(file,".eps",sep="")
      
      postscript(file,onefile=TRUE,horizontal=TRUE,paper="special",width=pwidth,height=pheight)
      cat(paste("# writing",file," \n"))
      
    }
    
    if (device=="ps" && split==TRUE){
      
      file1 <- paste(file,".mod.eps",sep="")
      
      postscript(file1,onefile=TRUE,horizontal=TRUE,paper="special",width=pwidth,height=pheight)
      cat(paste("# writing",file1," \n"))
      
    }
    
    if (phase==TRUE && split==FALSE) nlo <- layout(matrix(c(1,2,3,4),2,2,byrow=TRUE),widths=c(4,1))
    else nlo <- layout(matrix(c(1,2),ncol=2,byrow=TRUE),widths=c(4,1))
    
    
    #-----------------
    # Plotting modulus
    #-----------------
    
    if (logscale){
      if (units=="")
        image(x=t,z=log10(wt),col = colors,axes=FALSE,xlab="Time",ylab="Scale",frame.plot=TRUE,cex.lab=labsc,mgp=mgpv)
      else
        image(x=t,z=log10(wt),col = colors,axes=FALSE,xlab=paste("Time ","[",units,"]",sep=""),ylab=paste("Scale ","[",units,"]",sep=""),frame.plot=TRUE,cex.lab=labsc,mgp=mgpv)
    }
    else{
      if (units=="")
        image(x=t,z=wt,col = colors,axes=FALSE,xlab="Time",ylab="Scale",frame.plot=TRUE,cex.lab=labsc,mgp=mgpv)
      else
        image(x=t,z=wt,col = colors,axes=FALSE,xlab=paste("Time ","[",units,"]",sep=""),ylab=paste("Scale ","[",units,"]",sep=""),frame.plot=TRUE,cex.lab=labsc,mgp=mgpv)
    }
    
    text(t[1+as.integer(length(t)*0.1)],0.8,labtext,cex=2)
    
    if (sigplot==1 | sigplot==3){
      if (is.null(cv)==FALSE){                #Critical values
        if (cv[1]!=0 & cv[1]!=-1) plotcv(t,wt,cv)
      }
    }
    
    if (sigplot>1) plotat(t,wt,at,sigplot)
    
    if (!cut) plotcoi(t,s0,noctave,w0)              #cone of influence
    plotmarks(t,s0,noctave,markt,marks)   #additional guiding lines
    
    if (is.null(xax))
      axis(side=1,at=axTicks(1),cex.axis=labsc)
    else
      if (is.null(xlab))
        axis(side=1,xax,labels=as.character(xax),cex.axis=labsc)
    else
      axis(side=1,xax,labels=xlab,cex.axis=labsc)
    
    if (is.null(yax))
      axis(side=2,sn,labels=as.character(s),cex.axis=labsc)
    else
      if (is.null(ylab))
        axis(side=2,yax,labels=as.character(yax),cex.axis=labsc)
    else
      axis(side=2,yax,labels=ylab,cex.axis=labsc)
    
    title(main=plottitle,cex=labsc)
    
    image(z=rangebar,axes=FALSE,col=colors,frame.plot=TRUE,cex.lab=labsc,mgp=mgpv)
    
    if (is.null(cv)==FALSE){
      if (length(dim(cv))==0){
        for (i in 1:length(cv))
          if (cv[i]!=0) lines(c(-1,2),c(cv[i],cv[i]))
      }
    }
    
    if (!logscale)
      axis(side=2,(0:5)/5,labels=c("0","","","","","1"),cex.axis=labsc)
    else{
      labelv <- substr(as.character(c(0:5)*(max(log10(wt),na.rm=TRUE)-min(log10(wt),na.rm=TRUE))/5),1,4)
      axis(side=2,(0:5)/5,labels=labelv,cex.axis=labsc)
    }
    
    
    #-----------------
    # Plotting phase
    #-----------------  
    if (phase==TRUE){
      
      if (device=="ps" && split==TRUE){
        
        dev.off()
        
        file2 <- paste(file,".phs.eps",sep="")
        
        postscript(file2,onefile=TRUE,horizontal=TRUE,paper="special",width=10,height=5)
        cat(paste("# writing",file2," \n"))
        
      }
      
      colors <- rainbow(ncolors)
      if (color==FALSE){
        dummy <- gray((0:ncolors)/ncolors)
        colors[1:(ncolors/2)] <- dummy[(ncolors/2+1):ncolors]
        colors[(ncolors/2+1):ncolors] <- dummy[1:(ncolors/2)]
      }
      
      if (units=="")
        image(x=t,z=phs,col=colors,axes=FALSE,xlab="Time",ylab="Scale",frame.plot=TRUE,cex.lab=labsc,mgp=mgpv)    
      else
        image(x=t,z=phs,col=colors,axes=FALSE,xlab=paste("Time ","[",units,"]",sep=""),ylab=paste("Scale ","[",units,"]",sep=""),frame.plot=TRUE,cex.lab=labsc,mgp=mgpv)
      
      if (is.null(cv)==FALSE) plotcv(t,wt,cv)
      if (sigplot>1) plotat(t,wt,at,sigplot)
      if (!cut) plotcoi(t,s0,noctave,w0)
      plotmarks(t,s0,noctave,markt,marks)
      
      if (is.null(xax))
        axis(side=1,at=axTicks(1),cex.axis=labsc)
      else
        if (is.null(xlab))
          axis(side=1,xax,labels=as.character(xax),cex.axis=labsc)
      else
        axis(side=1,xax,labels=xlab,cex.axis=labsc)
      
      if (is.null(yax))
        axis(side=2,sn,labels=as.character(s),cex.axis=labsc)
      else
        if (is.null(ylab))
          axis(side=2,yax,labels=as.character(yax),cex.axis=labsc)
      else
        axis(side=2,yax,labels=ylab,cex.axis=labsc)
      
      
      title(main="Phase")
      
      image(z=rangebar,axes=FALSE,col=colors,frame.plot=TRUE,cex.lab=labsc,mgp=mgpv)
      axis(side=2,(0:4)/4,labels=c("-PI","","0","","PI"),cex.axis=labsc)
      
    }
    
    if (device=="ps") dev.off()
    
  }

##############################
# Surrogates
##############################

createwavesurrogates <- function(nsurr=1,wt=1,n,dt=1,s0=1,noctave=5,nvoice=10,w0=2*pi){
  
  surrmat <- matrix(rep(0,n*nsurr),ncol=nsurr)
  
  for (i in 1:nsurr){
    
    cat(paste("# Creating realization ",i,"\n",sep=""))
    
    x <- rnorm(n)
    xts <- ts(x,start=0,deltat=dt)
    
    xwt <- cwt.ts(xts,s0,noctave,nvoice,w0)
    wtsurr <- wt*xwt
    
    surri <- as.vector(invmorlet(wtsurr,0,dt,s0,noctave,nvoice,w0))
    
    surrmat[,i] <- Re(surri)
    
  }
  
  surrmat
  
}

surrwave <- function(x,...)
  UseMethod("surrwave")

surrwave.wt <- function(wt,nsurr=1,spec=TRUE){
  
  n <- length(wt$time)
  t0 <- wt$time[1]
  dt <- (wt$time[13]-t0)/12
  s0 <- wt$s0
  noctave <- wt$noctave
  nvoice <- wt$nvoice
  w0 <- wt$w0
  
  wt <- wt$modulus
  if (spec==TRUE)
    wt <- sqrt(wt)
  
  surrmat <- createwavesurrogates(nsurr,wt,n,dt,s0,noctave,nvoice,w0)
  
  surrts <- ts(surrmat,start=t0,deltat=dt)
  
  surrts
  
}

surrwave.matrix <- function(mat,nsurr=1,t0=0,dt=1,s0=1,noctave=5,nvoice=10,w0=2*pi,sw=0,tw=0,swabs=0,spec=FALSE){
  
  if (sw!=0 & swabs==0)
    swabs <- as.integer(sw*nvoice)
  
  scalevector <- 2^(0:(noctave*nvoice)/nvoice)*s0
  
  if ((noctave*nvoice+1)!=dim(mat)[2])
    cat("# ERROR! nscales unequal noctave*nvoice+1 ! \n")
  
  n <- dim(mat)[1]
  
  if (spec==FALSE)
    mat <- Mod(mat)
  else
    mat <- sqrt(Mod(mat))
  
  wtsm <- smooth.matrix(mat,swabs)
  wt <- smooth.time(wtsm,tw,dt,scalevector)
  
  surrmat <- createwavesurrogates(nsurr,wt,n,dt,s0,noctave,nvoice,w0)
  
  surrts <- ts(surrmat,start=t0,deltat=dt)
  
  surrts
  
}

surrwave.character <-
  function(file,nsurr=1,t0=0,dt=1,s0=1,noctave=5,nvoice=10,w0=2*pi,sw=0,tw=0,swabs=0,transpose=TRUE,spec=FALSE){
    
    if (sw!=0 & swabs==0)
      swabs <- as.integer(sw*nvoice)
    
    scalevector <- 2^(0:(noctave*nvoice)/nvoice)*s0
    
    if (transpose==FALSE)
      mat <- matrix(scan(file,comment.char="#"),ncol=nvoice*noctave+1,byrow=TRUE)    
    else 
      mat <- matrix(scan(file,comment.char="#"),ncol=nvoice*noctave+1,byrow=FALSE)
    
    if ((noctave*nvoice+1)!=dim(mat)[2])
      cat("# ERROR! nscales unequal noctave*nvoice+1 ! \n")
    
    n <- dim(mat)[1]
    
    if (spec==FALSE)
      mat <- Mod(mat)
    else
      mat <- sqrt(Mod(mat))
    
    wtsm <- smooth.matrix(mat,swabs)
    wt <- smooth.time(wtsm,tw,dt,scalevector)
    
    surrmat <- createwavesurrogates(nsurr,wt,n,dt,s0,noctave,nvoice,w0)
    
    surrts <- ts(surrmat,start=t0,deltat=dt)
    
    surrts
    
  }

surrwave.ts <- function(ts,nsurr=1,s0=1,noctave=5,nvoice=10,w0=2*pi,sw=0,tw=0,swabs=0){
  
  n <- length(ts)
  t0 <- time(ts)[1]
  dt <- deltat(ts)
  if (sw!=0 & swabs==0)
    swabs <- as.integer(sw*nvoice)
  
  scalevector <- 2^(0:(noctave*nvoice)/nvoice)*s0
  
  wt <- Mod(cwt.ts(ts,s0,noctave,nvoice,w0))
  
  wtsm <- smooth.matrix(wt,swabs)
  wt <- smooth.time(wtsm,tw,dt,scalevector)
  
  surrmat <- createwavesurrogates(nsurr,wt,n,dt,s0,noctave,nvoice,w0)
  
  surrts <- ts(surrmat,start=t0,deltat=dt)
  
  surrts
  
}

invmorlet <- function(wt,t0=0,dt=1,s0=1,noctave=5,nvoice=10,w0=2*pi){
  
  if ((noctave*nvoice+1)!=dim(wt)[2])
    cat("# ERROR! nscales unequal noctave*nvoice+1 ! \n")
  
  n <- dim(wt)[1]
  
  wt[is.na(wt)] <- 0
  
  tsRe <- rep(0,n)
  tsIm <- rep(0,n)
  
  wtRe <- t(Re(wt))
  wtIm <- t(Im(wt))
  
  z <- .C("invmorlet",
          as.double(as.vector(wtRe)),
          as.double(as.vector(wtIm)),
          as.integer(n),
          as.double(dt),
          as.double(s0),
          as.integer(noctave),
          as.integer(nvoice),     
          as.double(w0),
          tsRetmp = as.double(tsRe),
          tsImtmp = as.double(tsIm),
          PACKAGE="sowas")
  
  invvec=complex(real=z$tsRetmp,imaginary=z$tsImtmp)
  invts <- ts(data=invvec,start=t0,deltat=dt)
  
  invts
  
}

#################################
# INPUT / OUTPUT
#################################

readmatrix <- function(file,M){
  
  A <- matrix(scan(file,comment.char="#"),ncol=M,byrow=TRUE)
  
  A
  
}


readts <- function(file){
  
  A <- matrix(scan(file,comment.char="#"),ncol=2,byrow=TRUE)
  
  Adum <- A
  
  Adum[is.na(Adum)] <- 0
  
  t <- Adum %*% c(1,0)
  x <- A %*% c(0,1)
  
  N=length(t)
  f=1/(t[13]-t[1])*12
  
  if ((f>11) && (f<13)) f <- 12
  
  timeseries<-ts(data=x,start=t[1],frequency=f)
  
  timeseries
  
}

writematrix <- function(file,data,header="# R Matrix"){
  
  write(header,file)
  write(t(data),file,ncol(data),append=TRUE)
  
}

############################
# HELP FUNCTIONS
############################

smooth.matrix <- function(wt,swabs){
  
  if (swabs != 0)
    smwt <- t(filter(t(wt),rep(1,2*swabs+1)/(2*swabs+1)))
  else
    smwt <- wt
  
  smwt
  
}

smooth.time <- function(wt,tw,dt,scalevector){
  
  smwt <- wt
  
  if (tw != 0){
    for (i in 1:length(scalevector)){
      
      twi <- as.integer(scalevector[i]*tw/dt)
      smwt[,i] <- filter(wt[,i],rep(1,2*twi+1)/(2*twi+1))
      
    }
  }
  
  smwt
  
}


adjust.noctave <- function(N,dt,s0,tw,noctave){
  
  if (tw>0){  
    dumno <- as.integer((log(N*dt)-log(2*tw*s0))/log(2))
    if (dumno<noctave){
      cat("# noctave adjusted to time smoothing window \n")
      noctave <- dumno
    }
  }
  
  noctave
  
}

adjust.s0 <- function(s0,dt){
  
  if (s0<2*dt){
    s0 <- 2*dt
    cat(paste("# s0 set to ",s0," \n"))
  }
  
  s0
  
}

adjust.timeseries <- function(ts1,ts2){
  
  if (length(ts1)!=length(ts2)){
    tsint <- ts.intersect(ts1,ts2)
    dt <- deltat(ts1)
    ts1 <- ts(data=tsint[,1],start=time(tsint)[1],frequency=1/dt)
    ts2 <- ts(data=tsint[,2],start=time(tsint)[1],frequency=1/dt)
    t <- time(ts1)
  }
  
  list(ts1=ts1,ts2=ts2)
  
}

checkarealsiglevel <- function(sw,tw,w0,arealsiglevel,siglevel,type){
  
  if (type==0){
    
    swv <- c(0,0.5,1)
    twv <- c(0,1.5,3)
    w0v <- c(pi,2*pi)
    
    if (length(swv[swv==sw])==0 || length(twv[twv==tw])==0 ||
          length(w0v[w0v==w0])==0){
      arealsiglevel <- 0
      cat("# areawise test for spectrum currently \n")
      cat("# only possible for \n")
      cat("# sw = 0 \n") 
      cat("# tw = 0 \n") 
      cat("# w0 = 2pi \n") 
      cat("# No areawise test performed \n") 
    }
  }
  
  if (type==1){
    
    swv <- c(0.5)
    twv <- c(1.5)
    w0v <- c(2*pi)
    
    if (length(swv[swv==sw])==0 || length(twv[twv==tw])==0 ||
          length(w0v[w0v==w0])==0){
      arealsiglevel <- 0
      cat("# areawise test for coherence currently \n")
      cat("# only possible for \n")
      cat("# sw = 0.5 \n") 
      cat("# tw = 1.5 \n") 
      cat("# w0 = 2pi \n") 
      cat("# No areawise test performed \n") 
    }
  }
  
  
  if (arealsiglevel!=0){
    arealsiglevel <- 0.9
    siglevel <- 0.95
    cat("# currently only siglevel=0.95 and arealsiglevel=0.9 possible for areawise test \n")
  }
  
  list(siglevel=siglevel,arealsiglevel=arealsiglevel)
  
}

########################

as.wt <- function(modulus,phase=NULL,s0=NULL,noctave=NULL,nvoice=NULL,w0=NULL,dt=1,time=NULL,scales=NULL,critval=NULL,at=NULL,kernel=NULL,N=NULL,t0=NULL){
  
  if (is.null(scales))
    gotscales <- FALSE
  else
    gotscales <- TRUE
  
  if ((!gotscales) & (!is.null(s0)) & (!is.null(noctave)) &(!is.null(nvoice))){
    gotscales <- TRUE
    scales=2^(0:(noctave*nvoice)/nvoice)*s0
  }
  
  if (gotscales & (is.null(s0) | is.null(noctave) |
                     is.null(nvoice))){
    s0 <- scales[1]
    noctave <- log(scales[length(scales)]/s0)/log(2)
    nvoice <- (length(scales)-1)/noctave
  }
  
  
  if (!gotscales)
    cat("# ERROR! No scales given! \n")
  
  if (is.null(time)) gottimes <- FALSE
  else gottimes <- TRUE
  
  if ((!gottimes) & (!is.null(dt)) & (!is.null(t0)) &(!is.null(N))){
    gottimes <- TRUE
    time=(0:(N-1))*dt+t0
  }
  
  if (!gottimes)
    cat("# ERROR! No time vector given! \n")
  
  wcolist <- list(modulus=modulus,phase=phase,s0=s0,noctave=noctave,nvoice=nvoice,w0=w0,time=time,scales=scales,critval=critval,at=at,kernel=kernel)
  
  class(wcolist) <- "wt"
  
  wcolist
  
}

########################

wtcolors <- function(ncolors){
  
  upside <- rainbow(ncolors,start=0,end=.7)
  #upside <- heat.colors(ncolors+5)
  #upside <- upside[1:ncolors]
  
  
  down <- upside
  
  for (i in 1:ncolors){
    down[i] <- upside[ncolors+1-i]    
  }#
  
  down
  
}

####################

createwgn <- function(N,sig,dt){
  
  timeseries<-ts(data=rnorm(N,0,sig),start=0,deltat=dt)
  
  timeseries
  
}


createar <- function(N,a,sig,dt){
  
  if (a==0) a <- 0.000000001
  
  se <- sqrt(sig*sig*(1-a*a))
  tsMC <- ts(data=rep(0,N),start=0,deltat=dt)
  
  tsMC[1:N] <- arima.sim(list(ar = a), N, sd = se)
  
}

######################

rk <- function(N=1000,s=8,noctave=5,nvoice=10,w0=2*pi,plot=TRUE){
  
  t <- 1:N
  
  sunit <- s*(w0+sqrt(2+w0*w0))/(4*pi)
  
  s0 <- 4
  #s0unit <- s0*(w0+sqrt(2+w0*w0))/(4*pi)
  s0unit=s0/dt*w0/(2*pi)                  #(CORRECT) 
  s0log <- as.integer((log2(s0unit)-1)*nvoice+1.5)
  
  totnoct <- noctave+as.integer(s0log/nvoice)+1
  
  x <- morlet(N,N/2,sunit,w0)
  
  totts.cwt <- Mod(cwt(x,totnoct,nvoice,w0,plot=0))
  wt=totts.cwt[,s0log:(s0log+noctave*nvoice)]
  
  wt <- wt/max(wt)
  
  if (plot==TRUE) plotwt(wt,0,t,s0,noctave,w0,units="",plottitle="Reproducing Kernel")  
  
  wt
  
}

###################

addvalues <- function(nvoice,dnoctave,x,value){
  
  nr <- dim(x)[1]
  nc <- dim(x)[2]
  dnc <- nvoice*dnoctave
  
  y <- matrix(rep(value,nr*(nc+dnc)),nrow=nr)
  
  y[,(dnc+1):(dnc+nc)] <- x
  
  y
  
}

####################

scalematrix <- function(wt){
  
  # scales matrix, such that the maximal value is one
  # wt: matrix to be scaled
  
  mval <- max(wt,na.rm=TRUE)
  
  wt <- wt/mval
  
  wt
  
}


foldKernel <- function(F,swabs,tw,s,dt){
  
  # folds a matrix (e.g. kernel with smoothing window
  # F:  matrix input
  # swabs: smooting window width
  
  smsF <- smooth.matrix(F,swabs)
  smsF[is.na(smsF)] <- 0
  
  smtsF <- smooth.time(smsF,tw,dt,s)
  
  smtsF[is.na(smtsF)] <- 0
  
  scF <- scalematrix(smtsF)
  
  scF
  
}

kernelBitmap <- function(c,s0=1,w0=6,noctave=6,nvoice=20,swabs=0,tw=0,dt=1){
  # produces kernel bitmap
  # c:       height of contour, that defines area
  # s0:      lowest scale
  # noctave: number of octaves
  # nvoice:  number of voices per octave
  # swabs:      smoothing window length in scale direction
  # dt:      sampling time
  
  s <- s0*2^noctave
  is <- noctave*nvoice
  
  x <- s0*2^(((1:(nvoice*(noctave+2)))-1)/nvoice)
  
  t <- max(2000,max(x)*tw*2/dt*1.1)
  
  y <- ((0:(2*t))-t)/2*dt
  
  X <- matrix(x,ncol=nvoice*(noctave+2),nrow=2*t+1,byrow=T)
  Y <- matrix(y,ncol=nvoice*(noctave+2),nrow=2*t+1)
  
  F <- sqrt(2*s*X/(s*s+X*X))*exp(-0.5*(Y*Y+w0*w0*(X-s)*(X-s))/(s*s+X*X))
  
  F <- foldKernel(F,swabs,tw,x,dt)
  
  F[F<c] <- 0
  F[F>=c] <- 1
  
  is1 <- 1
  is2 <- nvoice*(noctave+1)
  it1 <- 1
  it2 <- 2*t+1
  
  L <- F[1:(2*t+1),is1]
  
  while (length(L[L!=0])==0) {
    is1 <- is1+1
    L <- F[1:(2*t+1),is1]
  }
  
  L <- F[1:(2*t+1),is2]
  
  while (length(L[L!=0])==0) {
    is2 <- is2-1
    L <- F[1:(2*t+1),is2]
  }
  
  
  L <- F[it1,1:(nvoice*(noctave+2))]
  
  while (length(L[L!=0])==0) {
    it1 <- it1+1
    L <- F[it1,1:(nvoice*(noctave+2))]
  }
  
  L <- F[it2,1:(nvoice*(noctave+2))]
  
  while (length(L[L!=0])==0) {
    it2 <- it2-1
    L <- F[it2,1:(nvoice*(noctave+2))]
  }
  
  kernel <- list(bitmap=F[(it1-1):(it2+1),(is1-1):(is2+1)],is=is-is1) 
  
  kernel
  
}

kernelRoot <- function(s0=1,w0=6,a=0,noctave=6,nvoice=20,swabs=0,tw=0,dt=1){
  
  tol <- 0.005
  cmin <- 0
  cmax <- 1
  cntr <- 0.5
  da <- a
  
  while (abs(da/a)>tol){
    
    da <- kernelArea(cntr,s0,w0,a,noctave,nvoice,swabs,tw,dt)
    if (da>0){
      cmin <- cntr
      cntr <- mean(c(cntr,cmax))
    }
    if (da<0){
      cmax <- cntr
      cntr <- mean(c(cntr,cmin))      
    }
  }
  
  cntr
  
}

kernelArea <- function(cntr,s0=1,w0=6,a=0,noctave=6,nvoice=20,swabs=0,tw=0,dt=1){
  
  # calulates area of reproducing kernel for smoothed data at scale s0*2^noctave
  # cntr:       height of contour line to define kernel area. This
  #          parameter is to be estimated!
  # s0:      lowest scale
  # w0:      parameter of Morlet Wavelet
  # a:       area offset (only needed, when finding root. Is set to
  #          desired area 
  # noctave: number of octaves
  # nvoice:  number of voices per octave
  # swabs:      smoothing window width in scale direction
  # dt:      sampling time
  
  s <- s0*2^noctave
  
  x <- s0*2^(((1:(nvoice*(noctave+2)))-1)/nvoice)
  
  t <- max(2000,max(x)*tw*2/dt*1.1)
  
  y <- ((0:(2*t))-t)/2*dt
  
  X <- matrix(x,ncol=nvoice*(noctave+2),nrow=2*t+1,byrow=T)
  Y <- matrix(y,ncol=nvoice*(noctave+2),nrow=2*t+1)
  
  F <- sqrt(2*s*X/(s*s+X*X))*exp(-0.5*(Y*Y+w0*w0*(X-s)*(X-s))/(s*s+X*X))
  
  F <- foldKernel(F,swabs,tw,x,dt)
  
  F[F>=cntr] <- 1
  F[F<cntr] <- 0
  
  area <- length(F[F==1])-a
  
  area
  
}


tobin <- function(x){
  
  # sets nonzero values to one
  
  y <- x/x
  y[is.na(y)] <- 0
  
  y
  
}

scaleKernel <- function(kernel,l){
  
  # scales kernel length in time direction proportional to scale
  # kernel: data bitmap of width n
  # l:      new width of kernel
  
  n <- nrow(kernel)
  m <- ncol(kernel)
  
  newKernel <- matrix(rep(0,m*l),nrow=l)
  
  d <- as.double(n)/as.double(l)
  
  for (i in 1:l){
    j <- as.integer((i-0.5)*d)
    if (j==0) j <- 1
    newKernel[i,1:m] <- kernel[j,1:m]
  }
  
  newKernel
  
}

slope <- function(w0,swabs,tw,nvoice,siglevel,arealsiglevel,type){
  
  sw <- swabs/nvoice
  
  if (type==0){ # wavelet spectrum
    
    if (tw == 0    & sw == 0 	 & w0 == 1 *pi) slp <-  5.82518 	 # w = 18.35831 
    if (tw == 1.5 	 & sw == 0 	 & w0 == 1 *pi) slp <-  24.69852 	 # w = 14.30709 
    if (tw == 3 	 & sw == 0 	 & w0 == 1 *pi) slp <-  35.48368 	 # w = 14.72354 
    if (tw == 0 	 & sw == 5 	 & w0 == 1 *pi) slp <-  7.347707 	 # w = 17.96942 
    if (tw == 1.5 	 & sw == 5 	 & w0 == 1 *pi) slp <-  28.24291 	 # w = 12.65993 
    if (tw == 3 	 & sw == 5 	 & w0 == 1 *pi) slp <-  51.13723 	 # w = 10.96359 
    if (tw == 0 	 & sw == 10 	 & w0 == 1 *pi) slp <-  10.47856 	 # w = 15.5941 
    if (tw == 1.5 	 & sw == 10 	 & w0 == 1 *pi) slp <-  45.07387 	 # w = 15.29793 
    if (tw == 3 	 & sw == 10 	 & w0 == 1 *pi) slp <-  52.82886 	 # w = 12.72361 
    
    if (tw == 0 	 & sw == 0 	 & w0 == 2 *pi) slp <-  8.718912 	 # w = 17.75933 
    if (tw == 1.5 	 & sw == 0 	 & w0 == 2 *pi) slp <-  11.88006 	 # w = 15.39648 
    if (tw == 3 	 & sw == 0 	 & w0 == 2 *pi) slp <-  26.55977 	 # w = 1.064384 
    if (tw == 0 	 & sw == 5 	 & w0 == 2 *pi) slp <-  14.64761 	 # w = 16.27518 
    if (tw == 1.5 	 & sw == 5 	 & w0 == 2 *pi) slp <-  28.27798 	 # w = 14.57059 
    if (tw == 3 	 & sw == 5 	 & w0 == 2 *pi) slp <-  63.54121 	 # w = 12.83778 
    if (tw == 0 	 & sw == 10 	 & w0 == 2 *pi) slp <-  27.78735 	 # w = 11.95813 
    if (tw == 1.5 	 & sw == 10 	 & w0 == 2 *pi) slp <-  41.27260 	 # w = 12.03379 
    if (tw == 3 	 & sw == 10 	 & w0 == 2 *pi) slp <-  67.37015 	 # w = 10.63935 
    
  }
  
  if (type==1){ # wavelet coherence
    
    if (tw==0   & sw==0   & w0==pi) slp <- 999 #not valid
    if (tw==1.5 & sw==0   & w0==pi) slp <- 1
    if (tw==3   & sw==0   & w0==pi) slp <- 1
    if (tw==0   & sw==0.5 & w0==pi) slp <- 1
    if (tw==1.5 & sw==0.5 & w0==pi) slp <- 1
    if (tw==3   & sw==0.5 & w0==pi) slp <- 1
    if (tw==0   & sw==1   & w0==pi) slp <- 1
    if (tw==1.5 & sw==1   & w0==pi) slp <- 1
    if (tw==3   & sw==1   & w0==pi) slp <- 1
    
    if (tw==0   & sw==0   & w0==2*pi) slp <- 999 #not valid
    if (tw==1.5 & sw==0   & w0==2*pi) slp <- 1
    if (tw==3   & sw==0   & w0==2*pi) slp <- 1
    if (tw==0   & sw==0.5 & w0==2*pi) slp <- 1
    if (tw==1.5 & sw==0.5 & w0==2*pi) slp <- 8.3
    if (tw==3   & sw==0.5 & w0==2*pi) slp <- 1
    if (tw==0   & sw==1   & w0==2*pi) slp <- 1
    if (tw==1.5 & sw==1   & w0==2*pi) slp <- 1
    if (tw==3   & sw==1   & w0==2*pi) slp <- 1
    
    if (tw==0   & sw==0   & w0==3*pi) slp <- 999 #not valid
    if (tw==1.5 & sw==0   & w0==3*pi) slp <- 1
    if (tw==3   & sw==0   & w0==3*pi) slp <- 1
    if (tw==0   & sw==0.5 & w0==3*pi) slp <- 1
    if (tw==1.5 & sw==0.5 & w0==3*pi) slp <- 1
    if (tw==3   & sw==0.5 & w0==3*pi) slp <- 1
    if (tw==0   & sw==1   & w0==3*pi) slp <- 1
    if (tw==1.5 & sw==1   & w0==3*pi) slp <- 1
    if (tw==3   & sw==1   & w0==3*pi) slp <- 1
    
    if (tw==0   & sw==0   & w0==4*pi) slp <- 999 #not valid
    if (tw==1.5 & sw==0   & w0==4*pi) slp <- 1
    if (tw==3   & sw==0   & w0==4*pi) slp <- 1
    if (tw==0   & sw==0.5 & w0==4*pi) slp <- 1
    if (tw==1.5 & sw==0.5 & w0==4*pi) slp <- 1
    if (tw==3   & sw==0.5 & w0==4*pi) slp <- 1
    if (tw==0   & sw==1   & w0==4*pi) slp <- 1
    if (tw==1.5 & sw==1   & w0==4*pi) slp <- 1
    if (tw==3   & sw==1   & w0==4*pi) slp <- 1
    
  }
  
  cat(paste("# slope ",slp,"\n",sep=""))
  
  slp
  
}


####################################################################
##
## File:        aux.r
##
## Description: Miscelaneous functions for clustering with kcca
##
## Modified:    june 2010
##
####################################################################


  #######################################################

  # Transforms a matrix of data (one observation by row)
  #     into an array where position[ , , i] gives
  #     the smoothed modulus of the i-th cwt observation 

  ########################################################


  toCWT  <- function(X, sw=  0,  tw=  0, swabs= 0,
                       nvoice= 12, noctave= 5, 
                       s0= 2, w0= 2*pi, lt= 24, dt= 0.5,
                       spectra = FALSE, smooth = TRUE,
                       scaled  = FALSE,
                     scalevector)
     { noctave  <- adjust.noctave(lt, dt, s0, tw, noctave)
       if(missing(scalevector)) 
          scalevector  <- 2^(0:(noctave * nvoice) / nvoice) * s0
       res <- lapply(1:nrow(X), function(n)
           { tsX         <- ts( X[n,] )
             tsCent      <- tsX - mean(tsX)
             if(scaled)  tsCent <- ts(scale(tsCent))           
             tsCent.cwt  <- cwt.ts(tsCent, s0, noctave, nvoice, w0)
             tsCent.cwt
           } )
	   if( spectra ) res <- lapply(res, function(l) Mod(l)^2 )
	   if( smooth  ) res <- lapply(res, smCWT, swabs = swabs,
	                               tw = tw, dt = dt, 
	                               scalevector = scalevector)
       resArray <- array(NA, c(nrow(res[[1]]), ncol(res[[1]]),
                               length(res)))
       for( l in 1:length(res) ) resArray[ , , l] <- res[[l]]
       resArray
     }


  # ===============================================================

  smCWT <- function(CWT, sw=  0,  tw=  0, swabs= 0,
                       nvoice= 12, noctave= 2, s0= 2, w0= 2*pi, 
					   lt= 24, dt= 0.5, scalevector )
		 {
#         noctave  <- adjust.noctave(lt, dt, s0, tw, noctave)
#         scalevector  <- 2^(0:(noctave * nvoice) / nvoice) * s0
         wsp     <- Mod(CWT)  
         smwsp   <- smooth.matrix(wsp, swabs)
         smsmwsp <- smooth.time(smwsp, tw, dt, scalevector)
         smsmwsp
       }


  # ===============================================================

  toDWT <- function(x, filter.number = 6, family = "DaubLeAsymm")
{ x2   <- spline(x, n = 2^ceiling( log(length(x), 2) ),
            method = 'natural')$y
  Dx2 <- wd(x2, family = family, filter.number = filter.number)$D
		Dx2
}

  # ===============================================================

  contrib <- function(x) 
      { J   <- log( length(x)+1, 2)
        nrj <- numeric(J)
        t0  <- 1
        t1  <- 0
        for( j in 1:J ) {
          t1     <- t1 + 2^(J-j)
          nrj[j] <- sqrt( sum( x[t0:t1]^2 ) )
          t0     <- t1 + 1
        }
        return(nrj)  
	  }


    # ========================================= distance for coh ===

  coherence <- function( x, y)
       { J <- log(length(x) + 1, 2)
	     t0 <- 1
		 sg2_x <- 0
		 sg2_y <- 0
		 sg_xy <- 0
	     for(j in 0:(J - 1))
		 {  t1 <- t0 + 2^(J - j)/2  - 1
		    tt <- t0:t1
		    sg2_x <- sg2_x + mean(x[t0:t1]^2)
			sg2_y <- sg2_y + mean(y[t0:t1]^2)
			sg_xy <- sg_xy + mean(x[t0:t1] * y[t0:t1])
            t0 <- t1 + 1
		 }
		res <- sg_xy^2 / sg2_x / sg2_y
		res
	   }


  vect2mat <- function(vect){
                 vect <- as.vector(vect)
                 matrix(vect[-(1:2)], delta) #, lscvect)
               }


  # =========================================  # myimg for graphics 
  jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", 
                                "cyan", "#7FFF7F", "yellow", 
                                "#FF7F00", "red", "#7F0000"))

  myimg <- function(MAT, x = 1:nrow(MAT), y = 1:col(MAT), ... )
              filled.contour(  x = x, y = y, z = MAT, 
			       xlab= 'Time', ylab= 'scale',
			       color.palette = jet.colors,
			       ... )

  # =========================================  # WER 
  

werDist <- function(data, ...) {
Xcwt4   <- toCWT(conso, noctave = noctave4, dt = 1,
                 scalevector = scalevector4,
                 lt = delta, smooth = FALSE, 
                 nvoice = nvoice)      # observations node with CWT

Xcwt2 <- matrix(0.0, nrow= n, ncol= 2 + delta * lscvect4)
Xcwt2 <- matrix(NA_complex_, nrow= n, ncol= 2 + length((c(Xcwt4[,,1]))))

for(i in 1:n) 
  Xcwt2[i,] <- c(delta, lscvect4, Xcwt4[,,i] / max(Mod(Xcwt4[,,i])) ) 

rm(Xcwt4)

## _.b WER^2 distances  ########
Xwer_dist    <- matrix(0.0, n, n)
for(i in 1:(n - 1)){
  cat(sprintf('\nIter: , %i', i))
  mat1   <- vect2mat(Xcwt2[i,])
  for(j in (i + 1):n){
    mat2 <- vect2mat(Xcwt2[j,])
    num     <- Mod(mat1 * Conj(mat2))
    WX      <- Mod(mat1 * Conj(mat1))
    WY      <- Mod(mat2 * Conj(mat2))
    smsmnum <- smCWT(num, scalevector = scalevector4)
    smsmWX  <- smCWT(WX,  scalevector = scalevector4)
    smsmWY  <- smCWT(WY,  scalevector = scalevector4)
    wer2    <- sum(colSums(smsmnum)^2)  /
      sum( sum(colSums(smsmWX) * colSums(smsmWY)) )
    Xwer_dist[i, j] <- sqrt(delta * lscvect4 * (1 - wer2))
    Xwer_dist[j, i] <- Xwer_dist[i, j]
  }
}
diag(Xwer_dist) <- numeric(n)

return(as.dist(Xwer_dist))
}
