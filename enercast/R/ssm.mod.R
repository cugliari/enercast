ssm.mod <- function(formulaXa,formulaXf=NULL,data,ini=NULL){
    n <- nrow(data)
    mf  <- model.frame(formulaXa,data)
    yt <- rbind(mf[,1])
    xa <- model.matrix(formulaXa,mf)

    ct <- matrix(0,1,1) 

    if (!is.null(formulaXf)){
        xf <-  model.matrix(formulaXf,data)
        p.f <- ncol(xf)
        ct <- as.matrix(xf)
    }

    #intercept en las 2 formulas?
    cond <- ('(Intercept)' %in% colnames(ct)) & ('(Intercept)' %in% colnames(xa))
    if(cond)
        ct <- ct[,!(colnames(ct) %in% '(Intercept)'),drop=FALSE] 

    ##    f1 <- update(f1,~ . - 1)

    #    warning('ss')
    p.a <- ncol(xa)
    dt  <-  matrix(rep(0,p.a),p.a,1)
    Tt <- diag(p.a)
    Zt <- array(NA,dim=c(1,p.a,n)); for (i in 1:n) Zt[1,,i] <- as.numeric(xa[i,])
    a0 <- coef(lm(yt~.-1,data=data.frame(yt=t(yt),xa)))
    P0 <- 1000000*diag(p.a)

    f1 <- function(x0,yt,a0,P0,ct,dt,Tt,Zt){
        k <- length(x0)
        hht <- diag(exp(x0[1:(k-1)]),ncol=k-1)
        ggt <- matrix(exp(x0[k]))
        -fkf(HHt=hht,GGt=ggt,yt=yt,a0=a0,P0=P0,dt=dt,ct=ct,Tt=Tt,Zt=Zt)$logLik
    }
    f1b <- function(x0,yt,a0,P0,ct,dt,Tt,Zt){
        k.a <- length(a0)
        k.f <- ncol(ct)
        ct <- t(ct%*%x0[1:k.f])
        hht <- diag(exp(x0[(k.f+1):(k.a+k.f)]))
        ggt <- matrix(exp(x0[k.a+k.f+1]))
        -fkf(HHt=hht,GGt=ggt,yt=yt,a0=a0,P0=P0,dt=dt,ct=ct,Tt=Tt,Zt=Zt)$logLik
    }
    
    if (is.null(ini)) {
        ini <- log(var(yt[1,])*rep(1/(p.a+1),(p.a+1)))
        if (!is.null(formulaXf)) {
            beta <- qr.solve(cbind(1,as.matrix(xf)),as.numeric(yt))
            if (any(apply(xa,2,var)==0)) beta <- beta[-1]
            ini <- c(beta,ini)
        }
    }

    if (is.null(formulaXf))   fkf.optim <- optim(par=ini,fn=f1,yt=yt,a0=a0,P0=P0,ct=ct,dt=dt,Zt=Zt,Tt=Tt,method='BFGS',control=list(trace=100,maxit=200))
    if (!is.null(formulaXf))  fkf.optim <- optim(par=ini,fn=f1b,yt=yt,a0=a0,P0=P0,ct=ct,dt=dt,Zt=Zt,Tt=Tt,method='BFGS',control=list(trace=100,maxit=200))


    return(list(optim=fkf.optim,fk=list(dt=dt,ct=ct,Tt=Tt,Zt=Zt,a0=a0,P0=P0),yt=yt))
}
