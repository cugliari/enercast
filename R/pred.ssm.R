#' predict.ssm
#'
#' Prediction in Univariate State Space Model
#' 
#' @param ssm.obj Estimated object from \code{ssm}
#' @param newdata data
#' @param for.each 'hour'
#'
#' @return list
#' @author Jairo Cugliari, Andres Castrillejo, Fernando Massa, Ignacio Ramirez
pred.ssm <- function(ssm.obj.i){
      yt <- ssm.obj.i$yt
      xa.pred <- ssm.obj.i$xa.pred
      xf.pred <- ssm.obj.i$xf.pred

          tita <- ssm.obj.i$optim$par
          if (ssm.obj.i$fk$ct!=0){
                    tita.xf <- tita[seq(ncol(ssm.obj.i$fk$ct))]
                            tita.xa <- tita[-seq(ncol(ssm.obj.i$fk$ct))]
                            ssm.obj.i$fk$ct <- t(ssm.obj.i$fk$ct%*%tita.xf)
                  } else {
                            tita.xa <- tita
                                    tita.xf <- NULL
                          }
          aux <- exp(tita.xa[1:(length(tita.xa)-1)])
          hht <- diag(length(aux))*aux
          ggt <- matrix(exp(tita.xa[length(tita.xa)]),1,1)
          fk1 <- fkf(a0=ssm.obj.i$fk$a0,P0=ssm.obj.i$fk$P0,dt=ssm.obj.i$fk$dt,ct=ssm.obj.i$fk$ct,Tt=ssm.obj.i$fk$Tt,Zt=ssm.obj.i$fk$Zt,yt=yt,HHt=hht,GGt=ggt)
          att <- fk1$att[,ncol(fk1$att)]
          Ptt <- fk1$Ptt[,,ncol(fk1$att)]

          if(class(xa.pred)!='matrix') xa.pred <- matrix(xa.pred,ncol=length(att),byrow=TRUE)

          pred <- xa.pred%*%att
          predcov <- xa.pred%*%Ptt%*%t(xa.pred)+ggt[1,1]
          if(!is.null(xf.pred)){
                    if(class(xf.pred)!='matrix') xf.pred <- as.matrix(xf.pred)
                            pred <- pred+xf.pred%*%tita.xf
                  }
          return(list(predt=pred,covt=predcov))
    }
