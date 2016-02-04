werdist <- function(data, ...) { #  series must be on rows 

n     <- nrow(data)
delta <- ncol(data)

## 2. Compute WER distance matrix ####

## _.a CWT -- Filtering the lowest freqs (>6m) ####
if(missing(nvoice)) nvoice <- 4
if(missing(noctave4)) noctave4 <- adjust.noctave(N = delta, dt = 1, s0 = 2,
                                                tw = 0, noctave = 13)
if(missing(scalevector4)) scalevector4  <- 2^(4:(noctave4 * nvoice) / nvoice) * 2

lscvect4      <- length(scalevector4)
lscvect <- lscvect4  

Xcwt4   <- toCWT(data, noctave = noctave4, dt = 1,
                 scalevector = scalevector4,
                 lt = delta, smooth = FALSE, 
                 nvoice = nvoice)      # observations node with CWT

Xcwt2 <- matrix(0.0, nrow= n, ncol= 2 + delta * lscvect)
Xcwt2 <- matrix(NA_complex_, nrow= n, ncol= 2 + length((c(Xcwt4[,,1]))))

for(i in 1:n) 
  Xcwt2[i,] <- c(delta, lscvect, Xcwt4[,,i] / max(Mod(Xcwt4[,,i])) ) 

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
    Xwer_dist[i, j] <- sqrt(delta * lscvect * (1 - wer2))
    Xwer_dist[j, i] <- Xwer_dist[i, j]
  }
}
diag(Xwer_dist) <- numeric(n)

return(Xwer_dist)
}

