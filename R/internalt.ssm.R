trigon <- function(data, frecH = 1, n = nrow(data), w = 24, q = 3){
  
  is.leapyear <- function(year){
  return(((year %% 4 == 0) & (year %% 100 != 0)) | (year %% 400 == 0))}

  if (!is.numeric(data$hour)) data$hour <- as.numeric(data$hour)
  if (!is.numeric(data$year)) data$year <- as.numeric(data$year)
    
  if (frecH == 1) xtime <- data$hour
  if (frecH == 7) xtime <- rep(1:(w * frecH), length = nrow(data))
  if (frecH == 365) {
    anio <- 365 + is.leapyear(unique(data$year))
    xtime <- unlist(sapply(anio, function(x) seq(1, w * x)))
    }

  RES <- NULL
  for (i in 1:q) {
    s1 <- sin(xtime * 2 * pi * i / (frecH * w))
    c1 <- cos(xtime * 2 * pi * i / (frecH * w))
    
    RES <- cbind(RES, s1, c1)
    
  }
  
  colnames(RES) <- paste(rep(c('sin', 'cos'), q), rep(1:q, each = 2), sep = '_')
  return(RES)
}


rcs <- function(x, m, nod){
  if (missing(nod))
    nod <- quantile(x, probs = seq(0, 1, 1 / (m + 1)))[-c(1, m + 2)]

  if (missing(m)) m <- length(nod)
  
  Cs <- matrix(NA, nrow = length(x), ncol = m - 2)
  
  for (i in 1:(m - 2))
    Cs[,i] <- apply(cbind(0, (x - nod[i])^3), 1, max) - 
    apply(cbind(0, (x - nod[m - 1])^3), 1, max) * (nod[m]     - nod[i]) / (nod[m] - nod[m - 1]) + 
    apply(cbind(0, (x - nod[m])^3),     1, max) * (nod[m - 1] - nod[i]) / (nod[m] - nod[m - 1])
  
  Cs <- cbind(x,Cs)
  colnames(Cs) <- paste('C', 0:(m - 2), sep = '')
  Cs <- sweep(Cs, MARGIN = 2, FUN = "/", 
              STATS = apply(Cs, 2, function(x) sqrt(sum(x^2))))
  return(Cs)
}

