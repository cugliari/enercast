plot.load.global <- function(yt, dataux, ...){
    if (missing(dataux)) {
        plot(yt, type = 'l', ylab = 'Load', xlab = 'Time', ...)
    } else {
      names(dataux)[1] <- 'year'
      aux1 <- by(dataux, dataux$year, nrow)
      aux2 <- cumsum(aux1) - (aux1 / 2)
        
      plot(yt, type = 'l', ylab = 'Load', xlab = 'Time', axes = FALSE, ...); box()

      axis(1, at = aux2,labels = names(aux1), tick = FALSE)
      axis(1, at = c(0, cumsum(aux1)), labels = FALSE, tick = TRUE)
      axis(2)      
    }
}

plot.load.wday <- function(data,...){
    auxwday <- factor(data$wday,levels = c('Mon', 'Tue', 'Wed', 'Thu',
                                           'Fri', 'Sat', 'Sun'),
                      labels = c('Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat', 'Sun'),
                      ordered = TRUE)
    auxsday <- factor(data$sdays, labels = c('regular', 'special'), ordered = TRUE)
    days <- unique(auxwday)
    mean.d <- NULL
    for (d in days) {
        DAT <- subset(data[(auxsday == 'regular' & auxwday == d),])
        mean.d <- rbind(mean.d, as.vector(by(DAT$load, DAT$hour, mean, na.rm = TRUE)))
    }
    row.names(mean.d) <- days
    plot(data$load, type = 'n', xlim = c(1, ncol(mean.d)), 
         ylim = c(min(mean.d), max(mean.d)),
         ylab = 'Load', xlab = 'Time', main = 'Mean Daily Load')
    aux <- seq(1,ncol(mean.d), by = ncol(mean.d) / length(days))
    for (d in 1:length(days)) {
        lines(mean.d[d, ], col = d, lwd = 2)
        mtext(row.names(mean.d)[d], col = d, cex = 0.7, side = 3, line = 0.5, 
              at = aux[d])
    }
}



plot.load.month <- function(data,...){
    auxmonth <- factor(data$month, labels = c('Jan', 'Feb', 'Mar', 'Apr',
                                              'May', 'Jun', 'Jul', 'Aug',
                                              'Sep', 'Oct', 'Nov', 'Dec'),
                       ordered = TRUE)
    auxsday <- factor(data$sdays, labels = c('regular', 'special'), ordered = TRUE)
    months <- unique(auxmonth)
    mean.m <- NULL
    for (m in months) {
        DAT <- subset(data[(auxsday == 'regular' & auxmonth == m),])
        mean.m <- rbind(mean.m, as.vector(by(DAT$load, DAT$hour, mean, na.rm = TRUE)))
    }

    row.names(mean.m) <- months
    
    plot(data$load, type = 'n', xlim = c(1, ncol(mean.m)), 
         ylim = c(min(mean.m), max(mean.m)), ylab = 'Load',
         xlab = 'Time', main = 'Mean Monthly Load')
    aux <- seq(min(mean.m), max(mean.m), 
               by = (max(mean.m) - min(mean.m)) / length(months))
    aux <- seq(1, ncol(mean.m), by = ncol(mean.m) / length(months))
    for (m in 1:length(months)) {
        lines(mean.m[m ,], col = m, lwd = 2)
        mtext(row.names(mean.m)[m], col = m, cex = 0.7, side = 3, line = 0.5, 
              at = aux[m])
    }
}



plot.load.corr <- function(yt, data = NULL, nw = 3, ...){
  if (is.null(data)) {
    s1 = 24
    s2 = 168
  } else {
    frec <- apply(data[ ,c('wday', 'hour')], 2, function(x) length(unique(x)))
    s1 <- frec[names(frec) == 'hour']
    s2 <- prod(frec)
  }
  lag.max = s2 * nw
  par(las = 1)
  
  layout(1:3)

  # ACF11
  acf(yt, lag.max = s2 * nw, axes = FALSE); box()
  axis(2, seq(-1, 1, 0.2))
  axis(1, seq(0, by = s2, length = nw + 1), 
       labels = c(NA, paste('week', seq(3))))
#  par(ask=TRUE,las=1)

  # ACF 2
  acf(diff(yt, 1), lag.max = s2 * nw, axes = FALSE); box()
  axis(2, seq(-1, 1, 0.2))
  axis(1, seq(0, by = s2, length = nw + 1), labels = c(NA, paste('week', seq(3))))
#  par(ask=TRUE,las=1)

  # ACF 3
  acf(diff(diff(yt), s1), lag.max = s2 * nw, axes = FALSE); box()
  axis(2, seq(-1, 1, 0.2))
  axis(1, seq(0, by = s2, length = nw + 1), labels = c(NA, paste('week', seq(3))))
}



plot.temp <- function(data, group, ...){
  
  if (!requireNamespace("hexbin", quietly = TRUE)) {
    stop("Package spams is needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if (missing(group)) group <- 1 # form <- data$load ~ data$temp  
  
  form <- data$load ~ data$temp | group

  rf <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, 'Spectral')))
  hexbin::hexbinplot(form, colramp = rf, xlab = "temperature",
                       ylab = "load", colorkey = FALSE) #, aspect = 1 )
}


