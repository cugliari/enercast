select.wkdata <-
  function(data,colnames=NULL){
    if ( (!is.null(class(data))) && (class((data)) != "wkdata"))
      stop("data is not of class wkdata")
    
    if (!is.null(colnames))
      data <- list( X=data$X[,colnames], S0=data$S0[colnames],
                    D0=data$D0[,colnames], p=data$p, J=data$J,
                    gr=data$gr[colnames])
    class(data) <- "wkdata"
    return(data)
  }

