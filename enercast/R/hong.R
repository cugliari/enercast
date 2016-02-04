hong <-
function(data){

    trend <- 1:nrow(data)
    
    if(!is.character(data$hour))   data$hour  <- as.character(data$hour)
    if(!is.character(data$month))  data$month <- as.character(data$month)
 
    hv <- lm(load ~ trend + wday * hour + month + month * temp + 
                    month * (temp^2) + month * (temp^3) +  hour * temp +  
                    hour  * (temp^2) +  hour * (temp^3),
             data = data)
    
    class(hv) <- c('lm', 'hv_class')
    return(hv)
}
