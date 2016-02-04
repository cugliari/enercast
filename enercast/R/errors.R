mape <- function(y, yhat){mean(abs(y - yhat)/abs(y))}
mae <- function(y, yhat){mean(abs(y - yhat))}
rmse <- function(y, yhat){sqrt(mean((y - yhat)^2))}
pdad <- function(y, yhat){mean(abs(y - yhat)/yhat) * 100}
