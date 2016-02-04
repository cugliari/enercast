\name{errors}
\alias{mape}
\alias{mae}
\alias{rmse}
\alias{pdad}
\title{Compute Predictions Errors}
\description{Calculate the predictions errors. Mean absolute
  percentage error(MAPE), mean absolute error (MAE), root mean square error (RMSE).}
\usage{
mape(y, yhat)
mae(y, yhat)
rmse(y, yhat)
pdad(y, yhat)
}
\arguments{
  \item{y}{vector of observed load}
  \item{yhat}{vector of predicted load}
}
\details{
  \itemize{
     \item Mean Absolute Percentage Error:  \deqn{mape(y,\hat{y})=\frac{1}{N}\sum\limits_{t=1}^N \frac{|y_t - \hat{y}_t|}{|y_t|}}
     \item Mean Absolute Error:             \deqn{mae(y,\hat{y}) =\frac{1}{N}\sum\limits_{t=1}^N |y_t - \hat{y}_t|}
     \item Root Mean Square Error:          \deqn{rmse(y,\hat{y})=\sqrt{\frac{1}{N}\sum\limits_{t=1}^N (y_t - \hat{y}_t)^2}}
     \item Promedio de Desviacion Absoluta de la Demanda:          \deqn{pdad(y,\hat{y})=\frac{1}{N}\sum\limits_{t=1}^N \frac{|y_t - \hat{y}_t|}{\hat{y}}*100}

   }
}
\references{
  Hyndman, R.J. and Athanasopoulos, G. 2014 \emph{Forecasting: principles and practice}. OTexts.
}
\author{Todos nosotros}


