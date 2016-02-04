\name{predict.hong}
\alias{predict.hong}
\title{Predict method for Hong Vanilla Model}
\description{Predicted values for T. Hong multiple regression model}
\usage{predict.hong(obj.hv, new.data = NULL)}
\arguments{
  \item{obj.hv}{Object of class inheriting from use of \code{hong} function.}
  \item{new.data}{Data frame with predictors. If omitted, the fitted values are used.}
}
\details{ \code{new.data} must contain the same predictors that data in \link{hong}. 

  Aditionally \code{new.data} must contain \code{trend} variable. This
  variable was computed like a sequence from 1 to number of rows of
  input data.frame in \link{hong} model function.

  In forecasting new consecutive periods, \code{trend} variable in
  new.data must be computed like a sequence from number of rows of
  original data.
}
\value{predicted values.}

\author{Todos nosotros}

\seealso{ \link[stats]{predict.lm} \link{hong}}
\examples{}
