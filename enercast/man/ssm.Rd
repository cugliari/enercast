\name{ssm}
\alias{ssm}
\title{Univariate State Space Model}
\description{VickyNando}
\usage{
ssm(formulaXa, formulaXf = NULL, data, ini = NULL, for.each = "hour", cores = detectCores() - 1)
}
\arguments{
  \item{formulaXa}{Formula with random regressors (\code{y~x1+x2+x3+x2:x3})}
  \item{formulaXf}{Formula with fixed regressors *without* response variable (\code{~x4+x5})}
  \item{data}{data.frame}
  \item{ini}{List of intial value parameters of length of 'for.each', 24 in
  hourly data}
  \item{for.each}{Variable of data.frame that indicate frequency}
  \item{cores}{For parallel}
}
\value{
  fkf.optim
   fk list(dt=dt,ct=ct,Tt=Tt,Zt=Zt,a0=a0,P0=P0)
}
\references{
Dordonnat, V. and Koopman, S. J. and Ooms, M. and Dessertaine, A. and
Collet, J. 2008 An hourly periodic state space model for modelling
French national electricity load. \emph{International Journal of
  Forecasting}, \bold{24(4)}:566--587.

Dordonnat, V. and Koopman, S. J. and Ooms, M. 2012 Dynamic factors in
periodic time-varying regressions with an application to hourly
electricity load modelling. \emph{Computational Statistics & Data
  Analysis}, \bold{56(11)}:3134--3152.
}
\author{Jairo Cugliari, Andres Castrillejo, Fernando Massa}


