\name{newvars.rf}
\alias{newvars.rf}
\title{Construct new variables for Random Forest Model}
\description{Construction of variables}
\usage{newvars.rf(data, type, t, ...)}
\arguments{
  \item{data}{\code{data.frame} contains 'load' and 'temp' columns}
  \item{type}{Character. 'load.lag', 'temp.lag', 'temp.sqcb',
    'temp.desv', 'load.desv'.}
  \item{t}{Lags}
  \item{\dots}{...}
}
\details{
           load.lag lags of load
           temp.lag lags of temperature
           temp.sqcb temperature square and cubic
           temp.desv deviation of temperature
           load.desv deviation of load
	 }
\references{
  Hong, T. and Wang, P. and Lee Willis, H. 2011 A naive multiple linear regression
  benchmark for short term load forecasting. In \emph{Power and Energy Society
    General Meeting, 2011 IEEE}, 1--6.
  
  Hong, T. and Pinson, P.and Fan, S. 2014 Global energy forecasting competition
  2012. \emph{International Journal of Forecasting}, \bold{30(2)}:357--363
}

\author{Jairo Cugliari, Andres Castrillejo, Fernando Massa}
