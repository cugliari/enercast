\name{internals.tsb}
\alias{comb.coef}
\alias{maInvert}
\alias{acft}
\alias{invento2}
\title{Internal functions for ARIMA estimation}
\description{blablabla}
\usage{
acft(ar=NULL,AR1=NULL,AR2=NULL,ma=NULL,MA1=NULL,MA2=NULL,s1=24,s2=168,k=500)
}
\arguments{
  \item{ar}{load data sequence}
  \item{ar}{ son los p parametros de la parte autoregresiva regular}
  \item{AR1}{ son los P1 parametros de la parte autoregresiva de la estacional s1}
  \item{AR2}{ son los P2 parametros de la parte autoregresiva de la estacional s2}
  \item{ma}{ son los q parametros de la parte de medias moviles regular}
  \item{MA1}{ son los Q1 parametros de la parte de medias moviles de la estacional s1}
  \item{MA2}{ son los Q2 parametros de la parte de medias moviles de la estacional s2}
  \item{s1}{ es el periodo de la primera estacionalidad}
  \item{s2}{ es el periodo de la segunda estacionalidad}
  \item{k}{ es el numero de rezagos}

}
\details{}
\value{}
\references{}
\author{Todos Nosotros}
\seealso{ \code{\link{tsb}}}
