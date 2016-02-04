\name{tsb}
\alias{tsb}
\title{Fit Time Series Benchmarck Model}
\description{Fit a double seasonal ARIMA type model used as benchmark in
STLF}
\usage{tsb(yt, p = 1, d = 1, q = 1, P1 = 1, D1 = 0, Q1 = 1, P2 = 0, D2 = 1, Q2 = 1, s1 = 24, s2 = 168, k = 200)}

\arguments{
  \item{yt}{load data sequence. type numeric.}
  \item{p, P1, P2}{AR order. type numeric.}
  \item{q, Q1, Q2}{MA order. type numeric.}
  \item{d, D1, D2}{Degree of differencing. type numeric.}
  \item{s1, s2}{Seasonal frequencies. type numeric.}
  \item{k}{number of lags}
}
\details{}
\value{ A list with ARIMA estimated model
 \item{ar}{(reg=ar,est1=AR1,est2=AR2)}
 \item{ma}{(reg=ma,est1=MA1,est2=MA2)}
 \item{d}{(d,D1,D2)}
 \item{frec}{}
 \item{dif}{}
 \item{iter}{number of iterations}
 \item{data}{load data}

}
\references{Brockwell, P.J.and Davis, R.A. 1991 \emph{Time Series: theory and methods} 2nd ed. Springer.}
\author{Todos Nosotros}
\seealso{\code{\link{predict.tsb}}}



 