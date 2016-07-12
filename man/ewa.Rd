\name{ewa}
\alias{ewa}
\title{Calibrate the EWA algorithm}
\description{\code{ewa} agregates expert predictions using exponential weights and a
 learning rate parameter \eqn{\eta}}
\usage{ewa(eta, y, experts, EPS = 1e-05)}
\arguments{
  \item{eta}{numeric value to start the initial grid for \eqn{\eta}}
  \item{y}{numeric vector containing the time series}
  \item{experts}{matrix of time points X individual predictions (experts) }
  \item{EPS}{Epsilon tolerance}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{A list containg: \code{w} the matrix of fitted weights and \code{eta} the obtain value for \eqn{\eta}.}
\author{Jairo Cugliari, Andres Castrillejo, Fernando Massa}

