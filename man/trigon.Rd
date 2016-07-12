\name{trigon}
\alias{trigon}
\title{crea vars ++}
\description{
mas
}
\usage{
trigon(data, frecH = 1, n = nrow(data), w = 24, q = 3)
}
\arguments{
  \item{data}{ \code{data.frame} contains year, month, day and hour columns}
  \item{frecH}{numero de dias por ciclo. diario=1, semanal=7 o anual=365}
  \item{n}{numero de filas del output }
  \item{w}{numero de observaciones por dia (24 horario; 48 demihorario...) }
  \item{q}{cantidad de variables generadas}
}

\value{matrix with \code{q} trigonometric variables}

\author{Jairo Cugliari, Andres Castrillejo, Fernando Massa}


