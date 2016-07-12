\name{hong}
\alias{hong}
\title{Fit Hong Vanilla multiple regression model}
\description{Compute Hong Vanilla multiple regression model}
\usage{
hong(data)
}
\arguments{
  \item{data}{Data frame with columns:
    \itemize{
     \item 'year'.  type character. (...2012, 2013, 2014...)
     \item 'month'. type character. (01, 02, ..., 12)
     \item 'wday'.  type character. (Mon Thu Wed Tue Fri Sat Sun)
     \item 'hour'.  type character. (01, 02, ..., 24)
     \item 'temp'.  type numeric.
     \item 'load'.  type numeric.
     }
}
}
\details{ The Hong Vanilla Model:
  \deqn{E(Load) = \beta_0 + \beta_1 Trend +\beta_2 Day \times Hour  +\beta_3 Month +\beta_4 Month  \times Temp}   
  \deqn{+\beta_5 Month  \times Temp^2  +\beta_6 Month  \times Temp^3 +\beta_7 Hour \times Temp}  
  \deqn{+\beta_8 Hour \times Temp^2  +\beta_9  Hour \times Temp^3 }
}



\value{An object of class \code{lm}}
\references{
  Hong, T. and Wang, P. and Lee Willis, H. 2011 A naive multiple linear regression
  benchmark for short term load forecasting. In \emph{Power and Energy Society
  General Meeting, 2011 IEEE}, 1--6.
  
  Hong, T. and Pinson, P.and Fan, S. 2014 Global energy forecasting competition
  2012. \emph{International Journal of Forecasting}, \bold{30(2)}:357--363
}
\author{Jairo Cugliari, Andres Castrillejo, Fernando Massa}

\seealso{ \link[stats]{lm} \link{predict.hong}}

