\name{plots.load}
\alias{plot.load.global}
\alias{plot.load.wday}
\alias{plot.load.month}
\alias{plot.temp}
\alias{plot.load.corr}
\title{Descriptive plots for load data}
\description{Plot load data: the complete time series, data by day of  week (mean) and by month (mean) excluding special days.}
\usage{
plot.load(yt,data=NULL, ...)
plot.load.wday(data,...)
plot.load.month(data,...)
plot.load.corr(yt,data=NULL,nw=3,...)
plot.temp(data, group=NULL,...)
}
\arguments{
  \item{yt}{load data sequence}
  \item{dataux}{Optional data.frame with 'year' variable}
  \item{data}{Data frame with columns:
    \itemize{
      \item 'year'.  type character. (...2012, 2013, 2014...)
      \item 'month'. type character. (01, 02, ..., 12)
      \item 'wday'.  type character. (Mon Thu Wed Tue Fri Sat Sun)
      \item 'hour'.  type character. (01, 02, ..., 24)
      \item 'temp'.  type numeric.
      \item 'load'.  type numeric.
      \item 'sdays'. type character. (0, 1) ('regular day','special day')
    }
  }
  \item{group}{type character}
  \item{nw}{Number of weeks to plot} 
  \item{\dots}{further arguments to be passed to 'plot'}
}
\author{Jairo Cugliari, Andres Castrillejo, Fernando Massa}

