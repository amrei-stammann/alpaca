% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R
\name{print.APEs}
\alias{print.APEs}
\title{Print \code{APEs}}
\usage{
\method{print}{APEs}(x, digits = max(3L, getOption("digits") - 3L), ...)
}
\arguments{
\item{x}{an object of class \code{"APEs"}.}

\item{digits}{unsigned integer indicating the number of decimal places. Default is 
\code{max(3L, getOption("digits") - 3L)}.}

\item{...}{other arguments.}
}
\description{
\code{\link{print.APEs}} is a generic function which displays some minimal information from 
objects returned by \code{\link{getAPEs}}.
}
\seealso{
\code{\link{getAPEs}}
}
