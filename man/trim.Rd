% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{trim}
\alias{trim}
\title{Trim an odd-length vector}
\usage{
trim(x, n = 5, LR = 1:2)
}
\arguments{
\item{x}{The vector to trim}

\item{n}{The number of elements to one or either side of the center element to keep}

\item{LR}{A vector indicating which side to keep: To keep \code{n} elements to the left, 
set \code{LR = 1}; to keep \code{n} elements to the right, set \code{LR = 2}; for
\code{n} elements to the left and right, set \code{LR = 1:2}.}
}
\value{
A trimmed \code{vector} according to \code{n} and \code{side}.
}
\description{
This function trims a vector with an odd length, leaving \code{n} elements to
one or either side of the center element.
}
