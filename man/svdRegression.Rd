% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{svdRegression}
\alias{svdRegression}
\title{Linear regression using a singular value decomposition}
\usage{
svdRegression(x, y)
}
\arguments{
\item{x}{A \code{matrix} or \code{vector} containing the predictors.}

\item{y}{A \code{vector} containing the response.}
}
\value{
A \code{list} containing the regression coefficients corresponding to the 
columns of x, the standard error estimates (*not* the variance) on those coefficients, 
and the eigenvalues of the decomposition (singular values squared).
}
\description{
Performs a multivariate linear regression using a singular value decomposition.
}
\details{
Uses the process described in Mandel (1982), "Use of the Singular value 
Decomposition in Regression Analysis".
}
