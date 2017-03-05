#' Section data into data blocks
#' 
#' Sections provided data using block size and overlap proportion
#' 
#' @param x A \code{data.frame} where each column is assumed to be a series.
#' @param blockSize An \code{integer} indicating the size of each section.
#' @param overlap An \code{numeric} between 0 and 1 indicating how much overlap should 
#' exist between sections.
#' 
#' @return A \code{list} whose elements contain a \code{data.frame} with the same number of 
#' columns as \code{x} and \code{blockSize} number of rows.
#' 
#' @export
sectionData <- function(x, blockSize, overlap = 0){
  increment <- ceiling(blockSize * (1-overlap))
  
  sectIdx <- seq(1, dim(as.data.frame(x))[1] - blockSize+1, by=increment)
  numSect <- length(sectIdx)
  
  blockIdxLst <- as.list(as.data.frame(mapply(":", sectIdx, sectIdx+blockSize-1)))
  
  sectionedData <- list()
  for (i in 1:numSect){
    sectionedData[[i]] <- as.data.frame(x)[blockIdxLst[[i]], , drop = FALSE]
  }
  
  attr(sectionedData, "class") <- c(class(sectionedData), "blocked")
  attr(sectionedData, "blockSize") <- blockSize
  attr(sectionedData, "overlap") <- overlap
  attr(sectionedData, "numSections") <- numSect
  attr(sectionedData, "n") <- dim(x)[1]
  
  sectionedData
}

#' Multiplies a vector by a data taper
#' 
#' Tapers the data using a discrete prolate spheroidal sequence (DPSS)
#' 
#' @param x A \code{list} of \code{data.frame}'s as a result of calling \code{sectionData()}.
#' @param dataTaper A \code{character} string indicating what type of taper to use.  
#' Currently only DPSS's tapers are available.
#' @param nw A \code{numeric} representing the time-bandwidth paramater for use in generating 
#' the DPSS's.
#' @param k A \code{numeric} indicating the number of tapers to be used.  Should be less 
#' than 2*nw.
#' 
#' @return A \code{list} of \code{lists}, each "sublist" a named matrix of length 
#' \code{blockSize} by \code{k} of tapered data for that series.
#' 
#' @details This is not implemented in an efficient manner - there needs to be a better way
#' written.  That will be in the future.
#' 
#' @export
taper <- function(x, dataTaper = "dpss", nw = 4, k = 7){
  if (!any(class(x) == "blocked")){
    stop("Call sectionData() on x before calling this function.")
  }
  
  slep <- dpss(n = attr(x, "blockSize"), nw = nw, k = k)$v
  
  tapered <- list()
  for (i in 1:attr(x, "numSections")){
    tmp <- list()
    
    for (j in 1:dim(x[[i]])[2]){
      tmp[[ colnames(x[[i]])[j] ]] <- apply(slep, 2, "*", x[[i]][, j])
    }
    tapered[[i]] <- tmp
  }
  
  mostattributes(tapered) <- attributes(x)
  attr(tapered, "dataTaper") <- dataTaper
  attr(tapered, "nw") <- nw
  attr(tapered, "k") <- k
  
  tapered
}

#' Linear regression using a singular value decomposition
#' 
#' Performs a multivariate linear regression using a singular value decomposition.
#' 
#' @param x A \code{matrix} or \code{vector} containing the predictors.
#' @param y A \code{vector} containing the response.
#' 
#' @return A \code{list} containing the regression coefficients corresponding to the 
#' columns of x, the standard error estimates (*not* the variance) on those coefficients, 
#' and the eigenvalues of the decomposition (singular values squared).
#' 
#' @details Uses the process described in Mandel (1982), "Use of the Singular value 
#' Decomposition in Regression Analysis".
#' 
#' @export
svdRegression <- function(x, y){
  # X = U D V'
  sng <- svd(x)
  
  beta <- sng$v %*% ((Conj(t(sng$u)) %*% y) / sng$d)
  
  stdErr <- apply(t(apply(sng$v, 1, "/", sng$d)), 1, sum) * sd.complex(y)
  
  list(coef = t(beta), stdErr = stdErr, ev = sng$d^2)
}

#' Standard deviation of complex values
#' 
#' Calculates the standard deviation of a complex-valued series.
#' 
#' @param x a \code{vector} containing real or complex values.
#' 
#' @details Calculates the standard deviation of either a real or complex series using
#' sd = E{(X-u)(X-u)'}
#' where the apostrophe denotes the complex conjugate.
#' 
#' @return A real valued scalar.
#' 
#' @export
sd.complex <- function(x){
  x <- as.vector(x)
  Re((1/(length(x) - 1)) * sum((x - mean(x)) * Conj(t(x - mean(x)))))
}