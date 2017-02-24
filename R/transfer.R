#' Estimate the frequency-domain transfer function
#' 
#' Estimates a transfer function between the columns of input \code{x} and the response
#' \code{y}.
#' 
#' @param x A \code{data.frame} whose columns are the time domain input series.
#' @param y A single column \code{data.frame} containing the response series.
#' @param blockSize A \code{numeric} indicating the block sizes into which the 
#' input and response series will be partitioned.
#' @param overlap A \code{numeric} between 0 and 1 indicating how much overlap should exist 
#' between adjacent blocks.
#' @param deltat A \code{numeric} indicating the sample rate.
#' @param nw A \code{numeric} indicating the time bandwidth parameter for estimating the 
#' Slepian data tapers.
#' @param k A \code{numeric} indicating the number of tapers to use - should be approximately
#' floor(2*nw - 1) and no larger than floor(2*nw).
#' @param nFFT A \code{numeric} indicating the number of frequency bins to use (i.e. setting 
#' the zeropadding amount).
#' @param freqRange NOT CURRENTLY IMPLEMENTED.
#' @param freqOffset NOT CURRENTLY IMPLEMENTED (don't chage this... ).
#' @param standardize Should the inputs and outputs be standardized to have mean = 0 and standard deviation = 1? 
#' 
#' @details Takes the times series inputs and response, divides these series into 
#' (optionally) overlapping blocks, tapers each block with Discrete 
#' Prolate Spheriodal Sequences (DPSS's or Slepian sequences), Fourier transforms each 
#' block, and then estimates the transfer function at each frequency between the Fourier 
#' transforms of the inputs and the response.
#' 
#' @return An object of class \code{transfer}, consisting of a complex matrix whose 
#' columns are the individual transfer function for each input, and several attributes
#' describing the transfer function estimate.
#' 
#' @export
tf <- function(x, y, blockSize = dim(x)[1], overlap = 0, deltat = 1, nw = 4, k = 7, nFFT = NULL
               , freqRange = NULL, freqOffset = NULL, standardize = FALSE){
  if( standardize ){
    stdPars <- vector( mode = "list" )
      stdPars$xmean <- sapply( x, mean )
      stdPars$ymean <- sapply( y, mean )
      stdPars$xsd <- sapply( x, sd )
      stdPars$ysd <- sapply( y, sd )
      stdPars <- data.frame( stdPars )
    std <- function( a ) (a - mean(a))/sd(a)
    x <- data.frame( lapply( x, std ) )
    y <- data.frame( lapply( y, std ) )
  }
  x2 <- sectionData(x, blockSize = blockSize, overlap = overlap)
  y2 <- sectionData(y, blockSize = blockSize, overlap = overlap)
  x3 <- taper(x2, nw = nw, k = k)
  y3 <- unlist(taper(y2, nw = nw, k = k), recursive = FALSE)
  
  time <- seq(0, (blockSize - 1)*deltat, by = deltat)
  
  if (is.null(nFFT)){
    nFFT <- 2^(floor(log2(blockSize)) + 3)
  }
  
  if (is.null(freqOffset)){
    freqOffset = 0
  }
  
  freq <- seq(0, 1/(2*deltat), by = 1/(nFFT*deltat))
  freqIdx <- 1:length(freq)
  
  H.mat <- olsTf(x = x3, y = y3, time = time, n = blockSize
             , npredictor = length(x3[[1]]), ntaper = k
             , freq = freq[freqIdx], fOffset = freqOffset)
  
  H <- data.frame(freq, H.mat)
  colnames(H) <- c("freq", colnames(x))
  
  attr(H, "class") <- c(attr(H, "class"), "transfer")
  attr(H, "blockSize") <- attr(x3, "blockSize")
  attr(H, "overlap") <- attr(x3, "overlap")
  attr(H, "numSections") <- attr(x3, "numSections")
  attr(H, "n") <- attr(x3, "n")
  attr(H, "dataTaper") <- attr(x3, "dataTaper")
  attr(H, "nw") <- attr(x3, "nw")
  attr(H, "k") <- attr(x3, "k")
  attr(H, "standardize") <- standardize
  if( standardize ) attr(H, "stdPars") <- stdPars
  
  H
}



#' Predict using an estimated frequency-domain transfer function
#' 
#' Using output from the \code{\link{tf}} function, and an input (multivariate) time series,
#' this function generates predicted values.
#' 
#' @param object An object of class \code{transfer}, from a call to the \code{\link{tf}} function.
#' @param newdata A \code{data.frame} whose columns are the time domain input series.
#' @param filterMod A \code{function} to be applied to the filter coefficients before convolution.
#' @param ... additional arguments passed to \code{filterMod}.
#' 
#' @details The transfer function estimate is used to calculate filter coefficients to be
#' applied in the time domain via convolution with the input series \code{newdata}.
#' Prior to the convolution, the filter coefficients can be modified using the
#' \code{filterMod} function.
#' 
#' @return A \code{data.frame} with the predicted values obtained by filtering 
#' the input series \code{newdata}.
#' 
#' @export

predict.transfer <- function( object, newdata, filterMod = trim, ... ){
  # Match up names of newdata to the object names
  tfNames <- names( object )
  nm <- match( names( newdata ), tfNames ) # The positions of newdata names in object names 
  nmNA <- which( is.na(nm) )
  nm2 <- which( !is.na(nm) ) # The newdata names that exist in object names
  if( length(nmNA) > 0 ) warning( "The names of newdata do not match the names of object" )
  if( length(nmNA) == length( nm ) ) stop( "None of the names of newdata match the names of object" )
  if( length(tfNames) - length(nm2) > 1 ) warning( "The newdata names do not match all of the names of object" )
  
  # Standardize inputs if required
  if( attr( object, "standardize" ) ){
    stdPars <- as.data.frame( t( attr( object, "stdPars" ) ) )[c("xmean","xsd"),]
    nm3 <- match( names( stdPars ), names( newdata[,nm2] ) ) # The positions of stdPars names in newdata names
    std <- function( x, sP ) (x - sP[1])/sP[2]
    newdata <- data.frame( mapply( FUN = std, x = newdata[,nm2], sP = stdPars[,nm3], SIMPLIFY = FALSE ) )
  }else{
    newdata <- newdata[,nm2]
  }
  
  # Rearrange transfer function coefficients
  objectC <- lapply( object, Conj )
  attributes( objectC ) <- attributes( object )
  objectC <- objectC[(nrow(object)-1):2, ]
  objectFull <- rbind( object, objectC )[,nm]
  
  # Inverse FFT of the coefficients
  fC <- lapply( objectFull, function(x,...) fft(x,...)/length(x), inverse = TRUE )
  
  # Rearrange the filter coefficients
  fC <- as.data.frame( lapply( fC, Re ) )
  
  # Only want n filter coefficients
  n <- attr( object, "n" )/2
  # The causal part of the filter
  ind <- n:1
  # The non-causal part of the filter
  ind <- c( ind, nrow( fC ):(nrow(fC)-n+2) )
  fC <- fC[ind,]
  
  # Apply the filterMod function to the coefficients
  fC <- data.frame( lapply( fC, filterMod, ... ) )
  
  # Compute prediction using filter
  out <- as.data.frame( mapply( filter, x = newdata, filter = fC ) )
  
  # Return prediction
  out <- rowSums(out)
out
}


#' Trim an odd-length vector
#' 
#' This function trims a vector with an odd length, leaving \code{n} elements to
#' one or either side of the center element.
#' 
#' @param x The vector to trim
#' @param n The number of elements to one or either side of the center element to keep
#' @param side A vector indicating which side to keep: To keep \code{n} elements to the left, 
#' set \code{side = 1}; to keep \code{n} elements to the right, set \code{side = 2}; for
#' \code{n} elements to the left and right, set \code{side = 1:2}.
#' 
#' @return A trimmed \code{vector} according to \code{n} and \code{side}.
#' 
#' @export
trim <- function(x, n = 5, side = 1:2){
  if( abs( length(side)-1 ) > 1 | !all( side %in% 1:2 ) ) stop( "Invalid side argument" )
  l <- length(x)
  if( l%%2 != 1 ) stop( "Length of x is not odd" )
  m <- ceiling( l/2 )
  i <- list( (m-n):m, m:(m+n) )
x[ unique( unlist(i[side]) ) ]
}