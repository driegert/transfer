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
#' @param method A \code{character} string indicating which method to use.  See details.
#' @param lOrd A vector with length == numColumns of x.  The order of the regression for svdBendat method.
#' @param adaptiveWeighting A \code{logical} indicating whether the eigencoefficients should 
#' be multiplied by their adaptive weights before performing the regression.
#' Note: Only implemented for \code{method = "svd"}.
#' @param freqRange NOT CURRENTLY IMPLEMENTED.
#' @param freqOffset NOT CURRENTLY IMPLEMENTED (don't chage this... ).
#' @param standardize Should the inputs and outputs be standardized to have mean = 0 and standard deviation = 1? 
#' @param prewhiten NOT CURRENTLY IMPLEMENTED.
#' @param removePeriodic NOT CURRENTLY IMPLEMENTED.
#' 
#' @details Takes the times series inputs and response, divides these series into 
#' (optionally) overlapping blocks, tapers each block with Discrete 
#' Prolate Spheriodal Sequences (DPSS's or Slepian sequences), Fourier transforms each 
#' block, uses the adaptive weights (\code{method = "svd"}) and then estimates the 
#' transfer function at each frequency between the Fourier 
#' transforms of the inputs and the response.
#' 
#' The \code{method} argument indicates how the transfer function should be estimated.  
#' \code{method = "sft"} does the following; 1) a slow Fourier transform is used, 2) no 
#' adaptive weights are used and 3) the matrix definition of the regression coefficients
#' are used:
#' $(X^{t}X)^{-1})X^{-1}y$
#' \code{method = "svd"} uses 1) the FFT, 2) adaptive weights, and 3) a singular value 
#' decomposition method for estimating the regression coefficients.
#' \code{method = "svdBendat"} decorrelates the inputs prior to estimating the transfer functions, then 
#' obtains the appropriate transfer function after (see Chapter 7.3 of Bendat & Piersol).
#' Overall - \code{method = "svdBendat"} is in all probability the better method to use.
#' 
#' @return An object of class \code{transfer}, consisting of a complex matrix whose 
#' columns are the individual transfer function for each input, and several attributes
#' describing the transfer function estimate.
#' 
#' @export
tf <- function(x, y, blockSize = dim(x)[1], overlap = 0, deltat = 1, nw = 4, k = 7, nFFT = NULL
               , method = c("svd", "sft", "robust", "svdBendat"), lOrd = NULL
               , adaptiveWeighting = TRUE
               , freqRange = NULL, freqOffset = NULL
               , standardize = TRUE, prewhiten = TRUE, removePeriodic = TRUE)
{
  x.nCol <- dim(x)[2]
  y.nCol <- dim(y)[2]
  
  # for a future implementation - a brighter tomorrow perhaps... (faster anyway... )
  # if (nodes == 1){
  #   warning("'nodes' is set to 1.  You may want to run this in parallel to gain some speed.")
  # }
  
  if (is.null(lOrd) & method[1] == "svdBendat"){
    lOrd <- 1:x.nCol
  } else if (method[1] == "svdBendat" & length(lOrd) != x.nCol){
    stop("Right now, length of 'lOrd' must have same length as number of columns of 'x' data.frame.")
  }
  
  # standardize the series by removing the mean and dividing by the standard deviation
  if( standardize ){
    stdPars <- vector( mode = "list" )
    stdPars$x <- data.frame( xmean = sapply( x, mean ), xsd = sapply( x, sd ) )
    stdPars$y <- data.frame( ymean = sapply( y, mean ), ysd = sapply( y, sd ) )
    x <- data.frame( lapply( x, std ) )
    y <- as.data.frame( lapply( y, std ) )
  }
  
  # number of frequencies bins to use (zero-pad length)
  if (is.null(nFFT)){
    nFFT <- 2^(floor(log2(blockSize)) + 3)
  }
  
  # This would be used if incorporating offset frequencies into the transfer function
  # regression.
  if (is.null(freqOffset)){
    freqOffset = 0
  }
  
  # determine the positive values of the frequencies.
  freq <- seq(0, 1/(2*deltat), by = 1/(nFFT*deltat))
  nfreq <- length(freq)
  
  # block the data (x2, y2 are a list of data.frames)
  x2 <- sectionData(x, blockSize = blockSize, overlap = overlap)
  y2 <- sectionData(y, blockSize = blockSize, overlap = overlap)
  
  if (method[1] != "sft"){
    numSections <- attr(x2, "numSections")
    
    x.wtEigenCoef <- blockedEigenCoef(x2, deltat = deltat, nw = nw, k = k
                                           , nFFT = nFFT, numSections = numSections
                                           , adaptiveWeighting = adaptiveWeighting)
    y.wtEigenCoef <- blockedEigenCoef(y2, deltat = deltat, nw = nw, k = k
                                           , nFFT = nFFT, numSections = numSections
                                           , adaptiveWeighting = adaptiveWeighting)
    
    # indexing helper function that grabs and stacks all the eigencoefficients at a frequency
    eigenByFreq <- function(obj, rowNum, numEl){
      matrix(unlist(lapply(obj, function(x, idx) x[idx, ], rowNum)), ncol = numEl)
    }
    
    x.design <- list()
    
    # stack the eigencoefficients by frequency from each block
    # form into a single list
    for (i in 1:nfreq){
      x.design[[i]] <- list(x = do.call(rbind, lapply(x.wtEigenCoef, eigenByFreq, rowNum = i, numEl = x.nCol))
                            , y = do.call(rbind, lapply(y.wtEigenCoef, eigenByFreq, rowNum = i, numEl = y.nCol)))
    }
    
    # this could be optimized... duplicate objects everywhere <sigh>
    if (method[1] == "robust"){
      H.all <- mapply(c, x.design, lapply(x.design, function(obj){ svdRegression(obj$x, obj$y) }), SIMPLIFY = F)
      #testing the timing
      H.tmp <- lapply(H.all[seq(1, 8193, length.out = 100)], robust.tf)
      H.mat <- matrix(unlist(H.tmp), ncol = x.nCol, byrow = T)
    } else if (method[1] == "svd"){
      H.tmp <- lapply(x.design, function(obj){ svdRegression(obj$x, obj$y) })
      H.mat <- do.call(rbind, lapply(H.tmp, "[[", "coef"))
    } else if (method[1] == "svdBendat"){
      H.tmp <- lapply(x.design, tfMiso, lOrd = lOrd)
      H.mat <- matrix(unlist(H.tmp), ncol = x.nCol, byrow = T)
    }
  } else if (method[1] == "sft") {
    # taper the blocked data
    x3 <- taper(x2, nw = nw, k = k)
    y3 <- unlist(taper(y2, nw = nw, k = k), recursive = FALSE)
    
    # determine an appropriate time vector (0 to T-dt)
    time <- seq(0, (blockSize - 1)*deltat, by = deltat)
    
    # freqIdx at which frequencies of the response the transfer function should be 
    # estimated (you could not use all frequencies if offsets were also used)
    freqIdx <- 1:length(freq)
    
    # perform the least squares complex-regresison.
    H.mat <- olsTf(x = x3, y = y3, time = time, n = blockSize
                   , npredictor = length(x3[[1]]), ntaper = k
                   , freq = freq[freqIdx], fOffset = freqOffset)
  }
  
  H <- data.frame(freq, H.mat)
  colnames(H) <- c("freq", colnames(x))
  
  attr(H, "class") <- c(attr(H, "class"), "transfer")
  attr(H, "blockSize") <- attr(x2, "blockSize")
  attr(H, "overlap") <- attr(x2, "overlap")
  attr(H, "numSections") <- attr(x2, "numSections")
  attr(H, "n") <- attr(x2, "n")
  attr(H, "dataTaper") <- attr(x2, "dataTaper")
  attr(H, "nw") <- attr(x2, "nw")
  attr(H, "k") <- attr(x2, "k")
  attr(H, "standardize") <- standardize
  if( standardize ) attr(H, "stdPars") <- stdPars
  attr(H, "adaptiveWeighting") <- adaptiveWeighting
  
  H
}

# helper function that performs the regressions based on residuals, then returns the proper 
# transfer functions.  This is mostly based on chapter 7.3 of Bendat and Piersol - Random Data
# Assumes the same time-domain block structure as what was passed to tf() above.
# Lord - the order that the regressions should take place before slamming everything back together 
# ... "SLAMMING" a la equation (7.100) in Bendat and Piersol.
tfMiso <- function(obj, lOrd){
  nReg <- length(lOrd) # number of intermediate regressions
  drow <- dim(obj$x)[1]
  dcol <- dim(obj$x)[2]
  
  X <- obj$x[, lOrd]
  
  L <- diag(as.complex(1), dcol, dcol)
  for (i in 2:nReg){
    L[1:(i-1), i] <- svdRegression(x = X[, 1:(i-1)], y = X[, i])$coef
    X[, i] <- X[, i] - (X[, 1:(i-1), drop = FALSE] %*% L[1:(i-1), i, drop = FALSE])
  }
  
  Ly <- t(svdRegression(x = X, y = obj$y)$coef)
  
  H <- svdRegression(x = L, y = Ly)$coef
  
  H[, order(lOrd)]
}

# obj is a list, containing the design matrix and 
# response vector, obj$x and obj$y respectively at a particular frequency (I guess.. .)
# H is the transfer function matrix (H.mat in the above)
robust.tf <- function(obj){
  # from page 195 of Chave and Jones:
  # 2a) compute the residuals
  r <- obj$y - (obj$x %*% t(obj$coef))
  # 2b) the scale using (5.40)
  d <- median(abs(r - median(r)))
  # 2c) residual sum of squares (r^{H} dot r)
  rss <- Re(Conj(t(r)) %*% r)[1]
  rss2 <- 3*rss
  
  # alpha for the Huber weights (usually 1.5?)
  alph <- 1.5
  
  #### START HERE WHEN CHECKING IT OUT
  # 2d) hat matrix diagonal (5.50)
  U <- diag(1, nrow = dim(obj$x)[1], ncol = dim(obj$x)[1])
  H <- sqrt(U) %*% obj$x %*% solve(hConj(obj$x) %*% U %*% obj$x) %*% hConj(obj$x) %*% sqrt(U)
  h <- diag(H)
  
  # number of predictors:
  p <- dim(obj$x)[2]
  
  Xi <- sum(diag(U))
  Xi.min <- Xi
  
  # Huber weights? (perhaps?)
  W <- diag(1, dim(obj$x)[1], dim(obj$x)[1])
  
  chi <- max( Re(h * p / sum(h)) )
  alpha <- pbeta(chi, p, Xi.min - p)
  alpha <- alpha - (1-alpha)/2
  chi <- qbeta(alpha, p, Xi.min - p)
  
  while(chi > qbeta(0.95, p, Xi.min - p)) {
    # print(paste("Outer: chi = ", chi))
    innerIterCount <- 0
    while (abs(rss2 - rss) / rss  > 0.01){
      if (innerIterCount > 5){ break; } else { innerIterCount <- innerIterCount + 1 }
      # print(paste("Inner - rss ratio: ", abs(rss2 - rss) / rss))
      rss <- rss2
      
      # Compute Huber weights - step 4
      V <- diag( ifelse(abs(r / d)[,,drop=T] <= alph, 1, alph / abs(r/d)[,,drop=T]),
                 length(r[,,drop=T]), length(r[,,drop=T]) )
      
      # y's are the "driving statistic" (5.48 in Chave & Jones)
      y <- Xi * h / p
      
      W <- diag(diag(W) * exp(exp(-chi^2)) * exp(-exp(chi*(y - chi))))
      w <- diag(W)
      
      elim <- which(round(abs(w*diag(V)), digits = 12) == 0)
      if (length(elim) > 0){
        print("This hit!")
        obj$x <- obj$x[-elim, ]
        obj$y <- obj$y[-elim]
      }
      
      # weight matrix
      U <- V %*% W
      
      # bounded influence estimator
      z <- solve(hConj(obj$x) %*% W %*% V %*% obj$x) %*% (hConj(obj$x) %*% W %*% V %*% obj$y)
      
      # residuals
      r <- obj$y - obj$x %*% z
      
      # new scale
      d <- median(abs(r - median(r))) / 0.44845
      
      # residual sum of squares
      rss2 <- Re(hConj(r) %*% U %*% r)
      
      H <- sqrt(U) %*% obj$x %*% solve(hConj(obj$x) %*% U %*% obj$x) %*% hConj(obj$x) %*% sqrt(U)
      h <- diag(H)
      
      Xi <- sum(diag(U))
    }
    # print(paste("Z: ", abs(z)))
    # chi <- chi - 0.0005 # I don't know what this should be ....
    alpha <- alpha - (1-alpha)/2
    chi <- qbeta(alpha, p, Xi.min - p)
    rss <- 3*rss2
  }
  
  z
}

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
#' Different from \link{tf} in that an FFT is used and additionally, the 
#' standard errors of the coefficient estimates are provided as well as the 
#' eigenvalues from the singular value decomposition.
#' 
#' @return An object of class \code{transfer}, consisting of a complex matrix whose 
#' columns are the individual transfer function for each input, and several attributes
#' describing the transfer function estimate.
#' 
#' @export
tf.svd <- function(x, y, blockSize = dim(x)[1], overlap = 0, deltat = 1, nw = 4
                    , k = 7, nFFT = NULL, freqRange = NULL, freqOffset = NULL
                    , standardize = TRUE, prewhiten = TRUE, removePeriodic = TRUE){
  # standardize the series by removing the mean and dividing by the standard deviation
  if( standardize ){
    stdPars <- vector( mode = "list" )
    stdPars$x <- data.frame( xmean = sapply( x, mean ), xsd = sapply( x, sd ) )
    stdPars$y <- data.frame( ymean = sapply( y, mean ), ysd = sapply( y, sd ) )
    std <- function( a ) (a - mean(a))/sd(a)
    x <- data.frame( lapply( x, std ) )
    y <- data.frame( lapply( y, std ) )
  }
  
  # I *think* we want to prewhiten based on the most data available?
  # (probably.. due to welch estimate, etc... )
  
  # block the data (x2, y2 are a list of data.frames)
  x2 <- sectionData(x, blockSize = blockSize, overlap = overlap)
  y2 <- sectionData(y, blockSize = blockSize, overlap = overlap)
  
  # determine an appropriate time vector (0 to T-dt)
  time <- seq(0, (blockSize - 1)*deltat, by = deltat)
  
  # number of frequencies bins to use (zero-pad length)
  if (is.null(nFFT)){
    nFFT <- 2^(floor(log2(blockSize)) + 3)
  }
  
  # determine the positive values of the frequencies.
  freq <- seq(0, 1/(2*deltat), by = 1/(nFFT*deltat))
  
  nfreq <- length(freq)
  
  # freqIdx at which frequencies of the response the transfer function should be 
  # estimated (you could not use all frequencies if offsets were also used)
  # NOT USED in this iteration of tf() ... 
  # freqIdx <- 1:length(freq)
  
  # taper the blocked data
  # x3 <- taper(x2, nw = nw, k = k)
  # y3 <- unlist(taper(y2, nw = nw, k = k), recursive = FALSE)
  
  # I don't think we actually need to store ALL of this (at least not the spec.mtm objects..)
  x.spec <- list()
  x.wtEigenCoef <- list()
  y.spec <- list()
  y.wtEigenCoef <- list()
  
  numSections <- attr(x2, "numSections")
  
  # used in the next for-loop ... 
  weightedEigen <- function(obj){ obj$mtm$eigenCoefs * obj$mtm$eigenCoefWt }
  
  for (i in 1:numSections){
    x.spec[[i]] <- lapply(x2[[i]], spec.mtm, deltat = deltat, dtUnits = "second", nw = nw
                          , k = k, nFFT = nFFT, plot = FALSE, returnInternals = TRUE)
    x.wtEigenCoef[[i]] <- lapply(x.spec[[i]], weightedEigen)
    y.spec[[i]] <- lapply(y2[[i]], spec.mtm, deltat = deltat, dtUnits = "second", nw = nw
                          , k = k, nFFT = nFFT, plot = FALSE, returnInternals = TRUE)
    y.wtEigenCoef[[i]] <- lapply(y.spec[[i]], weightedEigen)
  }
  
  
  # indexing helper function that grabs and stacks all the eigencoefficients at a frequency
  eigenByFreq <- function(obj, rowNum, numEl){
    matrix(unlist(lapply(obj, function(x, idx) x[idx, ], rowNum)), ncol = numEl)
  }
  
  
  x.design <- list()
  # y.design <- list()
  
  # stack the eigencoefficients by frequency from each block
  # form into a single list
  for (i in 1:nfreq){
    x.design[[i]] <- list(x = do.call(rbind, lapply(x.wtEigenCoef, eigenByFreq, rowNum = i, numEl = 3))
                          , y = do.call(rbind, lapply(y.wtEigenCoef, eigenByFreq, rowNum = i, numEl = 1)))
    # y.design[[i]] <- do.call(rbind, lapply(y.wtEigenCoef, eigenByFreq, rowNum = i, numEl = 1))
  }
  
  # This would be used if incorporating offset frequencies into the transfer function
  # regression.
  # NOT IMPLEMENTED YET (and not needed for the calling function)
  # if (is.null(freqOffset)){
  #   freqOffset = 0
  # }
  
  # perform the least squares complex-regresison.
  # H.mat <- olsTf(x = x3, y = y3, time = time, n = blockSize
  #                , npredictor = length(x3[[1]]), ntaper = k
  #                , freq = freq[freqIdx], fOffset = freqOffset)
  
  H.tmp <- lapply(x.design, function(obj){ svdRegression(obj$x, obj$y) })
  
  H.mat <- do.call(rbind, lapply(H.tmp, "[[", "coef"))
  # stdErr <- do.call(rbind, lapply(H.tmp, "[[", "stdErr"))
  # svdEigenval <- do.call(rbind, lapply(H.tmp, "[[", "ev"))
  
  H <- data.frame(freq, H.mat)
  colnames(H) <- c("freq", colnames(x))
  
  
  attr(H, "class") <- c(attr(H, "class"), "transfer")
  attr(H, "stdErr") <- do.call(rbind, lapply(H.tmp, "[[", "stdErr"))
  attr(H, "svdEv") <- do.call(rbind, lapply(H.tmp, "[[", "ev"))
  attr(H, "blockSize") <- attr(x2, "blockSize")
  attr(H, "overlap") <- attr(x2, "overlap")
  attr(H, "numSections") <- attr(x2, "numSections")
  attr(H, "n") <- attr(x2, "n")
  attr(H, "dataTaper") <- "Slepian"
  attr(H, "nw") <- nw
  attr(H, "k") <- k
  attr(H, "standardize") <- standardize
  if( standardize ) attr(H, "stdPars") <- stdPars
  
  # I made the standard errors and eigenvalues attributes on the transfer function
  # so that predict method wouldn't need any altering.  I don't know if this is the 
  # best way to do it though.
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
#' Defaults to \code{\link{trim}}.
#' @param sides Argument to \code{\link{filter}}: \code{sides = 1} is for filter coefficients
#' that apply to past values only (a causal filter); \code{sides = 2} is for filter coefficients
#' centered aroung lag 0.
#' @param returnInternals whether to return things computed during the prediction.
#' Currently, this returns the full vector(s) of filter coefficients as an attribute.
#' @param ... additional arguments passed to \code{filterMod}.
#' 
#' @details The transfer function estimate is used to calculate filter coefficients to be
#' applied in the time domain via convolution with the input series \code{newdata}.
#' Prior to the convolution, the filter coefficients can be modified using the
#' \code{filterMod} function. If \code{filterMod} produces a causal filter, ensure that
#' \code{sides = 1}; in any case, be sure that the output of \code{filterMod} conforms
#' to the requirements of \code{\link{filter}}.
#' 
#' The filter coefficients are put into vectors in the orientation expected by
#' \code{\link{filter}}. If \code{N} is the total length of a block, the causal
#' coefficients from \eqn{lag = 0} to \eqn{lag = (N-1)/2} are placed on
#' the "right", and the non-causal coefficients from \eqn{lag = -(N-1)/(2)}
#' to \eqn{lag = -1} are placed on the "left". This means that for a causal filter,
#' you would exclude filter coefficients with an index \emph{less} than N/2,
#' because there are an odd number of coefficients, with \eqn{lag = 0} in the middle,
#' and an equal number of coefficients on each side of zero.
#' 
#' @return A \code{data.frame} with the predicted values obtained by filtering 
#' the input series \code{newdata}.
#' 
#' @export

predict.transfer <- function( object, newdata, filterMod = trim, sides = 2, returnInternals = FALSE, ... ){
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
    stdPars <- as.data.frame( t( attr( object, "stdPars" )$x ) )[c("xmean","xsd"),]
    nm3 <- match( names( stdPars ), names( newdata[,nm2] ) ) # The positions of stdPars names in newdata names
    std <- function( x, sP ) (x - sP[1])/sP[2]
    newdata <- data.frame( mapply( FUN = std, x = newdata[,nm2], sP = stdPars[,nm3], SIMPLIFY = FALSE ) )
  }else{
    newdata <- newdata[,nm2]
  }
  
  ################################
  ### this is in impulseResponse()
  # Rearrange transfer function coefficients
  # objectC <- lapply( object, Conj )
  # attributes( objectC ) <- attributes( object )
  # objectC <- objectC[(nrow(object)-1):2, ]
  # objectFull <- rbind( object, objectC )[,nm]
  # 
  # # Inverse FFT of the transfer function coefficients
  # fC <- lapply( objectFull, function(x,...) fft(x,...)/length(x), inverse = TRUE )
  # 
  # # Rearrange the filter coefficients
  # fC <- as.data.frame( lapply( fC, Re ) )
  # 
  # # Only want n filter coefficients if n is odd, n-1 if n is even?
  # # We want an odd number of coefficients so we know that the middle one is lag 0
  # n2 <- ceiling( attr( object, "blockSize" )/2 )
  # # The causal part of the filter
  # ind <- 1:n2
  # # The non-causal part of the filter (has length n2-1)
  # ind <- c( (nrow(fC)-n2+2):nrow(fC), ind )
  # fC <- fC[ind,]
  ###
  ##############################
  
  fC <- impulseResponse(object) # if you uncomment the above block, delete this.
  
  # Apply the filterMod function to the coefficients
  fCMod <- data.frame( lapply( fC, filterMod, ... ) )
  
  # Compute prediction using filter
  out <- as.data.frame( mapply( filter, x = newdata, filter = fCMod, sides = sides ) )
  
  # Return prediction
  out <- data.frame( predict = rowSums(out) )
  if( attr( object, "standardize" ) ){
    yStdPars <- attr( object, "stdPars" )$y
    out$predict <- out$predict*yStdPars$ysd + yStdPars$ymean
  }
  if( returnInternals ){
    attr( out, "filterCoefs" ) <- fC
  }

out
}


#' Obtain the impusle response from a transfer object
#' 
#' Inverse Fourier transforms the transfer function to obtain the impulse response
#' 
#' @param object An object of class \code{transfer} obtained from \link{tf.svd} or \link{tf}.
#' 
#' @details Inverts the transfer function object returning causal and non-causal filter 
#' coefficients.  The non-causal filter coefficients are the left "half" of the vector with 
#' the causal parts being on the right (this is the form required for \link{filter}).
#' Aaron maybe want to make sure I didn't say this incorrectly O_O.
#' 
#' @return A vector containing the impulse response of the provided transfer function 
#' from \link{tf}.
#' 
#' @export

impulseResponse <- function(object, frequencyName = "freq"){
  tfNames <- names( object ) # don't actually need this I don't think
  nm <- tfNames[ tfNames != frequencyName ]
  
  # Rearrange transfer function coefficients
  objectC <- lapply( object, Conj )
  attributes( objectC ) <- attributes( object )
  objectC <- objectC[(nrow(object)-1):2, ]
  objectFull <- rbind( object, objectC )[,nm]
  
  # Inverse FFT of the transfer function coefficients
  if (is.data.frame(objectFull)){
    fC <- lapply( objectFull, function(x,...) fft(x,...)/length(x), inverse = TRUE )
  } else {
    fC <- fft(objectFull, inverse = TRUE) / length(objectFull)
  }
  
  # Rearrange the filter coefficients
  if (is.list(fC)){
    fC <- as.data.frame( lapply( fC, Re ) )
  } else {
    fC <- data.frame(Re(fC))
  }
  
  # Only want n filter coefficients if n is odd, n-1 if n is even?
  # We want an odd number of coefficients so we know that the middle one is lag 0
  n2 <- ceiling( attr( object, "blockSize" )/2 )
  # The causal part of the filter
  ind <- 1:n2
  # The non-causal part of the filter (has length n2-1)
  ind <- c( (nrow(fC)-n2+2):nrow(fC), ind )
  fC[ind,]
}


#' Trim an odd-length vector
#' 
#' This function trims a vector with an odd length, leaving \code{n} elements to
#' one or either side of the center element.
#' 
#' @param x The vector to trim
#' @param n The number of elements to one or either side of the center element to keep
#' @param LR A vector indicating which side to keep: To keep \code{n} elements to the left, 
#' set \code{LR = 1}; to keep \code{n} elements to the right, set \code{LR = 2}; for
#' \code{n} elements to the left and right, set \code{LR = 1:2}.
#' 
#' @return A trimmed \code{vector} according to \code{n} and \code{side}.
#' 
#' @export
trim <- function(x, n = 5, LR = 1:2){
  if( abs( length(LR)-1 ) > 1 | !all( LR %in% 1:2 ) ) stop( "Invalid LR argument" )
  l <- length(x)
  if( l%%2 != 1 ) stop( "Length of x is not odd" )
  m <- ceiling( l/2 )
  i <- list( (m-n):m, m:(m+n) )
x[ unique( unlist(i[LR]) ) ]
}
