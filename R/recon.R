## Frequency domain reconstruction:

#' Testing the fit with offset coherences
#' 
#' Frequency domain exploration of fit.
#' 
#' @param newData a \code{data.frame} containing the same columns that were used in the 
#' original transfer function fitting.
#' @param H an object of type \code{transfer} from tf().
#' @param nw time-bandwidth parameter (needed for window function estimation) - same value
#' as used to estimate \code{H}.
#' @param k an \code{integer} indicating the number of tapers used in estimation of 
#' \code{H} - although maybe they don't need to be the same for this part... 
#' 
#' @export
specOffRecon <- function(newData, H, nw, k, nFFT, deltat){
  dnames <- names(newData)
  dnames <- dnames[!(dnames == "freq")]
  spec <- list()
  for (i in 1:length(dnames)){
    spec[[ dnames[i] ]] <- spec.mtm(newData[, dnames[i] ], nw = nw, k = k
                                    , deltat = deltat, dtUnits = "second"
                                    , nFFT = nFFT, plot = FALSE)$spec
  }
  
  res <- rep(0, nFFT/2+1)
  for (i in 1:length(H)){
    cnames <- colnames(H[[i]]$coef)
    
    for (j in 1:length(cnames)){
      offNames <- unlist( strsplit( cnames[j], split = "..", fixed = TRUE ))
      if (length(offNames == 1)){
        res[i] <- res[i] + (abs(H[[i]]$coef[, offNames[1] ])^2)*spec[[ offNames[1] ]][i]
      } else {
        if (offNames[1] == "m"){ offSign <- -1 } else { offSign <- 1 }
        offIdx <- offSign * as.numeric(substr(offNames[2], 2, nchar(offNames[2])))
        if (offidx <= 0){ offIdx <- abs(offIdx) + 2 }
        res[i] <- res[i] + (abs(H[[i]]$coef[, offNames[1] ])^2) * spec[[ offNames[1] ]][offIdx]
      }
    }
  }
  
  freq <- seq(0, 1/(2*deltat), by = 1/(nFFT*deltat))
  cbind(freq, spec = res)
}


#' Reconstructs a response based on input tranfser functions
#' 
#' Based on multiple inputs and offset frequency coherence relationships
#' 
#' @param H A \code{transfer} object as returned by tf().
#' @param newData A \code{data.frame} containing the columns with the same name as those of \code{H}.
#' @param hTrim An \code{integer} indicating how points to keep on either side of the center (0 lag) point.
#' @param sides A value of 1 (causal) or 2 (non-causal)
#' 
#' @details This should be incorporated into predict.transfer() at some point.  But for now, ... 
#' 
#' @export
offsetRecon <- function(H, newData, hTrim = 5, sides = 2){
  
  parts <- list()
  M <- attr(H, "nFFT") #2^(floor(log2(dim(x)[1])) + 3)
  N <- dim(newData)[1]
  h.tmp <- impulseResponse(H, realPart = FALSE)
  
  if (sides == 2){
    h <- as.data.frame(lapply(h.tmp, trim, n = hTrim))
  } else if (sides == 1){
    h <- as.data.frame(lapply(h.tmp, trim, n = hTrim, LR = 2))
  } else {
    stop("Not a possible value of 'sides' argument.  1 or 2 are the only allowed values.")
  }
  
  nFilt <- dim(h)[1]
  for (i in 1:length(h)){
    offNames <- unlist( strsplit( names(h)[i], split = "..", fixed = TRUE ))
    if (length(offNames) == 1){
      parts[[i]] <- zFilter( newData[, offNames[1] ], filter = h[, i] )
    } else {
      pm <- substr(offNames[2], 1, 1)
      offIdx <- as.numeric( substr( offNames[2], 2, nchar(offNames[2]) ))
      eSgn <- 1
      if (pm == "m"){
        offIdx <- -offIdx
        eSgn <- -1
      }
      
      parts[[i]] <- zFilter( newData[, offNames[1]] * exp(-eSgn*2*pi*(offIdx / M))
                             , filter = h[, i] )
    }
  }
  
  recon <- rep(0, N)
  for (i in 1:length(parts)){
    recon <- recon + Re(parts[[i]])
  }
  
  
  ### FIX THIS!  FFT convolution works differently in terms of lags than filter() function. .. 
  if (sides == 2){
    recon[hTrim+1:N]
  } else {
    recon ### I don't know if this one is correct ... 
  }
  
}