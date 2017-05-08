## All code required for coherence

# x - dataset 1
# y - dataset 2
# xySectioned - data has already been blocked and treated.
### - x and y need to be the same length right now... otherwise this will not work.
# freqOffset - +/- freqOffset will be used to calc coherency
# slideFreq - how far to slide (in frequency) for each calculation
# slideIdx - how many frequency bins to move over to the next offset
## - set to 0 with to use W as the sliding amount
## - ignored if slideFreq is not NULL
# blockSize - how large the blocks should be
# offset for adjacent blocks (proportion = 0 < offset < 1)
# adaptive - whether to use adaptive weighting or not
# jackknife - use jackknife for variance estimates
# prewhiten - Prewhiten the blocks
# removePeriodic - removes line components (frequency domain approach)
# sigCutoff - what level of F-test to ues as the cutoff in the periodic removal
# nodes - how many cores for cluster use should this have
# clusterType - SOCK will be local, MPI would be for HPCVL

#' Calculates Frequency Domain Coherence
#' 
#' Estimates and returns the Loeve spectrum given series x and y.
#' 
#' @param x A \code{vector} of values for the first dataset.
#' @param y A \code{vector} of values for the second dataset - if null, assumes ... 
#' I'm not sure actually.
#' @param xySectioned A \code{logical} denoting whether the data has already been sectioned.
#' \code{FLASE} implies that the code will be sectioned according to \code{blockSize} and 
#' \code{overlap}.
#' @param n A \code{numeric} value giving the length of the \code{x} and \code{y} data sets. 
#' If \code{NULL}, will be determined in the function.
#' @param forwardCoh A \code{logical} value indicating if the forward coherence (\code{TRUE}) 
#' or reverse coherence \code{FALSE} should be computed.
#' @param freqRange A \code{vector} of two values containing the range of frequencies to use.
#' This range should be contained within [0, 1/(2*\code{dt})] (0 to Nyquist)
#' 
#' @return Returns a matrix containing the coherency.
#' 
#' @export
coherency.mtm <- function(x, y = NULL, xySectioned = FALSE, n = NULL, forwardCoh = TRUE
                          , freqRange = NULL, freqOffset = NULL, slideFreq = NULL, slideIdx = 1
                          , nw = NULL, k = NULL, deltat = 1, nFFT = NULL, nFFT.add = 3
                          , blockSize = NULL, overlap = NULL
                          , adaptive = TRUE, jackknife = FALSE, prewhiten = FALSE
                          , removePeriodic = TRUE, sigCutoff = NULL
                          , nodes = 1, clusterType = c("SOCK", "MPI")){
  ## Basic Setup ##
  
  if (!xySectioned){
    # not sure quite how different lengths will affect things
    nx <- length(x)
    ny <- length(y)
    
    # remove means, just in case.
    x <- x - mean(x)
    y <- y - mean(y)
  } else {
    blockSize <- length(x[[1]][[1]])
    if (is.null(n)){stop("Set n when xySectioned is TRUE.")}
    nx <- n
  }
  
  if (is.null(blockSize)){
    blockSize <- min(nx, ny)
  }
  
  if (is.null(nw)){
    nw <- 4
  }
  if (is.null(k)){
    k <- 7
  }
  if (is.null(nFFT)){
    nFFT <- 2^(floor(log2(blockSize)) + nFFT.add)
  }
  
  # how many indices to slide at each step.
  if (!is.null(slideFreq)){
    slideIdx <- round(slideFreq * deltat * nFFT)
  } else if (is.null(slideFreq) && slideIdx == 0){
    slideIdx <- floor((nw / (max(nx, ny) * deltat)) * (deltat * nFFT))
  }
  
  ## Calculate the Slepians and lambda's ##
  dpssIn <- dpss(n = blockSize, nw = nw, k = k)
  slep <- dpssIn$v
  ev <- dpssIn$eigen
  
  if (removePeriodic){
    Vk <- mvfft(rbind(slep, matrix(0, ncol = k, nrow = nFFT - blockSize)))
  }
  
  ## Don't need this right now, maybe in the future? Different blocksizes for each series? ... 
  # tmpY <- dpss(ny, nw = nw, k = k)
  # slepY <- tmpY$v
  # evY <- tmpY$eigen
  
  ## Sanity checks ##
  # Warns the user if linecomponents should be removed, but no F-test significance cutoff provided.
  if (removePeriodic && is.null(sigCutoff)){
    warning("Remove periodic components requested but no cutoff significance.  Set to (1 - 1/N).")
    sigCutoff <- 1 - 1/(min(c(length(x), length(y))))
  }
  
  # No frequency range means just a "simple" coherence at zero offset between x and y.
  if (is.null(freqRange)){
    cross <- crossSpectrum(x, y, slepX = slep, evX = ev, slepY = slep, evY = ev
                           , nFFT = nFFT, nw = nw, k = k, nFFT.add = nFFT.add
                           , adaptiveWeights = adaptive
                           , removePeriodic = removePeriodic, sigCutoff = sigCutoff
                           , returnInternals = TRUE)
    
    return(cross$Sxy / sqrt(cross$Sx * cross$Sy))
  }
  
  ## Setup Frequency Domain Indexing ##
  freq <- seq(0, 1/(2*deltat), by = 1/(deltat*nFFT))
  
  freqIdxStart <- head(which(freq >= freqRange[1]), 1)
  freqIdxEnd <- tail(which(freq <= freqRange[2]), 1)
  freqRangeIdx <- freqIdxStart:freqIdxEnd
  
  # number of indices to offset by
  numOffsetIdx <- floor(freqOffset * deltat * nFFT)
  if (slideIdx == 0){
    slideTmp <- 0
  } else {
    slideTmp <- c(rev(seq(-slideIdx, -numOffsetIdx, by = -slideIdx)), 0, seq(slideIdx, numOffsetIdx, by = slideIdx))
  }

  if (freqIdxStart - numOffsetIdx < 1 || freqIdxEnd + numOffsetIdx > (nFFT/2+1)){
    stop("Check freqRange and freqOffset - not enough frequency bins to make this combination work.")
  }
  
  ## Setup Time Domain Indexing ##
  incr <- ceiling(blockSize * (1-overlap))
  sectInd <- seq(1, nx-blockSize+1, by=incr)
  numSect <- length(sectInd)
  
  # create a list of vectors of indices - each element of the list contains the indices of the 
  ## block that will be used
  blockIdxLst <- as.list(as.data.frame(mapply(":", sectInd, sectInd+blockSize-1)))
  
  # might already be sectioned, probably not though...
  if (xySectioned){
    dat <- list()
    for (i in 1:length(x[[1]])){
      dat[[i]] <- list(x = x[[1]][[i]], y = x[[2]][[i]])
    }
  } else {
    # create a list of the data to use for parrallel
    dat <- list()
    for (i in 1:length(blockIdxLst)){
      dat[[i]] <- list(x = x[blockIdxLst[[i]]], y = y[blockIdxLst[[i]]])
    }
  }
  
  range2Start <- freqIdxStart + slideTmp
  range2End <- range2Start + length(freqRangeIdx) - 1
  # yIdxRange <- as.list(as.data.frame(mapply(":", range2Start, range2End))) # remove this - use C++ code...
  
  if (nodes <= 1){
    ## do something here - no parallelization required
    cohere <- list()
    for (i in 1:length(blockIdxLst)){
      cohere[[i]] <- coherencyHelper(dat[[i]], dpssIn = dpssIn, nw = nw, k = k, nFFT = nFFT
                                     , freqRangeIdx = freqRangeIdx, range2Start = range2Start
                                     , forwardCoh = forwardCoh
                                     , prewhiten = prewhiten, adaptive = adaptive
                                     , removePeriodic = removePeriodic, sigCutoff = sigCutoff
                                     , Vk = Vk)
    }
  } else {
    if (clusterType == "SOCK"){
      cl <- snow::makeCluster(rep("localhost", nodes), type = "SOCK")
      snow::clusterSetupRNG(cl)
      cohere <- snow::clusterApply(cl, dat, coherencyHelper
                                   , dpssIn = dpssIn, nw = nw, k = k, nFFT = nFFT
                                   , freqRangeIdx = freqRangeIdx, range2Start = range2Start
                                   , forwardCoh = forwardCoh
                                   , prewhiten = prewhiten, adaptive = adaptive,
                                   , removePeriodic = removePeriodic, sigCutoff = sigCutoff
                                   , Vk = Vk)
      snow::stopCluster(cl)
    } else if (clusterType == "MPI"){
      
    }
  }
  
  # msc <- matrix(0, nrow = dim(cohere[[1]])[1], ncol = dim(cohere[[1]])[2])
  # for (i in 1:length(cohere)){
  #   msc = msc + abs(cohere[[i]])^2
  # }
  # 
  # msc2 <- apply(msc, 1, mean)
  
  list(coherency = cohere, forwardCoh = forwardCoh
       , freqRange = freqRange, freqRangeIdx = freqRangeIdx, freqOffset = freqOffset
       , freqOffsetIdx = numOffsetIdx, slideFreq = slideFreq, slideIdx = slideIdx
       , nw = nw, k = k, deltat = deltat, nFFT = nFFT, nFFT.add = nFFT.add
       , blockSize = blockSize, overlap = overlap
       , adaptive = adaptive, jackknife = jackknife, prewhiten = prewhiten
       , removePeriodic = removePeriodic, sigCutoff = sigCutoff)
}

# should not be called by user.
coherencyHelper <- function(datLst, dpssIn, nw, k, nFFT, freqRangeIdx, range2Start
                            , forwardCoh
                            , prewhiten, adaptive, removePeriodic, sigCutoff, Vk){
  n <- length(datLst$x)
  
  if (removePeriodic){
    yk.x <- removeLineCoh(spec.mtm(datLst$x, nw = nw, k = k, nFFT = nFFT, Ftest = TRUE, plot = FALSE
                                   , dpssIn = dpssIn, returnInternals = TRUE)
                   , sigCutoff = sigCutoff, Vk = Vk, k = k, nFFT = nFFT)
    yk.y <- removeLineCoh(spec.mtm(datLst$y, nw = nw, k = k, nFFT = nFFT, Ftest = TRUE, plot = FALSE
                                    , dpssIn = dpssIn, returnInternals = TRUE)
                           , sigCutoff = sigCutoff, Vk = Vk, k = k, nFFT = nFFT)
  } else {
    yk.x <- spec.mtm(datLst$x, nw = nw, k = k, nFFT = nFFT, Ftest = TRUE, plot = FALSE
                     , dpssIn = dpssIn, returnInternals = TRUE)$mtm$eigenCoefs
    yk.y <- spec.mtm(datLst$y, nw = nw, k = k, nFFT = nFFT, Ftest = TRUE, plot = FALSE
                     , dpssIn = dpssIn, returnInternals = TRUE)$mtm$eigenCoefs
  }
  
  if (prewhiten){
    # implement this... for god's sakes man!
  }
  
  if (adaptive){
    dx <- adaptiveWeightsCpp(eigenSpec = abs(yk.x)^2, ev = dpssIn$eigen, k = k, nFFT = nFFT)
    dy <- adaptiveWeightsCpp(eigenSpec = abs(yk.y)^2, ev = dpssIn$eigen, k = k, nFFT = nFFT)
  }
  
  coherencyOffsetCpp(ykx = yk.x, yky = yk.y, dx = dx, dy = dy, nTaper = k
                     , nFreqRange = length(freqRangeIdx), nOffset = length(range2Start)
                     , range2Start = range2Start, freqRangeIdx = freqRangeIdx
                     , forwardCoh = forwardCoh)
  
}

# Don't think I need this - .coherency helper will do this..
# .crossSpectrum <- function(x, y, slep, ev, nFFT, nw, k, adaptive = TRUE){
#   
# }

removeLineCoh <- function(spec, sigCutoff, Vk, k, nFFT){
  fSigIdx <- findLocalFMax(spec, cutoff = sigCutoff)
  infl <- matrix(0, nrow = nFFT/2+1, ncol = k)
  
  if (length(fSigIdx) > 0){
    for (i in 1:k){
      for (j in 1:length(fSigIdx)){
        # store the 
        ### Conjugate this because I need to use the correct side of the tapered data
        infl[, i] <- infl[, i] + taRifx::shift(spec$mtm$cmv[fSigIdx[j]] * Vk[, i], -fSigIdx[j]+1)[1:(nFFT/2+1)]
      }
    }
  }
  
  eigencoef <- matrix(0, nrow = nFFT/2+1, ncol = k)
  for (i in 1:k){
    eigencoef[, i] <- spec$mtm$eigenCoefs[, i] - infl[, i]
  }
  
  eigencoef
}

# cohLst : result from coherency.mtm()
### /// fix this, should be for a specific offset.. silly David.
cohPhaseByFreq <- function(coh, offsetIdx){
  phse <- matrix(0, ncol = dim(coh$coherency[[1]])[2], nrow = length(offsetIdx))
  
  for (i in 1:length(coh$coherency)){
    phse = phse + atan2(Im(coh$coherency[[i]][offsetIdx, ]), Re(coh$coherency[[i]][offsetIdx, ]))
  }
  
  phse <- phse / length(coh$coherency)
}

#' Obtain marginal MSC
#' 
#' Maybe ... 
#' 
#' @export
mscFromCoh <- function(coh, jackknife = FALSE){
  msc <- rep(0, dim(coh$coherency[[1]])[1])
  jkVar <- rep(NA, length(msc))
    
  for (i in 1:length(coh$coherency)){
    msc <- msc + apply(abs(coh$coherency[[i]])^2, 1, mean)
  }
  
  if (jackknife){
    jkloo <- matrix(0, nrow = dim(coh$coherency[[1]])[1], ncol = dim(coh$coherency[[1]])[2])
    for (i in 1:dim(coh$coherency[[1]])[2]){
      jkloo[, i] <- apply(abs(coh$coherency[[1]][, -i])^2, 1, mean)
    }
    
    jkSlashDot <- apply(jkloo, 1, mean)
    jkVar <- (dim(jkloo)[2] - 1) * apply(matrix(mapply("-", jkloo, jkSlashDot)^2, ncol = dim(jkloo)[2])
                   , 1, mean)
  }
  
  list(msc = msc / length(coh$coherency), jkVar = jkVar)
}


#' Frequency axis calculation
#' 
#' Uses number of frequency bins and deltat values to create frequency axis.
#' 
#' @param nFFT A \code{numeric} giving total number of frequency bins.
#' @param freqOffset A \code{numeric} giving the maximum offset used.
#' @param deltat A \code{numeric} giving the sampling rate of the time series.
#' 
#' @return A \code{vector} of frequency values with range (0, nyquist) inclusive.
#' @export
freqOffsetAxis <- function(nFFT, freqOffset, deltat){
  c(rev(seq(0, -freqOffset, by = -1/(deltat * nFFT))), seq(0, freqOffset, by = 1/(deltat*nFFT))[-1])
}

offsetIdxFromFreq <- function(f, deltat, nFFT){
  f * deltat * nFFT
}

crossSpectrum <- function(x, y, slepX=NULL, evX = NULL, slepY=NULL, evY = NULL
                          , nFFT=NULL, nw=NULL, k=NULL, nFFT.add=3
                          , adaptiveWeights = TRUE, removePeriodic = TRUE, sigCutoff = 0.99
                          , returnInternals = TRUE){
  if (is.null(nFFT)){
    # set the number of frequency bins using the longer series
    nFFT <- 2^(floor(log2(max(length(x), length(y))) + nFFT.add))
  }
  if (is.null(nw) || is.null(k)){
    warning("Either nw or k not set, using nw=4, k=7.")
    nw <- 4; k <- 7;
  }
  
  ## Calculate the Slepians for x and y series independently if required.
  # assuems you will pass in slepians for X, but not Y, OR both at once.
  # X == null and Y != null should not occur. (... fix this )
  if (is.null(slepX) || is.null(evX)){
    require('multitaper')
    if (length(x) != length(y)){
      tmpX <- dpss(n=length(x), k=k, nw=nw)
      slepX <- tmpX$v
      evX <- tmpX$eigen
      tmpY <- dpss(n=length(y), k=k, nw=nw)
      slepY <- tmpY$v
      evY <- tmpY$eigen
    } else {
      tmp <- dpss(n=length(x), k=k, nw=nw)
      slepX <- tmp$v
      evX <- tmp$eigen
      slepY <- slepX
      evY <- evX
    }
  } else if (is.null(slepY) || is.null(evY)){
    if (length(x) != length(y)){
      tmpY <- dpss(n=length(y), k=k, nw=nw)
      slepY <- tmpY$v
      evY <- tmpY$eigen
    } else {
      slepY <- slepX
      evY <- evX
    }
  }
  
  vx <- matrix(0, ncol=k, nrow=nFFT)
  vy <- matrix(0, ncol=k, nrow=nFFT)
  
  vx[1:length(x), ] <- apply(slepX, MARGIN = 2, "*", x)
  vy[1:length(y), ] <- apply(slepY, MARGIN = 2, "*", y)
  
  if (removePeriodic){
    yk.x <- removeLine(x, nw = nw, k = k, nFFT = nFFT, sigCutoff = sigCutoff)$yk
    yk.y <- removeLine(y, nw = nw, k = k, nFFT = nFFT, sigCutoff = sigCutoff)$yk
  } else {
    yk.x <- mvfft(vx)
    yk.y <- mvfft(vy)
  }
  
  if (adaptiveWeights){
    # have to recalculate the weights due to taking out the periodic line components
    ## Bias should drop (you'd hope... )
    if (returnInternals){
      tmpX <- adaptiveWeights(eigenSpec = abs(yk.x)^2, ev = evX, weightsOnly = FALSE)
      d.x <- tmpX$weights
      Sx <- tmpX$adaptSpec
      tmpY <- adaptiveWeights(eigenSpec = abs(yk.y)^2, ev = evY, weightsOnly = FALSE)
      d.y <- tmpY$weights
      Sy <- tmpY$adaptSpec
    } else {
      d.x <- adaptiveWeights(eigenSpec = abs(yk.x)^2, ev = evX, weightsOnly = TRUE)$weights
      d.y <- adaptiveWeights(eigenSpec = abs(yk.y)^2, ev = evY, weightsOnly = TRUE)$weights
    }
    
    sk.xy <- yk.x * Conj(yk.y)
    s.xy <- apply(abs(d.x*d.y) * sk.xy, 1, sum) / apply(abs(d.x*d.y), 1, sum)
  } else {
    # Equally weighted estimates of spectra (why the hell would you use this... )
    d.x <- 1
    d.y <- 1
    Sx <- (1/k) * apply(abs(yk.x)^2, MARGIN = 1, sum)
    Sy <- (1/k) * apply(abs(yk.y)^2, MARGIN = 1, sum)
    # cross eigenspectra:
    sk.xy <- yk.x * Conj(yk.y)
    s.xy <- (1/k) * apply(sk.xy, 1, sum)
  }
  
  if (returnInternals){
    list(Sxy = s.xy, Sx = Sx, Sy = Sy, yk.x = yk.x, yk.y = yk.y
         , dk.x = d.x, dk.y = d.y, nFFT = nFFT, nw = nw, k = k)
  } else {
    s.xy
  }
  
  # Re(s.xy) # I don't think you want the real part of this... (Past David wrote this.. )
  ## This gets used in MSC where |s.xy|^2 is the numerator... soooo
}


## DTJ's formula (2nd at top of 1088 in 1982 paper)
# x - data or time-series - assumes ZERO-MEAN
# maybe add this argument, xIsEigencoef = FALSE, later?
removeLine <- function(x, nw = 4, k = 7, nFFT = NULL, nFFT.add = 5, sigCutoff = 0.99
                       , returnInfluenceMat = TRUE, reshape = FALSE){
  n <- length(x)
  if (is.null(nFFT)){
    nFFT <- 2^(floor(log2(length(x))) + nFFT.add)
  }
  
  Vk.tmp <- matrix(0, nrow = nFFT, ncol = k)
  dpss.tmp <- dpss(n = n, nw = nw, k = k)
  Vk.tmp[1:n, ] <- dpss.tmp$v
  ev <- dpss.tmp$eigen
  Vk <- mvfft(Vk.tmp)
  
  # Vk.centered <- Vk[c((nFFT/2+2):nFFT, 1:(nFFT/2+1)), ]
  
  spec <- spec.mtm(x, nw = nw, k = k, nFFT = nFFT, plot=FALSE, returnInternals = TRUE, Ftest=TRUE)
  
  fSigIdx <- findLocalFMax(spec, cutoff = sigCutoff)
  infl <- matrix(0, nrow = nFFT/2+1, ncol = k)
  
  if (length(fSigIdx) > 0){
    inflReshape <- list() # i-th element corresponds to i-th taper
    for (i in 1:k){
      inflReshape[[i]] <- matrix(0, nrow = nFFT/2+1, ncol = length(fSigIdx))
      for (j in 1:length(fSigIdx)){
        # store the 
        ### Conjugate this because I need to use the correct side of the tapered data
        # inflReshape[[i]][, j] <- shift(spec$mtm$cmv[fSigIdx[j]] * Vk.centered[, i], -fSigIdx[j]+1)[(nFFT/2):nFFT]
        inflReshape[[i]][, j] <- shift(spec$mtm$cmv[fSigIdx[j]] * Vk[, i], -fSigIdx[j]+1)[1:(nFFT/2+1)]
        infl[, i] <- infl[, i] + inflReshape[[i]][, j]
      }
    }
  }
  
  eigencoef <- matrix(0, nrow = nFFT/2+1, ncol = k)
  eigencoefFull <- matrix(0, nrow = nFFT, ncol = k)
  
  
  for (i in 1:k){
    eigencoef[, i] <- spec$mtm$eigenCoefs[, i] - infl[, i]
    eigencoefFull[, i] <- c(eigencoef[, i], Conj(rev(eigencoef[-c(1,nFFT/2+1), i])))
  }
  
  if (reshape){
    s.noLine <- adaptiveWeights(abs(eigencoefFull)^2, ev = spec$mtm$dpss$eigen) # noise spectrum
    
    # calculate CMV's to use in SD(f) calculation
    # mu <- cmvSimple(eigencoef, nw = nw, k = k, n = n)
    
    # I need help with this... will ask DJT tomorrow.
    # 1 / (2*pi*N) * sqrt(6 / ( (0.25) * mu^2 * n^3))
    
    
    
    # sum up the inflReshape elements to get the variance of the line component
    # spread this across a df bandwidth of ... (6*s.noLine$adaptSpec[sigFreqIdx[i]])/(pi*A^2*T^3)
    
    # varF <- VARIANCE OF THE FREQUENCIES
  }
  
  if (returnInfluenceMat){
    list(yk = eigencoefFull, sigFreqIdx = fSigIdx, nw = spec$mtm$nw, k = spec$mtm$k
         , nFFT = spec$mtm$nFFT, ev = ev)
  } else {
    eigencoefFull
  }
}

# checks for whether the significant point in the F-test is 
# larger than the previous and next point.  If yes, we have 
# a local max.
# - returns these frequencies
# obj needs to be either a) a spec.mtm object (Ftest class) or 
# b) a driegert.cmv object
findLocalFMax <- function(obj, cutoff){
  # Check whether this is a spec.mtm() object, or from my own CMV code.
  if (any(class(obj) == "Ftest")){
    Fval <- obj$mtm$Ftest
    k <- obj$mtm$k
  } else if (any(class(obj) == "driegert.cmv")){
    Fval <- obj$Ftest$Fval
    k <- obj$k
  } else {
    stop("obj needs to be of class 'driegert.cmv' or 'Ftest'.")
  }
  
  # Thanks to Frank Marshall for catching the missing "-2"
  fMaxInd <- which(Fval > qf(cutoff, 2, 2*k-2))
  maxes <- c()
  
  if (length(fMaxInd) == 0){
    return(maxes)
  }
  
  for (i in 1:length(fMaxInd)){
    if (fMaxInd[i] == 1 || fMaxInd[i] == length(Fval)){
      next
    }
    
    if (Fval[fMaxInd[i]] > Fval[fMaxInd[i]-1] && 
        Fval[fMaxInd[i]] > Fval[fMaxInd[i]+1]){
      maxes <- c(maxes, fMaxInd[i])
    }
  }
  
  maxes
}


# Adaptive MTM
## Full frequency array for eigenSpec unless halfFreqArray = TRUE
adaptiveWeights <- function(eigenSpec, ev, halfFreqArray = FALSE, weightsOnly = FALSE){
  s.hat2 <- 0.5 * (eigenSpec[, 1] + eigenSpec[, 2])
  s.hat <- 2*s.hat2
  
  # first estimate of the variance
  sigmasq <- sum(s.hat2) * (1/length(s.hat2))
  
  if (halfFreqArray){
    sigmasq <- 2 * sigmasq
  }
  
  B <- (1 - ev) * sigmasq
  
  d <- matrix(0, nrow = dim(eigenSpec)[1], ncol = dim(eigenSpec)[2])
  
  newVar = TRUE
  
  while (mean(abs(s.hat2 - s.hat)/s.hat2) > 0.0001){
    s.hat <- s.hat2
    
    for (i in 1:dim(eigenSpec)[2]){
      d[, i] <- sqrt(ev[i]) * s.hat / (ev[i] * s.hat + B[i])
    }
    
    s.hat2 <- apply(abs(d)^2 * eigenSpec, 1, sum) / apply(abs(d)^2, 1, sum)
    
    if (newVar){
      sigmasq <- sum(s.hat2) * (1/length(s.hat2))
      if (halfFreqArray){
        sigmasq <- 2 * sigmasq
      }
      B <- (1 - ev) * sigmasq
      newVar <- FALSE
    }
  }
  
  if (weightsOnly){
    list(adaptSpec = NULL, weights = d)
  } else {
    list(adaptSpec = s.hat2, weights = d)
  }
}

###### ^^^ This is all old code written ... Summer 2016-ish ^^^ ######
# Starting over now... see if we can do this a smidge more efficiently

#' Calculates the coherence between two series
#' 
#' Estimates the frequency domain coherence using the multitaper method.
#' @param x A \code{data.frame} whose columns are the time domain input series.
#' @param y A \code{numeric} vector containing the response series.
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
#' @param forward Indicates whether the forward (TRUE) or reverse (FALSE)
#' coherence should be calculated.
#' @param average An \code{integer} representing how the average across blocks 
#' should be calculated;
#' 0 - no average, return all the individual block information; 
#' 1 - average the cross and auto-spectra across blocks, then calculate the coherency
#' 2 - estimate the coherency for each block, average the coherencey across blocks
#' 3 - estimate the MSC for each block, average the MSC across blocks.
#' @param freqRange A \code{numeric} vector containing two elements with the start 
#' and end location for the band on which to estimate the coherence.
#' @param maxFreqOffset A \code{numeric} value indicating the maximum offset coherence to 
#' calculate in the specified band.
#' @param prewhiten NOT YET IMPLEMENTED
#' @param removePeriodic NOT YET IMPLEMENTED
#' @param sigCutoff NOT YET IMPLEMENTED
#' 
#' @details MSC stands for Magnitude Squared Coherence.
#' 
#' @export
coherence <- function(x, y = NULL, blockSize = length(x), overlap = 0, deltat = 1
                      , nw = 4, k = 7, nFFT = NULL, forward = TRUE
                      , average = 1, msc = FALSE
                      , freqRange = NULL, maxFreqOffset = NULL
                      , prewhiten = FALSE, removePeriodic = TRUE, sigCutoff = NULL)
{
  if (!any(average == 0:3)){ 
    stop("average must have an integer value between 0 and 3.")
  }
  
  if (is.null(freqRange) && !is.null(maxFreqOffset)){
    stop("maxFreqOffset implies that freqRange should be assigned.")
  }
  
  # number of frequencies bins to use (zero-pad length)
  if (is.null(nFFT)){
    nFFT <- 2^(floor(log2(blockSize)) + 3)
  }
  
  # determine the positive values of the frequencies.
  freq <- seq(0, 1/(2*deltat), by = 1/(nFFT*deltat))
  nfreq <- length(freq)
  
  # block the data (x2, y2 are a list of data.frames)
  if (is.null(y)){
    x2 <- sectionData(data.frame(x = x[, 1], y = x[, 1]), blockSize = blockSize, overlap = overlap)
  } else {
    x2 <- sectionData(data.frame(x = x[, 1], y = y[, 1]), blockSize = blockSize, overlap = overlap)
  }
  
  numSections <- attr(x2, "numSections")
  
  # wtEigenCoef <- blockedEigenCoef(x2, deltat = deltat, nw = nw, k = k
  #                                 , nFFT = nFFT, numSections = numSections
  #                                 , adaptiveWeighting = TRUE, returnWeights = FALSE)
  # 
  # if (is.null(freqRange)){
  #   spec <- lapply(wtEigenCoef, calculateSpec, forward = forward)
  #   freqIdx <- NULL
  # } else {
  #   freqIdx <- head(which(freq >= (freqRange[1] - maxFreqOffset) ), 1):tail(which(freq <= (freqRange[2] + maxFreqOffset)), 1)
  #   spec <- lapply(wtEigenCoef, calculateSpec, forward = forward, idx = freqIdx)
  # }
  
  
  if (is.null(freqRange)){
    freqIdx <- NULL
  } else {
    freqIdx <- head(which(freq >= (freqRange[1] - maxFreqOffset) ), 1):tail(which(freq <= (freqRange[2] + maxFreqOffset)), 1)
  }
  
  # spec <- lapply(wtEigenCoef, calculateSpec, forward = forward, idx = freqIdx)
  spec <- blockedSpec(x2, deltat = deltat, nw = nw, k = k, nFFT = nFFT
                      , numSections = numSections, adaptiveWeighting = TRUE, forward = TRUE
                      , idx = freqIdx)
  
  
  subFreq <- freq[freqIdx]
  info = list(allFreq = freq, blockSize = blockSize, overlap = overlap
              , numSections = numSections
              , deltat = deltat, nw = nw, k = k, nFFT = nFFT
              , forward = forward, average = average, msc = msc
              , freqRange = freqRange
              , freqRangeIdx = c(head(which(freq >= freqRange[1]), 1), tail(which(freq <= freqRange[2]), 1))
              , maxFreqOffset = maxFreqOffset
              , maxFreqOffsetIdx = tail(which(freq <= maxFreqOffset), 1) - 1 #-1 due to 0 frequency
              , freqIdx = freqIdx, prewhiten = prewhiten
              , removePeriodic = removePeriodic, sigCutoff = sigCutoff)
  
  if (average == 0){
    info$msc <- FALSE
    return(list(freq = subFreq, coh = spec, info = info))
  } else if (average == 1) {
    Sxy.ave <- Reduce('+', lapply(spec, "[[", "Sxy")) / numSections
    Sxx.ave <- Reduce('+', lapply(spec, "[[", "Sxx")) / numSections
    Syy.ave <- Reduce('+', lapply(spec, "[[", "Syy")) / numSections
    coh <- Sxy.ave / sqrt(Sxx.ave %*% t(Syy.ave))
  } else if (average == 2) {
    coh <- Reduce('+', lapply(spec, function(obj){ obj$Sxy / sqrt(obj$Sxx %*% t(obj$Syy)) })) / numSections
  } else if (average == 3) {
    coh <- Reduce('+', lapply(spec, function(obj){ abs(obj$Sxy)^2 / (obj$Sxx %*% t(obj$Syy)) })) / numSections
    info$msc <- TRUE
    return(list(freq = subFreq, coh = coh, info = info))
  }
  
  if (msc){
    list(freq = subFreq, coh = abs(coh)^2, info = info)
  } else {
    list(freq = subFreq, coh = coh, info = info)
  }
}

#' Convert frequency vs. frequency coherence to offset vs frequency
#' 
#' Converts from a "diagonal" coherence matrix to a "horizontal" coherence matrix.
#' 
#' @param obj An object as returned by \link{coherence}.
#' 
#' @details Takes the rows of the diagonal matrix centered at the frequencies 
#' in \code{freqRange} with \code{maxFreqOffset} values on either side.  These rows then 
#' become the columns of the returned offset coherence matrix.
#' 
#' The offset is from the frequency in x
#' 
#' @return A list containing the centre frequencies, offset values, offset coherences, and 
#' all the extra information returned by \link{coherence}.
#' 
#' @export
coherenceToOffset <- function(obj){
  startIdx <- obj$info$maxFreqOffsetIdx + 1
  endIdx <- length(obj$freq) - obj$info$maxFreqOffsetIdx
  
  sectIdx <- seq(startIdx - obj$info$maxFreqOffsetIdx, endIdx - obj$info$maxFreqOffsetIdx, by=1)
  
  blockIdxLst <- as.list(as.data.frame(mapply(":", sectIdx, sectIdx+(2*obj$info$maxFreqOffsetIdx))))
  
  offsetCoh <- unname(simplify2array(lapply(blockIdxLst, function(x) { obj$coh[x[1] + obj$info$maxFreqOffsetIdx, x] })))
  
  list(freq = obj$freq[startIdx:endIdx]
       , offset = freqOffsetAxis(obj$info$nFFT, obj$info$maxFreqOffset, obj$info$deltat)
       , offsetCoh = offsetCoh, info = obj$info)
}

#' Averages across the rows of the offset coherence matrix
#' 
#' Estimates the average coherence for each offset within a frequency band.
#' 
#' @param obj A \code{list} as returned by \link{coherenceToOffset}.
#' 
#' @detail None yet.. 
#' 
#' @export
marginalCoh <- function(obj, msc = TRUE){
  if (obj$info$msc){
    mcoh <- apply(obj$offsetCoh, 1, mean)
  } else if (msc == TRUE) {
    mcoh <- apply(abs(obj$offsetCoh)^2, 1, mean)
    obj$info$msc = TRUE
  } else {
    mcoh <- apply(obj$offsetCoh, 1, mean)
  }
  
  list(offset = obj$offset, marginalCoh = mcoh, info = obj$info)
}