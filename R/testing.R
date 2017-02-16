# testing
load("~/school_lab/contracts/health_canada2017/assets/data/edaData.RData")
x <- data.frame(no2 = no2.dly, o3 = o3.dly, pm25 = pm25.dly)
y <- data.frame(mort = mort$Mort.CP.b.A0)
x2 <- sectionData(x, blockSize = 750, overlap = 0.25)
y2 <- sectionData(y, blockSize = 750, overlap = 0.25)
x3 <- taper(x2)
y3 <- taper(y2)

deltat <- 24*3600
time <- seq(0, (attr(x2, "blockSize") - 1)*deltat, by = deltat)
ntaper <- attr(x3, "k")
blockSize <- attr(x3, "blockSize")
M <- 2^(floor(log2(blockSize)) + 3)
freq <- seq(0, 1/(2*deltat), by = 1/(M*deltat))
freqIdx <- 1:length(freq)

y4 <- unlist(y3, recursive = FALSE)

H.old <- olsTf(x = x3, y = y4, time = time, n = blockSize
            , npredictor = length(x3[[1]]), ntaper = ntaper
            , freq = freq[freqIdx], fOffset = 0)

H <- tf(x, y, blockSize = 750, overlap = 0.25, deltat = 24*3600, nw = 6, k = 10)
