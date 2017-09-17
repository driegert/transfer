## code to testerate ...

### Is the design matrix code giving me what I think it's giving me??

design.test <- function(xd, x.wt){
  problems <- FALSE
  for (i in 1:length(xd)){
    xx <- xd[[i]]$x
    cnames <- colnames(xx)
    splt <- strsplit(cnames, split = "..", fixed = TRUE)
    
    for (j in 1:length(splt)){
      if (length(splt[[j]]) == 1){
        diffie <- xx[, splt[[j]]] - c(x.wt[[1]][[ splt[[j]][1] ]][i, ]
                            , x.wt[[2]][[ splt[[j]][1] ]][i, ])
      } else {
        pm <- substr(splt[[j]][2], 1, 1)
        off <- as.numeric(substr(splt[[j]][2], 2, nchar(splt[[j]][2])))
        
        if (pm == "m"){ idx <- i - off } else {idx <- i + off }
        
        if (idx <= 0){
          idx <- abs(idx) + 2
          diffie <- xx[, paste(splt[[j]], collapse = "..")] - Conj(c(x.wt[[1]][[ splt[[j]][1] ]][idx, ]
                                                           , x.wt[[2]][[ splt[[j]][1] ]][idx, ]))
        } else {
          diffie <- xx[, paste(splt[[j]], collapse = "..")] - c(x.wt[[1]][[ splt[[j]][1] ]][idx, ]
                                                      , x.wt[[2]][[ splt[[j]][1] ]][idx, ])
        }
      }
      
      if (any(abs(diffie) != 0)){
        print(paste0("Frequency: ", i, " -- Column: ", j, " -- Name: ", cnames[j]))
        problems <- TRUE
      }
    }
  }
  if (problems){
    print("Something doesn't match.")
  } else {
    print("Everything matches!!")
  }
}

