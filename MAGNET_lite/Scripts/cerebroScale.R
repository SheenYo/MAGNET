cerebroScale2<-
  function (x, clamp, divData, center_zero) 
  {
    xmed <- median(x, na.rm = TRUE)
    xmad <- mad(x, constant = 1, na.rm = TRUE)
    xmin <- min(x, na.rm = TRUE)
    xmax <- max(x, na.rm = TRUE)
    avoidClamp <- max(abs(xmed - xmin), abs(xmed - xmax))/xmad
    fillMatrix <- x
    if (is.null(clamp)) {
      clamp <- avoidClamp + 1
    }
    outlrs <- clamp * xmad
    if (clamp <= 0) 
      stop("clamp must be >0")
    pctOL <- round(length(which(x[!is.na(x)] <= (xmed - (outlrs)) | 
                                  x[!is.na(x)] >= (xmed + (outlrs))))/length(x[!is.na(x)]) * 
                     100, 2)
    if (pctOL > 0) {
      warning(paste("The clamp value of ", clamp, " will clamp ", 
                    pctOL, "% of input values (outliers) to the min or max of the scaled range.", 
                    sep = ""))
    }
    if (divData == TRUE & center_zero == FALSE) {
      abvMed <- x[x >= xmed & x <= (xmed + outlrs) & !is.na(x)]
      belMed <- x[x <= xmed & x >= (xmed - outlrs) & !is.na(x)]
      if (length(which(!is.na(x)))%%2 == 0) {
        rightsc <- rescale(c(xmed, abvMed), c(0.5, 1))[-1]
        fillMatrix[x >= xmed & x <= (xmed + outlrs) & !is.na(x)] <- rightsc
        leftsc <- rescale(c(xmed, belMed), c(-1, 0))[-1]
        fillMatrix[x <= xmed & x >= (xmed - outlrs) & !is.na(x)] <- leftsc
        fillMatrix[x < (xmed - outlrs) & !is.na(x)] <- -1
        fillMatrix[x > (xmed + outlrs) & !is.na(x)] <- 1
        xScaled <- fillMatrix
      }
      if ((length(which(!is.na(x))))%%2 == 1) {
        rightsc <- rescale(abvMed, c(-1, 1))
        fillMatrix[x >= xmed & x <= (xmed + outlrs) & !is.na(x)] <- rightsc
        leftsc <- rescale(belMed, c(-1, 0))
        fillMatrix[x <= xmed & x >= (xmed - outlrs) & !is.na(x)] <- leftsc
        fillMatrix[x < (xmed - outlrs) & !is.na(x)] <- -1
        fillMatrix[x > (xmed + outlrs) & !is.na(x)] <- 1
        xScaled <- fillMatrix
      }
    }
    if (divData == TRUE & center_zero == TRUE) {
      abvMed <- x[x >= -1 & x <= (-1 + outlrs) & !is.na(x)]
      belMed <- x[x <= -1 & x >= (-1 - outlrs) & !is.na(x)]
      if (length(which(!is.na(x)))%%2 == 0) {
        rightsc <- rescale(c(-1, abvMed), c(-1, 1))[-1]
        fillMatrix[x >= -1 & x <= (-1 + outlrs) & !is.na(x)] <- rightsc
        leftsc <- rescale(c(-1, belMed), c(-1, 0))[-1]
        fillMatrix[x <= -1 & x >= (-1 - outlrs) & !is.na(x)] <- leftsc
        fillMatrix[x < (-1 - outlrs) & !is.na(x)] <- -1
        fillMatrix[x > (-1 + outlrs) & !is.na(x)] <- 1
        xScaled <- fillMatrix
      }
      if ((length(which(!is.na(x))))%%2 == 1) {
        rightsc <- rescale(abvMed, c(-1, 1))
        fillMatrix[x >= -1 & x <= (-1 + outlrs) & !is.na(x)] <- rightsc
        leftsc <- rescale(belMed, c(-1, 1))
        fillMatrix[x <= -1 & x >= (-1 - outlrs) & !is.na(x)] <- leftsc
        fillMatrix[x < (-1 - outlrs) & !is.na(x)] <- -1
        fillMatrix[x > (-1 + outlrs) & !is.na(x)] <- 1
        xScaled <- fillMatrix
      }
    }
    if (divData == FALSE) {
      nonoutlrs <- x[x >= (xmed - outlrs) & x <= (xmed + outlrs) & 
                       !is.na(x)]
      xsc <- rescale(nonoutlrs, c(-1, 1))
      fillMatrix[x >= (xmed - outlrs) & x <= (xmed + outlrs) & 
                   !is.na(x)] <- xsc
      fillMatrix[x < (xmed - outlrs) & !is.na(x)] <- -1
      fillMatrix[x > (xmed + outlrs) & !is.na(x)] <- 1
      xScaled <- fillMatrix
    }
    return(xScaled)
  }
