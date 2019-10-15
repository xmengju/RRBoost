temp1 <- function(y, wt, parms) {
  wmedian <- median(y*wt)
  rabs <- sum((y - wmedian)^2*wt)
  list(label= wmedian, deviance=rabs)
}

temp2 <- function(y, wt, x, parms, continuous) {
  
  n <- length(y)
  goodness <- NULL
  
  if(continuous){
    for(i in 1:(n-1)){
      err_l <- sum(abs(y[1:i] - median(y[1:i])))
      err_r <- sum(abs(y[(i+1):n] - median(y[(i+1):n])))
      goodness <- c(goodness, -err_l - err_r)
    }
    goodness <- goodness + max(abs(goodness))
    list(goodness= goodness, direction= sign(goodness))
    
  }else{
    ux <- sort(unique(x))
    n <- length(ux)
    
    goodness <- NULL
    for(i in 1:(n-1)){
      err_l <- sum(abs(y[x == ux[i]] - median(y[x == ux[i]])))
      err_r <- sum(abs(y[x != ux[i]] - median(y[x != ux[i]])))
      goodness <- c(goodness, -err_l - err_r)
    }
    
    list(goodness= goodness, direction= ux)
    
  }
}


temp3 <- function(y, offset, parms, wt) {
  if (!is.null(offset)) y <- y-offset
  list(y=y, parms=0, numresp=1, numy=1,
       summary= function(yval, dev, wt, ylevel, digits ) {
         paste("  median=", format(signif(yval, digits)),
               ", MAD=" , format(signif(dev/wt, digits)),
               sep='')
       })
}

alist <- list(eval=temp1, split=temp2, init=temp3)


temp4 <- function(y, wt, x, parms, continuous) {
  # Center y
  n <- length(y)
  y <- y- sum(y*wt)/sum(wt)
  
  if (continuous) {
    # continuous x variable
    temp <- cumsum(y*wt)[-n]
    
    left.wt  <- cumsum(wt)[-n]
    right.wt <- sum(wt) - left.wt
    lmean <- temp/left.wt
    rmean <- -temp/right.wt
    goodness <- (left.wt*lmean^2 + right.wt*rmean^2)/sum(wt*y^2)
    list(goodness= goodness, direction=sign(lmean))
  }
  else {
    # Categorical X variable
    ux <- sort(unique(x))
    wtsum <- tapply(wt, x, sum)
    ysum  <- tapply(y*wt, x, sum)
    means <- ysum/wtsum
    
    # For anova splits, we can order the categories by their means
    #  then use the same code as for a non-categorical
    ord <- order(means)
    n <- length(ord)
    temp <- cumsum(ysum[ord])[-n]
    left.wt  <- cumsum(wtsum[ord])[-n]
    right.wt <- sum(wt) - left.wt
    lmean <- temp/left.wt
    rmean <- -temp/right.wt
    list(goodness= (left.wt*lmean^2 + right.wt*rmean^2)/sum(wt*y^2),
         direction = ux[ord])
  }
}

# The init function:
#   fix up y to deal with offsets
#   return a dummy parms list
#   numresp is the number of values produced by the eval routine's "label"
#   numy is the number of columns for y
#   summary is a function used to print one line in summary.rpart
# In general, this function would also check for bad data, see rpart.poisson
#   for instace.
temp5 <- function(y, offset, parms, wt) {
  if (!is.null(offset)) y <- y-offset
  list(y=y, parms=0, numresp=1, numy=1,
       summary= function(yval, dev, wt, ylevel, digits ) {
         paste("  mean=", format(signif(yval, digits)),
               ", MSE=" , format(signif(dev/wt, digits)),
               sep='')
       })
}

