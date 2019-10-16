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
