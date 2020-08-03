### this file saves loss functions
### square loss, tukey's loss, huber's loss
## hampel's three part loss function
## add visualization of these loss functions

# square loss
# func.square <- function(x, cc = NULL) {
#   return(sum((x)^2)/2)
# }

func.square <- function(x, cc = NULL) {
  return((x)^2/2)
}


func.square.grad<- function(x, cc = NULL) {
  return(x)
}

func.square.grad.prime<- function(x, cc = NULL) {
  return(1)
}


# lad loss
func.lad <- function(x, cc = NULL) {
  return(abs(x))
}



func.lad.grad<- function(x, cc = NULL) {
  return(sign(x))
}


# tukey's bisquare loss
func.tukey <- function(r, cc= 4.685) {
  w <- as.numeric(abs(r) <= cc)
  v <- w*(1 - (1 - (r/cc)^2)^3)  +(1-w)*1
  return(v)
}


func.tukey.grad <- function(r, cc = 4.685) {
  w <- as.numeric(abs(r) <= cc )
  gg <- w*6*r*(1 - (r/cc)^2)^2/(cc^2)  +(1-w)*0
  return(gg)
}

func.tukey.grad.prime <- function(r, cc = 4.685) {
  w <- as.numeric(abs(r) <= cc )
  tmp <- (1 - (r/cc)^2)^2 - r*2*(1 - (r/cc)^2)*2*r/(cc^2)
  tmp <- (tmp*w + (1-w)*0)*6/cc^2
  return(tmp)
}

# Huber's loss
func.huber <- function(r, cc= 0.98) {
  res <- r^2
  res[abs(r) >cc] <- 2*cc*abs(r)[abs(r) >cc] - cc^2
  return (res)
}

func.huber.grad <- function(r, cc = 0.98) {
  res <- r
  res[abs(r) > cc] = sign(r)[abs(r) > cc]*cc
  return(res)
}

func.huber.grad.prime <- function(r, cc = 0.98) {
  return(0 + (abs(r) <=cc)*1)
}

# ##
# func.andrew <- function(r, cc= 1.339) {
#   res <- sum(sin(r[abs(r) <= pi*cc]/(2*cc))^2) + sum((abs(r) > pi*cc))
#   return (res )
# }
#
# func.andrew.grad <- function(r, cc = 1.339) {
#   res <- 0 + 0.5*sin(r*(abs(r) <= cc*pi)/cc)/cc
#   return(res)
# }
#
# func.andrew.grad.prime <- function(r, cc = 1.339) {
#   return(0 + 0.5*(abs(r)<=cc*pi)*cos(r/cc)*(1/cc^2))
# }
#
#
# func.welsch <- function(r, cc= 2.11) {
#
#   res <- sum(1 - exp(- (r/cc)^2/2))
#   return (res )
# }
#
# func.welsch.grad <- function(r, cc = 2.11) {
#   res <- r*exp(-(r/cc)^2/2)/(cc^2)
#   return(res)
# }
#
# func.welsch.grad.prime <- function(r, cc = 2.11) {
#
#   return((1 - (r/cc)^2)*exp(-(r/cc)^2/2)/cc^2)
# }
#

