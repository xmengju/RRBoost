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
# func.lad <- function(x, cc = NULL) {
#   return(mean(abs(x)))
# }

# lad loss
func.lad <- function(x, cc = NULL) {
  return(abs(x))
}



func.lad.grad<- function(x, cc = NULL) {
  return(sign(x))
}


# tukey's bisquare loss
# func.tukey <- function(r, cc= 4.685) {
#   w <- as.numeric(abs(r) <= cc)
#   v <- w*(1 - (1 - (r/cc)^2)^3)  +(1-w)*1
#   return(sum(v))
# }

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
# func.huber <- function(r, cc= 0.98) {
#   res <- r^2
#   res[abs(r) >cc] <- 2*cc*abs(r)[abs(r) >cc] - cc^2
#   return (sum(res) )
# }

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


# calculate efficiency  (for tukey's loss and Huber's loss)
cal.efficiency <- function(e, psi,psi.prime) {
  tmp_1 <- integrate(function(a, cc) (psi(a, cc)^2)*dnorm(a), cc=e, lower=-Inf, upper=+Inf)$value
  tmp_2 <- integrate(function(a, cc) psi.prime(a, cc)*dnorm(a), cc=e, lower=-Inf, upper=+Inf)$value

  return( 1/(tmp_1/tmp_2^2) )
}

# cc.huber <- uniroot( function(e) (cal.efficiency(e, func.huber.grad,func.huber.grad.prime)-.95), lower=0.1, upper=5)$root
# cc.tukey <- uniroot( function(e) (cal.efficiency(e, func.tukey.grad,func.tukey.grad.prime)-.95), lower=3, upper=5)$root
# cc.andrew <- uniroot( function(e) (cal.efficiency(e, func.andrew.grad,func.andrew.grad.prime)-.95), lower=1, upper=3.14)$root
# cc.welsch <- uniroot( function(e) (cal.efficiency(e, func.welsch.grad,func.welsch.grad.prime)-.95), lower=1, upper=3.14)$root


# calculate weights
cal.w <- function(x, psi, psi.prime, cc)
{
  res <- rep(psi.prime(0, cc= cc), length(x))
  res[x!=0] <- psi(x[x!=0],cc = cc)/x[x!=0]
  return (res)
}


Tukeys <- function(x, func.grad = func.tukey.grad, func.grad.prime = func.tukey.grad.prime, cc = 3.88, tol = 1e-8, max.it=50)
{
  mu.pre <- Inf
  mu.cur <- median(x)
  s.hat <- mad(x)
  it <- 0
  while( (abs(mu.pre -mu.cur) > tol) & (it < max.it) )
  {
    if(s.hat == 0) {s.hat = s.hat +10^{-7}}
    w <- cal.w((x - mu.cur)/s.hat, func.grad, func.grad.prime, cc)
    mu.pre <- mu.cur
    mu.cur <- sum(w*x)/sum(w)
    it <- it + 1
  }
  return(mu.cur)
}


Hubers <- function(x, func.huber.grad = func.huber.grad, func.huber.grad.prime = func.huber.grad.prime, cc = 0.98, tol = 1e-8, max.it=50, sigmam = "NA")
{
  mu.pre <- Inf
  mu.cur <- median(x)
  if(sigmam == "NA")
  { s.hat <- mad(x)
  }else{
    s.hat <- sigmam
  }

  it <- 0
  while( (abs(mu.pre -mu.cur) > tol) & (it < max.it) )
  {
    if(s.hat == 0) {s.hat = s.hat +10^{-7}}
    w <- cal.w((x - mu.cur)/(s.hat), func.huber.grad, func.huber.grad.prime, cc)
    mu.pre <- mu.cur
    mu.cur <- sum(w*x)/sum(w)
    it <- it + 1
  }
  return(mu.cur)
}


Andrews <- function(x, func.andrew.grad = func.andrew.grad, func.andrew.grad.prime = func.andrew.grad.prime, cc = 0.98, tol = 1e-8, max.it=50, sigmam = "NA")
{
  mu.pre <- Inf
  mu.cur <- median(x)
  if(sigmam == "NA")
  { s.hat <- mad(x)
  }else{
    s.hat <- sigmam
  }

  it <- 0
  while( (abs(mu.pre -mu.cur) > tol) & (it < max.it) )
  {
    if(s.hat == 0) {s.hat = s.hat +10^{-7}}
    w <- cal.w((x - mu.cur)/(s.hat), func.andrew.grad, func.andrew.grad.prime, cc)
    mu.pre <- mu.cur
    mu.cur <- sum(w*x)/sum(w)
    it <- it + 1
  }
  return(mu.cur)
}

Welschs <- function(x, func.welsch.grad = func.welsch.grad, func.welsch.grad.prime = func.welsch.grad.prime, cc = 0.98, tol = 1e-8, max.it=50, sigmam = "NA")
{
  mu.pre <- Inf
  mu.cur <- median(x)
  if(sigmam == "NA")
  { s.hat <- mad(x)
  }else{
    s.hat <- sigmam
  }

  it <- 0
  while( (abs(mu.pre -mu.cur) > tol) & (it < max.it) )
  {
    if(s.hat == 0) {s.hat = s.hat +10^{-7}}
    w <- cal.w((x - mu.cur)/s.hat, func.welsch.grad, func.welsch.grad.prime, cc)
    mu.pre <- mu.cur
    mu.cur <- sum(w*x)/sum(w)
    it <- it + 1
  }
  return(mu.cur)
}


## scale estimators
mscale <- function(u, delta=0.5, tuning.chi=1.547645, max.it=100, tol=1e-6) {
  # M-scale of a sample u
  # tol: accuracy
  # delta: breakdown point (right side)
  # Initial
  s0 <- median(abs(u))/.6745
  err <- tol + 1
  it <- 0
  while( (err > tol) && ( it < max.it) ) {
    it <- it+1
    s1 <- sqrt( s0^2 * mean(robustbase::Mchi(x=u/s0, cc = tuning.chi, psi='tukey')) / delta )
    err <- abs(s1-s0)/s0
    s0 <- s1
  }
  return(s0)
}

