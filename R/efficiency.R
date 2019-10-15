## calculate the efficiency for the tau-estimator  
func.tukey.tmp <- function(r, cc= 4.685) {
  w <- as.numeric(abs(r) <= cc)
  v <- w*(1 - (1 - (r/cc)^2)^3)  +(1-w)*1
  return(v)
}


cal_efficiency<- function(cc_2){
  cc_1 <- RobStatTM::lmrobdet.control(bb=.5, family='bisquare')$tuning.chi
  bb <- 0.5
  
  sigma0 <- 1
  term_1 <- 2*integrate(function(a) {func.tukey.tmp(a/sigma0, cc = cc_2)*dnorm(a)},lower=-Inf, upper=+Inf)$value
  term_2 <-  integrate(function(a) {func.tukey.grad(a/sigma0, cc = cc_2)* (a/sigma0)*dnorm(a)},lower=-Inf, upper=+Inf)$value
  term_3 <-  integrate(function(a) {func.tukey.grad(a/sigma0, cc = cc_1)* (a/sigma0)*dnorm(a)},lower=-Inf, upper=+Inf)$value
  
  W0 <- (term_1 - term_2)/term_3
  
  psi0<- function(u, cc_2){
    W0 * func.tukey.grad(u/sigma0, cc = cc_1) + func.tukey.grad(u/sigma0, cc = cc_2)
  }
  psi0.prime <- function(u, cc_2){
    W0 * func.tukey.grad.prime(u/sigma0, cc = cc_1) /sigma0 + func.tukey.grad.prime(u/sigma0, cc = cc_2)/sigma0
  }
  
  tmp_1 <- integrate(function(a, cc) (psi0(a, cc)^2)*dnorm(a), cc = cc_2, lower=-Inf, upper=+Inf)$value
  tmp_2 <- integrate(function(a, cc) psi0.prime(a,cc)*dnorm(a),  cc = cc_2, lower=-Inf, upper=+Inf)$value
  
  return(1/(tmp_1/tmp_2^2))
}


# ###  calculate the scale_M estimator 
# cal_scale_M <- function (u, delta = 0.5, func, tuning.chi = 1.547645, max.it = 100, tol = 1e-06) 
# {
#   s0 <- median(abs(u))/0.6745
#   err <- tol + 1
#   it <- 0
#   while ((err > tol) && (it < max.it)) {
#     it <- it + 1
#     s1 <- sqrt(s0^2 * mean(unlist(lapply(u/s0, FUN = function(x) func(x, cc = tuning.chi))))/delta)
#     err <- abs(s1 - s0)/s0
#     s0 <- s1
#   }
#   return(s0)
# }

# u = c(1:100)
# cal_scale_M(u, delta = 0.5, func.tukey, tuning.chi = 1.547645, max.it = 100, tol = 1e-06) 
# mscale(u)
# 
# cal_scale_M(u, delta = 0.203, func.rho, tuning.chi = 1.56, max.it = 100, tol = 1e-06) 
# 
# 


# func.rho <- function(t, cc){
#     tmp = (abs(t)< cc)
#     return(tmp*((t^2/2)*(1 - t^2/cc^2 + t^4/(3*cc^4))) + (1-tmp)*cc^2/6)
# 
# }
# 
# func.psi <- function(t, cc){
#   tmp = (abs(t)< cc)
#   return(tmp*t*(1 - t^2/cc^2)^2 + (1-tmp)*0)
# }
# 
# 
# func.psi.prime <- function(t, cc){
#   func.tukey.grad.prime(t,cc) * (cc^2/6)
# }
# # 
# 
# cal.efficiency <- function(e, psi,psi.prime) {
#   tmp_1 <- integrate(function(a, cc) (psi(a, cc)^2)*dnorm(a), cc=e, lower=-Inf, upper=+Inf)$value
#   tmp_2 <- integrate(function(a, cc) psi.prime(a, cc)*dnorm(a), cc=e, lower=-Inf, upper=+Inf)$value
#   
#   return( 1/(tmp_1/tmp_2^2) )
# }
# 
# cc <- uniroot( function(e) (cal.efficiency(e, func.psi,func.psi.prime)-.5), lower=0.1, upper=5)$root


# cal_efficiency_original <- function(cc_2){
#   cc_1 <- 1.56
#   #cc_2 <- 6.08
#   sigma0 <- 1
#   term_1 <- 2*integrate(function(a) {func.rho(a/sigma0, cc = cc_2)*dnorm(a)},lower=-Inf, upper=+Inf)$value
#   term_2 <-  integrate(function(a) {func.psi(a/sigma0, cc = cc_2)* (a/sigma0)*dnorm(a)},lower=-Inf, upper=+Inf)$value
#   term_3 <-  integrate(function(a) {func.psi(a/sigma0, cc = cc_1)* (a/sigma0)*dnorm(a)},lower=-Inf, upper=+Inf)$value
#   
#   W0 <- (term_1 - term_2)/term_3
#   psi0<- function(u, cc_2){
#     W0 * func.psi(u/sigma0, cc = cc_1) + func.psi(u/sigma0, cc = cc_2)
#   }
#   psi0.prime <- function(u, cc_2){
#     W0 * func.psi.prime(u/sigma0, cc = cc_1) /sigma0 + func.psi.prime(u/sigma0, cc = cc_2)/sigma0
#   }
#   tmp_1 <- integrate(function(a, cc) (psi0(a, cc)^2)*dnorm(a), cc = cc_2, lower=-Inf, upper=+Inf)$value
#   tmp_2 <- integrate(function(a, cc) psi0.prime(a,cc)*dnorm(a),  cc = cc_2, lower=-Inf, upper=+Inf)$value
#   
#   return(1/(tmp_1/tmp_2^2))
# }
#   
# tmp <- cal_efficiency_original(u)
# psi0 <- tmp[[1]]
# psi0.prime <- tmp[[2]]
# 
# tmp_1 <- integrate(function(a) (psi0(a)^2)*dnorm(a),  lower=-Inf, upper=+Inf)$value
# tmp_2 <- integrate(function(a) psi0.prime(a)*dnorm(a),  lower=-Inf, upper=+Inf)$value
# 1/(tmp_1/tmp_2^2) 



# cc_2 <- uniroot(function(e) (cal_efficiency_mine(e)-.95), lower=0.1, upper=10)$root
# tmp <- cal_efficiency_original(u)
# psi0 <- tmp[[1]]
# psi0.prime <- tmp[[2]]
# 
# cal_efficiency function()
# 
# 1/(tmp_1/tmp_2^2) 

