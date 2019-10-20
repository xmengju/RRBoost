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

