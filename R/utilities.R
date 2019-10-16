# initialize boosting parameters
init.boost <- function(type)
{
  switch (type,
          square = {
            func <- func.square
            func.grad <- func.square.grad
            func.grad.prime <- func.square.grad.prime
          },
          lad = {
            func <- func.lad
            func.grad <- func.lad.grad
            func.grad.prime <- 0
          },
          huber = {
            func <- func.huber
            func.grad <- func.huber.grad
            func.grad.prime <- func.huber.grad.prime
          },
          tukey = {
              func <- func.tukey
              func.grad <- func.tukey.grad
              func.grad.prime <- func.tukey.grad.prime
          }
        )

  return (list(func = func, func.grad = func.grad, func.grad.prime = func.grad.prime))
}



newton.search <- function( f_t_train, h_train, y_train, func, func.grad, func.grad.prime, ss, cc, min_sigma = FALSE,  alpha_init = c(0), Tol = 10^{-10}, T_newton = 10, bb = 0.5)
{
  tmp = c(NA, length(alpha_init))
  for(k in 1:length(alpha_init)) {
    tryCatch({
      alpha_t <- alpha_init[k]
      for(i in 1:T_newton) {
        r_t <- f_t_train + alpha_t*h_train - y_train
        if(min_sigma == FALSE) {
          f_prime_alpha <- sum(h_train*func.grad(r_t/ss, cc = cc))/ss
          f_prime_prime_alpha <- sum(h_train^2*func.grad.prime(r_t/ss, cc = cc))/(ss^2)
        }else{
          ss = mscale(r_t, tuning.chi = cc, delta = bb)
          f_prime_alpha = sum(func.grad(r_t/ss, cc = cc) *h_train * ss)/sum(func.grad(r_t/ss, cc = cc)*r_t)
          D  =  func.tukey.grad(r_t/ss, cc = cc)
          C = func.tukey.grad.prime(r_t/ss, cc = cc) *(h_train * ss - r_t*f_prime_alpha)/(ss^2)
          A <- sum(C*h_train*ss + D*h_train*f_prime_alpha)
          B <- sum(C*r_t + D*h_train)
          f_prime_prime_alpha <- (A*(sum(D*r_t)) - B*(sum(D*h_train*ss)))/(sum(D*r_t)^2)
        }

        if(f_prime_alpha == 0){
          alpha_t <- alpha_t
        }else{
          alpha_t <- alpha_t - f_prime_alpha/(f_prime_prime_alpha + 10^{-10})
        }

      }

      if(min_sigma == TRUE) {
        if (f_prime_prime_alpha >= 0){
        tmp[k] <-alpha_t
        }
      }else{
        if (f_prime_prime_alpha >= 0){
        tmp[k] <- alpha_t
        }
      }
    })

  }
    return(tmp[min(which(abs(tmp) == min(abs(tmp), na.rm = TRUE)))])
}


# bisection search, this is a recursive algorithm
bisection.search  <- function(f_t_train, h_train, y_train, func_line, func_line.grad, a1, a2, Tol = 10^{-6},step_num = NULL,  k = 0, min_sigma = FALSE)
{
  # [a1, a2] is the interval to search form
  c <-  (a1 + a2)/2
  if(min_sigma == FALSE){
   grad <-  func_line.grad((f_t_train + c*h_train - y_train))*(h_train)
   #print(head(grad))
  }else{
   grad <- -func_line.grad(y_train-f_t_train - c*h_train)*(h_train)
  }


  if(abs(a1 - a2) < Tol) {
    return(c)
  }

  if(length(step_num) == 0){
      if(sum(grad) < 0) {
        return(bisection.search(f_t_train, h_train, y_train,func_line, func_line.grad, c, a2,step_num = step_num,  min_sigma =  min_sigma))
      }
      if(sum(grad) >= 0) {
        return(bisection.search(f_t_train, h_train, y_train,func_line, func_line.grad, a1, c,step_num = step_num,  min_sigma =  min_sigma))
      }
  }else{
    if(k >= step_num){
      return(c)
    }
    if(sum(grad) < 0) {
      return(bisection.search(f_t_train, h_train, y_train,func_line, func_line.grad, c, a2,step_num = step_num,k = k+1,  min_sigma =  min_sigma))
    }
    if(sum(grad) >= 0) {
      return(bisection.search(f_t_train, h_train, y_train,func_line, func_line.grad, a1, c,step_num = step_num,k = k+1,  min_sigma =  min_sigma))
    }
  }
}

secant.method <- function(f, x0, x1, tol = 1e-7, n = 50) {
  if(x1 == Inf| x1 == -Inf){
    return (0)
  }
  for (i in 1:n) {
    x2 <- x1 - f(x1) / ((f(x1) - f(x0)) / (x1 - x0)) # Calculate the new x value
    #print(c(x2,x1,x0,f(x1), f(x1) - f(x0)))
    if(((f(x1) - f(x0)) == 0) | ((x1 - x0) == 0)){
      return(x1)
    }
    if (abs(x2 - x1) < tol) { # If the difference between the new value and the previous value is small enough, end iteration and output root.
      #print("terminate")
      return(x2)
    }
    # If the root was not determined in the previous iteration, update the values and proceed to the next iteration.
    x0 <- x1
    x1 <- x2
  }
  return(x2)
}


cal.tau <- function(u, cc_2){
  cc_1 <- RobStatTM::lmrobdet.control(bb=.5, family='bisquare')$tuning.chi
  s =  mscale(u,  tuning.chi= cc_1, delta = 0.5)
  return((s^2/length(u))* sum(unlist(lapply(u, function(x){func.tukey(x/s, cc = cc_2)}))))
}

rmse <- function(x){

  return(sqrt(mean(x^2)))
}

trmse <- function(trim_prop = NULL, trim_c = NULL, x){
  if(length(trim_c) != 0){
    idx <- (x < (median(x) + trim_c*mad(x))) & x > (median(x) - trim_c*mad(x))
  }
  if(length(trim_prop) != 0){
    idx <- (x < quantile(x, 1- trim_prop/2)) & (x > quantile(x, trim_prop/2))
  }
  return(list(trmse = rmse(x[idx]), idx = idx))
}


