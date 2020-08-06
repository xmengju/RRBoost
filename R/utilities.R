#' @import rpart stats
#'
# #' @importFrom utils head
# #'
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


rmse <- function(x){

  return(sqrt(mean(x^2)))
}

trmse <- function(trim_prop = NULL, trim_c = NULL, x){
  if(length(trim_c) != 0){
    idx <- (x < (median(x) + trim_c*mad(x))) & x > (median(x) - trim_c*mad(x))
  }else{
    if(length(trim_prop) != 0){
      idx <-  (abs(x) < quantile(abs(x), 1-trim_prop))
    }
  }
  return(list(trmse = rmse(x[idx]), idx = idx))
}
