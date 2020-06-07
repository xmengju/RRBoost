#' Boost.control
#'
#' Control for Boost fits
#'
#' Various parameters that control aspects of the `Boost` fit.
#'
#' @param n_init number of iterations for the 1st stage of RRBoost ($T_{1,max}$) (int)
#' @param cc_s  tuning constant of tukey's loss in SBoost (numeric)
#' @param eff_s  normal efficiency of tukey's loss in RRBoost (2nd stage) (numeric)
#' @param bb  breakdown point of the SBoost estimator (numeric)
#' one per explanatory variable.
#' @param trim_prop  trimming proportion if `trmse` is used as the performance metric (numeric)
#' @param trim_c the trimming constant if `trmse` is used as the performance metric (numeric)
#' @param max_depth_init the maximum depth of the initial LADTtree  (numeric, default 3)
#' @param min_leaf_size_init the minimum number of observations per node of the initial LADTtree (numeric)
#' @param cal_imp calculate variable importance  (TRUE or FALSE)
#' @param save_f save the function estimates at all iterations (TRUE or FALSE)
#' @param make_prediction make predictions using x_test  (TRUE or FALSE)
#' @param save_tree save trees at all iterations  (TRUE or FALSE)
#' @param save_all_err_rr save validation and test error of RRBoost trained with all combination of LADTree parameters (TRUE or FALSE)
#' @param shrinkage shrinkage parameter in boosting (numeric)
#' @return A list of all input parameters
#'
#' @author Xiaomeng Ju, \email{xmengju@stat.ubc.ca}
#'
#' @export
#' 
Boost.control <- function(n_init = 100,  cc_s  = NULL,  eff_m= NULL, bb = 0.5, trim_prop = NULL, trim_c = 3, max_depth_init = 3, min_leaf_size_init = 10, cal_imp = TRUE, save_f = FALSE, make_prediction = TRUE, save_tree = FALSE, precision = 4, save_all_err_rr = TRUE, shrinkage = 1, save_time= TRUE){

  if(length(cc_s) == 0){
    cc_s <- as.numeric(RobStatTM::lmrobdet.control(bb=.5, family='bisquare')$tuning.chi)
  }

  if(length(eff_m) == 0){
    eff_m <- 0.95
  }
  cc_m <-  as.numeric(RobStatTM::lmrobdet.control(efficiency=eff_m, family='bisquare')$tuning.psi)

  if(missing(save_all_err_rr)){
    save_all_err_rr <- TRUE
  }
  if(missing(save_time)){
    save_time <- TRUE
  }
  return(list(n_init = n_init,  cc_s = cc_s, cc_m = cc_m, bb = bb, trim_prop = trim_prop, trim_c = trim_c, max_depth_init = max_depth_init, min_leaf_size_init = min_leaf_size_init, cal_imp = cal_imp,  save_f = save_f, make_prediction = make_prediction, save_tree = save_tree, precision = precision,  save_all_err_rr =  save_all_err_rr, shrinkage = shrinkage, save_time = save_time))
}

init.boosting <- function(type)
{
  switch (type,
          L2Boost = {
          init_functions <- init.boost("square")
          },
          Robloss = {
            init_functions <- init.boost("huber")
          },
          LADBoost = {
            init_functions <- init.boost("lad")
          },
          MBoost = {
            init_functions <- init.boost("huber")
          },
          SBoost = {
            init_functions <- init.boost("tukey")
          },
          RRBoost = {
            init_functions <- init.boost("tukey")
          }
  )

  return(init_functions)
}


cal.neggrad <- function(type, x_train, y_train, f_t_train, init_status, ss, func, func.grad, cc)
{
  dat_tmp <- data.frame(x_train);

  if(type == "Robloss"){
    dat_tmp$neg_grad <-  - func.grad(f_t_train - y_train, cc = ss)
  }

  if(type == "MBoost") {
    dat_tmp$neg_grad <-  - func.grad(f_t_train - y_train, cc = ss)
  }

  if(type %in% c("LADBoost","L2Boost")) {
    dat_tmp$neg_grad <- - func.grad(f_t_train - y_train)
  }
  if(type == "SBoost" | (type == "RRBoost" & init_status == 0)) {
    tmp_numerator <- -func.grad((f_t_train - y_train)/ss, cc = cc)/ss
    tmp_denominator <- -sum(func.grad((f_t_train - y_train)/ss, cc = cc)*(y_train - f_t_train))/(ss^2)
    dat_tmp$neg_grad <- tmp_numerator/(tmp_denominator)
  }
  if(type == "RRBoost" & init_status == 1)
  {
    dat_tmp$neg_grad <- - func.grad((f_t_train - y_train)/ss, cc = cc)/ss
  }

  return(dat_tmp)
}


cal.ss <- function(type, f_t_train, y_train,  cc, bb) {

  ss <-1
  if(type == "Robloss"){
    ss <-  mad(f_t_train -y_train)
  }

  if(type %in% c("SBoost", "RRBoost")) {
      ss <- mscale(f_t_train - y_train,  tuning.chi=cc, delta = bb)
  }

  if(type == "MBoost") {
    ss  <- quantile(abs(f_t_train - y_train),0.9)
  }

  return(ss)
}



cal.alpha <- function(type,  f_t_train, h_train, y_train, func, ss, init_status, cc = 1.547) {
  switch(type,
         LADBoost = {
           ff = function(a,r,h){
             return(mean(func(r - a*h)))
           }
           return(optimize(ff, lower = -1, upper = 300, r = y_train - f_t_train, h = h_train)$minimum)
         },
         
         SBoost = {
           ff <- function(a, r, h) return(mscale(r - a*h))
           upper_region = c(0.5,10,100,300)
           tmp <-  tmp_val <- rep(NA, length(upper_region))
           for(i in 1:length(upper_region)){
             val = optimize(ff, lower = -1, upper = upper_region[i], r = y_train - f_t_train, h = h_train)
             tmp[i] <- val$minimum
             tmp_val[i] <- val$objective
            } 
          
           idx <- min(which(tmp_val == min(tmp_val)))
           order_val <- order(head(tmp_val, idx)) 
           if( sum(order_val == idx:1) == idx){
             return(tmp[idx])
           }else{
               tmp_order <- order_val  - c(max(order_val), head(order_val, length(order_val)-1))
               if(sum(tmp_order > 0) > 0){
                 tmp_idx <- min(which(tmp_order>0))-1
                 return(tmp[tmp_idx])
               }else{
                 return(tmp[1])
               }
           }
         },
         RRBoost = {
           if(init_status == 0) {
             ff <- function(a, r, h) return(mscale(r - a*h))
             upper_region = c(0.5,10,100,300)
             tmp <-  tmp_val <- rep(NA, length(upper_region))
             for(i in 1:length(upper_region)){
               val = optimize(ff, lower = -1, upper = upper_region[i], r = y_train - f_t_train, h = h_train)
               tmp[i] <- val$minimum
               tmp_val[i] <- val$objective
             } 
             
             idx <- min(which(tmp_val == min(tmp_val)))
             order_val <- order(head(tmp_val, idx)) 

             if( sum(order_val == idx:1) == idx){
               return(tmp[idx])
             }else{
               tmp_order <- order_val  - c(max(order_val), head(order_val, length(order_val)-1))
               if(sum(tmp_order > 0) > 0){
                 tmp_idx <- min(which(tmp_order>0))-1
                 return(tmp[tmp_idx])
               }else{
                 return(tmp[1])
               }
             }
           }else{
             ff <- function(a, r, h, c, s) return(mean(func( (r - a*h)/s,  c)))
             upper_region = c(0.5,10,100,300)
             tmp <- rep(NA, length(upper_region))
             tmp <-  tmp_val <- rep(NA, length(upper_region))
             for(i in 1:length(upper_region)){
               val = optimize(ff, lower = -1, upper = upper_region[i], r = y_train - f_t_train, h = h_train, c = cc, s = ss)
               tmp[i] <- val$minimum
               tmp_val[i] <- val$objective
             } 
             
             idx <- min(which(tmp_val == min(tmp_val)))
             order_val <- order(head(tmp_val, idx)) 
             
             if( sum(order_val == idx:1) == idx){ #continue going down
               return(tmp[idx])
             }else{
               tmp_order <- order_val  - c(max(order_val), head(order_val, length(order_val)-1))
               if(sum(tmp_order > 0) > 0){
                 tmp_idx <- min(which(tmp_order>0))-1
                 return(tmp[tmp_idx])
               }else{
                 return(tmp[1])
               }
             }
           }
         },
         L2Boost = {
           ff = function(a,r,h){
             return(mean(func(r - a*h)))
           }
           return(optimize(ff, lower = -1, upper = 10, r = y_train - f_t_train, h = h_train)$minimum)
         },
         MBoost = {
           ff = function(a,r,h, c){
             return(mean(func(r - a*h, c)))
           }
           return(optimize(ff, lower = -1, upper = 300, r = y_train - f_t_train, h = h_train, c = ss)$minimum)
         },
         Robloss = {
           ff = function(a,r,h, c){
             return(mean(func(r - a*h, c)))
           }
           return(optimize(ff, lower = -1, upper = 300, r = y_train - f_t_train, h = h_train, c = ss)$minimum)
         })
}



#' Boost
#'
#' A function to fit RRBoost which includes options to fit L2Boost, LADBoost, MBoost, Robloss, and SBoost
#'
#' A function to fit RRBoost which includes options to fit L2Boost, LADBoost, MBoost, Robloss, and SBoost
#'
#'@param x_train predictor matrix for training data (matrix/dataframe)
#'@param y_train response vector for training data (vector/dataframe)
#'@param x_val predictor matrix for validation data (matrix/dataframe)
#'@param y_train response vector for validation data (vector/dataframe)
#'@param x_test predictor matrix for test data (matrix/dataframe, optional, required when make_prediction in control = TRUE)
#'@param y_test response vector for test data (vector/dataframe,  optional, required when make_prediction in control = TRUE)
#'@param type type of the boosting method: "L2Boost", "MBoost", "Robloss", "SBoost", "RRBoost". (string)
#'@param error types of the error metric on the test set: "rmse","aad"(average absulute deviation), or "trmse" (trimmed rmse) (array)
#'@param y_init the initial estimator, "median" or "LADTree" (string)
#'@param max_depth the maximum depth of the tree learners (numeric)
#'@param niter number of iterations (for RRBoost T_{1,max} + T_{2,max}) (numeric)
#'@param control control parameters specified with Boost.control()
#'@return A list with the following components:
#'
#' \item{type}{type of the boosting estimator (e.g. 'RRBoost',"L2Boost")}
#' \item{control}{the input control parameters}
#' \item{error}{a vector of error values evaluated on the test set at early stopping time. The length of the vector depends on the `error` argument in the input.  (returned if make_prediction = TRUE in control).}
#' \item{tree_init}{the initial tree (rpart object, returned if y_init = "LADTree")}
#' \item{tree_list}{a list of trees fitted at each iteration (returned if save_tree = TRUE in control) }
#' \item{f_train_init}{a vector of initialized estimator of the training data}
#' \item{alpha}{a vector of base learners' coefficients}
#' \item{early_stop_idx}{early stopping iteration}
#' \item{when_init}{the early stopping time of the first stage of RRBoost (returned if type = "RRBoost)}
#' \item{loss_train}{a vector of training loss values}
#' \item{loss_val}{a vector of validation loss values}
#' \item{err_val}{a vector of validation aad error}
#' \item{err_train}{a vector of training aad error}
#' \item{err_test}{a matrix of test errors (returned if make_prediction = TRUE in control)}
#' \item{f_train}{matrix of training function estimates at all iterations (returned if save_f = TRUE in control)}
#' \item{f_val}{matrix of validation function estimates at all iterations (returned if save_f = TRUE in control)}
#' \item{f_test}{matrix of test function estimates at all iterations (returned if save_f = TRUE and make_prediction = TRUE in control)}
#' \item{var_importance}{vector of permutation importance at early stopping time (returned if cal_imp = TRUE in control)}
#' @author Xiaomeng Ju, \email{xmengju@stat.ubc.ca}
#'
#' @export
#'
Boost <- function(x_train, y_train, x_val, y_val, x_test, y_test, type = "L2Boost", error = c("rmse","aad"),   niter = 200, y_init = "median",  max_depth = 1, control = Boost.control(), newton = 0) {
  
  print(type)

  cc <- NA; n_init <- NA;

  save_f <- control$save_f
  save_tree <- control$save_tree
  save_time <-control$save_time
  make_prediction <- control$make_prediction
  bb <- control$bb
  precision <- control$precision
  shrinkage <- control$shrinkage
  var_select <- rep(1,ncol(x_train))  # will be updated when calling predict 
  
  if(type == "RRBoost"){
    n_init <- control$n_init
    cc <- control$cc_s
    cc_m <- control$cc_m
  }

  if(type == "SBoost"){
    cc <- control$cc_s
  }

  if(type == "MBoost"){
    cc <- control$cc_m
  }

  cal_imp <- control$cal_imp

  if(class(x_train) == "numeric") {
    x_train <- data.frame(x = x_train)
  }else{
    x_train <- data.frame(x_train)
  }
  if(class(x_val) == "numeric") {
    x_val <- data.frame(x = x_val)
  }else{
    x_val <- data.frame(x_val)
  }

  # initialize functions
  init_functions <- init.boosting(type)

  func  <- init_functions$func
  func.grad <- init_functions$func.grad
  func.grad.prime <-  init_functions$func.grad.prime

  if(save_time){
    time_vec <- rep(NA, niter+1)
  }
  # save initialized value
  f_train_init <- NULL

  # save the function value at early stopping time (flag outliers for variable importance calculation)
  f_train_early <- NULL
  f_val_early <- NULL

  # to save the predictions
  if(save_f == TRUE){
    f_train <-  matrix(NA, nrow(x_train), niter)
    f_val <-  matrix(NA,  nrow(x_val), niter)
  }

  # to save the trees
  tree_init <- NULL
  tree_list <- list()

  start_time_a <- Sys.time()
  # initialization
  if(y_init == "LADTree"){
      if(is.null(control$max_depth_init)) {
        max_depth_init <- 3
      }else{
        max_depth_init <- control$max_depth_init
      }

      if(is.null(control$min_leaf_size_init)) {
        min_leaf_size_init = 10
      }else{
        min_leaf_size_init <- control$min_leaf_size_init
      }

    #print(paste("max_depth_init = ", max_depth_init,"min_leaf_size_init = ", min_leaf_size_init))
    dat_tmp <- data.frame(x_train, y_train = y_train)
    tree_init <- rpart(y_train~ ., data = dat_tmp,control = rpart.control(maxdepth = max_depth_init, minbucket = min_leaf_size_init, xval = 0, cp = -Inf), method = alist)
    f_train_early <- f_train_init <- f_t_train <- predict(tree_init, newdata = x_train)
    f_val_early <- f_t_val <-  predict(tree_init, newdata = x_val)
  }else{
    # initialization either with median
    f_train_early <- f_train_init <- f_t_train <- rep(median(y_train), length(y_train))
    f_val_early <- f_t_val <- rep(median(y_train), length(y_val))
  }
  
  start_time_b <- Sys.time()
  time_vec[1] <- start_time_b - start_time_a
  
  # to save alpha
  alpha <- rep(NA, niter)

  # to save the error
  err_train <-   rep(NA, niter)
  err_val <- rep(NA, niter)

  # save the loss for SBoost and RRBoost (for early stop)
  loss_train <- rep(NA, niter)
  loss_val <- rep(NA, niter)

  # denote the transition from S-type to M-type
  init_status <- 0;

  # when S-type is finished
  when_init = NA


  for(i in 1:niter) {

    if(i%%200 ==0) {
     print(c("iteration", i))
    }

    if(init_status == 0) {
      ss <- cal.ss(type, f_t_train, y_train,  cc, bb)
    }

    dat_tmp <- cal.neggrad(type, x_train, y_train, f_t_train, init_status, ss, func, func.grad, cc)
    
    if(newton == 1 & type == "RRBoost" & init_status == 1){
      tmp_newton <- func.tukey.grad.prime((f_t_train - y_train)/ss ,cc = cc)
      idx_newton <- which(dat_tmp$neg_grad == 0)
      dat_tmp$neg_grad <- dat_tmp$neg_grad/func.tukey.grad.prime(  (f_t_train - y_train)/ss ,cc = cc)
      dat_tmp$neg_grad[idx_newton] <- 0
      #print(cbind(dat_tmp$neg_grad, (f_t_train - y_train)/ss, func.tukey.grad.prime( (f_t_train - y_train)/ss ,cc = cc), dat_tmp$neg_grad/func.tukey.grad.prime( (f_t_train - y_train)/ss ,cc = cc)))
    }
    
    tree.model <- rpart(neg_grad~ ., data = dat_tmp, control = rpart.control(maxdepth = max_depth, cp = 0))

    h_train <- predict(tree.model, newdata = data.frame(x_train))
    h_val <- predict(tree.model, newdata = data.frame(x_val))

    alpha[i] <- cal.alpha(type,  f_t_train, h_train, y_train, func, ss = ss, init_status, cc = cc)
    f_t_train <- f_t_train + shrinkage*alpha[i]* h_train
    f_t_val <- f_t_val +  shrinkage*alpha[i]*h_val

    tree_list[[i]] <- tree.model
    err_train[i] <- mean(abs(f_t_train - y_train))
    err_val[i] <- mean(abs(f_t_val - y_val))

    # record loss values for early stopping
    if(type == "SBoost" | (type == "RRBoost" & init_status == 0)){
      loss_val[i] <-cal.ss(type, f_t_val, y_val,  cc, bb)
      loss_train[i] <-ss
    }

    if(type == "RRBoost" & init_status == 1){
      loss_train[i] <- mean(func((f_t_train - y_train)/ss, cc = cc_m))
      loss_val[i] <- mean(func((f_t_val - y_val)/ss, cc = cc_m))
    }

    if(type == "LADBoost"){
      loss_train[i] <- mean(abs(f_t_train - y_train))
      loss_val[i] <-  mean(abs(f_t_val - y_val))
    }

    if(type %in% c("MBoost" ,"Robloss","L2Boost")){
      loss_train[i] <-  mean(func(f_t_train - y_train, cc = ss))
      loss_val[i] <- mean(func(f_t_val - y_val, cc = ss))
    }

    if(i == 1){
      # initiailze the early stopping f
      if(type == "RRBoost"){
        when_init <- 1
      }
      f_train_early  <- f_t_train
      f_val_early <- f_t_val
      early_stop_idx <- 1
    }else{
      
      if(type != "RRBoost"){
        if(type %in% c("MBoost", "Robloss")){
          if(round(err_val[i], precision) < min(round(err_val[1:(i-1)], precision))){
            early_stop_idx  <- i
            f_train_early  <- f_t_train
            f_val_early <- f_t_val
          }
        }else{
          if(round(loss_val[i], precision) < min(round(loss_val[1:(i-1)], precision)) ){
            early_stop_idx  <- i
            f_train_early  <- f_t_train
            f_val_early <- f_t_val
          }
        }
      }
      
      if(type == "RRBoost" & i <= n_init){
        if(round(loss_val[i], precision) < min(round(loss_val[1:(i-1)], precision))){
          when_init <- i
          early_stop_idx <- i
          f_train_early  <- f_t_train
          f_val_early <- f_t_val
        }
      }

      if(type == "RRBoost" &  (init_status == 1)){
          if(round(loss_val[i], precision) < min(round(loss_val[(n_init):(i-1)], precision))){
            early_stop_idx <- i
            f_train_early  <- f_t_train
            f_val_early <- f_t_val
          }
      }
      
      if(type == "RRBoost" & i == n_init){
        init_status <- 1
        f_t_train <- f_train_early  # rest the current one
        f_t_val <- f_val_early
        ss <-  mscale(f_t_train - y_train,  tuning.chi= cc, delta = bb) 
        cc <- cc_m
        loss_val[i] <- mean(func((f_t_val - y_val)/ss, cc = cc_m))
       }
    }
    if(save_f == TRUE){
      f_train[,i] <- f_t_train; f_val[,i] <- f_t_val;
    }
    if(save_time){
      time_vec[i+1] <- Sys.time() - start_time_a
    }
  }

  f_t_train <-   f_train_early
  f_t_val <- f_val_early
  tree_list <- tree_list[1:early_stop_idx]
  model <- list(type = type, control = control, error = error, y_init = y_init,  tree_init = tree_init, tree_list = tree_list, f_t_train = f_t_train, f_t_val = f_t_val,  f_train_init = f_train_init, alpha = alpha,  early_stop_idx = early_stop_idx, when_init = when_init, loss_train = loss_train, loss_val = loss_val,  err_val = err_val, err_train = err_train)
  
  if(make_prediction == TRUE){
    
    if(class(x_test) == "numeric") {
      x_test <- data.frame(x = x_test)
    }else{
      x_test <- data.frame(x_test)
    }
    
    res <- cal_predict(model, x_test, y_test)
    model$f_t_test <- res$f_t_test
    model$err_test <- res$err_test
    model$value <- res$value
    model$var_select <- res$var_select
    if(save_tree == TRUE){
      model$f_test <- res$f_test
    }
  }

  val_trmse <- trmse(control$trim_prop,control$trim_c, f_val_early - y_val)
  model$val_trmse <- val_trmse

  if(save_time){
    model$time_vec<- time_vec
  }
  if(cal_imp == TRUE){
    model$var_importance <- cal_imp_func(model, x_val, y_val)
  }

  if(save_tree == FALSE){
    model$tree_init = NULL
    model$tree_list = NULL
  }

  if(save_f == TRUE){
    model$f_train  <- f_train
    model$f_val  <- f_val
  }

  return(model)
}

#' Boost.validation
#'
#' A function to fit RRBoost where the initialization parameters are chosen by validation
#'
#' A function to fit RRBoost where the initialization parameters are chosen by validation
#'
#'@param x_train predictor matrix for training data (matrix/dataframe)
#'@param y_train response vector for training data (vector/dataframe)
#'@param x_val predictor matrix for validation data (matrix/dataframe)
#'@param y_train response vector for validation data (vector/dataframe)
#'@param x_test predictor matrix for test data (matrix/dataframe, optional, required when make_prediction in control = TRUE)
#'@param y_test response vector for test data (vector/dataframe,  optional, required when make_prediction in control = TRUE)
#'@param error types of the error metric on the test set: "rmse","aad"(average absulute deviation), or "trmse" (trimmed rmse) (array)
#'@param y_init the initial estimator, "median" or "LADTree" (string)
#'@param max_depth the maximum depth of the tree learners (numeric)
#'@param niter number of iterations (for RRBoost T_{1,max} + T_{2,max}) (numeric)
#'@param control control parameters specified with Boost.control()
#'@param max_depth_init_set a vector of possible values of the maximum depth of the initial LADTree that the algorithm choses from
#'@param min_leaf_size_init_set a vector of possible values of the minimum observations per node of the initial LADTree that the algorithm choses from
#'
#'@return A list with components
#' \item{the components of model}{an object returned by Boost that is trained with selected initialization parameters}
#' \item{param}{a vector of selected initialization parameters (return (0,0) if selected initialization is the median of the training responses)}
#' @author Xiaomeng Ju, \email{xmengju@stat.ubc.ca}
#'
#' @export
#'
Boost.validation <- function(x_train, y_train, x_val, y_val, x_test, y_test, type = "RRBoost", error = c("rmse","aad"),  niter = 1000, max_depth = 1, y_init = "LADTree", max_depth_init_set = c(1,2,3,4), min_leaf_size_init_set = c(10,20,30), control = Boost.control()){


  control_tmp <- control
  if(control$cal_imp == TRUE){
    control_tmp$cal_imp <- FALSE
    control_tmp$save_tree <- TRUE
  }

  model_best <- Boost(x_train, y_train, x_val, y_val, x_test, y_test, type = type, error = error,  niter = niter, y_init = "median", max_depth = max_depth, control =  control_tmp)
  #print(paste("median_sboost_val", model_best$err_test[model_best$when_init,], "median_rrboost_val", model_best$value))
  flagger_outlier <- which(abs(model_best$f_t_val - y_val)>3*mad(model_best$f_t_val - y_val))
  if(length(flagger_outlier)>=1){
    best_err <- mean(abs(model_best$f_t_val[-flagger_outlier] - y_val[-flagger_outlier]))  #test with tau-scale
  }else{
    best_err <- mean(abs(model_best$f_t_val - y_val))
  }
  
  params = c(0,0)
  errs_val <-  rep(NA, 1+ length(min_leaf_size_init_set)*length(max_depth_init_set))
  errs_test <- matrix(NA, 1+ length(min_leaf_size_init_set)*length(max_depth_init_set), length(error))
  errs_val[1] <- best_err
  
  if(control$make_prediction){
    errs_test[1,] <- as.numeric(model_best$value)
  }
  
  if(y_init == "LADTree") {
    model_pre_tree <- NA
    combs <- expand.grid(min_leafs= sort(min_leaf_size_init_set,TRUE), max_depths= max_depth_init_set)
    j_tmp <- rep(1, nrow(combs))
    
    tree_init <- list()
    for(j in 1:nrow(combs)) {
      min_leaf_size <- combs[j, 1]
      max_depths <- combs[j, 2]
      dat_tmp <- data.frame(x_train, y_train = y_train)
      tree_init[[j]] <- rpart(y_train~ ., data = dat_tmp,control = rpart.control(maxdepth = max_depths, minbucket = min_leaf_size, xval = 0, cp = -Inf), method = alist)
    }
    
    for(j in 1:length(max_depth_init_set)){
      for(k in 1:(length(min_leaf_size_init_set)-1)){
        idx_jk <- which(combs[,1] == sort(min_leaf_size_init_set, TRUE)[k] & combs[,2] == max_depth_init_set[j])
        idx_jk_plus<- which(combs[,1] == sort(min_leaf_size_init_set, TRUE)[k+1] & combs[,2] == max_depth_init_set[j])
        equal_tmp<- all.equal(tree_init[[idx_jk]],tree_init[[idx_jk_plus]]) == TRUE
        if(length(equal_tmp)==2 & sum(equal_tmp) == 0){
          j_tmp[ idx_jk_plus] <- 0
        }
      }
    }
     
    for(k in 1:length(min_leaf_size_init_set)){
      for(j in 1:(length(max_depth_init_set)-1)){
        idx_kj <- which(combs[,1] == sort(min_leaf_size_init_set, TRUE)[k] & combs[,2] == max_depth_init_set[j])
        idx_kj_plus<- which(combs[,1] == sort(min_leaf_size_init_set, TRUE)[k] & combs[,2] == max_depth_init_set[j+1])
        equal_tmp<- all.equal(tree_init[[idx_kj]],tree_init[[idx_kj_plus]]) == TRUE
        if(length(equal_tmp)==2 & sum(equal_tmp) == 0){
          j_tmp[idx_kj_plus] <- 0
        }
      }
    }
    
    print(j_tmp)
    
      for(j in 1:nrow(combs)) {
        
          if(j_tmp[j] == 1){
            
             min_leaf_size <- combs[j, 1]
             max_depths <- combs[j, 2]
             control_tmp$max_depth_init <- max_depths
             control_tmp$min_leaf_size_init  <- min_leaf_size
             
             model_tmp <- Boost(x_train, y_train, x_val, y_val, x_test, y_test, type = type, error= error,
                                niter = niter, y_init =  "LADTree", max_depth = max_depth,
                                control= control_tmp)
            
      
            if(length(flagger_outlier)>=1){
              err_tmp <- mean(abs(model_tmp$f_t_val[-flagger_outlier] - y_val[-flagger_outlier]))
            }else{
              err_tmp <- mean(abs(model_tmp$f_t_val - y_val))
            }
             
            errs_val[j+1] <- err_tmp
      
          if(control$make_prediction){
            errs_test[j+1,] <- as.numeric(model_tmp$value)
          }
      
        print(paste("leaf size:", min_leaf_size, " depths:", max_depths, " err(val):", round(err_tmp,4), " best err(val) :", round(best_err,4) ,sep = ""))
        if(err_tmp < best_err) {
          model_best <- model_tmp
          params <- combs[j, ]
          best_err <- err_tmp
          rm(model_tmp)
        }else{
          rm(model_tmp)
        }
      }
    }
  }

  if(control$cal_imp == TRUE){
    model_best$var_importance = cal_imp_func(model_best, x_val, y_val)
  }

  if(control$save_tree == FALSE){
      model_best$tree_list = NULL
      model_best$tree_init = NULL
  }
  model_best$params = params
  if(control$save_all_err_rr){
    rownames(errs_test) <- names(errs_val) <- c("median", paste(combs[,1], combs[,2], sep = " "))
    colnames(errs_test) <- error
    model_best$save_all_err_rr <- list(errs_val = errs_val, errs_test = errs_test)
  }
  
  #print(paste("best_val", model_best$value, "best_parameters", params))
  return(model_best)
}

cal_error <- function(control, error_type, f_t_test, y_test){
  if(error_type == "robrmse"){
    x_tmp <- rep(1, length(y_test))
    tauest <- FastTau(x=x_tmp, y=f_t_test - y_test, N=50, kk=2, tt=5, rr=2, approximate=0, seed=456)
    
    return(tauest$scale + median(f_t_test - y_test)^2)
  }
  
  if(error_type == "rmse"){
    return(rmse(f_t_test - y_test))
  }
  if(error_type == "tau"){
    x_tmp <- rep(1, length(y_test))
    tauest <- FastTau(x=x_tmp, y=f_t_test - y_test, N=50, kk=2, tt=5, rr=2, approximate=0, seed=456)
    return(as.numeric(tauest$scale))
  }
  if(error_type == "trmse"){
    return(trmse(control$trim_prop, control$trim_c, f_t_test - y_test)$trmse)
  }
  if(error_type == "trmse_mad"){
    return(trmse(NULL, control$trim_c, f_t_test - y_test)$trmse)
  }
  if(error_type == "trmse_prop_1"){
    return(trmse(0.1, NULL, f_t_test - y_test)$trmse)
  }
  if(error_type == "trmse_prop_2"){
    return(trmse(0.2, NULL, f_t_test - y_test)$trmse)
  }
  if(error_type == "aad"){
    return(mean(abs(f_t_test - y_test)))
  }
  if(error_type == "rrmse"){
    return((median(f_t_test - y_test)^2) + mad(f_t_test - y_test)^2)
  }
}


#' cal_predict
#'
#' A function to make predictions and calculate test error given an object returned by Boost and test data
#'
#' A function to make predictions and calculate test error given an object returned by Boost and test data
#'
#'@param model an object returned by Boost
#'@param x_test predictor matrix for test data (matrix/dataframe)
#'@param y_test response vector for test data (vector/dataframe)
#'@return A list with with the following components:
#'
#' \item{f_t_test}{predicted values with model using x_test as the predictors}
#' \item{err_test}{a matrix of test errors (returned if make_prediction = TRUE in control)}
#' \item{f_test}{matrix of test function estimates at all iterations (returned if save_f = TRUE in control)}
#' \item{value}{a vector of test error evaluated at early stopping time}

#' @author Xiaomeng Ju, \email{xmengju@stat.ubc.ca}
#'
#' @export
#' 
cal_predict <- function(model, x_test, y_test){

  if(class(x_test) == "numeric") {
    x_test <- data.frame(x = x_test)
  }else{
    x_test <- data.frame(x_test)
  }

  type <- model$type
  save_f <- model$control$save_f
  shrinkage <- model$control$shrinkage
  early_stop_idx <- model$early_stop_idx
  when_init <- model$when_init
  n_init <- model$control$n_init
  error <- model$error
  err_test <- data.frame(matrix(NA, nrow = early_stop_idx, ncol = length(error)))
  colnames(err_test) <- error
  control <- model$control
  
  var_select <- rep(0, ncol(x_test)) #1 means selected
  names(var_select) <- colnames(x_test)
  res <- list()

  if(save_f == TRUE){
    f_test <- matrix(NA, nrow(x_test), niter)
  }

  if(model$y_init == "median"){
    f_t_test <- f_t_test_init <- model$f_train_init[1]
  }
  if(model$y_init  == "LADTree"){
    f_t_test <- f_t_test_init <- predict(model$tree_init, newdata = x_test)
    frame <-model$tree_init$frame
    leaves <- frame$var == "<leaf>"
    used <- as.character(unique(frame$var[!leaves]))
    var_select[used] <- 1
  }

  if(type == "RRBoost" & (when_init < early_stop_idx)){
    for(i in c(1:when_init, ((n_init+1):early_stop_idx))){
      f_t_test  <-  f_t_test +    shrinkage*model$alpha[i] *predict(model$tree_list[[i]], newdata = x_test)
      frame <-model$tree_list[[i]]$frame
      leaves <- frame$var == "<leaf>"
      used <- as.character(unique(frame$var[!leaves]))
      var_select[used] <- 1
      if(save_f == TRUE){
        f_test[,i] <- f_t_test
      }
      for(error_type in error){
        err_test[i,error_type] <- cal_error(control, error_type, f_t_test, y_test)
      }
    }
  }else{
    for(i in 1:early_stop_idx){  #stopped at stage 1
      f_t_test  <-  f_t_test +   shrinkage*model$alpha[i] *predict(model$tree_list[[i]], newdata = x_test)
      frame <-model$tree_list[[i]]$frame
      leaves <- frame$var == "<leaf>"
      used <- as.character(unique(frame$var[!leaves]))
      var_select[used] <- 1
      if(save_f == TRUE){
        f_test[,i] <- f_t_test
      }
      for(error_type in error){
        err_test[i,error_type] <- cal_error(control, error_type, f_t_test, y_test)
      }
    }
  }

  
  if(save_f == TRUE){
    res$f_test <- f_test
  }
  res$f_t_test <- f_t_test
  res$err_test <- err_test
  res$value <- err_test[early_stop_idx,]
  res$var_select <- var_select
  return(res)
}

#' cal_imp_func
#'
#' A function to calculate variable importance given an object returned by Boost and validation data
#'
#' A function to calculate variable importance given an object returned by Boost and validation data
#'
#'@param model an object returned by Boost
#'@param x_train predictor matrix for test data (matrix/dataframe)
#'@param y_train response vector for test data (vector/dataframe)
#'@return
#' \item{var_importance}{a vector of permutation variable importance}
#' @author Xiaomeng Ju, \email{xmengju@stat.ubc.ca}
#'
#' @export
#' 
cal_imp_func <- function(model,  x_val, y_val){
  
  if(class(x_val) == "numeric") {
    x_val <- data.frame(x = x_val)
  }else{
    x_val <- data.frame(x_val)
  }
  
  var_imp <-  data.frame(t(rep(0, ncol(x_val))))
  names(var_imp)  <- colnames(x_val)
  when_init <- model$when_init
  early_stop_idx <- model$early_stop_idx
  shrinkage <- model$control$shrinkage
  val_trmse <- model$val_trmse 
  idx <- val_trmse$idx
  alpha <- model$alpha
  type <- model$type
  y_init <- model$y_init
  n_init <- model$control$n_init
  var_select <- model$var_select


  cal.imp.shuffle.j <- function(j){
    
    print(paste("calculating importance for", j, "th variable"))
    x_val_j <- x_val
    x_val_j[,j] <- sample(x_val_j[,j],length(x_val_j[,j]))
    
    if(y_init == "median"){
      f_t_val_j =  model$f_train_init[1]
    }
    if(y_init  == "LADTree"){
      f_t_val_j <- predict(model$tree_init, newdata = data.frame(x_val_j))
    }
    
    if(type == "RRBoost" & (n_init < early_stop_idx)){
      for(i in 1:when_init){
        f_t_val_j  <-  f_t_val_j  + shrinkage*alpha[i] *predict(model$tree_list[[i]], newdata = data.frame(x_val_j))
      }
      for(i in (n_init+1):early_stop_idx){
        f_t_val_j  <-  f_t_val_j  + shrinkage*alpha[i] *predict(model$tree_list[[i]], newdata = data.frame(x_val_j))
      }
    }else{
      for(i in 1:early_stop_idx){
        f_t_val_j  <-  f_t_val_j  + shrinkage*alpha[i] *predict(model$tree_list[[i]], newdata = data.frame(x_val_j))
      }
    }
    return(rmse(f_t_val_j[idx] - y_val[idx]) - val_trmse$trmse)
  }
  
  var_imp <- rep(0, ncol(x_val))
  var_idx <- which(var_select > 0)
  var_imp[var_idx] <- sapply(var_idx, cal.imp.shuffle.j)
  names(var_imp)<- colnames(x_val)
  return(var_imp)
}


