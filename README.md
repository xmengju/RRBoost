RRBoost: a robust boosting algorithm for regression problems
================
Xiaomeng Ju
2019-10-08

This repository contains `R`code implementing robust boosting algorithms to solve regression problems. Before illustrating its usage, we load needed R packages and install the rrboost package.

``` r
library(devtools)
```

    ## Loading required package: usethis

``` r
library(rpart)
#devtools::install_github("xmengju/RRBoost",auth_token ="09d286ef07a9e052c6aa0f02729a137250b91a76",force = TRUE)
library(rrboost)
```

To demonstrate the usage of the package, we will use a standard build-in dataset airfoil, which consists of 6 variables and 1503 observations. Let's print out the first 5 observations here.

``` r
data("airfoil")
head(cbind(x,y))
```

    ##   frequency angle chord_length velocity  thickness       y
    ## 1       800     0       0.3048     71.3 0.00266337 126.201
    ## 2      1000     0       0.3048     71.3 0.00266337 125.201
    ## 3      1250     0       0.3048     71.3 0.00266337 125.951
    ## 4      1600     0       0.3048     71.3 0.00266337 127.591
    ## 5      2000     0       0.3048     71.3 0.00266337 127.461
    ## 6      2500     0       0.3048     71.3 0.00266337 125.571

The goal here is to predict the response `y` using the 5 predictor variables. We split the dataset into 60% training, 20% validation, and 20% test sets. The indices of training, validation, and test sets are sampled as below.

``` r
set.seed(0)
idx_test <- sample.int(n = nrow(x), size = floor(.2*nrow(x)), replace = F)
idx_train <- sample(setdiff(1:nrow(x),idx_test), size = floor(.6*nrow(x)))
idx_val <- setdiff(1:nrow(x),c(idx_test, idx_train)) 
```

The aim of RRBoost is to protect the estimator from being affected by outliers. We added 20% oasymmetric outliers to the training and validation set to resemble a contaminated setting.

``` r
outliers <- sample(c(idx_train, idx_val), round(0.2*length(c(idx_train, idx_val))))
y[outliers] <- y[outliers] + rnorm(length(outliers), -20, 0.1)
```

The `Boost` function provides calculation of the `RRBoost` estimator, as well as previously proposed `L2Boost`,`LADBoost`,`MBoost`,and `Robloss`. For `RRBoost`, it allows user to define the initialization of the estimator: `median` or `LADTree`. We also provide a `Boost.validation` function to select the parameters of `LADTree` based on a the performance on the validation set.

Let us see how to use the functions to calculate the RRBoost with different initializations: initialized with the median of training ys (`model_RRBoost_median`), a LADTree with user specified parameters (`model_RRBoost_LADTree`), and a LADTree whose paramsters are selected by validation (`model_RRBoost_cv_LADTree`). We chose `rmse` as the type error what we measure on the test set. The other options are `aad`, `trmse`, and `rrmse` (see package descriotion), which can be used as a reference when the the user suspects the test set also has outliers.

The boosting estimators are calculated using decision stumps (`max_depth = 1`) as weak learners over 1000 iterations (`niter = 1000`). Setting `make_prediction =  TRUE, cal_imp = TRUE` in `control` allows the function to return predictions and permutation variable importance.

For `model_RRBoost_LADTree`, we specifiy the maximum depth of the LADTree `max_depth_init = 2` and the minimum number of observations per node `min_leaf_size_init = 20`.

For `model_RRBoost_cv_LADTree`, we set the possible combinations of parameters `max_depth_init_set = c(1,2,3,4,5)` and`min_leaf_size_init_set = c(10,20,30)`.

``` r
model_RRBoost_median = Boost(x_train = x[idx_train,], y_train = y[idx_train], x_val = x[idx_val,], y_val = y[idx_val], x_test = x[idx_test,], y_test = y[idx_test], type = "RRBoost",error = "rmse", y_init = "median", max_depth = 1, niter = 1000, control = Boost.control(make_prediction =  TRUE, cal_imp = TRUE))
```

    ## [1] "RRBoost"
    ## [1] "iteration" "200"      
    ## [1] "iteration" "400"      
    ## [1] "iteration" "600"      
    ## [1] "iteration" "800"      
    ## [1] "iteration" "1000"

``` r
 model_RRBoost_LADTree = Boost(x_train = x[idx_train,], y_train = y[idx_train], x_val = x[idx_val,], y_val = y[idx_val], x_test = x[idx_test,], y_test = y[idx_test], type = "RRBoost",error = "rmse", y_init = "LADTree", max_depth = 1, niter = 1000, control = Boost.control(max_depth_init = 2,min_leaf_size_init = 20, make_prediction =  TRUE, cal_imp = TRUE))
```

    ## [1] "RRBoost"
    ## [1] "max_depth_init =  2 min_leaf_size_init =  20"
    ## [1] "iteration" "200"      
    ## [1] "iteration" "400"      
    ## [1] "iteration" "600"      
    ## [1] "iteration" "800"      
    ## [1] "iteration" "1000"

``` r
model_RRBoost_cv_LADTree = Boost.validation(x_train = x[idx_train,], y_train = y[idx_train], x_val = x[idx_val,], y_val = y[idx_val], x_test = x[idx_test,], y_test = y[idx_test], type = "RRBoost", error = "rmse", y_init = "LADTree", max_depth = 1, niter = 1000, max_depth_init_set = c(1,2,3,4,5), min_leaf_size_init_set = c(10,20,30), control = Boost.control(make_prediction =  TRUE, cal_imp = TRUE))
```

    ## [1] "RRBoost"
    ## [1] "iteration" "200"      
    ## [1] "iteration" "400"      
    ## [1] "iteration" "600"      
    ## [1] "iteration" "800"      
    ## [1] "iteration" "1000"     
    ## [1] "RRBoost"
    ## [1] "max_depth_init =  1 min_leaf_size_init =  10"
    ## [1] "iteration" "200"      
    ## [1] "iteration" "400"      
    ## [1] "iteration" "600"      
    ## [1] "iteration" "800"      
    ## [1] "iteration" "1000"     
    ## [1] "RRBoost"
    ## [1] "max_depth_init =  1 min_leaf_size_init =  20"
    ## [1] "iteration" "200"      
    ## [1] "iteration" "400"      
    ## [1] "iteration" "600"      
    ## [1] "iteration" "800"      
    ## [1] "iteration" "1000"     
    ## [1] "RRBoost"
    ## [1] "max_depth_init =  1 min_leaf_size_init =  30"
    ## [1] "iteration" "200"      
    ## [1] "iteration" "400"      
    ## [1] "iteration" "600"      
    ## [1] "iteration" "800"      
    ## [1] "iteration" "1000"     
    ## [1] "RRBoost"
    ## [1] "max_depth_init =  2 min_leaf_size_init =  10"
    ## [1] "iteration" "200"      
    ## [1] "iteration" "400"      
    ## [1] "iteration" "600"      
    ## [1] "iteration" "800"      
    ## [1] "iteration" "1000"     
    ## [1] "RRBoost"
    ## [1] "max_depth_init =  2 min_leaf_size_init =  20"
    ## [1] "iteration" "200"      
    ## [1] "iteration" "400"      
    ## [1] "iteration" "600"      
    ## [1] "iteration" "800"      
    ## [1] "iteration" "1000"     
    ## [1] "RRBoost"
    ## [1] "max_depth_init =  2 min_leaf_size_init =  30"
    ## [1] "iteration" "200"      
    ## [1] "iteration" "400"      
    ## [1] "iteration" "600"      
    ## [1] "iteration" "800"      
    ## [1] "iteration" "1000"     
    ## [1] "RRBoost"
    ## [1] "max_depth_init =  3 min_leaf_size_init =  10"
    ## [1] "iteration" "200"      
    ## [1] "iteration" "400"      
    ## [1] "iteration" "600"      
    ## [1] "iteration" "800"      
    ## [1] "iteration" "1000"     
    ## [1] "RRBoost"
    ## [1] "max_depth_init =  3 min_leaf_size_init =  20"
    ## [1] "iteration" "200"      
    ## [1] "iteration" "400"      
    ## [1] "iteration" "600"      
    ## [1] "iteration" "800"      
    ## [1] "iteration" "1000"     
    ## [1] "RRBoost"
    ## [1] "max_depth_init =  3 min_leaf_size_init =  30"
    ## [1] "iteration" "200"      
    ## [1] "iteration" "400"      
    ## [1] "iteration" "600"      
    ## [1] "iteration" "800"      
    ## [1] "iteration" "1000"     
    ## [1] "RRBoost"
    ## [1] "max_depth_init =  4 min_leaf_size_init =  10"
    ## [1] "iteration" "200"      
    ## [1] "iteration" "400"      
    ## [1] "iteration" "600"      
    ## [1] "iteration" "800"      
    ## [1] "iteration" "1000"     
    ## [1] "RRBoost"
    ## [1] "max_depth_init =  4 min_leaf_size_init =  20"
    ## [1] "iteration" "200"      
    ## [1] "iteration" "400"      
    ## [1] "iteration" "600"      
    ## [1] "iteration" "800"      
    ## [1] "iteration" "1000"     
    ## [1] "RRBoost"
    ## [1] "max_depth_init =  4 min_leaf_size_init =  30"
    ## [1] "iteration" "200"      
    ## [1] "iteration" "400"      
    ## [1] "iteration" "600"      
    ## [1] "iteration" "800"      
    ## [1] "iteration" "1000"     
    ## [1] "RRBoost"
    ## [1] "max_depth_init =  5 min_leaf_size_init =  10"
    ## [1] "iteration" "200"      
    ## [1] "iteration" "400"      
    ## [1] "iteration" "600"      
    ## [1] "iteration" "800"      
    ## [1] "iteration" "1000"     
    ## [1] "RRBoost"
    ## [1] "max_depth_init =  5 min_leaf_size_init =  20"
    ## [1] "iteration" "200"      
    ## [1] "iteration" "400"      
    ## [1] "iteration" "600"      
    ## [1] "iteration" "800"      
    ## [1] "iteration" "1000"     
    ## [1] "RRBoost"
    ## [1] "max_depth_init =  5 min_leaf_size_init =  30"
    ## [1] "iteration" "200"      
    ## [1] "iteration" "400"      
    ## [1] "iteration" "600"      
    ## [1] "iteration" "800"      
    ## [1] "iteration" "1000"

In order to compare the predictive performance of the models, we access the the `value` field of the returned objects, which is the test error evaluated at the early stopping time. The early stopping time is determined by the validation set to avoid overfitting (see paper for more details). In this example, we specified the type of error to be `rmse`, the returned `value` will be the test `rmse` evaluated at the early stopping time. If we let the function return multiple types of errors, `value` will be a vector of those errors.

We compare `value` of different initializations

``` r
 print(model_RRBoost_median$value)
```

    ## [1] 9.051671

``` r
 print(model_RRBoost_LADTree$value)
```

    ## [1] 4.937753

``` r
 print(model_RRBoost_cv_LADTree$value)
```

    ## [1] 3.917782

The validation selected LADTree parameters are

``` r
 print(model_RRBoost_cv_LADTree$params)
```

    ##    min_leafs max_depths
    ## 13        10          5

When setting `cal_imp = TRUE` in `control`, we have access to the variable importance that is calculated with a permutation procedure described in the paper. Our best prediction model `model_RRBoost_cv_LADTree` shows that `frequency` and `thickness` are the top 2 most important variables.

``` r
 print(model_RRBoost_median$var_importance)
```

    ##   frequency angle chord_length   velocity thickness
    ## 1  1.537253     0            0 0.03372489 0.4712679

``` r
 print(model_RRBoost_LADTree$var_importance)
```

    ##   frequency      angle chord_length  velocity thickness
    ## 1  3.536077 0.06884003    0.7122387 0.2858896  3.106163

``` r
 print(model_RRBoost_cv_LADTree$var_importance)
```

    ##   frequency     angle chord_length  velocity thickness
    ## 1   4.27229 0.9042215     1.000275 0.2470756  3.689421

In the package, we also provide a way that separates training, making predictions, and calculating variable importance. It allows the flexiblity to make predictions on multiple test sets and calculate variable importance when needed. Note that to allow this separation, `save_tree = TRUE` is required in `control`. This saves the trainer weak learners, which are needed for calculating prediction and variable importance, as a model component.

``` r
 model = Boost(x_train = x[idx_train,], y_train = y[idx_train], x_val = x[idx_val,], y_val = y[idx_val], x_test = x[idx_test,], y_test = y[idx_test], type = "RRBoost",error = "rmse", y_init = "LADTree", max_depth = 1, niter = 1000, control = Boost.control(max_depth_init = 2,min_leaf_size_init = 20,save_tree = TRUE, make_prediction =  FALSE, cal_imp = FALSE))
```

    ## [1] "RRBoost"
    ## [1] "max_depth_init =  2 min_leaf_size_init =  20"
    ## [1] "iteration" "200"      
    ## [1] "iteration" "400"      
    ## [1] "iteration" "600"      
    ## [1] "iteration" "800"      
    ## [1] "iteration" "1000"

``` r
 prediction <- cal_predict(model, x_test = x[idx_test,], y_test = y[idx_test])
 var_importance <-  cal_imp(model, x_train = x[idx_train,], y_train = y[idx_train])
 print(prediction$value)
```

    ## [1] 4.937753

``` r
 print(var_importance)
```

    ##   frequency      angle chord_length  velocity thickness
    ## 1  3.027903 0.08711746    0.7676057 0.3144959  2.981596
