RRBoost: a robust boosting algorithm for regression problems
================
Xiaomeng Ju
2019-10-08

This repository contains `R` (in folder src) code implementing robust boosting algorithms to solve regression problems. Below is an example illustrating its use in `R`. We first install the package

``` r
library(devtools)
```

    ## Loading required package: usethis

``` r
library(rpart)
#devtools::install_github("xmengju/RRBoost",auth_token = "09d286ef07a9e052c6aa0f02729a137250b91a76",force = TRUE)
library(rrboost)
```

We apply it to the airfoil data set with added outliers. The data set is split into 60% training, 20% validation, and 20% test. In addition, 20% oasymmetric outliers are added to the training and validation set.

``` r
set.seed(0)
load("data/airfoil.RData")
idx_test <- sample.int(n = nrow(x), size = floor(.2*nrow(x)), replace = F)
idx_train <- sample(setdiff(1:nrow(x),idx_test), size = floor(.6*nrow(x)))
idx_val <- setdiff(1:nrow(x),c(idx_test, idx_train)) 
outliers <- sample(c(idx_train, idx_val), round(0.2*length(c(idx_train, idx_val))))
y[outliers] <- y[outliers] + rnorm(length(outliers), -20, 0.1)
```

We first calculate the RRBoost estimator initialized with the median of training ys.

``` r
model_RRBoost_median = Boost(x_train = x[idx_train,], y_train = y[idx_train], x_val = x[idx_val,], y_val = y[idx_val], x_test = x[idx_test,], y_test = y[idx_test], type = "RRBoost",error = c("rmse"), y_init = "median", shrinkage = 1, max_depth = 1, niter = 1000, control = Boost.control())
```

    ## [1] "RRBoost"
    ## [1] "iteration" "100"      
    ## [1] "iteration" "200"      
    ## [1] "iteration" "300"      
    ## [1] "iteration" "400"      
    ## [1] "iteration" "500"      
    ## [1] "iteration" "600"      
    ## [1] "iteration" "700"      
    ## [1] "iteration" "800"      
    ## [1] "iteration" "900"      
    ## [1] "iteration" "1000"

``` r
print(model_RRBoost_median$value)
```

    ## [1] 9.051671

Then we calculate the RRBoost estimator initialized with an LADTree, where the parameters of LADTree are selected by a validation procedure. The possibly values of the combination of parameters are specified with parameters 'max\_depth\_init\_set' and 'min\_leaf\_size\_init\_set'.

``` r
  model_RRBoost_cv_LADTree = Boost.validation(x_train = x[idx_train,], y_train = y[idx_train], x_val = x[idx_val,], y_val = y[idx_val], x_test = x[idx_test,], y_test = y[idx_test], type = "RRBoost", error = "rmse", y_init = "LADTree", shrinkage = 1, max_depth = 1, niter = 1000, max_depth_init_set = c(1,2,3,4,5), min_leaf_size_init_set = c(10,20,30), control = Boost.control())
```

    ## [1] "RRBoost"
    ## [1] "iteration" "100"      
    ## [1] "iteration" "200"      
    ## [1] "iteration" "300"      
    ## [1] "iteration" "400"      
    ## [1] "iteration" "500"      
    ## [1] "iteration" "600"      
    ## [1] "iteration" "700"      
    ## [1] "iteration" "800"      
    ## [1] "iteration" "900"      
    ## [1] "iteration" "1000"     
    ## [1] "RRBoost"
    ## [1] "max_depth_init"     "1"                  "min_leaf_size_init"
    ## [4] "10"                
    ## [1] "iteration" "100"      
    ## [1] "iteration" "200"      
    ## [1] "iteration" "300"      
    ## [1] "iteration" "400"      
    ## [1] "iteration" "500"      
    ## [1] "iteration" "600"      
    ## [1] "iteration" "700"      
    ## [1] "iteration" "800"      
    ## [1] "iteration" "900"      
    ## [1] "iteration" "1000"     
    ## [1]  6.314544  8.205009 10.000000  1.000000
    ## [1] "RRBoost"
    ## [1] "max_depth_init"     "1"                  "min_leaf_size_init"
    ## [4] "20"                
    ## [1] "iteration" "100"      
    ## [1] "iteration" "200"      
    ## [1] "iteration" "300"      
    ## [1] "iteration" "400"      
    ## [1] "iteration" "500"      
    ## [1] "iteration" "600"      
    ## [1] "iteration" "700"      
    ## [1] "iteration" "800"      
    ## [1] "iteration" "900"      
    ## [1] "iteration" "1000"     
    ## [1]  6.314544  6.314544 20.000000  1.000000
    ## [1] "RRBoost"
    ## [1] "max_depth_init"     "1"                  "min_leaf_size_init"
    ## [4] "30"                
    ## [1] "iteration" "100"      
    ## [1] "iteration" "200"      
    ## [1] "iteration" "300"      
    ## [1] "iteration" "400"      
    ## [1] "iteration" "500"      
    ## [1] "iteration" "600"      
    ## [1] "iteration" "700"      
    ## [1] "iteration" "800"      
    ## [1] "iteration" "900"      
    ## [1] "iteration" "1000"     
    ## [1]  6.314544  6.314544 30.000000  1.000000
    ## [1] "RRBoost"
    ## [1] "max_depth_init"     "2"                  "min_leaf_size_init"
    ## [4] "10"                
    ## [1] "iteration" "100"      
    ## [1] "iteration" "200"      
    ## [1] "iteration" "300"      
    ## [1] "iteration" "400"      
    ## [1] "iteration" "500"      
    ## [1] "iteration" "600"      
    ## [1] "iteration" "700"      
    ## [1] "iteration" "800"      
    ## [1] "iteration" "900"      
    ## [1] "iteration" "1000"     
    ## [1]  6.256681  6.314544 10.000000  2.000000
    ## [1] "RRBoost"
    ## [1] "max_depth_init"     "2"                  "min_leaf_size_init"
    ## [4] "20"                
    ## [1] "iteration" "100"      
    ## [1] "iteration" "200"      
    ## [1] "iteration" "300"      
    ## [1] "iteration" "400"      
    ## [1] "iteration" "500"      
    ## [1] "iteration" "600"      
    ## [1] "iteration" "700"      
    ## [1] "iteration" "800"      
    ## [1] "iteration" "900"      
    ## [1] "iteration" "1000"     
    ## [1]  6.256681  6.256681 20.000000  2.000000
    ## [1] "RRBoost"
    ## [1] "max_depth_init"     "2"                  "min_leaf_size_init"
    ## [4] "30"                
    ## [1] "iteration" "100"      
    ## [1] "iteration" "200"      
    ## [1] "iteration" "300"      
    ## [1] "iteration" "400"      
    ## [1] "iteration" "500"      
    ## [1] "iteration" "600"      
    ## [1] "iteration" "700"      
    ## [1] "iteration" "800"      
    ## [1] "iteration" "900"      
    ## [1] "iteration" "1000"     
    ## [1]  6.256681  6.256681 30.000000  2.000000
    ## [1] "RRBoost"
    ## [1] "max_depth_init"     "3"                  "min_leaf_size_init"
    ## [4] "10"                
    ## [1] "iteration" "100"      
    ## [1] "iteration" "200"      
    ## [1] "iteration" "300"      
    ## [1] "iteration" "400"      
    ## [1] "iteration" "500"      
    ## [1] "iteration" "600"      
    ## [1] "iteration" "700"      
    ## [1] "iteration" "800"      
    ## [1] "iteration" "900"      
    ## [1] "iteration" "1000"     
    ## [1]  5.994374  6.256681 10.000000  3.000000
    ## [1] "RRBoost"
    ## [1] "max_depth_init"     "3"                  "min_leaf_size_init"
    ## [4] "20"                
    ## [1] "iteration" "100"      
    ## [1] "iteration" "200"      
    ## [1] "iteration" "300"      
    ## [1] "iteration" "400"      
    ## [1] "iteration" "500"      
    ## [1] "iteration" "600"      
    ## [1] "iteration" "700"      
    ## [1] "iteration" "800"      
    ## [1] "iteration" "900"      
    ## [1] "iteration" "1000"     
    ## [1]  6.001818  5.994374 20.000000  3.000000
    ## [1] "RRBoost"
    ## [1] "max_depth_init"     "3"                  "min_leaf_size_init"
    ## [4] "30"                
    ## [1] "iteration" "100"      
    ## [1] "iteration" "200"      
    ## [1] "iteration" "300"      
    ## [1] "iteration" "400"      
    ## [1] "iteration" "500"      
    ## [1] "iteration" "600"      
    ## [1] "iteration" "700"      
    ## [1] "iteration" "800"      
    ## [1] "iteration" "900"      
    ## [1] "iteration" "1000"     
    ## [1]  6.068829  5.994374 30.000000  3.000000
    ## [1] "RRBoost"
    ## [1] "max_depth_init"     "4"                  "min_leaf_size_init"
    ## [4] "10"                
    ## [1] "iteration" "100"      
    ## [1] "iteration" "200"      
    ## [1] "iteration" "300"      
    ## [1] "iteration" "400"      
    ## [1] "iteration" "500"      
    ## [1] "iteration" "600"      
    ## [1] "iteration" "700"      
    ## [1] "iteration" "800"      
    ## [1] "iteration" "900"      
    ## [1] "iteration" "1000"     
    ## [1]  5.844458  5.994374 10.000000  4.000000
    ## [1] "RRBoost"
    ## [1] "max_depth_init"     "4"                  "min_leaf_size_init"
    ## [4] "20"                
    ## [1] "iteration" "100"      
    ## [1] "iteration" "200"      
    ## [1] "iteration" "300"      
    ## [1] "iteration" "400"      
    ## [1] "iteration" "500"      
    ## [1] "iteration" "600"      
    ## [1] "iteration" "700"      
    ## [1] "iteration" "800"      
    ## [1] "iteration" "900"      
    ## [1] "iteration" "1000"     
    ## [1]  5.906432  5.844458 20.000000  4.000000
    ## [1] "RRBoost"
    ## [1] "max_depth_init"     "4"                  "min_leaf_size_init"
    ## [4] "30"                
    ## [1] "iteration" "100"      
    ## [1] "iteration" "200"      
    ## [1] "iteration" "300"      
    ## [1] "iteration" "400"      
    ## [1] "iteration" "500"      
    ## [1] "iteration" "600"      
    ## [1] "iteration" "700"      
    ## [1] "iteration" "800"      
    ## [1] "iteration" "900"      
    ## [1] "iteration" "1000"     
    ## [1]  5.900764  5.844458 30.000000  4.000000
    ## [1] "RRBoost"
    ## [1] "max_depth_init"     "5"                  "min_leaf_size_init"
    ## [4] "10"                
    ## [1] "iteration" "100"      
    ## [1] "iteration" "200"      
    ## [1] "iteration" "300"      
    ## [1] "iteration" "400"      
    ## [1] "iteration" "500"      
    ## [1] "iteration" "600"      
    ## [1] "iteration" "700"      
    ## [1] "iteration" "800"      
    ## [1] "iteration" "900"      
    ## [1] "iteration" "1000"     
    ## [1]  5.632800  5.844458 10.000000  5.000000
    ## [1] "RRBoost"
    ## [1] "max_depth_init"     "5"                  "min_leaf_size_init"
    ## [4] "20"                
    ## [1] "iteration" "100"      
    ## [1] "iteration" "200"      
    ## [1] "iteration" "300"      
    ## [1] "iteration" "400"      
    ## [1] "iteration" "500"      
    ## [1] "iteration" "600"      
    ## [1] "iteration" "700"      
    ## [1] "iteration" "800"      
    ## [1] "iteration" "900"      
    ## [1] "iteration" "1000"     
    ## [1]  5.815917  5.632800 20.000000  5.000000
    ## [1] "RRBoost"
    ## [1] "max_depth_init"     "5"                  "min_leaf_size_init"
    ## [4] "30"                
    ## [1] "iteration" "100"      
    ## [1] "iteration" "200"      
    ## [1] "iteration" "300"      
    ## [1] "iteration" "400"      
    ## [1] "iteration" "500"      
    ## [1] "iteration" "600"      
    ## [1] "iteration" "700"      
    ## [1] "iteration" "800"      
    ## [1] "iteration" "900"      
    ## [1] "iteration" "1000"     
    ## [1]  5.971944  5.632800 30.000000  5.000000

We compare error at early stopping time: ("rmse" as we previously specified)

``` r
 print(model_RRBoost_median$value)
```

    ## [1] 9.051671

``` r
 print(model_RRBoost_cv_LADTree$model$value)
```

    ## [1] 3.917782

and permutation variable importance:

``` r
 print(model_RRBoost_median$var_importance)
```

    ##   frequency angle chord_length   velocity thickness
    ## 1  1.537253     0            0 0.03372489 0.4712679

``` r
 print(model_RRBoost_cv_LADTree$model$var_importance)
```

    ##   frequency     angle chord_length  velocity thickness
    ## 1  4.598502 0.6875571    0.9338656 0.2610984  3.823986

The values of the parameters of LADTree selected via validation are

``` r
 print(model_RRBoost_cv_LADTree$params)
```

    ##    min_leafs max_depths
    ## 13        10          5

To train RRBoost with user-specified parameters of the initial tree, run

``` r
 model_RRBoost_LADTree = Boost(x_train = x[idx_train,], y_train = y[idx_train], x_val = x[idx_val,], y_val = y[idx_val], x_test = x[idx_test,], y_test = y[idx_test], type = "RRBoost",error = "rmse", y_init = "LADTree", shrinkage = 1, max_depth = 1, niter = 1000, control = Boost.control(max_depth_init = 2,min_leaf_size_init = 20 ))
```

    ## [1] "RRBoost"
    ## [1] "max_depth_init"     "2"                  "min_leaf_size_init"
    ## [4] "20"                
    ## [1] "iteration" "100"      
    ## [1] "iteration" "200"      
    ## [1] "iteration" "300"      
    ## [1] "iteration" "400"      
    ## [1] "iteration" "500"      
    ## [1] "iteration" "600"      
    ## [1] "iteration" "700"      
    ## [1] "iteration" "800"      
    ## [1] "iteration" "900"      
    ## [1] "iteration" "1000"

We can also separate training, predicting, and calculating variable importance

``` r
 model = Boost(x_train = x[idx_train,], y_train = y[idx_train], x_val = x[idx_val,], y_val = y[idx_val], x_test = NULL, y_test = NULL, type = "RRBoost",error = "rmse", y_init = "LADTree", shrinkage = 1, max_depth = 1, niter = 1000, control = Boost.control(cal_imp = FALSE, save_tree = TRUE, make_prediction = FALSE, max_depth_init = 2,min_leaf_size_init = 20))
```

    ## [1] "RRBoost"
    ## [1] "max_depth_init"     "2"                  "min_leaf_size_init"
    ## [4] "20"                
    ## [1] "iteration" "100"      
    ## [1] "iteration" "200"      
    ## [1] "iteration" "300"      
    ## [1] "iteration" "400"      
    ## [1] "iteration" "500"      
    ## [1] "iteration" "600"      
    ## [1] "iteration" "700"      
    ## [1] "iteration" "800"      
    ## [1] "iteration" "900"      
    ## [1] "iteration" "1000"

``` r
 model <- cal_predict(model, x_test = x[idx_test,], y_test = y[idx_test])
 model <-  cal_imp(model, x_train = x[idx_train,], y_train = y[idx_train])
 print(model$value)
```

    ## [1] 4.937753

``` r
 print(model$var_importance)
```

    ##   frequency      angle chord_length  velocity thickness
    ## 1  3.027903 0.08711746    0.7676057 0.3144959  2.981596
