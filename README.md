RRBoost: a robust boosting algorithm for regression problems
================
Xiaomeng Ju and Matias Salibian Barrera
2020-07-31

This repository contains `R` code implementing a robust boosting
algorithm to fit non-parametric regression models. We refer the
interested reader to [paper URL](https://github.com).

This package can be installed as follows:

``` r
devtools::install_github("xmengju/RRBoost")
```

To demonstrate the usage of the package, we will use the `airfoil`
dataset. It consists of `n = 1503` observations with `p = 5` explanatory
variables, and the goal is to predict the response variable `y`. We load
the data and look at the first few observations:

``` r
data(airfoil)
head(airfoil)
```

    ##   frequency angle chord_length velocity  thickness       y
    ## 1       800     0       0.3048     71.3 0.00266337 126.201
    ## 2      1000     0       0.3048     71.3 0.00266337 125.201
    ## 3      1250     0       0.3048     71.3 0.00266337 125.951
    ## 4      1600     0       0.3048     71.3 0.00266337 127.591
    ## 5      2000     0       0.3048     71.3 0.00266337 127.461
    ## 6      2500     0       0.3048     71.3 0.00266337 125.571

In order to train our predictor we split the data set into a training
set (with 60% of the available data), a validation set and a test set
(both with 20% of the data). The following code creates these sets:

``` r
library(rrboost)
n <- nrow(airfoil) 
n0 <- floor( 0.2 * n ) 
set.seed(0)
idx_test <- sample(n, n0)
idx_train <- sample((1:n)[-idx_test], floor( 0.6 * n ) )
idx_val <- (1:n)[ -c(idx_test, idx_train) ] 
```

To illustrate the robustness of RRBoost, we added 20% of asymmetric
outliers to the training and validation set to resemble a contaminated
setting.

``` r
aircont <- airfoil
# indices of observations that may be contaminated
tmp <- c(idx_train, idx_val) 
# randomly sample 20% of them
outliers <- sample(tmp, round( 0.2 * length(tmp) ) )
# shift the response variable with a N(-20, 0.1^2) r.v.
aircont$y[outliers] <- aircont$y[outliers] + rnorm(length(outliers), -20, 0.1)
```

We now create the training, validation and test data sets:

``` r
xcont <- aircont[, -6]
ycont <- aircont$y
xtrain <- xcont[ idx_train, ]
ytrain <- ycont[ idx_train ]
xval <- xcont[ idx_val, ]
yval <- ycont[ idx_val ]
xtest <- xcont[ idx_test, ]
ytest <- ycont[ idx_test ]
```

The `Boost` function provides calculation of the `RRBoost` estimator, as
well as previously proposed `L2Boost`,`LADBoost`,`MBoost`,and
`Robloss`.  
For `RRBoost`, it allows user to define the initialization of the
estimator: `median` or `LADTree`. We also provide a `Boost.validation`
function to select the parameters of `LADTree` based on a the
performance on the validation set.

Let us see how to use the functions to calculate the RRBoost with
different initializations: initialized with the median of training ys
(`model_RRBoost_median`), a LADTree with user specified parameters
(`model_RRBoost_LADTree`), and a LADTree whose paramsters are selected
by validation (`model_RRBoost_cv_LADTree`). We chose `rmse` as the type
error what we measure on the test set. The other options are `aad`,
`trmse`, and `rrmse` (see package descriotion), which can be used as a
reference when the the user suspects the test set also has outliers.

The boosting estimators are calculated using decision stumps (`max_depth
= 1`) as weak learners over 1000 iterations (`niter = 1000`). Setting
`make_prediction = TRUE, cal_imp = TRUE` in `control` allows the
function to return predictions and permutation variable importance.

For `model_RRBoost_LADTree`, we specifiy the maximum depth of the
LADTree `max_depth_init = 2` and the minimum number of observations per
node `min_leaf_size_init = 20`.

For `model_RRBoost_cv_LADTree`, we set the possible combinations of
parameters `max_depth_init_set = c(1,2,3,4,5)`
and`min_leaf_size_init_set = c(10,20,30)`.

``` r
model_RRBoost_median = Boost(x_train = xtrain, y_train = ytrain, 
                             x_val = xval, y_val = yval, 
                             x_test = xtest, y_test = ytest, 
                             type = "RRBoost", error = "rmse", 
                             y_init = "median", max_depth = 1, 
                             niter = 1000, 
                             control = Boost.control(make_prediction =  TRUE,
                                                     cal_imp = TRUE))

model_RRBoost_LADTree = Boost(x_train = xtrain, y_train = ytrain, 
                              x_val = xval, y_val = yval, 
                              x_test = xtest, y_test = ytest, 
                              type = "RRBoost", error = "rmse", 
                              y_init = "LADTree", max_depth = 1, 
                              niter = 1000, 
                              control = Boost.control(max_depth_init = 2,
                                                      min_leaf_size_init = 20,
                                                      make_prediction =  TRUE,
                                                      cal_imp = TRUE))

model_RRBoost_cv_LADTree = Boost.validation(x_train = xtrain, 
                                            y_train = ytrain, 
                                            x_val = xval, y_val = yval, 
                                            x_test = xtest, 
                                            y_test = ytest, 
                                            type = "RRBoost", 
                                            error = "rmse", 
                                            y_init = "LADTree", 
                                            max_depth = 1, niter = 1000,
                                            max_depth_init_set = c(1,2,3,4,5),
                                            min_leaf_size_init_set = c(10,20,30),
                                            control = Boost.control(make_prediction =  TRUE, cal_imp = TRUE))
```

The validation selected LADTree parameters are

``` r
 print(model_RRBoost_cv_LADTree$params)
```

    ##    min_leafs max_depths
    ## 15        10          5

In order to compare the predictive performance of the models, we access
the `value` field of the returned objects, which is the test error
evaluated at the early stopping time. The early stopping time is
determined by the validation set to avoid overfitting (see paper for
more details). In this example, we specified the type of error to be
`rmse`, the returned `value` will be the test `rmse` evaluated at the
early stopping time. If we let the function return multiple types of
errors, `value` will be a vector of those errors.

We compare `value` of different initializations

``` r
print(c(median = model_RRBoost_median$value, 
        LADTree = model_RRBoost_LADTree$value,
        cv_LADTree = model_RRBoost_cv_LADTree$value))
```

    ##     median    LADTree cv_LADTree 
    ##   5.350651   5.472190   3.944847

When setting `cal_imp = TRUE` in `control`, we have access to the
variable importance that is calculated with a permutation procedure
described in the paper. Our best prediction model
`model_RRBoost_cv_LADTree` shows that `frequency` and `thickness` are
the top 2 most important variables.

``` r
 print(cbind(median = model_RRBoost_median$var_importance,
         LADTree = model_RRBoost_LADTree$var_importance,
         cv_LADTree = model_RRBoost_cv_LADTree$var_importance))
```

    ##                    median    LADTree cv_LADTree
    ## frequency     2.812442313 2.61368050  3.9228494
    ## angle        -0.008505854 0.06265083  0.1673166
    ## chord_length  0.166805204 0.22086603  0.6838844
    ## velocity      0.328773995 0.18913548  0.3955189
    ## thickness     2.121111371 2.17302992  3.0588768

In the package, we also provide a way that separates training, making
predictions, and calculating variable importance. It allows the
flexiblity to make predictions on multiple test sets and calculate
variable importance when needed. Note that to allow this separation,
`save_tree = TRUE` is required in `control`. This saves the trainer weak
learners, which are needed for calculating prediction and variable
importance, as a model component.

``` r
model = Boost(x_train = xtrain, y_train = ytrain, 
              x_val = xval, y_val = yval, 
              x_test = xtest, y_test = ytest, 
              type = "RRBoost", error = "rmse", 
              y_init = "LADTree", max_depth = 1, niter = 1000, 
              control = Boost.control(max_depth_init = 2, 
                                      min_leaf_size_init = 20, 
                                      save_tree = TRUE, 
                                      make_prediction =  FALSE, 
                                      cal_imp = FALSE))

prediction <- cal_predict(model, x_test = xtest, y_test = ytest)
# This does not work
# var_importance <-  cal_imp_func(model, x_train = xtrain, y_train = ytrain)
```

Sanity check: the results are the same

``` r
print(c(prediction$value, model_RRBoost_LADTree$value))
```

    ## [1] 5.47219 5.47219

``` r
# var_importance could not be calculated above
# print(c(var_importance, model_RRBoost_LADTree$var_importance))
```
