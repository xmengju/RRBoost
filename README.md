RRBoost: a robust boosting algorithm for regression problems
================
Xiaomeng Ju and Matias Salibian Barrera
2020-08-21

This repository contains `R` code implementing the robust boosting
algorithm described in [Ju X, Salibian-Barrera M. Robust Boosting for
Regression Problems](https://arxiv.org/abs/2002.02054](https://www.sciencedirect.com/science/article/pii/S0167947320301560). This method provides an
outlier-robust fit for non-parametric regression models.

The package can be installed as follows:

``` r
devtools::install_github("xmengju/RRBoost")
```

Below we illustrate the use of the package with the `airfoil` dataset.
It consists of `n = 1503` observations with `p = 5` explanatory
variables, and the goal is to predict the response variable `y`. We load
the data and look at the first few observations:

``` r
library(RRBoost)
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

In order to train our predictor, we split the data set into a `training`
set (with 60% of the available data), a `validation` set and a `test`
set (both with 20% of the data). We first randomly select the
observations for each of these three sets:

``` r
n <- nrow(airfoil) 
n0 <- floor( 0.2 * n ) 
set.seed(123)
idx_test <- sample(n, n0)
idx_train <- sample((1:n)[-idx_test], floor( 0.6 * n ) )
idx_val <- (1:n)[ -c(idx_test, idx_train) ] 
```

<!-- To illustrate the robustness of the `RRBoost` fit, we  -->

<!-- add 20% of asymmetric outliers to the `training` and  -->

<!-- `validation`. We randomly select 20% of the observations -->

<!-- among these two sets combined, and perturb the response -->

<!-- variable by adding a `N( -20, 0.1^2)` random variable.  -->

<!-- ```{r addoutlier} -->

<!-- aircont <- airfoil -->

<!-- # indices of observations that may be contaminated -->

<!-- tmp <- c(idx_train, idx_val)  -->

<!-- # randomly sample 20% of them -->

<!-- nnc <- round( 0.2 * length(tmp) ) -->

<!-- outliers <- sample(tmp, nnc) -->

<!-- # shift the response variable with a N(-20, 0.1^2) r.v. -->

<!-- # aircont$y[outliers] <- aircont$y[outliers] + rnorm(nnc, -20, 0.1) -->

<!-- ``` -->

We now create the matrices of explanatory variables (`x`) and vectors of
responses (`y`) corresponding to this partition.
<!-- Note that `ytrain` and `yval` may contain outliers. -->

``` r
xx <- airfoil[, -6]
yy <- airfoil$y
xtrain <- xx[ idx_train, ]
ytrain <- yy[ idx_train ]
xval <- xx[ idx_val, ]
yval <- yy[ idx_val ]
xtest <- xx[ idx_test, ]
ytest <- yy[ idx_test ]
```

The `Boost` function implements the `RRBoost` estimator, as well as the
following previously proposed boosting algorithms: `L2Boost`,
`LADBoost`, `MBoost` ([Friedman, J. H.
(2001)](https://doi.org/10.1214/aos/1013203451)) and `Robloss` ([Lutz et
al. (2008)](https://doi.org/10.1016/j.csda.2007.11.006)).

There are two choices for the initial fit used in the `RRBoost`
algorithm: `median` (corresponds to a constant initial fit equal to the
median of the response variable in the training set), and `LADTree` (the
initial fit is an L1 Tree as proposed in [Breiman
(1984)](https://doi.org/10.1201/9781315139470)). To construct the L1
initial fit the user can select the desired maximum depth and minimum
leaf size, or use the function `Boost.validation` to select these
parameters based on the predictive performance on the validation set.

We now train the `RRBoost` predictor using the following three initial
fits: `median`, an `LADTree` with parameters chosen a priori, and an
`LADTree` with parameteres selected using the validation set. The
performance of the resulting three predictors on the test set will be
measured using the root mean squared (prediction) error (setting the
argument `error = "rmse"`). Other possible options for this argument
are: `aad` (average absolute deviation), `trmse` (trimmed root mean
squared (prediction) error) and `rrmse` (a robust root mean squared
(prediction) error). For more information see the help page of the
function `Boost`.

The depth of the base learners in the boosting algorithm is set with the
argument `max_depth` (below we use decision stumps: `max_depth = 1`).
The argument `niter` specifies the number of iterations (epochs) to be
used (we set `niter = 1000`). Predictions (on the supplied test set) and
variable importance scores can be computed by setting `make_prediction =
TRUE` and `cal_imp = TRUE` in the argument `control`, as we do below.

For the fit computed with an initial `LADTree` with parameters chosen a
priori, we set its maximum depth to `max_depth_init = 2` and the minimum
number of observations per node at `min_leaf_size_init = 20`.

Finally, we also use as initial fit the best `LADTree` among those
constructed with all the possible combinations of `max_depth_init_set`
between 1 and 5 (`max_depth_init_set = 1:5`) and
`min_leaf_size_init_set` equal to 10, 20 or 30 (`min_leaf_size_init_set
= c(10,20,30)`).

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
                                            max_depth_init_set = 1:5,
                                            min_leaf_size_init_set = c(10,20,30),
                                            control = Boost.control(
                                              make_prediction =  TRUE,
                                              cal_imp = TRUE))
```

The parameters of the selected `LADTree` are returned in the `params`
entry of the returned object:

``` r
 print(model_RRBoost_cv_LADTree$params)
```

    ##    min_leafs max_depths
    ## 14        20          5

The predictive performance of each of the fits on the test set is stored
in the `value` entry of the returned objects. This is the test error
(using the criterion specified with the `error` argument of the `Boost`
call) evaluated at the early stopping time, which is determined using
the validation set to avoid overfitting (see paper for more details). In
this example the returned `value` is the test `rmse` evaluated at the
early stopping time. The user can evaluate more than one prediction
performance by passing a vector of strings to the argument `error` (e.g.
`Boost(..., error=c('rmse', 'rrmse'), ...)`). In that case the returned
`value` will be a vector corresponding to those prediction error
measures.

The best prediction performance on the test set was obtained by
selecting a deeper initial `LADTree`:

``` r
print(c(median = model_RRBoost_median$value, 
        LADTree = model_RRBoost_LADTree$value,
        cv_LADTree = model_RRBoost_cv_LADTree$value))
```

    ##     median    LADTree cv_LADTree 
    ##   4.622946   4.410887   3.443912

The variable selection scores are returned in the `var_importance`
entry.

``` r
 print(cbind(median = model_RRBoost_median$var_importance,
         LADTree = model_RRBoost_LADTree$var_importance,
         cv_LADTree = model_RRBoost_cv_LADTree$var_importance))
```

    ##                 median   LADTree cv_LADTree
    ## frequency    3.7844891 3.7814350 5.99593796
    ## angle        0.0651800 0.3201209 0.08082488
    ## chord_length 0.9628652 1.0416198 1.87940434
    ## velocity     0.5609014 0.5983525 0.65163300
    ## thickness    2.9582819 2.5745679 4.53945203

We note that for the top 2 explanatory variables are consistently found
to be `frequency` and `thickness`.

We can also separate the process of training the predictor, evaluating
it on a test set, and calculating the variable importance scores. In
this way we have the flexiblity to compute predictions on multiple test
sets and calculate variable importance when needed. In this case the
user needs to use `save_tree = TRUE` in the argument `control`. This
option includes the trainer weak learners in the returned object, which
are needed for calculating predictions and variable importance scores.
The following code illustrates this procedure:

``` r
model = Boost(x_train = xtrain, y_train = ytrain, 
              x_val = xval, y_val = yval, 
              # x_test = xtest, y_test = ytest,
              type = "RRBoost", error = "rmse", 
              y_init = "LADTree", max_depth = 1, niter = 1000, 
              control = Boost.control(max_depth_init = 2, 
                                      min_leaf_size_init = 20, 
                                      save_tree = TRUE, 
                                      make_prediction =  FALSE, #TRUE, 
                                      cal_imp = FALSE))

prediction <- cal_predict(model, x_test = xtest, y_test = ytest)
var_importance <-  cal_imp_func(model, x_val = xval, y_val= yval)
```

Sanity check: the results are the same as the ones obtained before:

``` r
print(all.equal(prediction$value, model_RRBoost_LADTree$value))
```

    ## [1] TRUE

``` r
print(all.equal(var_importance, model_RRBoost_LADTree$var_importance))
```

    ## [1] TRUE

Finally, we compare the predictive performance of `RRBoost` on the test
set with those of the other boosting algorithms implemented in the
package, namely: `L2Boost`, `LADBoost`, `MBoost`, and `Robloss`.

``` r
model_L2Boost = Boost(x_train = xtrain, y_train = ytrain, 
                              x_val = xval, y_val = yval, 
                              x_test = xtest, y_test = ytest, 
                              type = "L2Boost", error = "rmse", 
                              y_init = "median", max_depth = 1, 
                              niter = 1000, 
                              control = Boost.control(make_prediction =  TRUE,
                                                      cal_imp = TRUE))
model_LADBoost = Boost(x_train = xtrain, y_train = ytrain, 
                              x_val = xval, y_val = yval, 
                              x_test = xtest, y_test = ytest, 
                              type = "LADBoost", error = "rmse", 
                              y_init = "median", max_depth = 1, 
                              niter = 1000, 
                              control = Boost.control(make_prediction =  TRUE,
                                                      cal_imp = TRUE))
model_MBoost = Boost(x_train = xtrain, y_train = ytrain, 
                              x_val = xval, y_val = yval, 
                              x_test = xtest, y_test = ytest, 
                              type = "MBoost", error = "rmse", 
                              y_init = "median", max_depth = 1, 
                              niter = 1000, 
                              control = Boost.control(make_prediction =  TRUE,
                                                      cal_imp = TRUE))
model_Robloss = Boost(x_train = xtrain, y_train = ytrain, 
                              x_val = xval, y_val = yval, 
                              x_test = xtest, y_test = ytest, 
                              type = "Robloss", error = "rmse", 
                              y_init = "median", max_depth = 1, 
                              niter = 1000, 
                              control = Boost.control(make_prediction =  TRUE,
                                                      cal_imp = TRUE))
```

Note that `RRBoost` yields the best predictions on the test set:

``` r
print(c( RRBoost = model_RRBoost_cv_LADTree$value,  
         L2Boost = model_L2Boost$value, 
         LADTree = model_LADBoost$value,
         MBoost = model_MBoost$value, 
         Robloss = model_Robloss$value))
```

    ##  RRBoost  L2Boost  LADTree   MBoost  Robloss 
    ## 3.443912 4.521834 4.520258 4.483478 4.648868
