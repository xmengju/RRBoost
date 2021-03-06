% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Boost.R
\name{cal_imp_func}
\alias{cal_imp_func}
\title{Variable importance scores for the robust boosting algorithm RRBoost}
\usage{
cal_imp_func(model, x_val, y_val, trace = FALSE)
}
\arguments{
\item{model}{an object returned by \code{\link{Boost}}}

\item{x_val}{predictor matrix for validation data (matrix/dataframe)}

\item{y_val}{response vector for validation data (vector/dataframe)}

\item{trace}{logical indicating whether to print the variable under calculation for monitoring progress (defaults to \code{FALSE})}
}
\value{
a vector of permutation variable importance scores (one per explanatory variable)
}
\description{
This function  calculates variable importance scores for a previously
computed \code{RRBoost} fit.
}
\details{
This function computes permutation variable importance scores
given an object returned by \code{\link{Boost}} and a validation data set.
}
\examples{
\dontrun{
data(airfoil)
n <- nrow(airfoil)
n0 <- floor( 0.2 * n )
set.seed(123)
idx_test <- sample(n, n0)
idx_train <- sample((1:n)[-idx_test], floor( 0.6 * n ) )
idx_val <- (1:n)[ -c(idx_test, idx_train) ]
xx <- airfoil[, -6]
yy <- airfoil$y
xtrain <- xx[ idx_train, ]
ytrain <- yy[ idx_train ]
xval <- xx[ idx_val, ]
yval <- yy[ idx_val ]
xtest <- xx[ idx_test, ]
ytest <- yy[ idx_test ]
model = Boost(x_train = xtrain, y_train = ytrain,
     x_val = xval, y_val = yval,
     type = "RRBoost", error = "rmse",
     y_init = "LADTree", max_depth = 1, niter = 1000,
     control = Boost.control(max_depth_init = 2,
           min_leaf_size_init = 20, save_tree = TRUE,
           make_prediction =  FALSE, cal_imp = FALSE))
var_importance <-  cal_imp_func(model, x_val = xval, y_val= yval)
}

}
\author{
Xiaomeng Ju, \email{xmengju@stat.ubc.ca}
}
