% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.scope.logistic.R
\name{predict.scope.logistic}
\alias{predict.scope.logistic}
\title{Computes SCOPE logistic predictions}
\usage{
\method{predict}{scope.logistic}(object, newdata, probs = TRUE, include_intercept = TRUE, ...)
}
\arguments{
\item{object}{SCOPE model as outputted by scope.logistic. Must have simply.the.best = TRUE}

\item{newdata}{New covariates on which to make predictions. Must be of the same form as the model was trained on}

\item{probs}{If TRUE returns probabilities, if FALSE returns binary predictions}

\item{include_intercept}{If TRUE, a column of 1s will be added to the (continuous) design matrix. Must match format of training data.}

\item{...}{Additional arguments to pass to other \code{predict} methods}
}
\value{
Returns n-length vector of predictions
}
\description{
Computes SCOPE logistic predictions on new data
}
\seealso{
\code{\link{scope.logistic}}
}
