
#' Computes SCOPE predictions
#'
#' @description Computes SCOPE predictions on new data.
#'
#' @param object SCOPE model as outputted by scope. Must have simply.the.best = TRUE
#' @param newdata New covariates on which to make predictions. Must be of the same form as the model was trained on
#' @param interceptxlinear If intercept is to be included in the model and is not in the column span of the continuous variables in x, set to FALSE (default).
#' Must match format of training data.
#' @param ... Additional arguments to pass to other \code{predict} methods
#'
#' @return Returns n-length vector of predictions
#'
#' @seealso \code{\link{scope}}
#'
#' @export
predict.scope = function( object, newdata, interceptxlinear = FALSE,...)
{
  model<-object  
  model = model$beta.best
  p = dim(newdata)[ 2 ]
  factor_ind = rep(F, p)
  for ( j in 1:p ) {
    if (( class(newdata[ , j ]) != "numeric" ) && ( class(newdata[ , j ]) != "numeric" )) factor_ind[ j ] = T
  }
  xshrink = data.frame(newdata[ , factor_ind, drop = F ])
  if ( interceptxlinear == T ) {
    xlinear = as.matrix(newdata[ , !factor_ind, drop = F ])
  } else {
    xlinear = as.matrix(cbind(1, newdata[ , !factor_ind, drop = F ]))
  }
  for ( j in 1:dim(xshrink)[ 2 ] ) xshrink[ , j ] = as.factor(xshrink[ , j ])
  predvec = xlinear %*% model[[ 1 ]]
  for ( j in 1:dim(xshrink)[ 2 ] ) {
    predvec = predvec + model[[ 2 ]][[ j ]][ xshrink[ , j ] ]
  }
  return(predvec)
}
