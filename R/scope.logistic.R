
#' Computes solution for SCOPE logistic models
#'
#' @description Computes solution for SCOPE logistic models, computing along a path and iterating local quadratic approximations at each point. Performs K-fold cross-validation
#' for regularisation parameter lambda and can incorporate both linear and categorical (including logical) variables. This function uses a Proximal Newton algorithm and is not guaranteed to converge. It is
#' recommended for developer use only. See Stokell, Shah and Tibshirani (2020).
#'
#' @param x Data frame of covariates. Can include a mix of continuous and categorical covariates (no scaling of continuous covariates is performed within the program). By default an
#' intercept will be added; see interceptxlinear
#' @param y Response vector of length n
#' @param gamma	Concavity parameter in MCP; see Zhang (2010) Nearly unbiased estimation with minimax concave penalty
#' @param sd Standard deviation of noise used for calibration of parameter lambda. This is recommended to be left alone
#' @param AIC Controls whether information criterion for early stopping is AIC (=TRUE) or mBIC (=FALSE)
#' @param mBICconstant If using mBIC, this is the parameter m
#' @param default.lambdaseq If using automatically generated sequence of lambda values set to TRUE. Do not set to FALSE without good reason
#' @param default.length Length of sequence of automatically generated lambda values
#' @param lambda.min.ratio Ratio of largest to smallest value on sequence of lambda values
#' @param lambdaseq If default.lambdaseq = FALSE then add path of lambda values here
#' @param TerminateEpsilon Epsilon for convergence criterion, is multiplied by null deviance to get terminate criterion for objective value
#' @param interceptxlinear If intercept is to be included in the model and is not in the column span of the continuous variables in x, set to FALSE (default).
#' @param max.iter Maximum number of iterations at each local quadratic approximation
#' @param max.out.iter Maximum number of local quadratic approximations at each value of lambda
#' @param BICearlystop	If information criterion is to be used to stop computation early, set to TRUE
#' @param BICterminate	Specifies how many values of lambda to be computed after the minimum value of the information criterion is reached
#' @param silent If FALSE then progress updates will be printed as solutions are computed
#' @param K Number of folds in cross-validation. If K = 1, no cross-validation is performed
#' @param return.full.beta If TRUE then beta.full will be returned, else just the cross-validation optimal beta.best will be returned
#' @param simply.cross.validated If TRUE then cross-validation scores for each value of lambda will be returned, but not the estimates themselves
#' @param grid.safe As the automatically generated sequence of lambda values is adjusted during the first fold but fixed thereafter. For subsequent folds, this sets computation to
#' begin at a larger value of lambda to ensure that the first solution along the path is zero so as to maintain the advantages of the pathwise approach. This specifies how many larger
#' values of lambda should be used
#' @param blockorder By default the order in block coordinate descent is randomly sampled. Alternatively a permutation vector can be included here
#' @param FoldAssignment By default the assignments for cross-validation are randomly sampled automatically. Alternatively assignments can be included here
#' @param zero.penalty Fits logistic regression model with zero penalty. Primarily used for debugging, do not set to TRUE
#'
#' @return A list of objects. Some may not be returned depending on value of arguments K, simply.cross.validated, return.full.beta.
#' \itemize{
#' \item lambdaseq - A matrix of the values of lambda used to compute the solution path. Columns correspond to different points on the path, rows correspond to the categorical
#' variables. Lambda is scaled depending on the number of categories present in the data.
#' \item cverrors - Provided K > 1 then the cross-validation error for each point on the grid will be returned
#' \item beta.full - Contains full solution path. If K > 1 then will only be returned if simply.cross.validated = FALSE and return.full.beta = TRUE. First object [[ 1 ]] is
#' coefficients of continuous variables, [[ 2 ]] is a list of coefficients for categorical variables
#' \item beta.best - Contains solution at CV-optimal point. Requires K > 1 to be returned. This must not be NULL in order to use predict.scope. First object [[ 1 ]] is coefficients of
#' continuous variables, [[ 2 ]] is a list of coefficients for categorical variables
#' \item fold.assign - Contains fold assignments for cross-validation
#' }
#' 
#' @examples
#' \dontrun{
#' x = UniformDesignMatrix(200, 5, 5)
#' y = (as.integer(x[ , 1 ]) < 3) + rnorm(200)
#' y = as.integer(y > 0.8)
#' scope_mod = scope.logistic(x, y)
#' x_new = UniformDesignMatrix(10, 5, 5)
#' predict(scope_mod, x_new)
#' }
#'
#' @references \insertRef{stokell2020modelling}{CatReg}
#'
#' @export
scope.logistic = function ( x, y, gamma = 100, sd = NULL, AIC = TRUE, mBICconstant = 5, default.lambdaseq = TRUE, default.length = 160, lambda.min.ratio = 0.005,
                            lambdaseq = NULL, TerminateEpsilon = 1e-7, interceptxlinear = FALSE, max.iter = 1000, max.out.iter = 1000, BICearlystop = TRUE, BICterminate = 25,
                            silent = TRUE, K = 5, return.full.beta = FALSE, simply.cross.validated = FALSE, grid.safe = 10, blockorder = NULL, FoldAssignment = NULL, zero.penalty = FALSE ) {
  # y is response vector, xlinear is numeric matrix of covariates that will just be used for the unpenalized linear part of the regression
  # xshrink is a factor data frame where each column corresponds to a categorical variable, categories are numbered 1, 2, 3, ...
  inv.logit = function ( x ) return(exp(x) / ( 1 + exp(x) ))

  logit = function( x ) return(log( x / ( 1 + x )))
  # Important note is that we need to have no empty categories in 1, ..., c_j; if there are any empty ones in advance then ensure that
  # relabelling is performed before calling scope otherwise we'll get errors (I'll code this exception up myself in due course 7AUG NOW DONE)
  scopemod = list()
  attr(scopemod,"class")<-"scope.logistic"  
  scopemod$cverrors = NULL
  scopemod$lambdaseq = NULL
  scopemod$beta.full = NULL
  scopemod$beta.best = NULL
  scopemod$fold.assign = NULL
  

  n = length(y)
  p = dim(x)[ 2 ]
  factor_ind = rep(F, p)
  for ( j in 1:p ) {
    if (( class(x[ , j ]) != "numeric" ) && ( class(x[ , j ]) != "numeric" )) factor_ind[ j ] = T
  }
  xshrink = data.frame(x[ , factor_ind, drop = F ])
  if ( interceptxlinear == T ) {
    xlinear = as.matrix(x[ , !factor_ind, drop = F ])
  } else {
    xlinear = as.matrix(cbind(1, x[ , !factor_ind, drop = F ]))
  }
  plinear = dim(xlinear)[ 2 ]
  pshrink = dim(xshrink)[ 2 ]
  for ( j in 1:pshrink ) xshrink[ , j ] = as.factor(xshrink[ , j ])
  interceptxlinear = TRUE
  
  P = solve(t(xlinear) %*% xlinear) %*% t(xlinear)
  catnumbers = rep(0, pshrink)
  mcpContribution = rep(0, pshrink)
  catnames = list()
  # Store names (and number of) levels in the categorical variables
  for ( j in 1:pshrink ) {
    catnames[[ j  ]] = levels(xshrink[ , j ])
    catnumbers[ j ] = length(catnames[[ j ]])
  }

  if ( default.lambdaseq == TRUE ) {
    pathlength = default.length
  } else {
    pathlength = dim(lambdaseq)[ 2 ]
  }

  coefficientshrink = list()
  # Sample (and then fix) uniformly random permutation of variables for order of block descent
  blockorder = sample(1:pshrink)
  for ( j in 1:pshrink ) {
    coefficientshrink[[ j ]] = matrix(0, catnumbers[ j ], pathlength)
    rownames(coefficientshrink[[ j ]]) = catnames[[ j ]]
  }
  coefficientlinear = matrix(0, plinear, pathlength)
  beta = rep(0, plinear)
  subaverages = list()
  weights = list()
  weightsbool = list()
  for ( j in 1:pshrink ) {
    # initialize subaverages object, compute the weights (recall they're fixed a priori so this needs doing once)
    weights[[ j ]] = rep(0, catnumbers[ j ])
    weightsbool[[ j ]] = rep(FALSE, catnumbers[ j ])
    subaverages[[ j ]] = rep(0, catnumbers[ j ])
    for ( k in 1:catnumbers[ j ] ) {
      weights[[ j ]][ k ] = sum(xshrink[ , j ] == catnames[[ j ]][ k ]) / n
    }
    # weightsbool stores when a level has no observations, so that we never end up dividing by zero
    weightsbool[[ j ]] = (weights[[ j ]] > 0)
  }
  partialresiduals = y

  # Now a short routine to decide on the scaling for default tuning parameters
  beta = P %*% partialresiduals
  partialresiduals = partialresiduals - xlinear %*% beta
  minstdev = Inf
  for ( j in 1:pshrink ) {
    subaverages[[ j ]] = tapply(partialresiduals, xshrink[ , j ], mean)[ catnames[[ j ]] ]
    minstdev = min(minstdev, sqrt(var(subaverages[[ j ]][ weightsbool[[ j ]] ] * sqrt(n * weights[[ j ]][ weightsbool[[ j  ]] ]))))
  }

  if ( is.null(sd) ) {
    sd = 0.0125 * minstdev # We want this to be an underestimate; if we don't fit null model at pathpoint 1, we start again with 2x value
  }


  if ( default.lambdaseq == TRUE ) {
    baseseq = as.matrix(exp(seq(log( 8 * sd / sqrt(n)), log(lambda.min.ratio * 8 * sd / sqrt(n)), length = default.length)))
    lambdaseq = t(baseseq %*% catnumbers^(0.5)) # rescales the sequence of lambda by this, such that if all variables have no signal (and have
    # equal size categories) then they'll have consistent scaling for recovery of true (intercept-only) solution
  }

  # (Leftover from testing the LQA routine -- testing zero-penalty against simple examples in glm)
  if ( zero.penalty == TRUE ) {
    lambdaseq = lambdaseq[ , 1:2 ]
    lambdaseq[ , 2 ] = rep(0, pshrink)
  }
  # Split up the observations into different folds
  if ( K > 1 ) {
    if ( is.null(FoldAssignment) ) {
      FoldAssignment = as.integer(sample(ceiling((1:n)*K/n)))
      scopemod$fold.assign = FoldAssignment #16SEP
    }
    if ( default.lambdaseq == TRUE ){
      cverrors = matrix(0, n, default.length)
    } else {
      cverrors = matrix(0, n, dim(lambdaseq)[ 2 ])
    }
    removecounter = 0
    counter = 0

    for ( k in 1:K ) {
      if ( silent == FALSE ) print(paste0("Fold ", k))
      yfold = y[ (FoldAssignment != k), drop = FALSE ]
      xlinearfold = xlinear[ (FoldAssignment != k), , drop = FALSE ]
      xshrinkfold = xshrink[ (FoldAssignment != k), , drop = FALSE ]
      yremove = y[ (FoldAssignment == k), drop = FALSE ]
      xlinearremove = xlinear[ (FoldAssignment == k), , drop = FALSE ]
      xshrinkremove = xshrink[ (FoldAssignment == k), , drop = FALSE ]
      nremove = length(yremove)
      keepidentifier = rep(TRUE, nremove)
      # If a category is not present in the training data, remove all observations including it from the test data
      for ( i in 1:nremove ) {

        for ( j in 1:pshrink ) {

          if ( xshrinkremove[ i, j ] %in% xshrinkfold[ , j ] == FALSE ) {
            keepidentifier[ i ] = FALSE
            removecounter = removecounter + 1

          }
        }
      }
      yremove = yremove[ keepidentifier, drop = FALSE ]
      xlinearremove = xlinearremove[ keepidentifier, , drop = FALSE ]
      xshrinkremove = xshrinkremove[ keepidentifier, , drop = FALSE ]

      # Compute a single solution path for this fold
      cvsolution = core_scope.logistic(yfold, xlinearfold, xshrinkfold, blockorder, k, gamma, AIC, mBICconstant, lambdaseq, TerminateEpsilon,
                                           interceptxlinear, max.iter, max.out.iter, BICearlystop, BICterminate,
                                           silent, grid.safe)
      # Retain the scaling of the tuning parameter from the first fold
      if ( k == 1 ) {
        lambdaseq = cvsolution[[ 4 ]]
      }
      cvtemp = xlinearremove %*% cvsolution[[ 1 ]]
      for ( j in 1:pshrink ) {
        cvtemp = cvtemp + cvsolution[[ 2 ]][[ j ]][ xshrinkremove[ , j], ]
      }
      cverrorstemp = abs(yremove - inv.logit(cvtemp))


      cverrors[ (counter + 1):(counter + length(yremove)), 1:dim(cverrorstemp)[ 2 ] ] = as.numeric(cverrorstemp)

      counter = counter + length(yremove)
    }

    cverrors = as.matrix(cverrors[ 1:(n - removecounter), 1:dim(cverrorstemp)[ 2 ] ])
    if ( removecounter > 0 ) {
      warning(paste0(removecounter, " observations removed from test sets; number of evaluated predictions is ", n - removecounter, "."))
    }
    # Computes cross-validation errors and selects optimal lambda
    cverrors = colMeans(cverrors > 0.5)
    cverrors[ is.na(cverrors) ] = Inf # If we get NA as result of fitting saturated model
    scopemod$cverrors = cverrors
    if ( simply.cross.validated == TRUE ) {
      # If we only want to return the results of the cross-validation. Used if, for example, one wants to also cross-validate over gamma
      
      return(scopemod)
      # break
    }
    pathlengthfinal = which.min(cverrors)
    lambdaseqused = lambdaseq
    lambdaseq = lambdaseq[ , 1:pathlengthfinal ] # Terminated the tuning parameter sequence at optimal value
    
    if ( silent == FALSE ) print(paste0("Minimal cross-validation error = ", min(cverrors), " at pathpoint ", pathlengthfinal))
    
    scopemod$lambdaseq = lambdaseqused
    
    # Compute final solution with cross-validated tuning parameter
    fullsolution = core_scope.logistic(y, xlinear, xshrink, blockorder, 2, gamma, AIC, mBICconstant, lambdaseq, TerminateEpsilon,
                                            interceptxlinear, max.iter, max.out.iter, BICearlystop, BICterminate,
                                            silent, grid.safe)

    if ( return.full.beta == TRUE ) {
      scopemod$beta.full = list()
      scopemod$beta.full[[ 1 ]] = fullsolution[[ 1 ]]
      scopemod$beta.full[[ 2 ]] = fullsolution[[ 2 ]]
    }
    
    
    
    # This is if we only want to print out the solution at CV-optimal lambda (default TRUE)
    fullsolution[[ 1 ]] = matrix(fullsolution[[ 1 ]], plinear, pathlengthfinal)[ , pathlengthfinal ]
    for ( j in 1:pshrink ) {
      
      fullsolution[[ 2 ]][[ j ]] = as.matrix(fullsolution[[ 2 ]][[ j ]])[ , pathlengthfinal ]
      
      
    }
    fullsolution = list(fullsolution[[ 1 ]], fullsolution[[ 2 ]])
    scopemod$beta.best = list()
    scopemod$beta.best[[ 1 ]] = fullsolution[[ 1 ]]
    scopemod$beta.best[[ 2 ]] = fullsolution[[ 2 ]]
    
    
    
    return(scopemod)
    
    
  } else {

    solution = core_scope.logistic(y, xlinear, xshrink, blockorder, 1, gamma, AIC, mBICconstant, lambdaseq, TerminateEpsilon,
                                       interceptxlinear, max.iter, max.out.iter, BICearlystop, BICterminate,
                                       silent, grid.safe)

    pathlengthfinal = dim(solution[[ 2 ]][[ 1 ]])[ 2 ]

    
    if ( plinear > 1 ) {
      solution[[ 1 ]] = solution[[ 1 ]][ , 1:pathlengthfinal ]
    } else {
      solution[[ 1 ]] = solution[[ 1 ]][ 1:pathlengthfinal ]
    }
    
    for ( j in 1:pshrink ) {
      if ( pathlengthfinal > 1 ) {
        solution[[ 2 ]][[ j ]] = solution[[ 2 ]][[ j ]][ , 1:pathlengthfinal ]
      }
      
    }
    
    
    fullsolution = list()
    fullsolution$beta.full[[ 1 ]] = solution[[ 1 ]]
    fullsolution$beta.full[[ 2 ]] = solution[[ 2 ]]
    scopemod$beta.full = fullsolution
    lambdaseq = lambdaseq[ , 1:pathlengthfinal ]
    scopemod$lambdaseq = lambdaseq
    
    return(scopemod)
  }

}
