#' get FML values
#'
getFML <- function(regmodel, use_unbiasedCov = FALSE){
  expCov <- mxGetExpected(regmodel, component = "covariance", subname = "BaseModel")
  expMean <- mxGetExpected(regmodel, component = "means", subname = "BaseModel")
  means <- regmodel$BaseModel$data$means
  unbiasedCov <- regmodel$BaseModel$data$observed
  #biasedCov <- unbiasedCov*((regmodel$BaseModel$data$numObs-1)/regmodel$BaseModel$data$numObs)
  if(use_unbiasedCov == FALSE)
  {
    obsCov <- unbiasedCov*((regmodel$BaseModel$data$numObs-1)/regmodel$BaseModel$data$numObs)
  }
  else{
    obsCov <- regmodel$BaseModel$data$observed
  }
  obs_manif <- length(regmodel$BaseModel$manifestVars)
  if(is.na(means)| is.null(means)){

    F_ML <- .5*(log(det(expCov)) + tr(obsCov %*% solve(expCov)) - obs_manif - log(det(obsCov)))
  }
  else{
    F_ML <- .5*(log(det(expCov)) + tr(obsCov %*% solve(expCov)) + t(means - expMean) %*% solve(expCov) %*% (means - expMean) - obs_manif - log(det(obsCov)))
  }
  penalty_comb <- regmodel$penalty_combined$result / regmodel$BaseModel$data$numObs
  FML_pen <- 2*F_ML + penalty_comb# equivalent to formula in Jacobucci (2016)
  FML_pen_05 <- F_ML + .5*penalty_comb # reason for .5: jacobucci reports .5*(penalty + FML)in out$values

  return_values <- matrix(c(F_ML,FML_pen,FML_pen_05), ncol = 3, nrow = 1)
  colnames(return_values) <- c(".5*FML", "FML + penalty", ".5*(FML + penalty)")

  return(return_values)}
