#' get FML values
#'
#'Note: laremm is based on the R package \pkg{regsem}. Because of the early status of laremm, it is recommended to use regsem instead!
#'Computes FML values corresponding to the values reported in \pkg{regsem} for covariance based models
#'
#' @param regmodel already run regularized model
#' @param use_unbiasedCov wheather to use the unbiased covariance; default is FALSE which is equivalent to regsem
#' @param zeroThresh threshold for setting regularized parameters to zero. Default is .001 similar to \pkg{regsem}
#' @param setZero should parameters below zeroThresh be set to zero in all fit calculations. Default is FALSE, similar to \pkg{regsem}
#' @author Jannik Orzek
#' @import OpenMx ctsem
#'@export
getFML <- function(regmodel, zeroThresh = .001, use_unbiasedCov = FALSE, setZero = FALSE){
  if(setZero){
  # set regularized parameters below zeroThresh to zero
  newModel <- getEstimatedParameters(regmodel, zeroThresh = zeroThresh, setZero = setZero)$new_model
  expCov <- mxGetExpected(newModel, component = "covariance")
  expMean <- mxGetExpected(newModel, component = "means")
  }else{
    newModel <- regmodel
    expCov <- mxGetExpected(newModel, component = "covariance", subname = "BaseModel")
    expMean <- mxGetExpected(newModel, component = "means", subname = "BaseModel")
  }

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
  penalty_comb <- (regmodel$regfit_algebra$result - regmodel$BaseModel$fitfunction$result[[1]]) / regmodel$BaseModel$data$numObs
  FML_pen <- 2*F_ML + penalty_comb# equivalent to formula in Jacobucci (2016)
  FML_pen_05 <- F_ML + .5*penalty_comb # reason for .5: jacobucci reports .5*(penalty + FML)in out$values

  return_values <- matrix(c(F_ML,FML_pen,FML_pen_05), ncol = 3, nrow = 1)
  colnames(return_values) <- c(".5*FML", "FML + penalty", ".5*(FML + penalty)")

  ret = list("return_values" = return_values, "return_Model" = newModel)
  return(ret)}
