#' compute fit measures for laremm
#'
#' Note: laremm is based on the R package \pkg{regsem}. Because of the early status of laremm, it is recommended to use regsem instead!
#'
#' @param regmodel already run regularized model
#' @param model_type specify the type of model provided: ctsem or mxModel
#' @param fitfun fitfunction to be used in the fitting procedure. Either FML or FIML
#' @param cvsample mxData object with test sample data. Has to be of same data_type as the training data set
#' @param zeroThresh threshold for setting regularized parameters to zero. Default is .001 similar to \pkg{regsem}
#' @param setZero should parameters below zeroThresh be set to zero in all fit calculations. Default is FALSE, similar to \pkg{regsem}
#'
#' @author Jannik Orzek
#' @import OpenMx ctsem
#'
#' @export



getFitMeasures <- function (regmodel, model_type = "mxModel", fitfun = "FIML", cvsample = NULL, zeroThresh = .001, setZero= FALSE){

  # define return value:
  if(fitfun == "FML"){

    return_value <- data.frame("estimated_params" = NA,
                               "half_FML"= NA,
                               "FML_plus_penalty"= NA,
                               "half_FML_plus_penalty"= NA,
                               "m2LL"= NA,
                               "AIC"= NA,
                               "BIC"= NA,
                               "CV m2LL" = NA,
                               "CV AIC"= NA,
                               "CV BIC"= NA
    )
    # get F_ML:
    Fit <- getFML(regmodel, zeroThresh = zeroThresh, setZero = setZero)
    FML <- Fit$return_value

    return_value$half_FML <- FML[1,1]
    return_value$FML_plus_penalty <- FML[1,2]
    return_value$half_FML_plus_penalty <- FML[1,3]

  } else{
    return_value <- data.frame("estimated_params" = NA,
                               "m2LL"= NA,
                               "AIC"= NA,
                               "BIC"= NA,
                               "CV m2LL" = NA,
                               "CV AIC"= NA,
                               "CV BIC"= NA
    )
  }

  # get the number of estimated parameters:
  EstimatedParam <- getEstimatedParameters(regmodel, zeroThresh = zeroThresh, setZero = setZero)
  return_value$estimated_params <- EstimatedParam$`estimated parameters`[[1]]

  ### compute Fit Indices:

  # get -2LogL:
  if(fitfun == "FML"){
    ExpCov <- mxGetExpected(EstimatedParam$new_model, "covariance")
    ExpMean <- mxGetExpected(EstimatedParam$new_model, "means")

    # get the data set:
    ObsCov <- regmodel$BaseModel$data$observed
    NObs <- regmodel$BaseModel$data$numObs
    NManif <- length(regmodel$BaseModel$manifestVars)

    # compute the m2LogL:

    return_value$m2LL <- getM2LogL_cov(ExpCov = ExpCov, ObsCov = ObsCov, NObs = NObs, NManif = NManif)

  }else if(fitfun == "FIML"){
    #get the expected covariance/means
    if(model_type == "mxModel"){
    ExpCov <- mxGetExpected(EstimatedParam$new_model, "covariance")
    ExpMean <- mxGetExpected(EstimatedParam$new_model, "means")

    # get the data set:
    ObsData <- regmodel$BaseModel$data$observed

    # compute the m2LogL:
    return_value$m2LL <- getM2LogL_FIML(ExpCov = ExpCov, ExpMean = ExpMean, ObsData = ObsData)
    } else if(model_type == "ctsem"){
      temp_ct_model <- mxRun(EstimatedParam$new_model, useOptimizer = F, silent = T)
      return_value$m2LL <- temp_ct_model$fitfunction$result[[1]]
    }else{print("Error: Unknown model_type")}
  }

  # AIC
  AIC <- return_value$m2LL + 2* return_value$estimated_params
  return_value$AIC <- AIC
  # BIC
  BIC <- return_value$m2LL + log(regmodel$BaseModel$data$numObs) * return_value$estimated_params
  return_value$BIC <- BIC

  #### if Cross -Validation is used ####

  ## CV m2LL
  if( !is.null(cvsample) ){ # if cvsample (as mxData) is provided:
    if(fitfun == "FML"){
      ExpCov <- mxGetExpected(EstimatedParam$new_model, "covariance")
      ExpMean <- mxGetExpected(EstimatedParam$new_model, "means")

      # get the data set:
      ObsCov <- cvsample$observed
      NObs <- cvsample$numObs
      NManif <- length(regmodel$BaseModel$manifestVars)

      # compute the m2LogL:

      return_value$CV.m2LL <- getM2LogL_cov(ExpCov = ExpCov, ObsCov = ObsCov, NObs = NObs, NManif = NManif)

    }else if(fitfun == "FIML"){
      if(model_type == "mxModel"){

        #get the expected covariance/means

        ExpCov <- mxGetExpected(EstimatedParam$new_model, "covariance")
        ExpMean <- mxGetExpected(EstimatedParam$new_model, "means")

        # get the data set:
        ObsData <- cvsample$observed

        # compute the m2LogL:
        return_value$CV.m2LL <- getM2LogL_FIML(ExpCov = ExpCov, ExpMean = ExpMean, ObsData = ObsData)
      }else if(model_type=="ctsem"){

        CVModel <- EstimatedParam$new_model
        CVModel$data <- cvsample
        fit_CVModel <- mxRun(CVModel, useOptimizer = F, silent = T)
        return_value$CV.m2LL <- fit_CVModel$fitfunction$result[[1]]
      }
    }

    # AIC
    CV_AIC <- return_value$CV.m2LL + 2* return_value$estimated_params
    return_value$CV.AIC <- CV_AIC
    # BIC
    CV_BIC <- return_value$CV.m2LL + log(cvsample$numObs)* return_value$estimated_params
    return_value$CV.BIC <- CV_BIC
  }
  ret <- list("return_value" = return_value, "return_Model" = EstimatedParam$new_model)

  return(ret)


}
