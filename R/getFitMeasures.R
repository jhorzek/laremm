#' compute fit measures for laremm
#'
#' Note: laremm is based on the R package \pkg{regsem}. Because of the early status of laremm, it is recommended to use regsem instead!
#'
#' @param regmodel already run regularized model
#' @param model_type specify the type of model provided: ctsem or mxModel
#' @param fitfun fitfunction to be used in the fitting procedure. Either FML or FIML
#' @param ncp_rmsea should rmsea and ncp be computed? Only possible for covariance based models
#' @param cvsample mxData object with test sample data. Has to be of same data_type as the training data set
#' @param satmod saturated model. necessary for computation of ncp and rmsea in FIML models. In many cases, the OpenMx mxRefModels(model, run =TURE) function can be used to build this saturated model. Make sure to only provide the fitted saturated model, not the indipendence model
#' @param cv_satmod saturated model for cross validation. This model has to be based on the cv sample
#'
#' @export



getFitMeasures <- function (regmodel, model_type = "mxModel", fitfun = "FIML", ncp_rmsea = FALSE, cvsample = NULL, satmod= NULL, cv_satmod =NULL){

  return_value <- data.frame("estimated_params" = NA, "mxAIC"= NA, "lavaan_AIC"= NA,
                             "mxBIC"= NA, "lavaan_BIC"= NA, "ncp" = NA, "rmsea"= NA,
                             "CV_m2LL"= NA,"CV_AIC"= NA,"CV_BIC"= NA,
                             "CV_rmsea"= NA , "CV_ncp" = NA, "half_FML"= NA, "FML_plus_penalty"= NA,
                             "half_FML_plus_penalty"= NA)

  # get FML
  if(fitfun == "FML"){
  FML <- getFML(regmodel)

  return_value$half_FML <- FML[1,1]
  return_value$FML_plus_penalty <- FML[1,2]
  return_value$half_FML_plus_penalty <- FML[1,3]
  } else{
    return_value$half_FML <- NA
    return_value$FML_plus_penalty <- NA
    return_value$half_FML_and_penalty <- NA}

  ### compute AIC and BIC:
  # free parameter:
  matrices <- regmodel$BaseModel$matrices
  temp_model <- regmodel$BaseModel
  new_model <- regmodel$BaseModel
  estimated_params <- 0
  # values are only set to 0 and set as free = FALSE, if ther absolute value is smaller than .001 and if they have been penalized.
  # this procedure is equivalent to the one used in regsem
  # if penalty is on A:
  if(any( regmodel$selectedA_values$values == 1)){
    zero_matrix <- which(abs(matrices$A$values) < .001 & matrices$A$free & regmodel$selectedA_values$values == 1)
    matrices$A$values[zero_matrix] = 0 # values smaller than .001 and freely estimated are set to 0
    matrices$A$free[zero_matrix] = FALSE # these are assumed as fixed to compute the BIC / AIC
    new_model <- mxModel(new_model,  matrices$A)
  }
  # if penalty is on S
  if (any( regmodel$selectedS_values$values == 1)){
    zero_matrix <- which(abs(matrices$S$values) < .001 & matrices$S$free & regmodel$selectedS_values$values == 1)
    matrices$S$values[zero_matrix] = 0 # values smaller than .001 and freely estimated are set to 0
    matrices$S$free[zero_matrix] = FALSE # these are assumed as fixed to compute the BIC / AIC
    new_model <- mxModel(new_model,  matrices$S)
  }
  # if penalty is on DRIFT
  if (any( regmodel$selectedDrifts_values$values == 1)){
    zero_matrix <- which(abs(matrices$DRIFT$values) < .001 & matrices$DRIFT$free & regmodel$selectedDrifts_values$values == 1)
    matrices$DRIFT$values[zero_matrix] = 0 # values smaller than .001 and freely estimated are set to 0
    matrices$DRIFT$free[zero_matrix] = FALSE # these are assumed as fixed to compute the BIC / AIC
    new_model <- mxModel(new_model,  matrices$DRIFT)
  }
  # if penalty is on Lambda
  if (any( regmodel$selectedLambda_values$values == 1)){
    zero_matrix <- which(abs(matrices$LAMBDA$values) < .001 & matrices$LAMBDA$free & regmodel$selectedLambda_values$values == 1)
    matrices$LAMBDA$values[zero_matrix] = 0 # values smaller than .001 and freely estimated are set to 0
    matrices$LAMBDA$free[zero_matrix] = FALSE # these are assumed as fixed to compute the BIC / AIC
    new_model <- mxModel(new_model,  matrices$LAMBDA)
  }
  # if penalty is on Cint
  if (any( regmodel$selectedCint_values$values == 1)){
    zero_matrix <- which(abs(matrices$CINT$values) < .001 & matrices$CINT$free & regmodel$selectedCint_values$values == 1)
    matrices$CINT$values[zero_matrix] = 0 # values smaller than .001 and freely estimated are set to 0
    matrices$CINT$free[zero_matrix] = FALSE # these are assumed as fixed to compute the BIC / AIC
    new_model <- mxModel(new_model,  matrices$CINT)
  }
  # if penalty is on TDPREDEFFECT
  if (any( regmodel$selectedTdpredeffect_values$values == 1)){
    print("not yet implemented")
    #zero_matrix <- which(abs(matrices$$values) < .001 & matrices$CINT$free & regmodel$selectedTdpredeffect_values$values == 1)
    #matrices$CINT$values[zero_matrix] = 0 # values smaller than .001 and freely estimated are set to 0
    #matrices$CINT$free[zero_matrix] = FALSE # these are assumed as fixed to compute the BIC / AIC
    #new_model <- mxModel(new_model,  matrices$CINT)
  }
  # if penalty is on TDPREDEFFECT
  if (any( regmodel$selectedDiffusion_values$values == 1)){
    print("not yet implemented")
    #zero_matrix <- which(abs(matrices$$values) < .001 & matrices$CINT$free & regmodel$selectedDiffusion_values$values == 1)
    #matrices$CINT$values[zero_matrix] = 0 # values smaller than .001 and freely estimated are set to 0
    #matrices$CINT$free[zero_matrix] = FALSE # these are assumed as fixed to compute the BIC / AIC
    #new_model <- mxModel(new_model,  matrices$CINT)
  }

  # get estimated parameters:
  for (matrix in matrices){
    estimated_params <- estimated_params + sum(matrix$free)
  }
  return_value$estimated_params <- estimated_params

  #m2LL_no_pen <- regmodel$regfit_algebra$result - regmodel$penalty_combined$result # m2LL without the penalty term

  #in sample AIC:
  if(model_type == "mxModel"){
  AIC_new <- AIC(mxRun(new_model, useOptimizer = F, silent = T))
  if(fitfun == "FML"){
  lavAIC <- AIC_new+regmodel$BaseModel$data$numObs*length(regmodel$BaseModel$manifestVars)*log(2*pi)}
  else if(fitfun == "FIML"){
    lavAIC <- AIC_new # identical if FIML is used!
  }
  } # equivalent to lavaans AIC result
  else if(model_type == "ctsem"){
    AIC_new <- AIC(mxRun(new_model, useOptimizer = F, silent = T))
    lavAIC <- NA
    }
  return_value$mxAIC <- AIC_new
  return_value$lavaan_AIC <- lavAIC

  # in sample BIC:
  if(model_type == "mxModel"){
    BIC_new <- BIC(mxRun(new_model, useOptimizer = F, silent = T))
    if(fitfun == "FML"){
    lavBIC <- BIC_new + regmodel$BaseModel$data$numObs*length(regmodel$BaseModel$manifestVars)*log(2*pi)}
    else if(fitfun == "FIML"){
      lavBIC <- BIC_new # identical if FIML is used!
    }
    } # equivalent to lavaans BIC result
  else if(model_type == "ctsem"){
    BIC_new <- BIC(mxRun(new_model, useOptimizer = F, silent = T))
    lavBIC <- NA
  }
  return_value$mxBIC <- BIC_new
  return_value$lavaan_BIC <-lavBIC

  # in sample RMSEA & ncp
  if(ncp_rmsea == TRUE){
    if(fitfun == "FIML"){
      print("rmsea and ncp not tested for FIML")
      ## step 1: get saturated  model
      satmod <- satmod
      #get minus 2 LL
      m2LL_reg <- regmodel$BaseModel$fitfunction$result[1] # m2LL of regularized model
      m2LL_sat <- satmod$fitfunction$result[1] # m2LL of saturated model
      chisq <- m2LL_reg - m2LL_sat
      ## step 2: get difference in df:

      # get new df of regularized model:
      temp <- mxRun(new_model, useOptimizer = F, silent = T)
      df_reg <- summary(temp$BaseModel)$observedStatistics - summary(temp$BaseModel)$estimatedParameters

      # get df of saturated model:
      df_sat <- summary(satmod)$degreesOfFreedom

      # compute df
      df <- df_reg- df_sat

      if(chisq < df){
        ncp_value <- 0
        rmsea_value <- 0
      } else{
        ncp_value <- (chisq - df)/(nrow(regmodel$BaseModel$data$observed)-1)
        rmsea_value <- sqrt(chisq - df)/sqrt(df * (nrow(regmodel$BaseModel$data$observed))-1)# Note: the formula is not equivalent to results from OpenMx: insted of N-1, N is used in OpenMx
      }
      return_value$ncp <- ncp_value
      return_value$rmsea <- rmsea_value}




    else if( fitfun == "FML"){

  # step 1: get FML equivalent to Jacobucci: .5*formula in Jacobucci_2016 without penalty
  half_FML <- FML[1,1]

  # step 2: get new df:
  temp <- mxRun(new_model, useOptimizer = F, silent = T)
  df <- summary(temp$BaseModel)$observedStatistics - summary(temp$BaseModel)$estimatedParameters

  # step 3: compute rmsea and ncp:

  if(half_FML*2*regmodel$BaseModel$data$numObs < df){
    ncp_value <- 0
    rmsea_value <- 0
  }
  else{
    ncp_value <- (half_FML*2*regmodel$BaseModel$data$numObs - df)/(regmodel$BaseModel$data$numObs-1)
    rmsea_value <- sqrt(half_FML*2*regmodel$BaseModel$data$numObs - df)/sqrt(df * (regmodel$BaseModel$data$numObs)-1)# Note: the formula is not equivalent to results from OpenMx: insted of N-1, N is used in OpenMx
  }
  return_value$ncp <- ncp_value
  return_value$rmsea <- rmsea_value}
}
  else{
    return_value$ncp <- NA
    return_value$rmsea <- NA}

  ## CV m2LL
  if( !is.null(cvsample) ){ # if cvsample (as mxData) is provided:
    CVModel <- new_model # uses model with parameters <.001 set to zero; cross validates this sparse model
    CVModel$data <- cvsample
    fit_CVModel <- mxRun(CVModel, useOptimizer = F, silent = T)
    return_value$CV_m2LL <-fit_CVModel$output$Minus2LogLikelihood
    return_value$CV_AIC <-AIC(fit_CVModel)
    return_value$CV_BIC <-BIC(fit_CVModel)

    if(ncp_rmsea==TRUE ){
      if(fitfun == "FML"){
      CV_FitM <- getFML(fit_CVModel)

      # step 1: get FML equivalent to Jacobucci: .5*formula in Jacobucci_2016 without penalty
      CV_half_FML <- CV_FitM[1,1]

      # step 2: get new df:
      df <- summary(fit_CVModel$BaseModel)$observedStatistics - summary(fit_CVModel$BaseModel)$estimatedParameters

      # step 3: compute rmsea and ncp:
      if(half_FML*2*fit_CVModel$BaseModel$data$numObs < df){
        CV_ncp <- 0
        CV_rmsea <- 0
      }
      else{
        CV_ncp <- (CV_half_FML*2*regmodel$BaseModel$data$numObs - df)/(regmodel$BaseModel$data$numObs-1)
        CV_rmsea <- sqrt(CV_half_FML*2*regmodel$BaseModel$data$numObs - df)/sqrt(df * (regmodel$BaseModel$data$numObs-1))
      }
      return_value$CV_ncp <- CV_ncp
      return_value$CV_rmsea <- CV_rmsea}
      else if(fitfun == "FIML"){
        print("rmsea and ncp not tested for FIML")
        ## step 1: get saturated  model
        cv_satmod <- cv_satmod
        # step 1: get minus 2 LL
        m2LL_cv <- fit_CVModel$fitfunction$result[1,1] # m2LL of CV model
        m2LL_sat <- cv_satmod$fitfunction$result[1] # m2LL of saturated model

        chisq <- m2LL_sat - m2LL_cv
        ## step 2: get difference in df:

        # get new df of regularized model:

        df_cv <- summary(fit_CVModel$BaseModel)$observedStatistics - summary(fit_CVModel$BaseModel)$estimatedParameters
        # get df of saturated model:
        df_sat <- summary(cv_satmod)$degreesOfFreedom
        # compute df
        df <- df_cv- df_sat

        if(chisq < df){
          ncp_value <- 0
          rmsea_value <- 0
        }
        else{
          ncp_value <- (chisq - df)/(nrow(fit_CVModel$data$observed)-1)
          rmsea_value <- sqrt(chisq - df)/sqrt(df * (nrow(fit_CVModel$data$observed)-1))
        }
        return_value$CV_ncp <- CV_ncp
        return_value$CV_rmsea <- CV_rmsea}
      }
    else if(ncp_rmsea!=TRUE){
      return_value$CV_ncp <- NA
      return_value$CV_rmsea <- NA}
  }


  return(return_value)


}
