#' fitregModels
#'
#' Note: laremm is based on the R package \pkg{regsem}. Because of the early status of laremm, it is recommended to use regsem instead!
#' fitregModels creates a regularized model from a mxModel or ctsem. It then runs this model with multiple penalty values.
#'
#' @param model mxModel or ctsem object
#' @param model_type specify the type of model provided: ctsem or mxModel
#' @param fitfun fitfunction to be used in the fitting procedure. Either FML or FIML
#' @param data_type type of data in the model. Either "cov" or "raw"
#' @param penalty_type so far only "lasso" implemented
#' @param pen_start lowest penalty value to evaluate. Recommended: 0
#' @param pen_end highest penalty value to evaluate
#' @param pen_stepsize increse of penalty with each iteration. e.g. if pen_start = 0, pen_end = 1, pen_stepsize = .1, fitRegModels will iterate over pen = 0, pen = .1, pen = .2, ...
#' @param pen_on string vector with matrices that should be regularized. Possible are combinations of "A", "S", "DRIFT"
#' @param selectedDrifts drift values to regularize. Possible are "all", "cross", "auto" or providing a matrix of the same size as the drift matrix with ones for every parameter to regularize and 0 for every non-regularized parameter
#' @param driftexpo specifiy if the regularization will be performed on the raw drift matrix or on the exponential of the drift matrix (discrete time parameters)
#' @param DRIFT_dt provide the discrete time points for which the drift will be regularized. A vector with multiple values is possible
#' @param selectedA A values to regularize. Possible are "all", or providing a matrix of the same size as the A matrix with ones for every parameter to regularize and 0 for every non-regularized parameter
#' @param selectedS S values to regularize. Possible are "all", or providing a matrix of the same size as the S matrix with ones for every parameter to regularize and 0 for every non-regularized parameter
#' @param fit_index which fit index should be used to find the best model? Possible are AIC and BIC. RMSEA and NCP for covariance based models
#' @param ncp_rmsea should rmsea and ncp be computed? Only possible for covariance based models
#' @param CV should a cross validation be computed? If TRUE, provide a Test_Sample
#' @param Test_Sample mxData object with test sample data. Has to be of same data_type as the training data set
#'
#' @export
#'
fitRegModels <- function(model, model_type = "ctsem", fitfun = "FIML", data_type= "raw",pen_type = "lasso",
                         pen_on = "none",selectedDrifts ="none", driftexpo = TRUE, selectedA = "none", selectedS = "none",
                         pen_start = 0, pen_end = 1, pen_stepsize = .01,
                         fit_index = "BIC",
                         ncp_rmsea =FALSE,
                         CV = FALSE,
                         Test_Sample = NULL){
  # in future:
  #, selectedLambda = "none",
  #selectedCint = "none",
  #selectedTdpredeffect = "none",
  #selectedDiffusion = "none",

  pen_values = seq(from = pen_start, to = pen_end, by = pen_stepsize) # iteration through these values


  if(CV == FALSE){
  # matrix for results:
  results <- matrix(NA, nrow = 8, ncol = length(pen_values))
  rownames(results) <- c("penalty", "estimated_Parameters", "mx_train_AIC","lav_train_AIC",
                         "mx_train_BIC", "lav_train_BIC","RMSEA", "NCP")
  counter <- 1

  # computation:
  pb <- txtProgressBar(min = pen_start, max = pen_end, style = 3) # progress-bar
  for(i in pen_values){
    results["penalty",counter] <- i
    reg_Model <- createRegModel_new(model, model_type = model_type, fitfun = fitfun, data_type= data_type, pen_type = pen_type, pen_value = i,
                                pen_on = pen_on, selectedDrifts =selectedDrifts, driftexpo = driftexpo, selectedA = selectedA, selectedS = selectedS)

    reg_Model <- mxOption(reg_Model, "Calculate Hessian", "No") # might cause errors; check
    reg_Model <- mxOption(reg_Model, "Standard Errors", "No") # might cause errors; check
    fit_reg_Model <- mxRun(reg_Model, silent = T) # run Model; starting values can be very critical as the model tends to get stuck in local minima either close to the model parameters without penalty or all parameters set to 0

    ### compute AIC and BIC:
    FitM <- getFitMeasures(regmodel = fit_reg_Model, fitfun = fitfun, ncp_rmsea = ncp_rmsea,  model_type = model_type, cvsample = NULL)

    results["estimated_Parameters",counter] <- FitM$estimated_params # estimated parameters
    results["mx_train_AIC",counter] <- FitM$mxAIC # mxAIC
    results["lav_train_AIC",counter] <- FitM$lavaan_AIC # lavaan equivalent AIC
    results["mx_train_BIC",counter] <- FitM$mxBIC # mxBIC
    results["lav_train_BIC",counter] <- FitM$lavaan_BIC # lavBIC
    results["RMSEA",counter] <- FitM$rmsea # regsem RMSEA
    results["NCP",counter] <- FitM$ncp # regsem ncp

    setTxtProgressBar(pb, i)
    counter <- counter+1
  }

  # Find Minima / best penalty value

  minimum_mxAIC <- results[1,which(results["mx_train_AIC",]==min(as.numeric(results["mx_train_AIC",])))]
  minimum_lavAIC <- results[1,which(results["lav_train_AIC",]==min(as.numeric(results["lav_train_AIC",])))] # ?quivalent zu minimalem AIC von regsem


  minimum_mxBIC <- results[1,which(results["mx_train_BIC",]==min(results["mx_train_BIC",]))]
  minimum_lavBIC <- results[1,which(results["lav_train_BIC",]==min(results["lav_train_BIC",]))] # ?quivalent zu minimalem BIC von regsem!!

  minimum_Rmsea <- results[1,which(results["RMSEA",]==min(results["RMSEA",]))] # ?quivalent zu minimalem BIC von regsem!!
  minimum_Ncp <- results[1,which(results["NCP",]==min(results["NCP",]))] # ?quivalent zu minimalem BIC von regsem!!


  # getting parameters:
  if(fit_index == "AIC"){
    best_penalty = minimum_mxAIC
  reg_Model_AIC <- createRegModel_new(model, model_type = model_type, fitfun = fitfun, data_type= data_type, pen_type = pen_type, pen_value = minimum_mxAIC,
                                  pen_on = pen_on, selectedDrifts =selectedDrifts, driftexpo = driftexpo, selectedA = selectedA,
                                  selectedS = selectedS)
  fit_reg_Model_BIC <- mxRun(reg_Model_BIC, silent = T)
  fit_reg_Model_AIC <- mxRun(reg_Model_AIC, silent = T)
  out <- list("best_penalty" = minimum_mxAIC, "bestmodel" = fit_reg_Model_AIC, "fit_measures" = results)
  }

  if(fit_index == "BIC"){
    best_penalty = minimum_mxBIC
  reg_Model_BIC <- createRegModel_new(model, model_type = model_type, fitfun = fitfun, data_type= data_type, pen_type = pen_type, pen_value = minimum_mxBIC,
                                  pen_on = pen_on, selectedDrifts =selectedDrifts, driftexpo = driftexpo, selectedA = selectedA,
                                  selectedS = selectedS)
  fit_reg_Model_BIC <- mxRun(reg_Model_BIC, silent = T)
  out <- list("best_penalty" = minimum_mxBIC, "bestmodel" = fit_reg_Model_BIC, "fit_measures" = results)
  }

  if(fit_index == "RMSEA"){
    best_penalty = minimum_Rmsea
  reg_Model_rmsea <- createRegModel_new(model, model_type = model_type, fitfun = fitfun, data_type= data_type, pen_type = pen_type, pen_value = minimum_Rmsea,
                                    pen_on = pen_on, selectedDrifts =selectedDrifts, driftexpo = driftexpo, selectedA = selectedA,
                                    selectedS = selectedS)
  fit_reg_Model_rmsea <- mxRun(reg_Model_rmsea, silent = T)
  out <- list("best_penalty" = minimum_Rmsea, "bestmodel" = fit_reg_Model_rmsea, "fit_measures" = results)
  }

  if(fit_index == "NCP"){
    best_penalty = minimum_Ncp
    reg_Model_ncp <- createRegModel_new(model, model_type = model_type, fitfun = fitfun, data_type= data_type, pen_type = pen_type, pen_value = minimum_Ncp,
                                          pen_on = pen_on, selectedDrifts =selectedDrifts, driftexpo = driftexpo, selectedA = selectedA,
                                          selectedS = selectedS)
    fit_reg_Model_ncp <- mxRun(reg_Model_ncp, silent = T)
    out <- list("best_penalty" = minimum_Rmsea, "bestmodel" = fit_reg_Model_ncp, "fit_measures" = results)
  }


  return(out)

  }

    if(CV == TRUE){
      if(is.null(Test_Sample)){
        print("Error: Please provide Test Samples as mxData sets.")
        break
      }
      # assumes that the sample given in model is the train data and tests it against the provided test data

        # matrix for results:
        results <- matrix(NA, nrow = 13, ncol = length(pen_values))
        rownames(results) <- c("penalty", "estimated_Parameters", "mx_train_AIC","lav_train_AIC",
                               "mx_train_BIC", "lav_train_BIC","RMSEA","NCP", "CV_m2LL", "CV_AIC","CV_BIC", "CV_RMSEA","CV_NCP")
        counter <- 1
        pb <- txtProgressBar(min = pen_start, max = pen_end, style = 3) # progress-bar
      for(i in pen_values){

        results["penalty",counter] <- i # save penalty value


        train_reg_Model <- createRegModel_new(model, model_type = model_type, fitfun = fitfun, data_type= data_type, pen_type = pen_type, pen_value = i,
                                    pen_on = pen_on, selectedDrifts =selectedDrifts, driftexpo = driftexpo, selectedA = selectedA,
                                    selectedS = selectedS)

        train_reg_Model <- mxOption(train_reg_Model, "Calculate Hessian", "No") # might cause errors; check
        train_reg_Model <- mxOption(train_reg_Model, "Standard Errors", "No") # might cause errors; check
        fit_train_reg_Model <- mxRun(train_reg_Model, silent = T) # run Model; starting values can be very critical as the model tends to get stuck in local minima either close to the model parameters without penalty or all parameters set to 0

        ### compute AIC and BIC:
        FitM <- getFitMeasures(regmodel = fit_train_reg_Model, model_type = model_type, fitfun = fitfun, ncp_rmsea = ncp_rmsea, cvsample = Test_Sample)

        results["estimated_Parameters",counter] <- FitM$estimated_params # estimated parameters
        results["mx_train_AIC",counter] <- FitM$mxAIC # mxAIC
        results["lav_train_AIC",counter] <- FitM$lavaan_AIC # lavaan equivalent AIC
        results["mx_train_BIC",counter] <- FitM$mxBIC # mxBIC
        results["lav_train_BIC",counter] <- FitM$lavaan_BIC # lavBIC
        results["RMSEA",counter] <- FitM$rmsea # regsem RMSEA
        results["NCP",counter] <- FitM$ncp # regsem RMSEA

        ## CV
        results["CV_m2LL",counter] <- FitM$CV_m2LL
        results["CV_AIC",counter] <- FitM$CV_AIC
        results["CV_BIC",counter] <- FitM$CV_BIC
        results["CV_RMSEA",counter] <- FitM$CV_rmsea
        results["CV_NCP",counter] <- FitM$CV_ncp

        setTxtProgressBar(pb, i)
        counter <- counter+1
      }

        ## search best penalty values:
        # train set
        minimum_train_mxAIC <- results[1,which(results["mx_train_AIC",]==min(as.numeric(results["mx_train_AIC",])))]
        minimum_train_lavAIC <- results[1,which(results["lav_train_AIC",]==min(as.numeric(results["lav_train_AIC",])))] # ?quivalent zu minimalem AIC von regsem?


        minimum_train_mxBIC <- results[1,which(results["mx_train_BIC",]==min(results["mx_train_BIC",]))]
        minimum_train_lavBIC <- results[1,which(results["lav_train_BIC",]==min(results["lav_train_BIC",]))] # ?quivalent zu minimalem BIC von regsem!!

        minimum_Rmsea <- results[1,which(results["RMSEA",]==min(results["RMSEA",]))] # ?quivalent zu minimalem BIC von regsem!!
        minimum_Ncp <- results[1,which(results["NCP",]==min(results["NCP",]))] # ?quivalent zu minimalem BIC von regsem!!

        # CV
        minimum_CVm2LL <- results[1,which(results["CV_m2LL",]==min(results["CV_m2LL",]))]
        minimum_CV_AIC <- results[1,which(results["CV_AIC",]==min(results["CV_AIC",]))]
        minimum_CV_BIC <- results[1,which(results["CV_BIC",]==min(results["CV_BIC",]))]
        minimum_CV_RMSEA <- results[1,which(results["CV_RMSEA",]==min(results["CV_RMSEA",]))]
        minimum_CV_NCP <- results[1,which(results["CV_NCP",]==min(results["CV_NCP",]))]


        # getting parameters:
        if(fit_index == "CV_m2LL"){
          best_penalty = minimum_CVm2LL
          reg_Model_CV_m2LL <- createRegModel_new(model, model_type = model_type, fitfun = fitfun, data_type= data_type, pen_type = pen_type, pen_value = minimum_CVm2LL,
                                                  pen_on = pen_on, selectedDrifts =selectedDrifts, driftexpo = driftexpo, selectedA = selectedA,
                                                  selectedS = selectedS)
          fit_reg_Model_CV_m2LL <- mxRun(reg_Model_CV_m2LL, silent = T)
          out <- list("best_penalty" = minimum_CVm2LL, "bestmodel" = fit_reg_Model_CV_m2LL, "fit_measures" = results)
        }

        if(fit_index == "AIC"){
          best_penalty = minimum_CV_AIC
          reg_Model_AIC <- createRegModel_new(model, model_type = model_type, fitfun = fitfun, data_type= data_type, pen_type = pen_type, pen_value = minimum_CV_AIC,
                                              pen_on = pen_on, selectedDrifts =selectedDrifts, driftexpo = driftexpo, selectedA = selectedA,
                                              selectedS = selectedS)
          fit_reg_Model_AIC <- mxRun(reg_Model_AIC, silent = T)
          out <- list("best_penalty" = minimum_CV_AIC, "bestmodel" = fit_reg_Model_AIC, "fit_measures" = results)
        }

        if(fit_index == "BIC"){
          best_penalty = minimum_CV_BIC
          reg_Model_BIC <- createRegModel_new(model, model_type = model_type, fitfun = fitfun, data_type= data_type, pen_type = pen_type, pen_value = minimum_CV_BIC,
                                              pen_on = pen_on, selectedDrifts =selectedDrifts, driftexpo = driftexpo, selectedA = selectedA,
                                              selectedS = selectedS)
          fit_reg_Model_BIC <- mxRun(reg_Model_BIC, silent = T)
          out <- list("best_penalty" = minimum_CV_BIC, "bestmodel" = fit_reg_Model_BIC, "fit_measures" = results)
        }

        if(fit_index == "RMSEA"){
          best_penalty = minimum_CV_RMSEA
          reg_Model_rmsea <- createRegModel_new(model, model_type = model_type, fitfun = fitfun, data_type= data_type, pen_type = pen_type, pen_value = minimum_CV_RMSEA,
                                                pen_on = pen_on, selectedDrifts =selectedDrifts, driftexpo = driftexpo, selectedA = selectedA,
                                                selectedS = selectedS)
          fit_reg_Model_rmsea <- mxRun(reg_Model_rmsea, silent = T)
          out <- list("best_penalty" = minimum_CV_RMSEA, "bestmodel" = fit_reg_Model_rmsea, "fit_measures" = results)
        }
        if(fit_index == "NCP"){
          best_penalty = minimum_CV_NCP
          reg_Model_ncp <- createRegModel_new(model, model_type = model_type, fitfun = fitfun, data_type= data_type, pen_type = pen_type, pen_value = minimum_CV_NCP,
                                                pen_on = pen_on, selectedDrifts =selectedDrifts, driftexpo = driftexpo, selectedA = selectedA,
                                                selectedS = selectedS)
          fit_reg_Model_ncp <- mxRun(reg_Model_ncp, silent = T)
          out <- list("best_penalty" = minimum_CV_NCP, "bestmodel" = fit_reg_Model_ncp, "fit_measures" = results)
        }



        return(out)



    }

}
