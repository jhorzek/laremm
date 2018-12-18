#' fitregModelsTryHard
#'
#' Note: laremm is based on the R package \pkg{regsem}. Because of the early status of laremm, it is recommended to use regsem instead!
#' fitregModelsTryHard is identical to fitRegModels, but uses OpenMx TryHard instead of mxRun
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
#' @param satmod saturated model. necessary for computation of ncp and rmsea in FIML models. In many cases, the OpenMx mxRefModels(model, run =TURE) function can be used to build this saturated model. Make sure to only provide the fitted saturated model, not the indipendence model
#' @param cv_satmod saturated model for cross validation. This model has to be based on the cv sample
#' @param CV should a cross validation be computed? If TRUE, provide a Test_Sample
#' @param Test_Sample mxData object with test sample data. Has to be of same data_type as the training data set
#' @examples
#' # The following example is taken from the regsem help to demonstrate the equivalence of both methods:
#'
#' library(lavaan)
#' library(OpenMx)
#' # put variables on same scale for regsem
#' HS <- data.frame(scale(HolzingerSwineford1939[,7:15]))
#'
#' # define variables:
#' latent = c("f1")
#' manifest = c("x1","x2","x3","x4","x5", "x6", "x7", "x8", "x9")
#'
#' # define paths:
#' loadings <- mxPath(from = latent, to = manifest, free = c(F,T,T,T,T,T,T,T,T), values = 1)
#' lcov <- mxPath(from = latent, arrows = 2, free = T, values = 1)
#' lmanif <- mxPath(from = manifest, arrows =2 , free =T, values = 1)
#'
#' # define model:
#' myModel <- mxModel(name = "myModel", latentVars = latent, manifestVars = manifest, type = "RAM",
#'                    mxData(observed = HS, type = "raw"), loadings, lcov, lmanif,
#'                    mxPath(from = "one", to = manifest, free = T)
#'                   )
#'
#' fit_myModel <- mxTryHard(myModel)
#' summary(fit_myModel)
#'
#' # create regularized model:
#'
#' selectedA <- matrix(0, ncol = ncol(fit_myModel$A$values), nrow = nrow(fit_myModel$A$values))
#' selectedA[c(2,3,7,8,9),10] <-1
#'
#'
#' reg_model <- fitRegModels(model = fit_myModel, model_type = "mxModel", fitfun = "FIML",
#'                           pen_on = "A", selectedA = selectedA,
#'                           pen_start = 0, pen_end = .05, pen_stepsize = .01,
#'                           ncp_rmsea = F
#'                           )
#' summary(reg_model$bestmodel)
#'
#' round(reg_model$fit_measures,5)
#'
#'### use laremm in ctsem ####
#'library(ctsem)
#'
#'set.seed(12)
#'
#'## define the population model:
#'
#'# set the drift matrix. Note that drift eta_1_eta2 is set to equal 0 in the population.
#'ct_drift <- matrix(c(-.3,.2,0,-.5), ncol = 2)
#'
#'generatingModel<-ctModel(Tpoints=10,n.latent=2,n.TDpred=0,n.TIpred=0,n.manifest=2,
#'                         MANIFESTVAR=diag(0,2),
#'                         LAMBDA=diag(1,2),
#'                         DRIFT=ct_drift,
#'                         DIFFUSION=matrix(c(.5,0,0,.5),2),
#'                         CINT=matrix(c(0,0),nrow=2),
#'                         T0MEANS=matrix(0,ncol=1,nrow=2),
#'                         T0VAR=diag(1,2))
#'
#'# simulate a training data set
#'traindata <- ctGenerate(generatingModel,n.subjects = 100)
#'
#'## Build the analysis model. Note that drift eta1_eta2 is freely estimated
#'# although it is 0 in the population.
#'
#'myModel <- ctModel(Tpoints=10,n.latent=2,n.TDpred=0,n.TIpred=0,n.manifest=2,
#'                   LAMBDA=diag(1,2),
#'                   MANIFESTVAR=diag(0,2),
#'                   CINT=matrix(c(0,0),nrow=2),
#'                   DIFFUSION=matrix(c('eta1_eta1',0,0,'eta2_eta2'),2),
#'                   T0MEANS=matrix(0,ncol=1,nrow=2),
#'                   T0VAR=diag(1,2))
#'
#'# fit the model using ctsem:
#'fit_myModel <- ctFit(traindata, myModel)
#'summary(fit_myModel)
#'
#'# regularize the model:
#'library(laremm)
#'
#'# start regularization:
#'reg_myModel <- fitRegModels(model = fit_myModel, model_type = "ctsem",
#'                            fitfun = "FIML",data_type = "raw",
#'                            pen_on = "DRIFT", selectedDrifts = "cross",
#'                            pen_start = 0, pen_end = 1, pen_stepsize = .1)
#'
#'# show the best value for penalty term (tuning parameter):
#'reg_myModel$best_penalty
#'
#'# show summary of parameters:
#'summary(reg_myModel$bestmodel)
#'
#'#### additional Cross - validation: #####
#'set.seed(15)
#'#simulate new dataset from same population:
#'testdata <- ctGenerate(generatingModel,n.subjects = 100)
#'# note: ctsem renames the rows and columns when fitting a model. To get
#'# the right names, we fit the ctsem model with the new dataset and then extract
#'# the dataset, where rows and columns have now been renamed
#'
#'fit_myModel <- ctFit(testdata, myModel, useOptimizer = F)
#'testdata <- fit_myModel$mxobj$data
#'
#'# fit models with cross-validation:
#'cv_reg_myModel <- fitRegModels(model = fit_myModel, model_type = "ctsem",
#'                               fitfun = "FIML",data_type = "raw",pen_on = "DRIFT",
#'                               selectedDrifts = "cross", pen_start = 0,
#'                               pen_end = 1, pen_stepsize = .1, CV = TRUE,
#'                               Test_Sample = testdata)
#'# show the best value for penalty term (tuning parameter):
#'cv_reg_myModel$best_penalty
#'# show summary of parameters:
#'summary(cv_reg_myModel$bestmodel)
#'
#'
#' @export
#'

fitRegModelsTryHard <- function (model, model_type = "ctsem", fitfun = "FIML", data_type = "raw",
                          pen_type = "lasso", pen_on = "none", selectedDrifts = "none",
                          driftexpo = TRUE, selectedA = "none", selectedS = "none",
                          pen_start = 0, pen_end = 1, pen_stepsize = 0.01, fit_index = "BIC",
                          ncp_rmsea = FALSE, satmod = NULL, CV = FALSE, cv_satmod = NULL,
                          Test_Sample = NULL)
{
  pen_values = seq(from = pen_start, to = pen_end, by = pen_stepsize)
  if (CV == FALSE) {
    results <- matrix(NA, nrow = 10, ncol = length(pen_values))
    rownames(results) <- c("penalty", "estimated_Parameters",
                           "mx_train_AIC", "lav_train_AIC", "mx_train_BIC",
                           "lav_train_BIC", "RMSEA", "NCP", "negative_Variances",
                           "convergence")
    counter <- 1
    pb <- txtProgressBar(min = pen_start, max = pen_end,
                         style = 3)
    for (i in pen_values) {
      results["penalty", counter] <- i
      reg_Model <- createRegModel(model, model_type = model_type,
                                  fitfun = fitfun, data_type = data_type, pen_type = pen_type,
                                  pen_value = i, pen_on = pen_on, selectedDrifts = selectedDrifts,
                                  driftexpo = driftexpo, selectedA = selectedA,
                                  selectedS = selectedS)
      reg_Model <- mxOption(reg_Model, "Calculate Hessian",
                            "No")
      reg_Model <- mxOption(reg_Model, "Standard Errors",
                            "No")
      fit_reg_Model <- mxTryHard(reg_Model, silent = T)
      results["convergence", counter] <- fit_reg_Model$output$status$code
      variances = diag(nrow(fit_reg_Model$BaseModel$S$values)) ==
        1
      if (any(fit_reg_Model$BaseModel$S$values[variances] <
              0)) {
        results["negative_Variances", counter] <- 1
      }
      else if (model_type == "ctsem" && any(fit_reg_Model$BaseModel$DIFFUSION$result <
                                            0)) {
        results["negative_Variances", counter] <- 1
      }
      else (results["negative_Variances", counter] <- 0)
      FitM <- getFitMeasures(regmodel = fit_reg_Model,
                             fitfun = fitfun, ncp_rmsea = ncp_rmsea, model_type = model_type,
                             cvsample = NULL, satmod = satmod, cv_satmod = cv_satmod)
      results["estimated_Parameters", counter] <- FitM$estimated_params
      results["mx_train_AIC", counter] <- FitM$mxAIC
      results["lav_train_AIC", counter] <- FitM$lavaan_AIC
      results["mx_train_BIC", counter] <- FitM$mxBIC
      results["lav_train_BIC", counter] <- FitM$lavaan_BIC
      results["RMSEA", counter] <- FitM$rmsea
      results["NCP", counter] <- FitM$ncp
      setTxtProgressBar(pb, i)
      counter <- counter + 1
    }
    convergeSubset <- (results["convergence", ] == 0) ==
      (results["negative_Variances", ] == 0)
    convergeSubset <- results[, convergeSubset]
    minimum_mxAIC <- convergeSubset[1, which(convergeSubset["mx_train_AIC",
                                                            ] == min(as.numeric(convergeSubset["mx_train_AIC",
                                                                                               ])))]
    minimum_lavAIC <- convergeSubset[1, which(convergeSubset["lav_train_AIC",
                                                             ] == min(as.numeric(convergeSubset["lav_train_AIC",
                                                                                                ])))]
    minimum_mxBIC <- convergeSubset[1, which(convergeSubset["mx_train_BIC",
                                                            ] == min(convergeSubset["mx_train_BIC", ]))]
    minimum_lavBIC <- convergeSubset[1, which(convergeSubset["lav_train_BIC",
                                                             ] == min(convergeSubset["lav_train_BIC", ]))]
    minimum_Rmsea <- convergeSubset[1, which(convergeSubset["RMSEA",
                                                            ] == min(convergeSubset["RMSEA", ]))]
    minimum_Ncp <- convergeSubset[1, which(convergeSubset["NCP",
                                                          ] == min(convergeSubset["NCP", ]))]
    if (fit_index == "AIC") {
      best_penalty = minimum_mxAIC
      reg_Model_AIC <- createRegModel(model, model_type = model_type,
                                      fitfun = fitfun, data_type = data_type, pen_type = pen_type,
                                      pen_value = minimum_mxAIC, pen_on = pen_on, selectedDrifts = selectedDrifts,
                                      driftexpo = driftexpo, selectedA = selectedA,
                                      selectedS = selectedS)
      reg_Model_AIC <- mxOption(reg_Model_AIC, "Calculate Hessian",
                                "No")
      reg_Model_AIC <- mxOption(reg_Model_AIC, "Standard Errors",
                                "No")
      fit_reg_Model_AIC <- mxTryHard(reg_Model_AIC, silent = T)
      out <- list(best_penalty = minimum_mxAIC, bestmodel = fit_reg_Model_AIC,
                  fit_measures = results)
    }
    if (fit_index == "BIC") {
      best_penalty = minimum_mxBIC
      reg_Model_BIC <- createRegModel(model, model_type = model_type,
                                      fitfun = fitfun, data_type = data_type, pen_type = pen_type,
                                      pen_value = minimum_mxBIC, pen_on = pen_on, selectedDrifts = selectedDrifts,
                                      driftexpo = driftexpo, selectedA = selectedA,
                                      selectedS = selectedS)
      reg_Model_BIC <- mxOption(reg_Model_BIC, "Calculate Hessian",
                                "No")
      reg_Model_BIC <- mxOption(reg_Model_BIC, "Standard Errors",
                                "No")
      fit_reg_Model_BIC <- mxTryHard(reg_Model_BIC, silent = T)
      out <- list(best_penalty = minimum_mxBIC, bestmodel = fit_reg_Model_BIC,
                  fit_measures = results)
    }
    if (fit_index == "RMSEA") {
      best_penalty = minimum_Rmsea
      reg_Model_rmsea <- createRegModel(model, model_type = model_type,
                                        fitfun = fitfun, data_type = data_type, pen_type = pen_type,
                                        pen_value = minimum_Rmsea, pen_on = pen_on, selectedDrifts = selectedDrifts,
                                        driftexpo = driftexpo, selectedA = selectedA,
                                        selectedS = selectedS)
      reg_Model_rmsea <- mxOption(reg_Model_rmsea, "Calculate Hessian",
                                  "No")
      reg_Model_rmsea <- mxOption(reg_Model_rmsea, "Standard Errors",
                                  "No")
      fit_reg_Model_rmsea <- mxTryHard(reg_Model_rmsea, silent = T)
      out <- list(best_penalty = minimum_Rmsea, bestmodel = fit_reg_Model_rmsea,
                  fit_measures = results)
    }
    if (fit_index == "NCP") {
      best_penalty = minimum_Ncp
      reg_Model_ncp <- createRegModel(model, model_type = model_type,
                                      fitfun = fitfun, data_type = data_type, pen_type = pen_type,
                                      pen_value = minimum_Ncp, pen_on = pen_on, selectedDrifts = selectedDrifts,
                                      driftexpo = driftexpo, selectedA = selectedA,
                                      selectedS = selectedS)
      reg_Model_ncp <- mxOption(reg_Model_ncp, "Calculate Hessian",
                                "No")
      reg_Model_ncp <- mxOption(reg_Model_ncp, "Standard Errors",
                                "No")
      fit_reg_Model_ncp <- mxTryHard(reg_Model_ncp, silent = T)
      out <- list(best_penalty = minimum_Rmsea, bestmodel = fit_reg_Model_ncp,
                  fit_measures = results)
    }
    return(out)
  }
  if (CV == TRUE) {
    if (is.null(Test_Sample)) {
      print("Error: Please provide Test Samples as mxData sets.")
    }
    results <- matrix(NA, nrow = 15, ncol = length(pen_values))
    rownames(results) <- c("penalty", "estimated_Parameters",
                           "mx_train_AIC", "lav_train_AIC", "mx_train_BIC",
                           "lav_train_BIC", "RMSEA", "NCP", "CV_m2LL", "CV_AIC",
                           "CV_BIC", "CV_RMSEA", "CV_NCP", "negative_Variances",
                           "convergence")
    counter <- 1
    pb <- txtProgressBar(min = pen_start, max = pen_end,
                         style = 3)
    for (i in pen_values) {
      results["penalty", counter] <- i
      train_reg_Model <- createRegModel(model, model_type = model_type,
                                        fitfun = fitfun, data_type = data_type, pen_type = pen_type,
                                        pen_value = i, pen_on = pen_on, selectedDrifts = selectedDrifts,
                                        driftexpo = driftexpo, selectedA = selectedA,
                                        selectedS = selectedS)
      train_reg_Model <- mxOption(train_reg_Model, "Calculate Hessian",
                                  "No")
      train_reg_Model <- mxOption(train_reg_Model, "Standard Errors",
                                  "No")
      fit_train_reg_Model <- mxTryHard(train_reg_Model, silent = T)
      results["convergence", counter] <- fit_train_reg_Model$output$status$code
      variances = diag(nrow(fit_train_reg_Model$BaseModel$S$values)) ==
        1
      if (any(fit_train_reg_Model$BaseModel$S$values[variances] <
              0)) {
        results["negative_Variances", counter] <- 1
      }
      else (results["negative_Variances", counter] <- 0)
      FitM <- getFitMeasures(regmodel = fit_train_reg_Model,
                             model_type = model_type, fitfun = fitfun, ncp_rmsea = ncp_rmsea,
                             cvsample = Test_Sample, satmod = satmod, cv_satmod = cv_satmod)
      results["estimated_Parameters", counter] <- FitM$estimated_params
      results["mx_train_AIC", counter] <- FitM$mxAIC
      results["lav_train_AIC", counter] <- FitM$lavaan_AIC
      results["mx_train_BIC", counter] <- FitM$mxBIC
      results["lav_train_BIC", counter] <- FitM$lavaan_BIC
      results["RMSEA", counter] <- FitM$rmsea
      results["NCP", counter] <- FitM$ncp
      results["CV_m2LL", counter] <- FitM$CV_m2LL
      results["CV_AIC", counter] <- FitM$CV_AIC
      results["CV_BIC", counter] <- FitM$CV_BIC
      results["CV_RMSEA", counter] <- FitM$CV_rmsea
      results["CV_NCP", counter] <- FitM$CV_ncp
      setTxtProgressBar(pb, i)
      counter <- counter + 1
    }
    convergeSubset <- (results["convergence", ] == 0) ==
      (results["negative_Variances", ] == 0)
    convergeSubset <- results[, convergeSubset]
    minimum_train_mxAIC <- convergeSubset[1, which(convergeSubset["mx_train_AIC",
                                                                  ] == min(as.numeric(convergeSubset["mx_train_AIC",
                                                                                                     ])))]
    minimum_train_lavAIC <- convergeSubset[1, which(convergeSubset["lav_train_AIC",
                                                                   ] == min(as.numeric(convergeSubset["lav_train_AIC",
                                                                                                      ])))]
    minimum_train_mxBIC <- convergeSubset[1, which(convergeSubset["mx_train_BIC",
                                                                  ] == min(convergeSubset["mx_train_BIC", ]))]
    minimum_train_lavBIC <- convergeSubset[1, which(convergeSubset["lav_train_BIC",
                                                                   ] == min(convergeSubset["lav_train_BIC", ]))]
    minimum_Rmsea <- convergeSubset[1, which(convergeSubset["RMSEA",
                                                            ] == min(convergeSubset["RMSEA", ]))]
    minimum_Ncp <- convergeSubset[1, which(convergeSubset["NCP",
                                                          ] == min(convergeSubset["NCP", ]))]
    minimum_CVm2LL <- convergeSubset[1, which(convergeSubset["CV_m2LL",
                                                             ] == min(convergeSubset["CV_m2LL", ]))]
    minimum_CV_AIC <- convergeSubset[1, which(convergeSubset["CV_AIC",
                                                             ] == min(convergeSubset["CV_AIC", ]))]
    minimum_CV_BIC <- convergeSubset[1, which(convergeSubset["CV_BIC",
                                                             ] == min(convergeSubset["CV_BIC", ]))]
    minimum_CV_RMSEA <- convergeSubset[1, which(convergeSubset["CV_RMSEA",
                                                               ] == min(convergeSubset["CV_RMSEA", ]))]
    minimum_CV_NCP <- convergeSubset[1, which(convergeSubset["CV_NCP",
                                                             ] == min(convergeSubset["CV_NCP", ]))]
    if (fit_index == "CV_m2LL") {
      best_penalty = minimum_CVm2LL
      reg_Model_CV_m2LL <- createRegModel(model, model_type = model_type,
                                          fitfun = fitfun, data_type = data_type, pen_type = pen_type,
                                          pen_value = minimum_CVm2LL, pen_on = pen_on,
                                          selectedDrifts = selectedDrifts, driftexpo = driftexpo,
                                          selectedA = selectedA, selectedS = selectedS)
      reg_Model_CV_m2LL <- mxOption(reg_Model_CV_m2LL,
                                    "Calculate Hessian", "No")
      reg_Model_CV_m2LL <- mxOption(reg_Model_CV_m2LL,
                                    "Standard Errors", "No")
      fit_reg_Model_CV_m2LL <- mxTryHard(reg_Model_CV_m2LL,
                                     silent = T)
      out <- list(best_penalty = minimum_CVm2LL, bestmodel = fit_reg_Model_CV_m2LL,
                  fit_measures = results)
    }
    if (fit_index == "AIC") {
      best_penalty = minimum_CV_AIC
      reg_Model_AIC <- createRegModel(model, model_type = model_type,
                                      fitfun = fitfun, data_type = data_type, pen_type = pen_type,
                                      pen_value = minimum_CV_AIC, pen_on = pen_on,
                                      selectedDrifts = selectedDrifts, driftexpo = driftexpo,
                                      selectedA = selectedA, selectedS = selectedS)
      reg_Model_AIC <- mxOption(reg_Model_AIC, "Calculate Hessian",
                                "No")
      reg_Model_AIC <- mxOption(reg_Model_AIC, "Standard Errors",
                                "No")
      fit_reg_Model_AIC <- mxTryHard(reg_Model_AIC, silent = T)
      out <- list(best_penalty = minimum_CV_AIC, bestmodel = fit_reg_Model_AIC,
                  fit_measures = results)
    }
    if (fit_index == "BIC") {
      best_penalty = minimum_CV_BIC
      reg_Model_BIC <- createRegModel(model, model_type = model_type,
                                      fitfun = fitfun, data_type = data_type, pen_type = pen_type,
                                      pen_value = minimum_CV_BIC, pen_on = pen_on,
                                      selectedDrifts = selectedDrifts, driftexpo = driftexpo,
                                      selectedA = selectedA, selectedS = selectedS)
      reg_Model_BIC <- mxOption(reg_Model_BIC, "Calculate Hessian",
                                "No")
      reg_Model_BIC <- mxOption(reg_Model_BIC, "Standard Errors",
                                "No")
      fit_reg_Model_BIC <- mxTryHard(reg_Model_BIC, silent = T)
      out <- list(best_penalty = minimum_CV_BIC, bestmodel = fit_reg_Model_BIC,
                  fit_measures = results)
    }
    if (fit_index == "RMSEA") {
      best_penalty = minimum_CV_RMSEA
      reg_Model_rmsea <- createRegModel(model, model_type = model_type,
                                        fitfun = fitfun, data_type = data_type, pen_type = pen_type,
                                        pen_value = minimum_CV_RMSEA, pen_on = pen_on,
                                        selectedDrifts = selectedDrifts, driftexpo = driftexpo,
                                        selectedA = selectedA, selectedS = selectedS)
      reg_Model_rmsea <- mxOption(reg_Model_rmsea, "Calculate Hessian",
                                  "No")
      reg_Model_rmsea <- mxOption(reg_Model_rmsea, "Standard Errors",
                                  "No")
      fit_reg_Model_rmsea <- mxTryHard(reg_Model_rmsea, silent = T)
      out <- list(best_penalty = minimum_CV_RMSEA, bestmodel = fit_reg_Model_rmsea,
                  fit_measures = results)
    }
    if (fit_index == "NCP") {
      best_penalty = minimum_CV_NCP
      reg_Model_ncp <- createRegModel(model, model_type = model_type,
                                      fitfun = fitfun, data_type = data_type, pen_type = pen_type,
                                      pen_value = minimum_CV_NCP, pen_on = pen_on,
                                      selectedDrifts = selectedDrifts, driftexpo = driftexpo,
                                      selectedA = selectedA, selectedS = selectedS)
      reg_Model_ncp <- mxOption(reg_Model_ncp, "Calculate Hessian",
                                "No")
      reg_Model_ncp <- mxOption(reg_Model_ncp, "Standard Errors",
                                "No")
      fit_reg_Model_ncp <- mxTryHard(reg_Model_ncp, silent = T)
      out <- list(best_penalty = minimum_CV_NCP, bestmodel = fit_reg_Model_ncp,
                  fit_measures = results)
    }
    return(out)
  }
}
