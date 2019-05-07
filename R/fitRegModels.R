#' fitRegModels
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
#' @param fit_index which fit index should be used to find the best model? Possible are AIC and BIC, CV_m2LL, CV_AIC, CV_BIC
#' @param CV should a cross validation be computed? If TRUE, provide a Test_Sample
#' @param Test_Sample mxData object with test sample data. Has to be of same data_type as the training data set
#' @param zeroThresh threshold for evaluating regularized parameters as zero. Default is .001 similar to \pkg{regsem}
#' @param setZero should parameters below zeroThresh be set to zero in all fit calculations. Default is FALSE, similar to \pkg{regsem}
#'
#'
#' @author Jannik Orzek
#' @import OpenMx ctsem
#'
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
#' )
#'
#' fit_myModel <- mxRun(myModel)
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
#'                           pen_start = 0, pen_end = .05, pen_stepsize = .01
#'                           )
#' summary(reg_model)
#' reg_model$`fit measures`
#'
#' ### use laremm in ctsem ####
#' library(ctsem)
#'
#' set.seed(12)
#'
#' ## define the population model:
#'
#' # set the drift matrix. Note that drift eta_1_eta2 is set to equal 0 in the population.
#' ct_drift <- matrix(c(-.3,.2,0,-.5), ncol = 2)
#'
#' generatingModel<-ctModel(Tpoints=10,n.latent=2,n.TDpred=0,n.TIpred=0,n.manifest=2,
#'                          MANIFESTVAR=diag(0,2),
#'                          LAMBDA=diag(1,2),
#'                          DRIFT=ct_drift,
#'                          DIFFUSION=matrix(c(.5,0,0,.5),2),
#'                          CINT=matrix(c(0,0),nrow=2),
#'                          T0MEANS=matrix(0,ncol=1,nrow=2),
#'                          T0VAR=diag(1,2))
#'
#' # simulate a training data set
#' traindata <- ctGenerate(generatingModel,n.subjects = 100)
#'
#' ## Build the analysis model. Note that drift eta1_eta2 is freely estimated
#' # although it is 0 in the population.
#'
#' myModel <- ctModel(Tpoints=10,n.latent=2,n.TDpred=0,n.TIpred=0,n.manifest=2,
#'                    LAMBDA=diag(1,2),
#'                    MANIFESTVAR=diag(0,2),
#'                    CINT=matrix(c(0,0),nrow=2),
#'                    DIFFUSION=matrix(c('eta1_eta1',0,0,'eta2_eta2'),2),
#'                    T0MEANS=matrix(0,ncol=1,nrow=2),
#'                    T0VAR=diag(1,2))
#'
#' # fit the model using ctsem:
#' fit_myModel <- ctFit(traindata, myModel)
#' fit_myModel$mxobj$DRIFT$values
#'
#' # regularize the model:
#' library(laremm)
#'
#' # start regularization:
#' reg_myModel <- fitRegModels(model = fit_myModel, model_type = "ctsem",
#'                             fitfun = "FIML",data_type = "raw",
#'                             pen_on = "DRIFT", selectedDrifts = "cross",
#'                             pen_start = 0, pen_end = 1, pen_stepsize = .1)
#'
#' # show the best value for penalty term (tuning parameter):
#' reg_myModel$`best penalty`
#'
#' # show summary of parameters:
#' summary(reg_myModel)
#'
#' #### additional Cross - validation: #####
#' set.seed(15)
#' #simulate new dataset from same population:
#' testdata <- ctGenerate(generatingModel,n.subjects = 100)
#' # note: ctsem renames the rows and columns when fitting a model. To get
#' # the right names, we fit the ctsem model with the new dataset and then extract
#' # the dataset, where rows and columns have now been renamed
#'
#' fit_myModel <- ctFit(testdata, myModel, useOptimizer = F)
#' testdata <- fit_myModel$mxobj$data
#'
#' # fit models with cross-validation:
#' cv_reg_myModel <- fitRegModels(model = fit_myModel, model_type = "ctsem",
#'                                fitfun = "FIML",data_type = "raw",pen_on = "DRIFT",
#'                                selectedDrifts = "cross", pen_start = 0,
#'                                pen_end = 1, pen_stepsize = .1, CV = TRUE,
#'                                Test_Sample = testdata, fit_index = "CV_BIC")
#' # show the summary:
#' summary(cv_reg_myModel)
#'
#' @export
#'
fitRegModels <- function(model, model_type = "ctsem", fitfun = "FIML", data_type= "raw",pen_type = "lasso",
                         pen_on = "none",selectedDrifts ="none", selectedA = "none", selectedS = "none",
                         pen_start = 0, pen_end = 1, pen_stepsize = .01,
                         fit_index = "BIC",
                         CV = FALSE,
                         Test_Sample = NULL,
                         zeroThresh = .001,
                         setZero = FALSE,
                         driftexpo = TRUE,
                         DRIFT_dt =1){
  # in future:
  #, selectedLambda = "none",
  #selectedCint = "none",
  #selectedTdpredeffect = "none",
  #selectedDiffusion = "none",

  call <- mget(names(formals()),sys.frame(sys.nframe()))

  pen_values = seq(from = pen_start, to = pen_end, by = pen_stepsize) # iteration through these values


  if(CV == FALSE){
    # matrix for results:
    results <- matrix(NA, nrow = 7, ncol = length(pen_values))
    rownames(results) <- c("penalty", "estimated Parameters", "m2LL","AIC",
                           "BIC", "negative variances","convergence")
    counter <- 1

    # computation:
    pb <- txtProgressBar(min = pen_start, max = pen_end, style = 3) # progress-bar
    for(i in pen_values){
      results["penalty",counter] <- i
      reg_Model <- createRegModel(model, model_type = model_type, fitfun = fitfun, data_type= data_type, pen_type = pen_type, pen_value = i,
                                  pen_on = pen_on, selectedDrifts =selectedDrifts, driftexpo = driftexpo, DRIFT_dt = DRIFT_dt, selectedA = selectedA, selectedS = selectedS)

      reg_Model <- mxOption(reg_Model, "Calculate Hessian", "No") # might cause errors; check
      reg_Model <- mxOption(reg_Model, "Standard Errors", "No") # might cause errors; check
      fit_reg_Model <- mxRun(reg_Model, silent = T) # run Model; starting values can be very critical as the model tends to get stuck in local minima either close to the model parameters without penalty or all parameters set to 0

      results["convergence",counter] <- fit_reg_Model$output$status$code# check convergence

      variances = diag(nrow(fit_reg_Model$BaseModel$S$values))==1

      if(any(fit_reg_Model$BaseModel$S$values[variances] <0)){
        results["negative variances",counter] <- 1 # check negative variances
      #}else if(model_type=="ctsem" && any(fit_reg_Model$BaseModel$DIFFUSION$result <0)){
      #  results["negative variances",counter] <- 1 # check negative diffusion in ctsem
      }else(
        results["negative variances",counter] <- 0
      )

      ### compute AIC and BIC:
      Fit <- getFitMeasures(regmodel = fit_reg_Model, fitfun = fitfun, model_type = model_type, cvsample = NULL, zeroThresh = zeroThresh, setZero = setZero)
      FitM <- Fit$return_value
      FitModel <- Fit$return_Model

      results["estimated Parameters",counter] <- FitM$estimated_params # estimated parameters
      results["m2LL",counter] <- FitM$m2LL # -2LogL
      results["AIC",counter] <- FitM$AIC # AIC
      results["BIC",counter] <- FitM$BIC # BIC

      setTxtProgressBar(pb, i)
      counter <- counter+1
    }

    # Find Minima / best penalty value
    convergeSubset <- (results["convergence",] == 0) == (results["negative variances",]==0)
    convergeSubset <- results[,convergeSubset]

    minimum_m2LL <- convergeSubset[1,which(convergeSubset["m2LL",]==min(as.numeric(convergeSubset["m2LL",])))]

    minimum_AIC <- convergeSubset[1,which(convergeSubset["AIC",]==min(as.numeric(convergeSubset["AIC",])))]

    minimum_BIC <- convergeSubset[1,which(convergeSubset["BIC",]==min(convergeSubset["BIC",]))]


      # getting parameters:
    if(fit_index == "m2LL"){
      best_penalty = minimum_m2LL
      reg_Model_m2LL <- createRegModel(model, model_type = model_type, fitfun = fitfun, data_type= data_type, pen_type = pen_type, pen_value = best_penalty,
                                      pen_on = pen_on, selectedDrifts =selectedDrifts,DRIFT_dt = DRIFT_dt, driftexpo = driftexpo, selectedA = selectedA,
                                      selectedS = selectedS)
      reg_Model_m2LL <- mxOption(reg_Model_m2LL, "Calculate Hessian", "No") # might cause errors; check
      reg_Model_m2LL <- mxOption(reg_Model_m2LL, "Standard Errors", "No") # might cause errors; check

      fit_reg_Model_m2LL <- mxRun(reg_Model_m2LL, silent = T)
      #set parameters below zeroThresh to zero
      fit_reg_Model_m2LL <- getFitMeasures(regmodel = reg_Model_m2LL, fitfun = fitfun, model_type = model_type, cvsample = NULL, zeroThresh = zeroThresh, setZero = setZero)
      fit_reg_Model_m2LL <- fit_reg_Model_m2LL$return_Model

      out <- list("best penalty" = minimum_m2LL, "bestmodel" = fit_reg_Model_m2LL, "fit measures" = t(results), "call" = call)
    }

    if(fit_index == "AIC"){
      best_penalty = minimum_AIC
    reg_Model_AIC <- createRegModel(model, model_type = model_type, fitfun = fitfun, data_type= data_type, pen_type = pen_type, pen_value = best_penalty,
                                    pen_on = pen_on, selectedDrifts =selectedDrifts, driftexpo = driftexpo, DRIFT_dt = DRIFT_dt,selectedA = selectedA,
                                    selectedS = selectedS)
    reg_Model_AIC <- mxOption(reg_Model_AIC, "Calculate Hessian", "No") # might cause errors; check
    reg_Model_AIC <- mxOption(reg_Model_AIC, "Standard Errors", "No") # might cause errors; check

    fit_reg_Model_AIC <- mxRun(reg_Model_AIC, silent = T)
    #set parameters below zeroThresh to zero
    fit_reg_Model_AIC <- getFitMeasures(regmodel = fit_reg_Model_AIC, fitfun = fitfun, model_type = model_type, cvsample = NULL, zeroThresh = zeroThresh, setZero = setZero)
    fit_reg_Model_AIC <- fit_reg_Model_AIC$return_Model

    out <- list("best penalty" = minimum_AIC, "bestmodel" = fit_reg_Model_AIC, "fit measures" = t(results), "call" = call)
    }

    if(fit_index == "BIC"){
      best_penalty = minimum_BIC
    reg_Model_BIC <- createRegModel(model, model_type = model_type, fitfun = fitfun, data_type= data_type, pen_type = pen_type, pen_value = best_penalty,
                                    pen_on = pen_on, selectedDrifts =selectedDrifts, driftexpo = driftexpo,DRIFT_dt = DRIFT_dt, selectedA = selectedA,
                                    selectedS = selectedS)
    reg_Model_BIC <- mxOption(reg_Model_BIC, "Calculate Hessian", "No") # might cause errors; check
    reg_Model_BIC <- mxOption(reg_Model_BIC, "Standard Errors", "No") # might cause errors; check

    fit_reg_Model_BIC <- mxRun(reg_Model_BIC, silent = T)
    #set parameters below zeroThresh to zero
    fit_reg_Model_BIC <- getFitMeasures(regmodel = fit_reg_Model_BIC, fitfun = fitfun, model_type = model_type, cvsample = NULL, zeroThresh = zeroThresh, setZero = setZero)
    fit_reg_Model_BIC <- fit_reg_Model_BIC$return_Model


    out <- list("best penalty" = minimum_BIC, "bestmodel" = fit_reg_Model_BIC, "fit measures" = t(results), "call" = call)
    }
    class(out) <- "FitLaremmObject"

    return(out)

  }

    if(CV == TRUE){
      if(is.null(Test_Sample)){
        print("Error: Please provide Test Samples as mxData sets.")

      }
      # assumes that the sample given in model is the train data and tests it against the provided test data

        # matrix for results:
        results <- matrix(NA, nrow = 10, ncol = length(pen_values))
        rownames(results) <- c("penalty", "estimated Parameters", "m2LL","AIC",
                               "BIC", "CV m2LL","CV AIC","CV BIC", "negative variances","convergence")
        counter <- 1
        pb <- txtProgressBar(min = pen_start, max = pen_end, style = 3) # progress-bar
      for(i in pen_values){

        results["penalty",counter] <- i # save penalty value


        train_reg_Model <- createRegModel(model, model_type = model_type, fitfun = fitfun, data_type= data_type, pen_type = pen_type, pen_value = i,
                                    pen_on = pen_on, selectedDrifts =selectedDrifts, driftexpo = driftexpo,DRIFT_dt = DRIFT_dt, selectedA = selectedA,
                                    selectedS = selectedS)

        train_reg_Model <- mxOption(train_reg_Model, "Calculate Hessian", "No") # might cause errors; check
        train_reg_Model <- mxOption(train_reg_Model, "Standard Errors", "No") # might cause errors; check
        fit_train_reg_Model <- mxRun(train_reg_Model, silent = T) # run Model; starting values can be very critical as the model tends to get stuck in local minima either close to the model parameters without penalty or all parameters set to 0

        results["convergence",counter] <- fit_train_reg_Model$output$status$code# check convergence

        variances = diag(nrow(fit_train_reg_Model$BaseModel$S$values))==1

        if(any(fit_train_reg_Model$BaseModel$S$values[variances] <0)){
          results["negative variances",counter] <- 1 # check negative variances
        }else(
          results["negative variances",counter] <- 0
        )

        ### compute AIC and BIC:
        Fit <- getFitMeasures(regmodel = fit_train_reg_Model, model_type = model_type, fitfun = fitfun, cvsample = Test_Sample, zeroThresh = zeroThresh, setZero = setZero)
        FitM <- Fit$return_value
        FitModel <- Fit$return_Model

        results["estimated Parameters",counter] <- FitM$estimated_params # estimated parameters
        results["m2LL",counter] <- FitM$m2LL # -2LogL
        results["AIC",counter] <- FitM$AIC # AIC
        results["BIC",counter] <- FitM$BIC # BIC
        results["CV m2LL",counter] <- FitM$CV.m2LL # Cross Validated -2LogL
        results["CV AIC",counter] <- FitM$CV.AIC # Cross Validated AIC
        results["CV BIC",counter] <- FitM$CV.BIC # Cross Validated BIC

        setTxtProgressBar(pb, i)
        counter <- counter+1
      }

        ## search best penalty values:
        convergeSubset <- (results["convergence",] == 0) == (results["negative variances",]==0)
        convergeSubset <- results[,convergeSubset]

        minimum_m2LL <- convergeSubset[1,which(convergeSubset["m2LL",]==min(as.numeric(convergeSubset["m2LL",])))]

        minimum_AIC <- convergeSubset[1,which(convergeSubset["AIC",]==min(as.numeric(convergeSubset["AIC",])))]

        minimum_BIC <- convergeSubset[1,which(convergeSubset["BIC",]==min(convergeSubset["BIC",]))]

        minimum_CVm2LL <- convergeSubset[1,which(convergeSubset["CV m2LL",]==min(convergeSubset["CV m2LL",]))]

        minimum_CVAIC <- convergeSubset[1,which(convergeSubset["CV AIC",]==min(convergeSubset["CV AIC",]))]

        minimum_CVBIC <- convergeSubset[1,which(convergeSubset["CV BIC",]==min(convergeSubset["CV BIC",]))]


        # getting parameters:
        if(fit_index == "m2LL"){
          best_penalty = minimum_m2LL
          reg_Model_m2LL <- createRegModel(model, model_type = model_type, fitfun = fitfun, data_type= data_type, pen_type = pen_type, pen_value = best_penalty,
                                           pen_on = pen_on, selectedDrifts =selectedDrifts, driftexpo = driftexpo, DRIFT_dt=DRIFT_dt, selectedA = selectedA,
                                           selectedS = selectedS)
          reg_Model_m2LL <- mxOption(reg_Model_m2LL, "Calculate Hessian", "No") # might cause errors; check
          reg_Model_m2LL <- mxOption(reg_Model_m2LL, "Standard Errors", "No") # might cause errors; check

          fit_reg_Model_m2LL <- mxRun(reg_Model_m2LL, silent = T)
          #set parameters below zeroThresh to zero
          fit_reg_Model_m2LL <- getFitMeasures(regmodel = fit_reg_Model_m2LL, fitfun = fitfun, model_type = model_type, cvsample = NULL, zeroThresh = zeroThresh, setZero = setZero)
          fit_reg_Model_m2LL <- fit_reg_Model_m2LL$return_Model

          out <- list("best penalty" = minimum_m2LL, "bestmodel" = fit_reg_Model_m2LL, "fit measures" = t(results), "call" = call)
        }

        if(fit_index == "AIC"){
          best_penalty = minimum_AIC
          reg_Model_AIC <- createRegModel(model, model_type = model_type, fitfun = fitfun, data_type= data_type, pen_type = pen_type, pen_value = best_penalty,
                                          pen_on = pen_on, selectedDrifts =selectedDrifts, driftexpo = driftexpo, DRIFT_dt = DRIFT_dt, selectedA = selectedA,
                                          selectedS = selectedS)
          reg_Model_AIC <- mxOption(reg_Model_AIC, "Calculate Hessian", "No") # might cause errors; check
          reg_Model_AIC <- mxOption(reg_Model_AIC, "Standard Errors", "No") # might cause errors; check

          fit_reg_Model_AIC <- mxRun(reg_Model_AIC, silent = T)
          #set parameters below zeroThresh to zero
          fit_reg_Model_AIC <- getFitMeasures(regmodel = fit_reg_Model_AIC, fitfun = fitfun, model_type = model_type, cvsample = NULL, zeroThresh = zeroThresh, setZero = setZero)
          fit_reg_Model_AIC <- fit_reg_Model_AIC$return_Model

          out <- list("best penalty" = minimum_AIC, "bestmodel" = fit_reg_Model_AIC, "fit measures" = t(results), "call" = call)
        }

        if(fit_index == "BIC"){
          best_penalty = minimum_BIC
          reg_Model_BIC <- createRegModel(model, model_type = model_type, fitfun = fitfun, data_type= data_type, pen_type = pen_type, pen_value = best_penalty,
                                          pen_on = pen_on, selectedDrifts =selectedDrifts, driftexpo = driftexpo, DRIFT_dt = DRIFT_dt, selectedA = selectedA,
                                          selectedS = selectedS)
          reg_Model_BIC <- mxOption(reg_Model_BIC, "Calculate Hessian", "No") # might cause errors; check
          reg_Model_BIC <- mxOption(reg_Model_BIC, "Standard Errors", "No") # might cause errors; check

          fit_reg_Model_BIC <- mxRun(reg_Model_BIC, silent = T)
          #set parameters below zeroThresh to zero
          fit_reg_Model_BIC <- getFitMeasures(regmodel = fit_reg_Model_BIC, fitfun = fitfun, model_type = model_type, cvsample = NULL, zeroThresh = zeroThresh, setZero = setZero)
          fit_reg_Model_BIC <- fit_reg_Model_BIC$return_Model

          out <- list("best penalty" = minimum_BIC, "bestmodel" = fit_reg_Model_BIC, "fit measures" = t(results), "call" = call)
        }

        if(fit_index == "CV_m2LL"){
          best_penalty = minimum_CVm2LL
          reg_Model_CVm2LL <- createRegModel(model, model_type = model_type, fitfun = fitfun, data_type= data_type, pen_type = pen_type, pen_value = best_penalty,
                                             pen_on = pen_on, selectedDrifts =selectedDrifts, driftexpo = driftexpo, DRIFT_dt = DRIFT_dt, selectedA = selectedA,
                                             selectedS = selectedS)
          reg_Model_CVm2LL <- mxOption(reg_Model_CVm2LL, "Calculate Hessian", "No") # might cause errors; check
          reg_Model_CVm2LL <- mxOption(reg_Model_CVm2LL, "Standard Errors", "No") # might cause errors; check

          fit_reg_Model_CVm2LL <- mxRun(reg_Model_CVm2LL, silent = T)
          #set parameters below zeroThresh to zero
          fit_reg_Model_CVm2LL <- getFitMeasures(regmodel = fit_reg_Model_CVm2LL, fitfun = fitfun, model_type = model_type, cvsample = NULL, zeroThresh = zeroThresh, setZero = setZero)
          fit_reg_Model_CVm2LL <- fit_reg_Model_CVm2LL$return_Model

          out <- list("best penalty" = minimum_CVm2LL, "bestmodel" = fit_reg_Model_CVm2LL, "fit measures" = t(results), "call" = call)
        }

        if(fit_index == "CV_AIC"){
          best_penalty = minimum_CVAIC
          reg_Model_CVAIC <- createRegModel(model, model_type = model_type, fitfun = fitfun, data_type= data_type, pen_type = pen_type, pen_value = best_penalty,
                                            pen_on = pen_on, selectedDrifts =selectedDrifts, driftexpo = driftexpo, DRIFT_dt = DRIFT_dt, selectedA = selectedA,
                                            selectedS = selectedS)
          reg_Model_CVAIC <- mxOption(reg_Model_CVAIC, "Calculate Hessian", "No") # might cause errors; check
          reg_Model_CVAIC <- mxOption(reg_Model_CVAIC, "Standard Errors", "No") # might cause errors; check

          fit_reg_Model_CVAIC <- mxRun(reg_Model_CVAIC, silent = T)
          #set parameters below zeroThresh to zero
          fit_reg_Model_CVAIC <- getFitMeasures(regmodel = fit_reg_Model_CVAIC, fitfun = fitfun, model_type = model_type, cvsample = NULL, zeroThresh = zeroThresh, setZero = setZero)
          fit_reg_Model_CVAIC <- fit_reg_Model_CVAIC$return_Model

          out <- list("best penalty" = minimum_CVAIC, "bestmodel" = fit_reg_Model_CVAIC, "fit measures" = t(results), "call" = call)
        }

        if(fit_index == "CV_BIC"){
          best_penalty = minimum_CVBIC
          reg_Model_CVBIC <- createRegModel(model, model_type = model_type, fitfun = fitfun, data_type= data_type, pen_type = pen_type, pen_value = best_penalty,
                                            pen_on = pen_on, selectedDrifts =selectedDrifts, driftexpo = driftexpo, DRIFT_dt = DRIFT_dt, selectedA = selectedA,
                                            selectedS = selectedS)
          reg_Model_CVBIC <- mxOption(reg_Model_CVBIC, "Calculate Hessian", "No") # might cause errors; check
          reg_Model_CVBIC <- mxOption(reg_Model_CVBIC, "Standard Errors", "No") # might cause errors; check

          fit_reg_Model_CVBIC <- mxRun(reg_Model_CVBIC, silent = T)
          #set parameters below zeroThresh to zero
          fit_reg_Model_CVBIC <- getFitMeasures(regmodel = fit_reg_Model_CVBIC, fitfun = fitfun, model_type = model_type, cvsample = NULL, zeroThresh = zeroThresh, setZero = setZero)
          fit_reg_Model_CVBIC <- fit_reg_Model_CVBIC$return_Model

          out <- list("best penalty" = minimum_CVBIC, "bestmodel" = fit_reg_Model_CVBIC, "fit measures" = t(results), "call" = call)
        }


        class(out) <- "FitLaremmObject"



        return(out)



    }

}
