#' kFoldCV
#'
#' Note: laremm is based on the R package \pkg{regsem}. Because of the early status of laremm, it is recommended to use regsem instead!
#' kFoldCV uses k fold cross-validation with fitRegModels
#'
#' @param k specifies the number of splits (e.g. k = 5 for 5-fold-CV)
#' @param model mxModel with the full data set in model$data. The data set will be split by kFoldCV
#' @param model_type specify the type of model provided: only mxModel supported
#' @param fitfun fitfunction to be used in the fitting procedure. Currently only FIML implemented
#' @param data_type type of data in the model. Only "raw" supported
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
#' @param zeroThresh threshold for evaluating regularized parameters as zero. Default is .001 similar to \pkg{regsem}
#' @param setZero should parameters below zeroThresh be set to zero in all fit calculations. Default is FALSE, similar to \pkg{regsem}
#'
#'
#' @author Jannik Orzek
#' @import OpenMx caret
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
#' kFolRreg_model <- kFoldCV(model = fit_myModel, model_type = "mxModel", fitfun = "FIML",
#'                           pen_on = "A", selectedA = selectedA,
#'                           pen_start = 0, pen_end = .05, pen_stepsize = .01,
#'                           k = 5
#'                           )
#' summary(kFolRreg_model)
#' # inspect results:
#' kFolRreg_model$`CV results`
#'
#' @export
#'

kFoldCV <- function(k = 5, model,
                    model_type = "mxModel",
                    fitfun = "FIML",
                    data_type= "raw",
                    pen_type = "lasso",
                    pen_on = "none",
                    selectedDrifts ="none",
                    driftexpo = TRUE,
                    selectedA = "none",
                    selectedS = "none",
                    pen_start = 0,
                    pen_end = 1,
                    pen_stepsize = .01,
                    zeroThresh = .001,
                    setZero = FALSE){

  if(!model_type == "mxModel"){
    disp("Currently only implemented for mxModels.")
    return(NA)
  }
  if(!data_type == "raw"){
    disp("Currently only implemented for raw data")
    return(NA)
  }
  if(!fitfun == "FIML"){
    disp("Currently only implemented for FIML")
    return(NA)
  }

  fit_index = "CV_m2LL"
  CV = TRUE

  # split dataset
  sampleSize <- model$data$numObs
  Folds <- createFolds(c(1:sampleSize), k)
  full_raw_data <- model$data$observed

  Res <- matrix(NA, nrow = length(seq(from = pen_start, to = pen_end, by = pen_stepsize)), ncol = k+4)
  colnames(Res) <- c("penalty", "mean CV_m2LL",paste("fold", 1:k), "negative variances", "convergence problems")
  Res[,"penalty"] <- seq(from = pen_start, to = pen_end, by = pen_stepsize)
  Res[,"negative variances"] <- 0
  Res[,"convergence problems"] <- 0

  fold <- 1
  for(Fold in Folds){
    Test_Sample <- full_raw_data[Fold,]
    Train_Sample <- full_raw_data[-Fold,]

    Train_Sample <- scale(Train_Sample)
    Test_Sample <- scale(Test_Sample)

    tempModel <- model
    tempModel$data <- mxData(observed = Train_Sample, type = "raw")
    Test_Sample <- mxData(observed = Test_Sample, type = "raw")

    tempFit <- fitRegModels(model = tempModel, model_type = model_type, fitfun = fitfun, data_type = data_type, pen_type, pen_on,
                 selectedDrifts = selectedDrifts, driftexpo = driftexpo, selectedA = selectedA,
                 selectedS = selectedS, pen_start = pen_start, pen_end = pen_end,
                 pen_stepsize = pen_stepsize, fit_index = fit_index, CV = CV,
                 Test_Sample = Test_Sample, zeroThresh = zeroThresh, setZero = setZero)

    Res[,paste("fold", fold)] <- tempFit$`fit measures`[,'CV m2LL']
    Res[,"negative variances"] <- Res[,"negative variances"] + tempFit$`fit measures`[,'negative variances']
    Res[,"convergence problems"] <- Res[,"convergence problems"] + tempFit$`fit measures`[,'convergence']

    cat(paste("\n Completed CV for fold", fold, "of", k, "\n"))

    fold <- fold + 1
  }

  # mean the m2LLs:
  Res[,"mean CV_m2LL"] <- apply(Res[, paste("fold", 1:k)], 1, mean)

  # only use runs without problems:
  Res_final <- Res[Res[,"negative variances"] == 0 && Res[,"convergence problems"] == 0,]

  # find best penalty value:
  best_penalty <- Res_final[which(Res_final[,"mean CV_m2LL"] == min(Res_final[,"mean CV_m2LL"])), "penalty"]

  # fit best penalty model with full data set:

  finalModel <- createRegModel(model, model_type, fitfun, data_type,pen_type, pen_value = best_penalty,
  pen_on,selectedDrifts, driftexpo, DRIFT_dt, selectedA, selectedS)

  ffinalModel <- mxRun(finalModel)

  ret <- list("CV results" = Res, "final Model" = ffinalModel, "best penalty" = best_penalty, "k" = k)
  class(ret) <- "CVlaremm"

  return(ret)
}
