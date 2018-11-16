#' createRegModel
#'
#' Note: laremm is based on the R package \pkg{regsem}. Because of the early status of laremm, it is recommended to use regsem instead!
#' createRegModel creates a regularized model from a mxModel or ctsem.
#'
#' @param model mxModel or ctsem object
#' @param model_type specify the type of model provided: ctsem or mxModel
#' @param fitfun fitfunction to be used in the fitting procedure. Either FML or FIML
#' @param data_type type of data in the model. Either "cov" or "raw"
#' @param penalty_type so far only "lasso" implemented
#' @param pen_value numeric value of penalty size
#' @param pen_on string vector with matrices that should be regularized. Possible are combinations of "A", "S", "DRIFT"
#' @param selectedDrifts drift values to regularize. Possible are "all", "cross", "auto" or providing a matrix of the same size as the drift matrix with ones for every parameter to regularize and 0 for every non-regularized parameter
#' @param driftexpo specifiy if the regularization will be performed on the raw drift matrix or on the exponential of the drift matrix (discrete time parameters)
#' @param DRIFT_dt provide the discrete time points for which the drift will be regularized. A vector with multiple values is possible
#' @param selectedA A values to regularize. Possible are "all", or providing a matrix of the same size as the A matrix with ones for every parameter to regularize and 0 for every non-regularized parameter
#' @param selectedS S values to regularize. Possible are "all", or providing a matrix of the same size as the S matrix with ones for every parameter to regularize and 0 for every non-regularized parameter
#'
#' @export


createRegModel <- function(model, model_type = "ctsem", fitfun = "FIML", data_type= "raw",pen_type = "lasso", pen_value = 0,
                               pen_on = "none",selectedDrifts ="none", driftexpo = TRUE, DRIFT_dt =1, selectedA = "none", selectedS = "none"){
  # to be implemented later:
  #, selectedLambda = "none",
  #selectedCint = "none",
  #selectedTdpredeffect = "none",
  #selectedDiffusion = "none",


  # get the mxModel object of the provided model:
  if(model_type=="ctsem"){
    modelobject = model$mxobj #read mxmodel-object and ctfit-object
    ctobject = model$ctmodelobj}
  else if(model_type =="mxModel"){
    modelobject = model
  }

  if(data_type == "raw")
    numObs <- nrow(modelobject$data$observed) # numObs is number of observed persons
  else if (data_type == "cov"){
    numObs <- modelobject$data$numObs
  }

  numObs <- mxMatrix(type= "Full", nrow= 1, ncol = 1, free = FALSE, values = numObs,name = "numObs") # define peanlty value
  pen_value <- mxMatrix(type= "Full", nrow= 1, ncol = 1, free = FALSE, values = pen_value,name = "pen_value") # define peanlty value

  # Basis for model building:
  BaseModel <- mxModel(modelobject, name = "BaseModel") # has all the parameters and the base fit function (FML or FIML)

  outModel <- mxModel(model= "outModel", # model that is returned
                      BaseModel, # BaseModel is a submodel of outModel. The elements of BaseModel can be accessed by the outModel
                      numObs,
                      pen_value)



  # What matrices should be penalized?

  # set base fitfunction:

  fitfun_string <- "BaseModel.fitfunction"


  ## Define the penalty on A
  if("a" %in% sapply(pen_on, tolower)){ # if a is in pen_on

    if (class(selectedA)=="matrix"){ # if selected A is a matrix
      selectedA_values = selectedA # use this matrix
    }

    else if (selectedA =="none"){ # if selectedA is none
      selectedA_values = matrix(0,nrow(modelobject$A$values),ncol(modelobject$A$values),byrow = T) # set all values to 0
    }

    else if (selectedA=="all"){ # if selected A is all
      selectedA_values = matrix(1,nrow(modelobject$A$values),ncol(modelobject$A$values),byrow = T) # set all values to 1
    }

    # define the penalty values
    selectedA_values <- mxMatrix(type = "Full", values = selectedA_values, free = F, name ="selectedA_values")
    penalty_A = mxAlgebra(numObs*(pen_value*(t(abs(cvectorize(BaseModel.A))) %*% cvectorize(selectedA_values))), name = "penalty_A") # N*penalty to keep equivalence to Jacobucci's penalty function

    # defining the regularized fitfunction:
    fitfun_string <- paste(fitfun_string,"penalty_A", sep = " + ")

    # expand Model
    outModel <- mxModel(outModel,
                        penalty_A,
                        selectedA_values
                        )
  }

  # Define the penalty on S
  if("s" %in% sapply(pen_on, tolower)){

    if (class(selectedS)=="matrix"){
      selectedS_values = selectedS
    }
    else if (selectedS =="none"){
      selectedS_values = matrix(0,nrow(modelobject$S$values),ncol(modelobject$S$values),byrow = T)
    }

    else if (selectedS=="all"){
      selectedS_values = matrix(1,nrow(modelobject$S$values),ncol(modelobject$S$values),byrow = T)
    }

    # define the penalty values
    selectedS_values <- mxMatrix(type = "Full", values = selectedS_values, free = F, name ="selectedS_values")
    penalty_S = mxAlgebra(numObs(pen_value*(t(abs(cvectorize(BaseModel.S))) %*% cvectorize(selectedS_values))), name = "penalty_S")

    # defining the regularized fitfunction:
    fitfun_string <- paste(fitfun_string,"penalty_S", sep = " + ")

    # expand Model
    outModel <- mxModel(outModel,
                        penalty_S,
                        selectedS_values
                        )


  }

  # Define penalty on Drift
  if("drift" %in% sapply(pen_on, tolower)){
    if (class(selectedDrifts)=="matrix"){
      selectedDrifts_values = selectedDrifts
    }
    else if (selectedDrifts=="none"){
      selectedDrifts_values = matrix(0,nrow(modelobject$DRIFT$values),ncol(modelobject$DRIFT$values),byrow = T)
    }
    else if (selectedDrifts=="all"){
      selectedDrifts_values = matrix(1,nrow(modelobject$DRIFT$values),ncol(modelobject$DRIFT$values),byrow = T)
    }
    else if (selectedDrifts=="auto"){
      selectedDrifts_values = diag(1,nrow(modelobject$DRIFT$values))
    }
    else if (selectedDrifts =="cross"){
      selectedDrifts_values = matrix(1,nrow(modelobject$DRIFT$values),ncol(modelobject$DRIFT$values),byrow = T) - diag(1,nrow(modelobject$DRIFT$values))
    }


    # define the penalty values
    selectedDrifts_values <- mxMatrix("Full", values = selectedDrifts_values, free = FALSE, name = "selectedDrifts_values")
    if (driftexpo == TRUE){
      # get delta t values for which to penalize:

      drift_count = 1
      drift_string <- rep(NA, length(DRIFT_dt))
      for (i in DRIFT_dt){
        drift_string[drift_count] <- paste("numObs*(pen_value*(t(abs(cvectorize(expm(BaseModel.DRIFT*", i, "))))%*%cvectorize(selectedDrifts_values)))", sep = "")
        drift_count <- drift_count+1
      }

      drift_string_combined <- paste(drift_string, collapse = "+")

      penalty_DRIFT = mxAlgebraFromString(drift_string_combined, name = "penalty_DRIFT")
    }
    else if (driftexpo == FALSE){
      penalty_DRIFT = mxAlgebra(numObs*(pen_value*(t(abs(cvectorize(BaseModel.DRIFT)))%*%cvectorize(selectedDrifts_values))), name = "penalty_DRIFT")

    }

    # defining the regularized fitfunction:
    fitfun_string <- paste(fitfun_string,"penalty_DRIFT", sep = " + ")


    # expand Model
    outModel <- mxModel(outModel,
                        penalty_DRIFT,
                        selectedDrifts_values
    )

  }

  # define fitfunction:
  regfit_algebra <- mxAlgebraFromString(fitfun_string, name = "regfit_algebra")
  regfit_fun <- mxFitFunctionAlgebra("regfit_algebra")

  # complete model
  outModel <- mxModel(outModel,
                      regfit_algebra,
                      regfit_fun
  )
  outModel <- mxOption(outModel, "Calculate Hessian", "No")
  outModel <- mxOption(outModel, "Standard Errors", "No")

  # return model
  return(outModel)
}
