getEstimatedParameters <- function(regmodel, zeroThresh = .001){

  # free parameter:
  matrices <- regmodel$BaseModel$matrices

  new_model <- regmodel$BaseModel

  estimated_params <- 0

  # values are only set to 0 and set as free = FALSE, if ther absolute value is smaller than zeroThresh and if they have been penalized.
  # this procedure is equivalent to the one used in regsem
  # if penalty is on A:
  if(regmodel$pen_value$values >0){
    if(any( regmodel$selectedA_values$values == 1)){
      zero_matrix <- which(abs(matrices$A$values) < zeroThresh & matrices$A$free & regmodel$selectedA_values$values == 1)
      matrices$A$values[zero_matrix] = 0 # values smaller than zeroThresh and freely estimated are set to 0
      matrices$A$free[zero_matrix] = FALSE # these are assumed as fixed to compute the BIC / AIC
      new_model <- mxModel(new_model,  matrices$A)
    }
    # if penalty is on S
    if (any( regmodel$selectedS_values$values == 1)){
      zero_matrix <- which(abs(matrices$S$values) < zeroThresh & matrices$S$free & regmodel$selectedS_values$values == 1)
      matrices$S$values[zero_matrix] = 0 # values smaller than zeroThresh and freely estimated are set to 0
      matrices$S$free[zero_matrix] = FALSE # these are assumed as fixed to compute the BIC / AIC
      new_model <- mxModel(new_model,  matrices$S)
    }
    # if penalty is on DRIFT
    if (any( regmodel$selectedDrifts_values$values == 1)){
      zero_matrix <- which(abs(matrices$DRIFT$values) < zeroThresh & matrices$DRIFT$free & regmodel$selectedDrifts_values$values == 1)
      matrices$DRIFT$values[zero_matrix] = 0 # values smaller than zeroThresh and freely estimated are set to 0
      matrices$DRIFT$free[zero_matrix] = FALSE # these are assumed as fixed to compute the BIC / AIC
      new_model <- mxModel(new_model,  matrices$DRIFT)
    }
    # if penalty is on Lambda
    if (any( regmodel$selectedLambda_values$values == 1)){
      zero_matrix <- which(abs(matrices$LAMBDA$values) < zeroThresh & matrices$LAMBDA$free & regmodel$selectedLambda_values$values == 1)
      matrices$LAMBDA$values[zero_matrix] = 0 # values smaller than zeroThresh and freely estimated are set to 0
      matrices$LAMBDA$free[zero_matrix] = FALSE # these are assumed as fixed to compute the BIC / AIC
      new_model <- mxModel(new_model,  matrices$LAMBDA)
    }
    # if penalty is on Cint
    if (any( regmodel$selectedCint_values$values == 1)){
      zero_matrix <- which(abs(matrices$CINT$values) < zeroThresh & matrices$CINT$free & regmodel$selectedCint_values$values == 1)
      matrices$CINT$values[zero_matrix] = 0 # values smaller than zeroThresh and freely estimated are set to 0
      matrices$CINT$free[zero_matrix] = FALSE # these are assumed as fixed to compute the BIC / AIC
      new_model <- mxModel(new_model,  matrices$CINT)
    }
    # if penalty is on TDPREDEFFECT
    if (any( regmodel$selectedTdpredeffect_values$values == 1)){
      stop("not yet implemented")
      #zero_matrix <- which(abs(matrices$$values) < zeroThresh & matrices$CINT$free & regmodel$selectedTdpredeffect_values$values == 1)
      #matrices$CINT$values[zero_matrix] = 0 # values smaller than zeroThresh and freely estimated are set to 0
      #matrices$CINT$free[zero_matrix] = FALSE # these are assumed as fixed to compute the BIC / AIC
      #new_model <- mxModel(new_model,  matrices$CINT)
    }
    # if penalty is on TDPREDEFFECT
    if (any( regmodel$selectedDiffusion_values$values == 1)){
      stop("not yet implemented")
      #zero_matrix <- which(abs(matrices$$values) < zeroThresh & matrices$CINT$free & regmodel$selectedDiffusion_values$values == 1)
      #matrices$CINT$values[zero_matrix] = 0 # values smaller than zeroThresh and freely estimated are set to 0
      #matrices$CINT$free[zero_matrix] = FALSE # these are assumed as fixed to compute the BIC / AIC
      #new_model <- mxModel(new_model,  matrices$CINT)
    }
  }

  # get estimated parameters:
  for(matrix in matrices){
    if(any(!is.na(matrix$labels))){
      # elements without labels that are free:
      sum1 <- sum(is.na(matrix$labels)&& matrix$free)
      # unique elements with labels that are free
      sum2 <- length(unique(matrix$labels[matrix$free]))
      estimated_params <- estimated_params + sum1 + sum2
    }else{
      estimated_params <- estimated_params + sum(matrix$free)
    }

  }

  retval = list("new_model" = new_model, "estimated parameters" = estimated_params)

  return(retval)

}
