#' Summary function for FitLaremmObject
#'
#' @param object Object from fitRegModels.
#' @param ... additional parameters.
#' @export
summary.FitLaremmObject <- function(object,...){

  if(object$call$model_type == "mxModel"){
    # get fit measures
    m2LL <- object$`fit measures`[which(object$`fit measures`[,"penalty"] == object$`best penalty`), "m2LL"]
    AIC <- object$`fit measures`[which(object$`fit measures`[,"penalty"] == object$`best penalty`), "AIC"]
    BIC <- object$`fit measures`[which(object$`fit measures`[,"penalty"] == object$`best penalty`), "BIC"]
    if(object$call$CV){
      CV_m2LL <- object$`fit measures`[which(object$`fit measures`[,"penalty"] == object$`best penalty`), "CV m2LL"]
    }else{CV_m2LL <- NA}

    estimatedParameters <- object$`fit measures`[which(object$`fit measures`[,"penalty"] == object$`best penalty`), "estimated Parameters"]
    convergence <- if(any(!object$`fit measures`[,"convergence"] ==0)){"Problems with convergence occured"}else{"No convergence problems"}
    negativeVariances <- if(any(!object$`fit measures`[,"negative variances"] ==0)){"Negative variances occured"}else{"No negative variances"}

    # get parameters
    # A
    parNames <- c(object$bestmodel$manifestVars, object$bestmodel$latentVars)
    Aparam <- round(object$bestmodel$A$values[object$bestmodel$A$free], 3)
    expLabel <- matrix(NA, ncol = ncol(object$bestmodel$A$values), nrow = nrow(object$bestmodel$A$values))
    for(row in 1:nrow(expLabel)){
      for(col in 1:ncol(expLabel)){
        expLabel[row,col] = paste(parNames[col] , "->", parNames[row])
      }
    }
    labels <- object$bestmodel$A$labels[object$bestmodel$A$free]
    if(any(!is.na(labels))){
      AparMat <- matrix(c(expLabel[object$bestmodel$A$free],labels, Aparam), nrow = length(expLabel[object$bestmodel$A$free]), ncol = 3)

    }else{
      AparMat <- matrix(c(expLabel[object$bestmodel$A$free], Aparam), nrow = length(expLabel[object$bestmodel$A$free]), ncol = 2)
    }

    # S
    Sparam <- round(object$bestmodel$S$values[object$bestmodel$S$free], 3)
    expLabel <- matrix(NA, ncol = ncol(object$bestmodel$S$values), nrow = nrow(object$bestmodel$S$values))
    for(row in 1:nrow(expLabel)){
      for(col in 1:ncol(expLabel)){
        expLabel[row,col] = paste(parNames[col] , "<->", parNames[row])
      }
    }
    Slabels <- object$bestmodel$S$labels[object$bestmodel$S$free]
    if(any(!is.na(Slabels))){
      SparMat <- matrix(c(expLabel[object$bestmodel$S$free],Slabels, Sparam), nrow = length(expLabel[object$bestmodel$S$free]), ncol = 3)

    }else{
      SparMat <- matrix(c(expLabel[object$bestmodel$S$free], Sparam), nrow = length(expLabel[object$bestmodel$S$free]), ncol = 2)
    }
    # M
    if(any(object$bestmodel$M$free)){
      Mparam <- round(object$bestmodel$M$values[object$bestmodel$M$free], 3)
      expLabel <- matrix(NA, ncol = ncol(object$bestmodel$M$values), nrow = nrow(object$bestmodel$M$values))
      for(row in 1:nrow(expLabel)){
        for(col in 1:ncol(expLabel)){
          expLabel[row,col] = paste("intercept", parNames[col])
        }
      }
      Mlabels <- object$bestmodel$M$labels[object$bestmodel$M$free]
      if(any(!is.na(Mlabels))){
        MparMat <- matrix(c(expLabel[object$bestmodel$M$free],Mlabels, Mparam), nrow = length(expLabel[object$bestmodel$M$free]), ncol = 3)

      }else{
        MparMat <- matrix(c(expLabel[object$bestmodel$M$free], Mparam), nrow = length(expLabel[object$bestmodel$M$free]), ncol = 2)
      }
    }
    if(any(object$bestmodel$M$free)){
    retList <- list(
      "object type" = object$call$model_type,
      "penalty value" = paste("Best penalty value based on", object$call$fit_index ,"was" , object$`best penalty`),
      "A-parameter" = AparMat,
      "S-pararameter" = SparMat,
      "M-parameters" = MparMat,
      "m2LL" = m2LL,
      "AIC" = AIC,
      "BIC" = BIC,
      "CV m2LL" =  CV_m2LL,
      "estimated parameters" = estimatedParameters,
      "negative variances" = negativeVariances,
      "convergence" = convergence


      )
    }else{
      retList <- list(
        "object type" = object$call$model_type,
        "penalty value" = paste("Best penalty value based on", object$call$fit_index ,"was" , object$`best penalty`),
        "A-parameter" = AparMat,
        "S-pararameter" = SparMat,
        #"M-parameters" = MparMat,
        "m2LL" = m2LL,
        "AIC" = AIC,
        "BIC" = BIC,
        "CV m2LL" =  CV_m2LL,
        "estimated parameters" = estimatedParameters,
        "negative variances" = negativeVariances,
        "convergence" = convergence

      )

    }

  }else if(object$call$model_type == "ctsem"){
    print("Not fully implemented. Elements can be accessed with object$bestmodel$BaseModel...")
    # get fit measures
    m2LL <- object$`fit measures`[which(object$`fit measures`[,"penalty"] == object$`best penalty`), "m2LL"]
    AIC <- object$`fit measures`[which(object$`fit measures`[,"penalty"] == object$`best penalty`), "AIC"]
    BIC <- object$`fit measures`[which(object$`fit measures`[,"penalty"] == object$`best penalty`), "BIC"]
    if(object$call$CV){
      CV_m2LL <- object$`fit measures`[which(object$`fit measures`[,"penalty"] == object$`best penalty`), "CV m2LL"]
    }else{CV_m2LL <- NA}

    estimatedParameters <- object$`fit measures`[which(object$`fit measures`[,"penalty"] == object$`best penalty`), "estimated Parameters"]
    convergence <- if(any(!object$`fit measures`[,"convergence"] ==0)){"Problems with convergence occured"}else{"No convergence problems"}
    negativeVariances <- if(any(!object$`fit measures`[,"negative variances"] ==0)){"Negative variances occured"}else{"No negative variances"}

    Drift <- object$bestmodel$BaseModel$DRIFT$values
    retList <- list(
      "object type" = object$call$model_type,
      "penalty value" = paste("Best penalty value based on", object$call$fit_index ,"was" , object$`best penalty`),
      "Drift" = Drift,
      "m2LL" = m2LL,
      "AIC" = AIC,
      "BIC" = BIC,
      "CV m2LL" =  CV_m2LL,
      "estimated parameters" = estimatedParameters,
      "negative variances" = negativeVariances,
      "convergence" = convergence

    )

  }



  class(retList) <- "summary.FitLaremmObject"
  print(retList)
  #return(list(retList, "Full Model" = object$bestmodel))
}

