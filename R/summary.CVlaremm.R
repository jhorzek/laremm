#' Summary function for kFoldCV
#'
#' @param object Object from kFoldCV
#' @param ... additional parameters.
#' @export
summary.CVlaremm <- function(object,...){
  k <- object$k
    # get fit measures
    CV_m2LL <- object$`CV results`[which(object$`CV results`[,"penalty"] == object$`best penalty`), "sum CV_m2LL"]

    convergence <- if(any(!object$`CV results`[,"convergence problems"] ==0)){"Problems with convergence occured"}else{"No convergence problems"}
    negativeVariances <- if(any(!object$`CV results`[,"negative variances"] ==0)){"Negative variances occured"}else{"No negative variances"}

    # get parameters
    # A
    parNames <- c(object$`final Model`$BaseModel$manifestVars, object$`final Model`$BaseModel$latentVars)
    Aparam <- round(object$`final Model`$BaseModel$A$values[object$`final Model`$BaseModel$A$free], 3)
    expLabel <- matrix(NA, ncol = ncol(object$`final Model`$BaseModel$A$values), nrow = nrow(object$`final Model`$BaseModel$A$values))
    for(row in 1:nrow(expLabel)){
      for(col in 1:ncol(expLabel)){
        expLabel[row,col] = paste(parNames[col] , "->", parNames[row])
      }
    }
    labels <- object$`final Model`$BaseModel$A$labels[object$`final Model`$BaseModel$A$free]
    if(any(!is.na(labels))){
      AparMat <- matrix(c(expLabel[object$`final Model`$BaseModel$A$free],labels, Aparam), nrow = length(expLabel[object$`final Model`$BaseModel$A$free]), ncol = 3)

    }else{
      AparMat <- matrix(c(expLabel[object$`final Model`$BaseModel$A$free], Aparam), nrow = length(expLabel[object$`final Model`$BaseModel$A$free]), ncol = 2)
    }

    # S
    Sparam <- round(object$`final Model`$BaseModel$S$values[object$`final Model`$BaseModel$S$free], 3)
    expLabel <- matrix(NA, ncol = ncol(object$`final Model`$BaseModel$S$values), nrow = nrow(object$`final Model`$BaseModel$S$values))
    for(row in 1:nrow(expLabel)){
      for(col in 1:ncol(expLabel)){
        expLabel[row,col] = paste(parNames[col] , "<->", parNames[row])
      }
    }
    Slabels <- object$`final Model`$BaseModel$S$labels[object$`final Model`$BaseModel$S$free]
    if(any(!is.na(Slabels))){
      SparMat <- matrix(c(expLabel[object$`final Model`$BaseModel$S$free],Slabels, Sparam), nrow = length(expLabel[object$`final Model`$BaseModel$S$free]), ncol = 3)

    }else{
      SparMat <- matrix(c(expLabel[object$`final Model`$BaseModel$S$free], Sparam), nrow = length(expLabel[object$`final Model`$BaseModel$S$free]), ncol = 2)
    }
    # M
    if(any(object$`final Model`$BaseModel$M$free)){
      Mparam <- round(object$`final Model`$BaseModel$M$values[object$`final Model`$BaseModel$M$free], 3)
      expLabel <- matrix(NA, ncol = ncol(object$`final Model`$BaseModel$M$values), nrow = nrow(object$`final Model`$BaseModel$M$values))
      for(row in 1:nrow(expLabel)){
        for(col in 1:ncol(expLabel)){
          expLabel[row,col] = paste("intercept", parNames[col])
        }
      }
      Mlabels <- object$`final Model`$BaseModel$M$labels[object$`final Model`$BaseModel$M$free]
      if(any(!is.na(Mlabels))){
        MparMat <- matrix(c(expLabel[object$`final Model`$BaseModel$M$free],Mlabels, Mparam), nrow = length(expLabel[object$`final Model`$BaseModel$M$free]), ncol = 3)

      }else{
        MparMat <- matrix(c(expLabel[object$`final Model`$BaseModel$M$free], Mparam), nrow = length(expLabel[object$`final Model`$BaseModel$M$free]), ncol = 2)
      }
    }
    if(any(object$`final Model`$BaseModel$M$free)){
      retList <- list(
        "object type" = class(object),
        "penalty value" = paste("Best penalty value based on", k ,"fold CV was" , object$`best penalty`),
        "A-parameter" = AparMat,
        "S-pararameter" = SparMat,
        "M-parameters" = MparMat,
        "sum CV m2LL" =  CV_m2LL,
        "negative variances" = negativeVariances,
        "convergence" = convergence


      )
    }else{
      retList <- list(
        "object type" = class(object),
        "penalty value" = paste("Best penalty value based on", k ,"fold CV was" , object$`best penalty`),
        "A-parameter" = AparMat,
        "S-pararameter" = SparMat,
        #"M-parameters" = MparMat,
        "CV m2LL" =  CV_m2LL,
        "negative variances" = negativeVariances,
        "convergence" = convergence

      )

    }

  class(retList) <- "summary.CVlaremm"
  print(retList)
  #return(list(retList, "Full Model" = object$bestmodel))
}

