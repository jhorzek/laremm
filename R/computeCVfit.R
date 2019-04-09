#' compute fit of cross valudation
#' @author Jannik Orzek
#' @import OpenMx ctsem
#' @export
computeCVfit <- function(model, ctmodelobj, testset, numObs, model_type = "ctsem", fitfunction = "FIML"){
  if(model_type == "ctsem" & fitfunction == "FIML"){
    if (is.null(ctmodelobj)){print("Please provide ctmodelobject")
      break}

  expCov_test <- mxGetExpected(model, component = 'covariance', subname = "BaseModel") # compute the expected covairance; is independent of the observes data
  expMean_test <- mxGetExpected(model, component = 'means', subname = "BaseModel") # compute the expected means; is independent of the observes data
  frequencyvector_test <- matrix(c(1), nrow = nrow(testset)) # m2LL function requests a frequency - vector; here: vector with 1s of length of dataset

  filtervector <- c(1:(ctmodelobj$n.manifest*ctmodelobj$Tpoints+ctmodelobj$n.TDpred*(ctmodelobj$Tpoints)+ctmodelobj$n.TIpred)) ### ACHTUNG: leichte Ver?nderung zur Vorgabe in ctsem; dort: model$ctmodelobj$n.TDpred*(model$ctmodelobj$Tpoints-1); das -1 f?hrt allerdings zu falschen Ergebnissen (auch tdpred f?ngt bei 0 an)
  #     #Filtervector, which separates the observations from the time intervalls

  fit_testset <- m2LL(data_unique = as.matrix(testset[,filtervector]),
                      data_frequency = frequencyvector_test,
                      Sigma = as.matrix(expCov_test),
                      Mu = as.matrix(expMean_test) ) # compute m2LL for the test-set
  return(fit_testset) # return the m2LL
  }
  if(model_type == "mxModel"){
    if(fitfunction == "FIML"){
      expCov_test <- mxGetExpected(model, component = 'covariance', subname = "BaseModel") # compute the expected covairance; is independent of the observes data
      expMean_test <- mxGetExpected(model, component = 'means', subname = "BaseModel") # compute the expected means; is independent of the observes data
      frequencyvector_test <- matrix(c(1), nrow = nrow(testset)) # m2LL function requests a frequency - vector; here: vector with 1s of length of dataset

      fit_testset <- m2LL(data_unique = as.matrix(testset),
                          data_frequency = frequencyvector_test,
                          Sigma = as.matrix(expCov_test),
                          Mu = as.matrix(expMean_test) ) # compute m2LL for the test-set
      return(fit_testset) # return the m2LL

    }
    if(fitfunction == "FML"){
      expCov_test <- mxGetExpected(model, component = 'covariance', subname = "BaseModel") # compute the expected covairance; is independent of the observes data
      expMean_test <- mxGetExpected(model, component = 'means', subname = "BaseModel") # compute the expected means; is independent of the observes data

      Num_Obs <- numObs
      biasCov <- testset*((Num_Obs-1)/(Num_Obs))
      if(is.null(model$BaseModel$data$means) | is.na(model$BaseModel$data$means)) {
        testset_FML <- Num_Obs*( log(det(expCov_test)) + tr(biasCov %*% solve(expCov_test)))
      }
      else(
        testset_FML <- Num_Obs*( log(det(expCov_test)) + tr(biasCov %*% solve(expCov_test)) + t(means - expMean_test) %*% solve(expCov_test) %*% (means - expMean_test))
      )
      return(testset_FML)

    }

  }
}
