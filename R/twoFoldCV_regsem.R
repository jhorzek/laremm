#' compute one 2 - fold - CV with regsem
#'
#' Note: This is an extension to the package \pkg{regsem}. \pkg{regsem} can no longer perform 2-fold-CV; this functionality is now added with this script. For multiple 2-Fold-CVs use multiToFoldCV_regsem.
#'
#' @param regsemmodel already run regsemmodel model object
#' @param cov_test_data covairance matrix of the test sample
#' @param numObs_test_data Number of observations in the test sample
#' @param use_unbiasedCov use unbiased covariance for computation of FML. Default is FALSE
#'
#' @export

twoFoldCV_regsem <- function(regsemmodel, cov_test_data, numObs_test_data, use_unbiasedCov = FALSE){

  # get expected covariance:
  expCov <- regsemmodel$Imp_Cov
  obsCov_unbiased <- cov_test_data
  obsCov_biased <- obsCov_unbiased*(numObs_test_data-1)/(numObs_test_data)
  obs_manif <- regsemmodel$nvar

  if(use_unbiasedCov == FALSE)
  {
  obsCov <- obsCov_biased
  }
  else{
  obsCov <- obsCov_unbiased
  }

  F_ML <- .5*(log(det(expCov)) + tr(obsCov %*% solve(expCov)) - obs_manif - log(det(obsCov)))


  return_values <- matrix(c(F_ML), ncol = 1, nrow = 1)
  colnames(return_values) <- c(".5*FML")

  return(return_values)
  
}
