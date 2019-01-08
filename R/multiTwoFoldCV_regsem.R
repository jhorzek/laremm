#' compute multiple 2 - fold - CVs with regsem
#'
#' Note: This is an extension to the package \pkg{regsem}. \pkg{regsem} can no longer perform 2-fold-CV; this functionality is now added with this script. 
#'
#' @param lavaanobject already run lavaan model object
#' @param type only use "lasso" 
#' @param pars_pen Vector with parameters that should be penalized
#' @param lambda_values Vector with penalty values that are to be used
#' @param cov_test_data covairance matrix of the test sample
#' @param numObs_test_data Number of observations in the test sample
#' @param use_unbiasedCov use unbiased covariance for computation of FML. Default is FALSE
#'
#' @export


multiTwoFoldCV_regsem <- function(lavaanobject, type = "lasso", pars_pen, lambda_values, cov_test_data, numObs_test_data, use_unbiasedCov = FALSE){
  iteration = 1
  results <- matrix(NA, nrow = length(lambda_values), ncol = 3)
  colnames(results) <- c("penalty", "TrainFML", "TestFML")
  for(pen in lambda_values){
    results[iteration, 1] = pen
    regsemmodel = regsem(lavaanobject,type = type, pars_pen= pars_pen, lambda = pen, use_unbiasedCov)
    results[iteration, 2] = regsemmodel$fit
    results[iteration, 3] = twoFoldCV_regsem(regsemmodel, cov_test_data, numObs_test_data = nrow(test_data))
    iteration = iteration + 1 
  }
  return(results)
}


