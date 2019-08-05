#' create cv folds for ctsem
#'
#' @author Jannik Orzek
#' @import caret
#' @export
createCVfolds <- function(data, k_folds){
  # used for creating folds in ctsem
  rowind <- c(1:nrow(data)) # create rowindicator vector; caret only works on wide format
  foldindicators <- createFolds(rowind, k = k_folds)
  return(foldindicators)
}
