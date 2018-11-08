#' create cv folds for ctsem
#'
createCVfolds <- function(data, k_folds){
  # used for creating folds in ctsem
  library(caret)
  rowind <- c(1:nrow(data)) # create rowindicator vector; caret only works on wide format
  foldindicators <- createFolds(rowind, k = k_folds)
  return(foldindicators)
}
