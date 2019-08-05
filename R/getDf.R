#' get degrees of freedom
getDf <- function(regmodel){
  print("Warning: so far only for FML with mxModels")
  estimatedParam <- 0 # number of free parameters in the model
  observedParam <- 0 # number of observed parameters
  matrices <- regmodel$BaseModel$matrices
  for (matrix in matrices){
    estimatedParam <- estimatedParam +sum(matrix$free==TRUE)
  }
  observedParam <- length(regmodel$BaseModel$manifestVars)*(length(regmodel$BaseModel$manifestVars)+1)/2
  if(!is.na(regmodel$BaseModel$expectation$M)){
    observedParam <- observedParam + length(regmodel$BaseModel$manifestVars) # add mean structure
  }
  df <- observedParam - estimatedParam # calculate df
  return(df)
}

