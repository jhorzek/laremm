% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/createRegModel.R
\name{createRegModel}
\alias{createRegModel}
\title{createRegModel}
\usage{
createRegModel(model, model_type = "ctsem", fitfun = "FIML",
  data_type = "raw", pen_type = "lasso", pen_value = 0,
  pen_on = "none", selectedDrifts = "none", driftexpo = TRUE,
  DRIFT_dt = 1, selectedA = "none", selectedS = "none")
}
\arguments{
\item{model}{mxModel or ctsem object}

\item{model_type}{specify the type of model provided: ctsem or mxModel}

\item{fitfun}{fitfunction to be used in the fitting procedure. Either FML or FIML}

\item{data_type}{type of data in the model. Either "cov" or "raw"}

\item{pen_value}{numeric value of penalty size}

\item{pen_on}{string vector with matrices that should be regularized. Possible are combinations of "A", "S", "DRIFT"}

\item{selectedDrifts}{drift values to regularize. Possible are "all", "cross", "auto" or providing a matrix of the same size as the drift matrix with ones for every parameter to regularize and 0 for every non-regularized parameter}

\item{driftexpo}{specifiy if the regularization will be performed on the raw drift matrix or on the exponential of the drift matrix (discrete time parameters)}

\item{DRIFT_dt}{provide the discrete time points for which the drift will be regularized. A vector with multiple values is possible}

\item{selectedA}{A values to regularize. Possible are "all", or providing a matrix of the same size as the A matrix with ones for every parameter to regularize and 0 for every non-regularized parameter}

\item{selectedS}{S values to regularize. Possible are "all", or providing a matrix of the same size as the S matrix with ones for every parameter to regularize and 0 for every non-regularized parameter}

\item{penalty_type}{so far only "lasso" implemented}
}
\description{
Note: laremm is based on the R package \pkg{regsem}. Because of the early status of laremm, it is recommended to use regsem instead!
createRegModel creates a regularized model from a mxModel or ctsem.
}
\examples{
# The following example is taken from the regsem help to demonstrate the equivalence of both methods:

library(lavaan)
library(OpenMx)
# put variables on same scale for regsem
HS <- data.frame(scale(HolzingerSwineford1939[,7:15]))

# define variables:
latent = c("f1")
manifest = c("x1","x2","x3","x4","x5", "x6", "x7", "x8", "x9")

# define paths:
loadings <- mxPath(from = latent, to = manifest, free = c(F,T,T,T,T,T,T,T,T), values = 1)
lcov <- mxPath(from = latent, arrows = 2, free = T, values = 1)
lmanif <- mxPath(from = manifest, arrows =2 , free =T, values = 1)

# define model:
myModel <- mxModel(name = "myModel", latentVars = latent, manifestVars = manifest, type = "RAM",
                   mxData(observed = HS, type = "raw"), loadings, lcov, lmanif,
                   mxPath(from = "one", to = manifest, free = T)
)

fit_myModel <- mxRun(myModel)
round(fit_myModel$A$values,5)

# create regularized model:

selectedA <- matrix(0, ncol = ncol(fit_myModel$A$values), nrow = nrow(fit_myModel$A$values))
selectedA[c(2,3,7,8,9),10] <-1


reg_model <- createRegModel(model = fit_myModel, model_type = "mxModel", fitfun = "FIML", data_type = "raw",
                            pen_on = "A", selectedA = selectedA, pen_value = .05
)
fit_reg_model <- mxRun(reg_model)

round(fit_reg_model$BaseModel$A$values,5)

}
\author{
Jannik Orzek
}
