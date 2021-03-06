% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getM2LogL_cov.R
\name{getM2LogL_cov}
\alias{getM2LogL_cov}
\title{getM2LogL_cov}
\usage{
getM2LogL_cov(ExpCov, ObsCov, NManif, NObs, useBiasedCov = T)
}
\arguments{
\item{ExpCov}{Expected covariance matrix}

\item{ObsCov}{Observed covariance matrix}

\item{NManif}{Number of manifest variables in the data set}

\item{NObs}{Sample size}

\item{useBiasedCov}{if TRUE, the biased covariance is used to compute the -2LogL}
}
\description{
Note: laremm is based on the R package \pkg{regsem}. Because of the early status of laremm, it is recommended to use regsem instead!
}
\details{
This function computes the -2log(Likelihood) based on the formula reported in Bollen 1989, p. 133
}
\examples{

}
\author{
Jannik Orzek
}
