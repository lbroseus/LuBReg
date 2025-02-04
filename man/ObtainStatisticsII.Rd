% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Utils.R
\name{ObtainStatisticsII}
\alias{ObtainStatisticsII}
\title{Obtain regression statistics for multiple effects}
\usage{
ObtainStatisticsII(fit, data, exposures, transformToMvalue)
}
\arguments{
\item{fit}{An object created by MASS:mlr robust linear regression function}

\item{data}{A data.frame containing a subset of methylation data of one CpG (y) and exposure-covariates data (x)}

\item{exposures}{A string with the names exposures of interest}
}
\value{
A data.frame containing results of the regression for one model for one exposure
}
\description{
Obtain regression statistics for multiple effects
}
