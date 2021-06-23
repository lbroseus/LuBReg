################################################################################
# Functions for fitting multivariate robust linear models (RLM) 
# And significance testing
# With binary phenotype (eg: case/control) and paired design
# DNAme ~ Status + Confounders
################################################################################

################################################################################
#' Run differential analysis using limma: CpG ~ Treatment + Confounders
#' Warning: does not allow missing data in confounders
#' 
#' @param meth_data A numeric matrix containing methylation data (CpGs as rows, samples as columns)
#' @param covariates A data.frame containing all variables needed for regression (exposures and confounders)
#' @param treatment A character string with the status name (must be a column in covariates)
#' @param clinical_confounders A character vector defining clinical confounders the regression will be adjusted for (must be columns in covariates)
#' @param technical_confounders A character vector defining technical confounders the regression will be adjusted for (must be columns in covariates)
#' @param transformToMvalue Boolean. Whether to transform CpG values into M-values. Default: FALSE.
#' @param paired Boolean. Whether samples are paired. Default: FALSE.
#' @param pair.lab String. Column name were to find pairing. Default: "pair".
#' @param robust Boolean. Whether to perform robust linear regression. Default: FALSE.
#' @param maxit The number of iterations for each robust regression.
#' @param verbose Boolean. Whether to output messages. Default: TRUE.
#'
#' @import limma
#' @import magrittr
#' 
#' @importFrom stats as.formula model.matrix
#'
#' @return An object of class MArrayLM as output by the limma::eBayes() function.
#' 
casectrlAnalysis <- function(meth_data,
                             covariates,
                             treatment,
                             clinical_confounders = 0,
                             technical_confounders = 0,
                             transformToMvalue = FALSE,
                             paired = FALSE,
                             pair.lab = "pair",
                             robust = FALSE,
                             verbose = TRUE){
  
  
  # Design specification
  formula <-
    stats::as.formula(paste("~", treatment, "+",
                            paste(clinical_confounders, collapse = " + "), "+",
                            paste(technical_confounders, collapse = " + ")
    )
    )
  
  if( verbose ) cat("[LIMMA] Linear model to be fitted: \n"); print(formula)
  
  design <- stats::model.matrix(object = formula, data = covariates)
  
  # Transform beta to M-values
  if( transformToMvalue ){
    if( verbose ) cat("[LIMMA] Transforming input methylation data into M-values...")
    meth_data <- apply(X = meth_data, MARGIN = 1:2, FUN = logit2)
    if( verbose ) cat("OK.\n")
  }
  
  if( paired ) cor <- limma::duplicateCorrelation(meth_data, 
                                     design = design, 
                                     block = covariates[,pair.lab])
  
  if( verbose ) cat("[LIMMA] Fitting linear models...")
  fit <- limma::lmFit(object = meth_data,
                      design = design,
                      block = ifelse(paired, covariates[,pair.lab], NULL),
                      correlation = ifelse(paired, cor$consensus, NULL))
  # Apply eBayes
  fit <- eBayes(fit)
  if( verbose ) cat("OK.\n")
  
  return( fit )
}
