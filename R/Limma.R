################################################################################
# Using limma/lmFit to perform EWAS
# Continuous phenotype/exposure
################################################################################

################################################################################
#' Run linear regression using limma: CpG ~ Exposure + Confounders
#' Warning: does not allow missing data in confounders/exposures
#' 
#'
#' @param meth_data A numeric matrix containing methylation data (CpGs as rows, samples as columns)
#' @param covariates A data.frame containing all variables needed for regression (exposures and confounders)
#' @param exposure A character string with the exposure name (must be a column in covariates)
#' @param clinical_confounders A character vector defining clinical confounders the regression will be adjusted for (must be columns in covariates)
#' @param technical_confounders A character vector defining technical confounders the regression will be adjusted for (must be columns in covariates)
#' @param transformToMvalue Boolean. Whether to transform CpG values into M-values. Default: FALSE.
#' @param robust Boolean. Whether to perform robust linear regression. Default: FALSE.
#' @param maxit The number of iterations for each robust regression.
#' @param verbose Boolean. Whether to output messages. Default: TRUE.
#' @param id.lab String. Column name were to find sample ids. Default: "id".
#'
#' @import limma
#' @import magrittr
#' 
#' @importFrom stats as.formula model.matrix
#'
#' @return An object of class MArrayLM as output by the limma::eBayes() function.

runLimma <- function(meth_data, 
                     covariates,
                     exposure,
                     clinical_confounders = 0,
                     technical_confounders = 0,
                     transformToMvalue = FALSE,
                     robust = FALSE,
                     verbose = TRUE,
                     id.lab = "id"){
  
  if(!(id.lab %in% colnames(covariates))) stop(paste("Input covariate matrix must have a column labelled", id.lab, ".\n"))
  
  if( !identical(covariates[,id.lab], colnames(meth_data)) ) stop("Covariates and methylation data do not match \n")
  
  if(!(exposure %in% colnames(covariates))) stop("Specified exposure not in input covariate matrix. \n")
  
  covariates[,exposure] <- as.numeric(covariates[,exposure])
  
  nIndiv <- nrow(covariates)
  
  covariates <- na.omit(covariates[, unique(c(id.lab, exposure, clinical_confounders, technical_confounders))])
  
  if( nrow(covariates) != nIndiv ){
    warning(paste(nIndiv-nrow(covariates), 
                  "individual(s) with missing covariate data had to be removed from the analysis."))
    meth_data <- meth_data[,intersect(colnames(meth_data), covariates[,id.lab])]
  }
  
  # Design specification
  formula <-
    stats::as.formula(paste("~", exposure, "+",
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
  
  # Fitting the linear model
  if( verbose ) cat("[LIMMA] Fitting the linear model...")
  fit <- limma::lmFit(object = meth_data, 
                      design = design,
                      method = ifelse(robust, "robust", "ls"), 
                      maxit = maxit)
  #apply eBayes
  fit <- limma::eBayes(fit)
  
  if( verbose ) cat("OK.\n")
  
  return( fit )

}