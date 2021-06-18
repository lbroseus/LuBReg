################################################################################
# Functions for running conditional lifecourse models
################################################################################

################################################################################
#' Run robust linear regression for methylation at specific loci ~ exposures
#'
#' @param CpG A numeric vector containing methylation data from 1 CpG
#' @param covariate_data A data.frame containing all variables needed for regression (exposures and confounders)
#' @param exposures A character string with the exposure name(s)
#' @param clinical_confounders A character vector containing names of the confounders
#' @param technical_confounders A character vector defining technical confounders the regression will be adjusted to
#' @param transformToMvalue Boolean. Whether CpG values should be turned to M-values.
#' @param maxit The number of iterations for each robust regression
#'
#' @importFrom stats as.formula p.adjust sd
#' 
#' @export
#'
#' @return A named vector of p values for regressions for each CpG

multiRobustLinearRegression <-
  function(CpG,
           covariate_data,
           exposures,
           clinical_confounders = 0,
           technical_confounders = 0,
           transformToMvalue = FALSE,
           maxit) {
    
    if( transformToMvalue & max(CpG, na.rm = T)<1 & min(CpG, na.rm = T)>0){
      CpG <- logit2(CpG)
    }
    
    # Keeping CpG as a matrix will preserve the row names (apply coerces it to numeric vector)
    CpG <- as.matrix(CpG)
    
    # Change the CpG colname name to "y" as it will be used as such in the regression formula
    colnames(CpG) <- "y"
    
    # Change ID to rownames so covariates can be merged with methylation dataset
    CpG <- data.frame(id = rownames(CpG), CpG)
    
    if(!("id" %in% colnames(covariate_data)))  covariate_data <- data.frame(id = rownames(covariate_data), covariate_data)
    
    # Create a subset of data containing methylation data of one CpG (y) and exposure-covariates data (x)
    data <- merge(x = covariate_data, y = CpG, by = "id")
    
    # Create regression formula
    formula <-
      stats::as.formula(paste(
        "y ~",
        paste(exposures, collapse = " + "),
        "+",
        paste(clinical_confounders, collapse = " + "),
        "+",
        paste(technical_confounders, collapse = " + ")
      ))
    
    # Fit the robust linear regression
    fit <- MASS::rlm(formula = formula, data = data, maxit = maxit)
    
    # Obtain regression statistics
    sfit <- ObtainStatisticsII(fit = fit, 
                               data = data, 
                               exposures = exposures,
                               transformToMvalue = transformToMvalue)
    rownames(sfit) <- NULL
    
    return(sfit)
}

#' Robust linear regression for mean methylation level in function of multiple
#' exposure: methylation ~ exposure1 + exposure2 + exposure3 + W + epsilon
#'
#' @param meth_data A matrix containing methylation data, CpG in rows.
#' @param covariate_data A data.frame containing all variables needed for regression (exposures and confounders)
#' @param clinical_confounders A character vector containing names of the confounders
#' @param technical_confounders A character vector defining technical confounders the regression will be adjusted to
#' @param exposures A character vector naming the exposures
#' @param expo_labels A character vector containing full names of exposures
#' @param transformToMvalue Boolean: whether input data should be transformed to Mvalue (ie: logit transformation)
#' 
#' @param maxit The number of iterations for each robust regression
#' @param path Path for saving the result file
#' @param file_name File name without extension
#' @param ncores The number of cores used for parallel computing, by default all available cores
#'
#' @return A data.frame with annotated statistics for each CpG: beta estimate for the association with exposure; uncorrected and corrected p values; average mean meth level with SE;
#' @export
#' @import dplyr
#' @import here
#' @import parallel
#' @import stringr
#' @importFrom stats p.adjust
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach registerDoSEQ
#' @importFrom tibble rownames_to_column
#' @importFrom data.table fwrite
#' @importFrom bigstatsr nb_cores
#'

LocusWiseMultiRLM <-
  function(meth_data,
           covariate_data,
           exposures,
           clinical_confounders,
           technical_confounders,
           exp_label,
           transformToMvalue = FALSE,
           maxit,
           ncores = bigstatsr::nb_cores(),
           path,
           file_name){

    
    # Create a set of copies of R running in parallel and communicating over sockets
    cl <- parallel::makeCluster(getOption("cl.cores", ncores))
    
    # Register the parallel backend
    doParallel::registerDoParallel(cl)
    
    # Explicitly register a sequential parallel backend. This will prevent a warning message
    # from being issued if the %dopar% function is called and no parallel backend has been registered.
    foreach::registerDoSEQ()
    
    # Stop the cluster
    on.exit(parallel::stopCluster(cl))
    
    # 1. Run robust linear regressions for each CpG in function of exposure(s), adjusted for covariates
    
    # Create an object where regression results will be stored
    
    res_all_expo <- list()
    
    cat("Regression model includes exposures",
        paste(exposures, collapse = " "), "\n")
      
    # Run robust linear regressions separately for each CpG
    # Adjust p values for multiple testing (BH)
    
    result_loci_regr <- parallel::parApply(
        cl = cl,
        X = meth_data,
        MARGIN = 1,
        FUN = multiRobustLinearRegression,
        covariate_data = covariate_data,
        exposures = exposures,
        clinical_confounders = clinical_confounders,
        technical_confounders = technical_confounders,
        transformToMvalue = transformToMvalue,
        maxit = maxit
      )
      
    # Transform results for an exposure into a data frame
    result_loci_regr_df <- result_loci_regr %>%
      do.call(rbind, .) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("CpG")
    
    cat("Saving results...")
    saveRDS(result_loci_regr_df, file = here::here(path, paste0(file_name, ".", exp_label, ".rds")))
    
    # 4.Save the confounders set and the technical confounders set to a file
    data.table::fwrite(list(clinical_confounders),
                       file = here::here(path, paste0(file_name, 
                                                      "_clinical_confounders", 
                                                      ".", exp_label,
                                                      ".txt")))
    data.table::fwrite(list(technical_confounders),
                       file = here::here(path, paste0(file_name, 
                                                      "_technical_confounders",
                                                      ".", exp_label,
                                                      ".txt")))
    cat("OK. \n")
    
    return(result_loci_regr_df)
  }