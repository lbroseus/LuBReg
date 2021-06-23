################################################################################
# Functions for fitting multivariate robust linear models (RLM)
# DNAme 
################################################################################
#' Run robust linear regression: CpG ~ Exposure + Confounders
#'
#' @param CpG A numeric vector containing methylation data from 1 CpG
#' @param covariates A data.frame containing all variables needed for regression (exposures and confounders)
#' @param exposure A character string with the exposure name (must be a column in covariates)
#' @param clinical_confounders A character vector defining clinical confounders the regression will be adjusted for (must be columns in covariates)
#' @param technical_confounders A character vector defining technical confounders the regression will be adjusted for (must be columns in covariates)
#' @param maxit The number of iterations for each robust regression
#' @param transformToMvalue Boolean. Whether to transform CpG values into M-values. Default: FALSE.
#'
#' @importFrom stats as.formula p.adjust sd
#' @importFrom tibble column_to_rownames
#'
#' @return A named vector of p values for regressions for each CpG

RLM <- function(CpG, 
                covariates,
                exposure,
                clinical_confounders = 0,
                technical_confounders = 0,
                maxit, 
                transformToMvalue = FALSE){
  
  # If required, transform beta- to M-values
  if( transformToMvalue & max(CpG, na.rm = T)<1 & min(CpG, na.rm = T)>0){
    CpG <- logit2(CpG)
  }
  # Keep CpG as a matrix will preserve the row names (apply coerces it to numeric vector)
  CpG <- as.matrix(CpG)
  
  # Change the CpG colname name to "y"; will be used as such in the regression formula
  colnames(CpG) <- "y"
  
  # Change ID to rownames so covariates can be merged with methylation dataset
  CpG <- data.frame(id = rownames(CpG), CpG)
  
  if(!("id" %in% colnames(covariates))){
    covariates <- data.frame(id = rownames(covariates), covariates)
  }
  
  # Create a subset of data containing methylation data of one CpG (y) 
  # and covariate data (x)
  data <- merge(x = covariates, y = CpG, by = "id")
  
  # Create regression formula
  formula <-
    stats::as.formula(paste(
      "y ~",
      exposure,
      "+",
      paste(clinical_confounders, collapse = " + "),
      "+",
      paste(technical_confounders, collapse = " + ")
    ))
  
  # Fit the robust linear regression
  fit <- MASS::rlm(formula = formula, data = data, maxit = maxit)
  
  #--------------------------------------------------------------------------#
  # Obtain regression statistics
  #--------------------------------------------------------------------------#
  
  # Calculate p values using Wald test
  wald <- survey::regTermTest(model = fit,
                              test.terms = exposure,
                              null = NULL,
                              df = Inf,
                              method = "Wald")
  # Calculate CIs
  CIs <- stats::confint.default(object = fit, parm = exposure, level = 0.95) %>%
    as.data.frame() %>%
    dplyr::rename(conf_low = `2.5 %`, conf_high = `97.5 %`)
  
  # Obtain standard error of the estimate
  SE <- summary(fit)$coefficients[exposure, "Std. Error"]
  
  # Obtain estimate
  Estimate <- fit$coefficients[exposure]
  
  sfit <- cbind(Estimate, SE, CIs, raw_p_value = as.numeric(wald$p))
  
  if( transformToMvalue ){
    
    Intercept <- fit$coefficients["(Intercept)"]
    
    MeanBeta <- m2beta(Intercept+Estimate) - m2beta(Intercept)
    
    # Combine estimates, CIs and p values
    sfit <- cbind(MeanBeta, sfit)
  }
  
  return(sfit)
}
################################################################################
#' Results of the RLMs for mean methylation level in function of exposures.
#' One exposure at a time.
#'
#' @param meth_data A matrix containing methylation data (samples in columns).
#' @param covariates A data.frame containing all variables needed for regression (exposures and confounders)
#' @param exposures A character vector naming the exposures
#' @param clinical_confounders A character vector containing names of the confounders
#' @param technical_confounders A character vector defining technical confounders the regression will be adjusted to
#' @param maxit The number of iterations for each robust regression
#' @param path Path for saving the result file
#' @param file_name File name without extension
#' @param ncores The number of cores used for parallel computing, by default all available cores
#' @param transformToMvalue Boolean: whether input data should be transformed to Mvalue
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

LocusWiseRLM <-
  function(meth_data,
           covariates,
           exposures,
           clinical_confounders,
           technical_confounders,
           maxit,
           transformToMvalue = FALSE,
           ncores = bigstatsr::nb_cores(),
           path,
           file_name) {
    # Set the parallel computing conditions
    
    # Create a set of copies of R running in parallel and communicating over sockets
    cl <- parallel::makeCluster(getOption("cl.cores", ncores))
    
    # Register the parallel backend
    doParallel::registerDoParallel(cl)
    
    # Explicitly register a sequential parallel backend. This will prevent a warning message
    # from being issued if the %dopar% function is called and no parallel backend has been registered.
    foreach::registerDoSEQ()
    
    # Stop the cluster
    on.exit(parallel::stopCluster(cl))
    
    # 1. Run robust linear regressions for each CpG in function of exposure, adjusted for covariates
    
    # Create an object where regression results for all exposures will be stored
    
    res_all_expo <- list()
    
    for (i in seq_along(exposures)){
      
      cat("Progress:", i, "/", length(exposures), "exposures\n")
      
      # Run robust linear regressions separately for each CpG and adjust p values for multiple testing
      result_loci_regr <- parallel::parApply(
        cl = cl,
        X = meth_data,
        MARGIN = 1,
        FUN = RLM,
        exposure = exposures[i],
        covariates = covariates,
        clinical_confounders = clinical_confounders,
        technical_confounders = technical_confounders,
        maxit = maxit, 
        transformToMvalue = transformToMvalue
      )
      
      # Transform results for an exposure into a data frame
      result_loci_regr_df <- result_loci_regr %>%
        do.call(rbind, .) %>%
        as.data.frame() %>%
        tibble::rownames_to_column("CpG") %>%
        dplyr::mutate(Exposure = exposures[i])
      
      # Adjust p values
      result_loci_regr_table <- adjustPvalues(result_loci_regr_df = result_loci_regr_df,
                                              method = "BH")
      
      # Append the result list
      res_all_expo[[i]] <- result_loci_regr_table
      names(res_all_expo)[i] <- exposures[i]
    }
    
    #Shape results:
    res_all_expo <- do.call(rbind.data.frame, res_all_expo)
    rownames(res_all_expo) <- NULL
    
    cat("Saving results...")
    saveRDS(res_all_expo, file = here::here(path, paste0(file_name, ".rds")))
    
    # 4. Save the confounders set and the technical confounders set to a file
    data.table::fwrite(list(clinical_confounders),
                       file = here::here(path, paste0(file_name, "_clinical_confounders.txt")))
    data.table::fwrite(list(technical_confounders),
                       file = here::here(path, paste0(file_name, "_technical_confounders.txt")))
    cat("OK. \n")
    
    return(res_all_expo)
  }