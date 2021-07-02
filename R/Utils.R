################################################################################
#' Transform beta-value to M-value
#'
#' @param p A proportion (eg: betavalue)
#'
#' @return logit2(p)

logit2 <- function(p) return( log2(p/(1-p)) )

################################################################################
#' Transform M-value to beta-value
#'
#' @param m A numeric value (eg: M-value)
#'
#' @return  2^m/(2^m + 1)

m2beta <- function(m) return( 2^m/(2^m + 1) )

################################################################################
#' Obtain regression statistics for one
#'
#'
#' @param fit An object created by MASS:mlr robust linear regression function
#' @param data A data.frame containing a subset of methylation data of one CpG (y) and exposure-covariates data (x)
#' @param exposure A string with the exposure name
#' @import survey
#' @importFrom stats confint.default
#' @importFrom broom tidy
#'
#' @return A data.frame containing results of the regression for one model for one exposure

myObtainStatistics <- function(fit, data, exposure) {
  
  # Calculate p values using Wald test
  wald <-
    survey::regTermTest(model = fit,
                        test.terms = exposure,
                        null = NULL,
                        df = Inf,
                        method = "Wald")
  # Calculate CIs
  CIs <-
    stats::confint.default(object = fit,
                           parm = exposure,
                           level = 0.95) %>%
    as.data.frame() %>%
    dplyr::rename(conf_low = `2.5 %`, conf_high = `97.5 %`)
  
  # Obtain estimate with it 95% CI
  Estimate_CI <- paste0(
    format(round(fit$coefficients[exposure], 4), nsmall = 4, scientific = FALSE),
    " (",
    format(round(CIs$conf_low, 4), nsmall = 4, scientific = FALSE),
    ";",
    format(round(CIs$conf_high, 4), nsmall = 4, scientific = FALSE),
    ")"
  )
  
  # Obtain estimate
  Estimate <- fit$coefficients[exposure]
  
  # Obtain standard error of the estimate
  SE <- summary(fit)$coefficients[exposure, "Std. Error"]
  
  #}
  
  # Combine estimates, CIs and p values
  sfit <-
    cbind(
      Estimate,
      SE,
      CIs,
      Estimate_CI,
      raw_p_value = as.numeric(wald$p)
    )
  
  return(sfit)
}

################################################################################
#' Obtain regression statistics for multiple effects
#'
#' @param fit An object created by MASS:mlr robust linear regression function
#' @param data A data.frame containing a subset of methylation data of one CpG (y) and exposure-covariates data (x)
#' @param exposures A string with the names exposures of interest
#' @import survey
#' @importFrom stats confint.default
#' @importFrom broom tidy
#' 
#' @export
#'
#' @return A data.frame containing results of the regression for one model for one exposure

ObtainStatisticsII <- function(fit, data, exposures, transformToMvalue){
  
  # Calculate p values using Wald test
  raw_p_value.effect <- tapply(exposures, 
                 1:length(exposures),
                 FUN = function(exposure) 
                        as.numeric(survey::regTermTest(model = fit,
                                            test.terms = exposure,
                                            null = NULL,
                                            df = Inf,
                                            method = "Wald")$p)
                 )
  raw_p_val.model <-  as.numeric(survey::regTermTest(model = fit,
                                                     test.terms = exposures,
                                                     null = NULL,
                                                     df = Inf,
                                                     method = "Wald")$p)
  # Calculate CIs
  CIs <-
    stats::confint.default(object = fit,
                           parm = exposures,
                           level = 0.95) %>%
    as.data.frame() %>%
    dplyr::rename(conf_low = `2.5 %`, conf_high = `97.5 %`)
  
  # Obtain estimate with it 95% CI
  Estimate_CI <- paste0(
    format(round(fit$coefficients[exposures], 4), nsmall = 4, scientific = FALSE),
    " (",
    format(round(CIs$conf_low, 4), nsmall = 4, scientific = FALSE),
    ";",
    format(round(CIs$conf_high, 4), nsmall = 4, scientific = FALSE),
    ")"
  )
  
  # Obtain estimates
  Estimates <- fit$coefficients[exposures]
  
  Intercept <- fit$coefficients["(Intercept)"]
  
  if( transformToMvalue ){
    MeanBeta <- tapply(Estimates, 
                     1:length(exposures),
                     function(est) m2beta(Intercept+est)-m2beta(Intercept))
  }
  
  # Obtain standard error of the estimate
  SE <- summary(fit)$coefficients[exposures, "Std. Error"]
  
  # Combine estimates, CIs and p values
  if(transformToMvalue){
    sfit <- cbind(exposures, Estimates, SE, CIs, Estimate_CI,
                  MeanBeta,
                  raw_p_value.effect, raw_p_val.model)
  }else{
    sfit <- cbind(exposures, Estimates, SE, CIs, Estimate_CI,
                  raw_p_value.effect, raw_p_val.model)
  }
  return(sfit)
}

#' Run linear regression for methylation at specific loci ~ 1 exposure at a time (ExWAS) adjusted for diferent set of confounders
#'
#' @param covariates A data.frame containing all variables needed for regression (exposures and confounders)
#' @param exposure A character string with the exposure name
#' @param technical_confounders A character vector defining technical confounders the regression will be adjusted to
#' @param maxit The number of iterations for each robust regression
#' @param CpG A numeric vector containing methylation data from 1 CpG
#' @param confounders A character vector containing names of the confounders
#'
#' @importFrom stats as.formula p.adjust sd
#' @importFrom tibble column_to_rownames
#'
#' @return A named vector of p values for regressions for each CpG

myLinearRegressionLoci <-
  function(CpG,
           exposure,
           covariates,
           confounders = 0,
           technical_confounders = 0,
           transformToMvalue = transformToMvalue) {
    
    if( transformToMvalue & max(CpG, na.rm = T)<1 & min(CpG, na.rm = T)>0){
      CpG <- logit2(CpG)
    }
    # Keeping CpG as a matrix will preserve the row names (apply coerces it to numeric vector)
    CpG <- as.matrix(CpG)
    
    # Change the CpG colname name to "y" as it will be used as such in the regression formula
    colnames(CpG) <- "y"
    
    # Change ID to rownames so covariates can be merged with methylation dataset
    CpG <- data.frame(id = rownames(CpG), CpG)
    
    if(!("id" %in% colnames(covariates)))  covariates <- data.frame(id = rownames(covariates), covariates)
    # Create a subset of data containing methylation data of one CpG (y) and exposure-covariates data (x)
    data <- merge(x = covariates, y = CpG, by = "id")
    
    # Create regression formula
    formula <-
      stats::as.formula(paste(
        "y ~ ",
        exposure,
        "+",
        paste(confounders, collapse = " + ")
      )
      )
    
    # Fit the robust linear regression for one exposure at a time
    fit <- lm(formula = formula, data = data)
    fit <- summary(fit)
    
    # Obtain regression statistics
    sfit <- 
      cbind.data.frame(
        Estimate = fit$coefficients[2,1],
        SE = fit$coefficients[2,2],
        raw_p_value = fit$coefficients[2,4]
      )
    
    return(sfit)
  }


################################################################################
#' Adjust p values for multiple testing
#'
#' @param result_loci_regr_df A data.frame of uncorrected p values for regressions for each CpG
#' @param method Correction method, a character string. Can be abbreviated.
#'
#' @return A data.frame of uncorrected and corrected p values for regressions for each CpG
#' @seealso [p.adjust()] for p values adjustment method
#'

adjustPvalues <- function(result_loci_regr_df, method) {
  
  result_loci_regr_table <- result_loci_regr_df %>%
    dplyr::mutate(
      # Add Benjamini-Hochberg corrected p value
      p_value_FDR = stats::p.adjust(
        raw_p_value,
        method = method
      )
    )
  
  return(result_loci_regr_table)
}

################################################################################
#' Synchronize methylation and meta datasets (eg: by trimester)
#'
#' @param mmList A 2-entry list where mmList[[1]] is a matrix containing full methylation data 
#' and mmList[[2]] is a data.frame containing all covariates
#' @param rm_samples IDs of samples to be removed from analyses (eg: spurious ones, outliers)
#'
#' @importFrom dplyr filter arrange
#' 
#' @export
#' 
#' @return A 2-entry list with synchronized data (ie: same samples, matching order)

synchronizeDatasets <- function(mmList, rm_samples = NULL){
  
  names(mmList) <- c("methylData", "metaData")
  
  #----------------------------------------#
  # Processing of mmList[["methylData"]]
  #----------------------------------------#
  
  # rm specified samples
  if(!is.null(rm_samples) ) mmList[["methylData"]] <- mmList[["methylData"]][, !(colnames(mmList[["methylData"]]) %in% rm_samples)]
  #rm samples not in the meta data set
  mmList[["methylData"]] <- mmList[["methylData"]][,colnames(mmList[["methylData"]]) %in% mmList[["metaData"]]$id]
  #order samples by name
  mmList[["methylData"]] <- mmList[["methylData"]][, order(colnames(mmList[["methylData"]]))]
  
  cat("Methylation dataset has", 
      ncol(mmList[["methylData"]]), "individuals and spans",
      nrow(mmList[["methylData"]]), "sites.\n")
  
  #----------------------------------------#
  # Processing of mmList[["metaData"]]
  #----------------------------------------#
  
  mmList[["metaData"]] <- mmList[["metaData"]] %>% 
    dplyr::filter(id %in% colnames(mmList[["methylData"]])) %>%
    dplyr::arrange(id)
  
  #----------------------------------------#
  
  #if(!identical(mmList[["metaData"]]$id, colnames(mmList[["methylData"]])) ) stop("Datasets do not match..!\n")
  
  #----------------------------------------#
  
  return( mmList )
  
}

################################################################################
#' Results of the linear regressions for mean methylation level in function of exposure 
#'
#' @param meth_data A matrix containing methylation data (samples in columns).
#' @param covariates A data.frame containing all variables needed for regression (exposures and confounders)
#' @param clinical_confounders A character vector containing names of the confounders
#' @param technical_confounders A character vector defining technical confounders the regression will be adjusted to
#' @param path Path for saving the result file
#' @param file_name File name without extension
#' @param exposures A character vector naming the exposures
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

LocusWiseLM <-
  function(meth_data,
           covariates,
           exposures,
           clinical_confounders,
           technical_confounders,
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
        FUN = myLinearRegressionLoci,
        exposure = exposures[i],
        covariates = covariates,
        clinical_confounders = clinical_confounders,
        technical_confounders = technical_confounders,
        transformToMvalue = transformToMvalue
      )
      
      # Transform results for an exposure into a data frame
      result_loci_regr_df <- do.call(rbind,result_loci_regr) %>%
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
    saveRDS(res_all_expo, file = here::here(path, paste0(file_name, ".rds")))
    
    # 4. Save the confounders set and the technical confounders set to a file
    data.table::fwrite(list(clinical_confounders),
                       file = here::here(path, paste0(file_name, "_confounders.txt")))
    if(length(technical_confounders)>0){
      data.table::fwrite(list(technical_confounders),
                         file = here::here(path, paste0(file_name, "_technical_confounders.txt")))
    }
    cat("OK \n")  
    
    return(res_all_expo)
  }

