################################################################################
#' Run Mixed Linear regression for methylation at a specific locus and repeated
#' exposure. Likelihood ratio test
#'
#' @param covariates A data.frame containing all variables needed for regression (exposures and confounders)
#' @param exposure A character string with the exposure name
#' @param clinical_confounders A character vector defining clinical confounders the regression will be adjusted to
#' @param technical_confounders A character vector defining technical confounders the regression will be adjusted to
#' @param CpG A numeric vector containing methylation data from one CpG position.
#'
#' @importFrom stats as.formula sd
#' @importFrom lmerTest lmer
#' @importFrom dplyr filter
#'
#' @return A named vector of p values for regressions for each CpG

myMLRegressionLoci.LRT <-
  function(CpG,
           exposure,
           covariates,
           clinical_confounders = 0,
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
    data <- dplyr::filter(data, !is.na(y))
    
    formula <-
      stats::as.formula(paste(
        exposure, " ~ ",
        " y ",
        "+",
        paste(clinical_confounders, collapse = " + "),
        "+",
        paste(technical_confounders, collapse = " + "),
        " + (1 | id)"
      ))
    fit <- lmerTest::lmer(formula = formula1, data = data)
    
    # Calculate CIs
    suppressMessages(
      CIs <- stats::confint(object = fit, parm = "y", level = 0.95) %>%
        as.data.frame() %>%
        dplyr::rename(conf_low = `2.5 %`, conf_high = `97.5 %`)
    )
    
    # Obtain estimate
    Estimate <- summary(fit)$coefficients["y", "Estimate"]
    
    # Obtain standard error of the estimate
    SE <- summary(fit)$coefficients["y", "Std. Error"]
    
    # pval using LRT+anova
    pval <- anova(fit)["y","Pr(>F)"]
    
    # Combine estimates, CIs and p values
    sfit <-
      cbind(
        Estimate,
        SE,
        CIs,
        raw_p_value = pval
      )
    
    if( transformToMvalue ){
      # Return adjusted regression coefficient in terms of beta values
      Intercept <- summary(fit1)$coefficients["(Intercept)", "Estimate"]
      MeanBeta <- as.numeric(m2beta(Intercept+Estimate)-m2beta(Intercept))
      sfit <- cbind(MeanBeta, sfit)
    }
    
    return(sfit)
  }
################################################################################

################################################################################
#' Run Mixed Linear regression for methylation at a specific locus and repeated
#' exposure. Wald test.
#'
#' @param covariates A data.frame containing all variables needed for regression (exposures and confounders)
#' @param exposure A character string with the exposure name
#' @param clinical_confounders A character vector defining clinical confounders the regression will be adjusted to
#' @param technical_confounders A character vector defining technical confounders the regression will be adjusted to
#' @param CpG A numeric vector containing methylation data from 1 CpG
#'
#' @importFrom stats as.formula p.adjust sd
#' @importFrom lmerTest lmer
#' @importFrom dplyr filter
#'
#' @return A named vector of p values for regressions for each CpG

myMLRegressionLoci.Wald <-
  function(CpG,
           exposure,
           covariates,
           clinical_confounders = 0,
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
    data <- dplyr::filter(data, !is.na(y))
    
    # Create regression formula
    formula <-
      stats::as.formula(paste(
        exposure, " ~ ",
        " y ",
        "+",
        paste(clinical_confounders, collapse = " + "),
        "+",
        paste(technical_confounders, collapse = " + "),
        " + (1 | id)"
      ))
    fit <- lmerTest::lmer(formula = formula, data = data)
    
    # Calculate CIs
    suppressMessages(
      CIs <- stats::confint(object = fit, parm = "y", level = 0.95) %>%
        as.data.frame() %>%
        dplyr::rename(conf_low = `2.5 %`, conf_high = `97.5 %`)
    )
    
    # Obtain estimates
    Estimate <- summary(fit)$coefficients["y", "Estimate"]
    
    # Obtain standard error of the estimate
    SE <- summary(fit)$coefficients["y", "Std. Error"]
    
    pval <- summary(fit)$coefficients["y", "Pr(>|t|)"]
    
    # Combine estimates, CIs and p values
    sfit <-
      cbind(
        Estimate,
        SE,
        CIs,
        raw_p_value = pval
      )
    
    if( transformToMvalue ){
      # Return adjusted regression coefficient in terms of beta values
      Intercept <- summary(fit)$coefficients["(Intercept)", "Estimate"]
      MeanBeta <- as.numeric(m2beta(Intercept+Estimate)-m2beta(Intercept))
      sfit <- cbind(MeanBeta, sfit)
    }
    
    return(sfit)
  }
################################################################################

################################################################################
#' Results of the mixed linear regressions for mean methylation level in function of repeated exposure measurements
#'
#' @param meth_data A matrix containing methylation data (samples in columns).
#' @param covariates A data.frame containing all variables needed for regression (exposures and confounders).
#' @param exposures A character vector naming the exposures.
#' @param clinical_confounders A character vector containing names of the confounders.
#' @param technical_confounders A character vector defining technical confounders the regression will be adjusted for.
#' @param transformToMvalue Boolean: whether input data should be transformed to Mvalue.
#' @param method Method to be used to test the association: "Wald" test or "LRT". Default: "Wald".
#' @param path Path for saving the result file.
#' @param file_name File name without extension.
#' @param ncores The number of cores used for parallel computing, by default all available cores.
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
#'

LocusWiseLME <-
  function(meth_data,
           covariates,
           exposures,
           clinical_confounders,
           technical_confounders,
           transformToMvalue = FALSE,
           method = "Wald",
           path,
           file_name,
           ncores = 1){
    
    # Check input parameters
    if(!method %in% c("Wald", "LRT")) stop("method must be either 'Wald' or 'LRT'")
    cat("Using", method, "method for testing association. \n")
    
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
      if( method == "Wald"){ 
        
        result_loci_regr <- parallel::parApply(
          cl = cl,
          X = meth_data,
          MARGIN = 1,
          FUN = myMLRegressionLoci.Wald,
          exposure = exposures[i],
          covariates = covariates,
          clinical_confounders = clinical_confounders,
          technical_confounders = technical_confounders,
          transformToMvalue = transformToMvalue)
        
      }else{
        
        result_loci_regr <- parallel::parApply(
          cl = cl,
          X = meth_data,
          MARGIN = 1,
          FUN = myMLRegressionLoci.LRT,
          exposure = exposures[i],
          covariates = covariates,
          clinical_confounders = clinical_confounders,
          technical_confounders = technical_confounders,
          transformToMvalue = transformToMvalue)
      }
      
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
    
    saveRDS(res_all_expo, file = here::here(path, paste0(file_name, ".rds")))
    
    # 4. Save the confounders set and the technical confounders set to a file
    data.table::fwrite(list(clinical_confounders),
                       file = here::here(path, paste0(file_name, "clinical_confounders.txt")))
    data.table::fwrite(list(technical_confounders),
                       file = here::here(path, paste0(file_name, "_technical_confounders.txt")))
    
    return(res_all_expo)
  }
################################################################################