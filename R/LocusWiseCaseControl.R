################################################################################
# `Functions for case-control studies
################################################################################

################################################################################
#' Student T-Test
#'
#' @param formula A formula defining the elements of the regression
#' @param data A data.frame containing global methylation variable, exposure variables, confounders and technical confounders
#'
#' @return An object created by lm linear regression function

myTTest <- function(formula, data, paired = FALSE) {
  ttest <- t.test(formula = formula, data = data, paired = paired)
  return(ttest)
}

################################################################################
#' Run T-Test for methylation at specific loci ~ 1 exposure at a time (ExWAS) adjusted 
#'
#' @param expo_cov A data.frame containing all variables needed for regression (exposures and confounders)
#' @param exposure A character string with the exposure name
#' @param technical_confounders A character vector defining technical confounders the regression will be adjusted to
#' @param CpG A numeric vector containing methylation data from 1 CpG
#' @param confounders A character vector containing names of the confounders
#'
#' @importFrom stats as.formula p.adjust sd
#' @importFrom lme4 lmer
#' @importFrom dplyr filter
#'
#' @return A named vector of p values for regressions for each CpG

myTTestLoci <-
  function(CpG,
           exposure,
           expo_cov,
           confounders = 0,
           technical_confounders = 0,
           transformToMvalue = transformToMvalue,
           paired = FALSE) {
    
    if( transformToMvalue & max(CpG, na.rm = T)<1 & min(CpG, na.rm = T)>0){
      CpG <- logit2(CpG)
    }
    # Keeping CpG as a matrix will preserve the row names (apply coerces it to numeric vector)
    CpG <- as.matrix(CpG)
    
    # Change the CpG colname name to "y" as it will be used as such in the regression formula
    colnames(CpG) <- "y"
    
    # Change ID to rownames so covariates can be merged with methylation dataset
    CpG <- data.frame(id = rownames(CpG), CpG)
    
    if(!("id" %in% colnames(expo_cov)))  expo_cov <- data.frame(id = rownames(expo_cov), expo_cov)
    # Create a subset of data containing methylation data of one CpG (y) and exposure-covariates data (x)
    data <- merge(x = expo_cov, y = CpG, by = "id")
    data <- dplyr::filter(data, !is.na(y))
    
    # Create model formula
    formula <-
      stats::as.formula(paste("y ~ ", exposure))
    
    test <- myTTest(formula = formula, data = data, paired = paired)
    
    CIs <- data.frame(conf_inf = test$conf.low[1], conf_high = test$conf.low[2])
    # Obtain estimate
    Estimate <- test$estimate
    
    # Obtain standard error of the estimate
    SE <- test$stderr
    
    pval <- as.numeric(test$pval)
    
    # Combine estimates, CIs and p values
    sfit <-
      cbind(
        Estimate,
        SE,
        CIs,
        raw_p_value = pval
      )
    
    return(sfit)
  }

################################################################################
#' Results of a locus-wise differential analysis between two groups
#'
#' @param expo_cov A data.frame containing all variables needed for regression (exposures and confounders)
#' @param confounders A character vector containing names of the confounders
#' @param technical_confounders A character vector defining technical confounders the regression will be adjusted to
#' @param maxit The number of iterations for each robust regression
#' @param path Path for saving the result file
#' @param file_name File name without extension
#' @param meth_data A matrix containing methylation data
#' @param exposures A character vector naming the exposures
#' @param ncores The number of cores used for parallel computing, by default all available cores
#' @param annotation_object An Illumina object used to annotate files
#' @param expo_labels A character vector containing full names of exposures
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

LocusWiseTTest <-
  function(meth_data,
           expo_cov,
           exposures,
           confounders,
           technical_confounders,
           paired,
           maxit,
           transformToMvalue = FALSE,
           ncores = bigstatsr::nb_cores(),
           annotation_object,
           expo_labels,
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
        MARGIN = 2,
        FUN = myTTestLoci,
        exposure = exposures[i],
        expo_cov = expo_cov,
        confounders = confounders,
        technical_confounders = technical_confounders,
        transformToMvalue = transformToMvalue,
        paired = paired
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
    data.table::fwrite(list(confounders),
                       file = here::here(path, paste0(file_name, "_confounders.txt")))
    data.table::fwrite(list(technical_confounders),
                       file = here::here(path, paste0(file_name, "_technical_confounders.txt")))
    cat("OK. \n")
    
    return(res_all_expo)
  }