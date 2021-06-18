
################################################################################
#' Results of the linear regressions for mean methylation level in function of exposure 
#'
#' @param meth_data A matrix containing methylation data
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
        MARGIN = 2,
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

################################################################################
#' Results of the mixed linear regressions for mean methylation level in function of repeated exposure measurements
#'
#' @param covariates A data.frame containing all variables needed for regression (exposures and confounders)
#' @param clinical_confounders A character vector containing names of the confounders
#' @param technical_confounders A character vector defining technical confounders the regression will be adjusted to
#' @param path Path for saving the result file
#' @param file_name File name without extension
#' @param meth_data A matrix containing methylation data
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

LocusWiseLME <-
  function(meth_data,
           covariates,
           exposures,
           clinical_confounders,
           technical_confounders,
           transformToMvalue = FALSE,
           ncores = 1,
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
        FUN = myMLRegressionLoci,
        exposure = exposures[i],
        covariates = covariates,
        clinical_confounders = clinical_confounders,
        technical_confounders = technical_confounders,
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
    
    saveRDS(res_all_expo, file = here::here(path, paste0(file_name, ".rds")))
    
    # 4. Save the confounders set and the technical confounders set to a file
    data.table::fwrite(list(clinical_confounders),
                       file = here::here(path, paste0(file_name, "clinical_confounders.txt")))
    data.table::fwrite(list(technical_confounders),
                       file = here::here(path, paste0(file_name, "_technical_confounders.txt")))
    
    return(res_all_expo)
  }



