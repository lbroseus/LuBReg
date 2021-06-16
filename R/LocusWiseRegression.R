# Adapted from Paulina J.'s
# quiets concerns of R CMD check re: the .'s that appear in pipelines
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(".", "var", "%dopar%"))
}

################################################################################
#' Results of the robust linear regressions for mean methylation level in function of exposures
#' One exposure at a time.
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

LocusWiseRLM <-
  function(meth_data,
           expo_cov,
           exposures,
           confounders,
           technical_confounders,
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
          FUN = myRobustLinearRegressionLoci,
          exposure = exposures[i],
          expo_cov = expo_cov,
          confounders = confounders,
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
    data.table::fwrite(list(confounders),
                       file = here::here(path, paste0(file_name, "_confounders.txt")))
    data.table::fwrite(list(technical_confounders),
                       file = here::here(path, paste0(file_name, "_technical_confounders.txt")))
    cat("OK. \n")
    
    return(res_all_expo)
}

################################################################################
#' Results of the linear regressions for mean methylation level in function of exposure 
#'
#' @param expo_cov A data.frame containing all variables needed for regression (exposures and confounders)
#' @param confounders A character vector containing names of the confounders
#' @param technical_confounders A character vector defining technical confounders the regression will be adjusted to
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

LocusWiseLM <-
  function(meth_data,
           expo_cov,
           exposures,
           confounders,
           technical_confounders,
           transformToMvalue = FALSE,
           ncores = bigstatsr::nb_cores(),
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
        FUN = myLinearRegressionLoci,
        exposure = exposures[i],
        expo_cov = expo_cov,
        confounders = confounders,
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
    data.table::fwrite(list(confounders),
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

LocusWiseLME <-
  function(meth_data,
           expo_cov,
           exposures,
           confounders,
           technical_confounders,
           maxit,
           transformToMvalue = FALSE,
           ncores = bigstatsr::nb_cores(),
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
        FUN = myMLRegressionLoci ,
        exposure = exposures[i],
        expo_cov = expo_cov,
        confounders = confounders,
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
    data.table::fwrite(list(confounders),
                       file = here::here(path, paste0(file_name, "_confounders.txt")))
    data.table::fwrite(list(technical_confounders),
                       file = here::here(path, paste0(file_name, "_technical_confounders.txt")))
    
    return(res_all_expo)
  }


################################################################################
#' Results of the mixed linear regressions for mean methylation level in function of repeated exposure measurements
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

LocusWiseTwoStageMLM <-
  function(meth_data,
           expo_cov,
           exposures,
           confounders,
           technical_confounders,
           maxit,
           transformToMvalue = FALSE,
           ncores = bigstatsr::nb_cores(),
           expo_labels,
           path,
           file_name) {

    # Create a set of copies of R running in parallel and communicating over sockets
    cl <- parallel::makeCluster(getOption("cl.cores", ncores))
    
    # Register the parallel backend
    doParallel::registerDoParallel(cl)
    
    # Explicitly register a sequential parallel backend. This will prevent a warning message
    # from being issued if the %dopar% function is called and no parallel backend has been registered.
    foreach::registerDoSEQ()
    
    # Stop the cluster
    on.exit(parallel::stopCluster(cl))
    
    # 1. Run models for each CpG in function of exposure, adjusted for covariates
    # 2. Create an object where regression results for all exposures will be stored
    
    res_all_expo <- list()
    
    for (i in seq_along(exposures)){
      
      cat("Progress:", i, "/", length(exposures), "exposures\n")
      
      # Run models separately for each CpG and adjust pvalues for multiple testing
      result_loci_regr <- parallel::parApply(
        cl = cl,
        X = meth_data,
        MARGIN = 2,
        FUN = twoStageMLM,
        exposure = exposures[i],
        expo_cov = expo_cov,
        confounders = confounders,
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
    data.table::fwrite(list(confounders),
                       file = here::here(path, paste0(file_name, "_confounders.txt")))
    data.table::fwrite(list(technical_confounders),
                       file = here::here(path, paste0(file_name, "_technical_confounders.txt")))
    
    return(res_all_expo)
  }
################################################################################



################################################################################
#' Results of the linear regressions for mean methylation level in function of exposure 
#'
#' @param expo_cov A data.frame containing all variables needed for regression (exposures and confounders)
#' @param confounders A character vector containing names of the confounders
#' @param technical_confounders A character vector defining technical confounders the regression will be adjusted to
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

LocusWiseTimeMLReg <-
  function(meth_data,
           covariates.df,
           exposures,
           time,
           clinical_confounders,
           technical_confounders,
           transformToMvalue = FALSE,
           ncores = bigstatsr::nb_cores(),
           path,
           file_name){
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
      
      # Run ML regressions separately for each CpG and adjust p values for multiple testing
      result_loci_regr <- parallel::parApply(
        cl = cl,
        X = meth_data,
        MARGIN = 2,
        FUN = TimeMLReg,
        covariates.df = covariates.df,
        mainEffect = exposures[i],
        time = time,
        testTimeInteraction = TRUE,
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
      
      res_all_expo[[i]] <- result_loci_regr_df
      
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
