################################################################################
#' Run Mixed Linear regression for methylation at specific loci ~ 1 exposure at a time (ExWAS) adjusted 
#'
#' @param CpG A numeric vector containing methylation data from one CpG
#' @param covariates A data.frame containing all variables needed for regression (exposures and confounders)
#' @param exposure A character string with the exposure name
#' @param clinical_confounders A character vector containing names of the confounders
#' @param technical_confounders A character vector defining technical confounders the regression will be adjusted for
#' @param transformToMvalue Boolean: whether input data should be transformed to M-value
#' @param robust Boolean. Whether to perform robust linear regressions in the second stage. Default: TRUE.
#' @param maxit The number of iterations for each robust regression

#' @importFrom stats as.formula p.adjust sd
#' @importFrom lmerTest lmer 
#' @importFrom lme4 ranef
#' @importFrom survey regTermTest
#' @importFrom dplyr filter distinct
#'
#' @return A named vector of p values for regressions for each CpG

twoStageMLM <- function(CpG,
                        covariates,
                        exposure,
                        clinical_confounders = 0,
                        technical_confounders = 0,
                        transformToMvalue = transformToMvalue,
                        robust = TRUE,
                        maxit = 25){
    
    #--------------------------------------------------------------------------#
    # Preparing data
    #--------------------------------------------------------------------------#
  
    ## If required, transform beta- to M-values
    if( transformToMvalue & max(CpG, na.rm = T)<1 & min(CpG, na.rm = T)>0){
      CpG <- logit2(CpG)
    }
    ## Keeping CpG as a matrix will preserve the row names 
    ## (apply coerces it to numeric vector)
    CpG <- as.matrix(CpG)
    
    ## Change the CpG colname to "y" as it will be used as such in the regression formula
    colnames(CpG) <- "y"
    
    ## Change ID to rownames so that covariates and methylation datasets can be merged
    CpG <- data.frame(id = rownames(CpG), CpG)
    
    if(!("id" %in% colnames(covariates)))  covariates <- data.frame(id = rownames(covariates), covariates)
    
    #--------------------------------------------------------------------------#
    # Two-stage model
    #--------------------------------------------------------------------------#
    ## Stage one: mixed model 
    ## (assuming 'constant' exposure, repeatedly measured)
    
    ### ML regression formula: exposure ~ 1 + (1 | id)
    formula.repeat <- stats::as.formula(paste(exposure, " ~ (1 | id)"))
    ### Fit mixed linear model
    fit <- lmerTest::lmer(formula = formula.repeat, data = covariates)
    ### Gather fixed and random effects
    fixed <- as.vector(summary(fit)$coefficients["(Intercept)","Estimate"])
    random <- lme4::ranef(fit)$id
    ### Compute individual coefficient (pop fixed + indiv random)
    exposure.estim <- data.frame(id = rownames(random),
                                 exposure = random$`(Intercept)`+fixed)
    
    #--------------------------------------------------------------------------#
    ## Stage two: multi-linear regression: 
    ## outcome ~ (exposure.fixed + exposure.random) + confounders
    
    ### Create the new data set of covariates
    covariates <- covariates %>% dplyr::distinct(id, .keep_all = TRUE)
    covariates <- merge(x = covariates, y = CpG, by = "id")
    covariates <- dplyr::filter(covariates, !is.na(y))
    covariates <- merge(covariates, exposure.estim, by ='id')
    
    ### Create regression formula
    formula.ewas <-
      stats::as.formula(paste(
        "y", " ~ exposure + ",
        paste(clinical_confounders, collapse = " + "),
        "+",
        paste(technical_confounders, collapse = " + ")
      ))
    
    ### Fit (robust) linear regression   
    if( robust ){
      
      fit <- MASS::rlm(formula = formula.ewas, data = covariates, maxit = maxit)
      
      Intercept <- fit$coefficients["(Intercept)"]
      
      Estimate <- fit$coefficients["exposure"] %>% as.numeric()
      
    }else{
      
      fit <- lm(formula = formula.ewas, data = covariates)
      
      Intercept <- summary(fit)$coefficients["(Intercept)", "Estimate"]
      
      Estimate <- summary(fit)$coefficients["exposure", "Estimate"] %>% as.numeric()
      
    }
    
    #--------------------------------------------------------------------------#
    ## Compute estimates
    
    #### Calculate p-value using Wald test 
    raw_p_value <- survey::regTermTest(model = fit,
                                       test.terms = "exposure",
                                       null = NULL,
                                       df = Inf,
                                       method = "Wald")$p
    
    ### Calculate 95% CIs
    suppressMessages(
      CIs <- stats::confint.default(object = fit, 
                            parm = "exposure", level = 0.95) %>%
        as.data.frame() %>%
        dplyr::rename(conf_low = `2.5 %`, conf_high = `97.5 %`)
    )
    
    ### Standard errors
    SE <- summary(fit)$coefficients["exposure", "Std. Error"] %>% as.numeric()
    
    sfit <- cbind(Estimate, SE, CIs, raw_p_value = raw_p_value)
    
    if( transformToMvalue ){
      # Return adjusted regression coefficient in terms of beta values
      MeanBeta <- as.numeric(m2beta(Intercept+Estimate)-m2beta(Intercept))
      sfit <- cbind(MeanBeta, sfit)
    }
    
    return(sfit)
  }

################################################################################
#' Results of the mixed linear regressions for mean methylation level in function of repeated exposure measurements
#'
#' @param meth_data A matrix containing methylation data (samples in columns).
#' @param exposures A character vector naming the exposures.
#' @param covariates A data.frame containing all variables needed for regression (exposures and confounders).
#' @param clinical_confounders A character vector containing names of the confounders.
#' @param technical_confounders A character vector defining technical confounders the regression will be adjusted for.
#' @param transformToMvalue Boolean. Whether input data should be transformed to M-value
#' @param robust Boolean. Whether to perform robust linear regressions in the second stage. Default: TRUE.
#' @param maxit The number of iterations for each robust regression.
#' @param path Path for saving the result file.
#' @param file_name File name without extension.
#' @param ncores The number of cores used for parallel computing, by default half of all available cores.
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

LocusWiseTwoStageMLM <- function(meth_data,
                                 covariates,
                                 exposures,
                                 clinical_confounders,
                                 technical_confounders,
                                 transformToMvalue = FALSE,
                                 robust = TRUE,
                                 maxit,
                                 ncores = ceiling(bigstatsr::nb_cores()/2),
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
    
    # 1. Run models for each CpG in function of exposure, adjusted for covariates
    # 2. Create an object where regression results for all exposures will be stored
    
    res_all_expo <- list()
    
    for (i in seq_along(exposures)){
      
      cat("Progress:", i, "/", length(exposures), "exposures\n")
      
      # Run models separately for each CpG and adjust pvalues for multiple testing
      result_loci_regr <- parallel::parApply(
        cl = cl,
        X = meth_data,
        MARGIN = 1,
        FUN = twoStageMLM,
        exposure = exposures[i],
        covariates = covariates,
        clinical_confounders = clinical_confounders,
        technical_confounders = technical_confounders,
        transformToMvalue = transformToMvalue,
        robust = robust
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
                       file = here::here(path, paste0(file_name, "_clinical_confounders.txt")))
    data.table::fwrite(list(technical_confounders),
                       file = here::here(path, paste0(file_name, "_technical_confounders.txt")))
    
    return(res_all_expo)
  }
################################################################################