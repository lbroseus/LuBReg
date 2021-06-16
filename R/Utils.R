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
#' Fit the robust linear regression
#'
#' @param formula A formula defining the elements of the regression
#' @param data A data.frame containing global methylation variable, exposure variables, confounders and technical confounders
#' @param maxit The number of iterations for each robust regression
#' @importFrom MASS rlm
#'
#' @return An object created by MASS:mlr robust linear regression function

myRobustLinearRegression <- function(formula, data, maxit) {
  fit <- MASS::rlm(formula = formula,
                   data = data,
                   maxit = maxit)
  return(fit)
}

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

ObtainStatisticsII <- function(fit, data, exposures){
  
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
  Estimate <- fit$coefficients[exposures]
  
  # Obtain standard error of the estimate
  SE <- summary(fit)$coefficients[exposures, "Std. Error"]
  
  # Combine estimates, CIs and p values
  sfit <-
    cbind(
      exposures,
      Estimate,
      SE,
      CIs,
      Estimate_CI,
      raw_p_value.effect,
      raw_p_val.model
    )
  
  return(sfit)
}

################################################################################
#' Run robust linear regression for methylation at specific loci ~ 1 exposure at a time (ExWAS) adjusted for diferent set of confounders
#'
#' @param expo_cov A data.frame containing all variables needed for regression (exposures and confounders)
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

myRobustLinearRegressionLoci <-
  function(CpG,
           exposure,
           expo_cov,
           confounders = 0,
           technical_confounders = 0,
           maxit, 
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

    if(!("id" %in% colnames(expo_cov)))  expo_cov <- data.frame(id = rownames(expo_cov), expo_cov)
    # Create a subset of data containing methylation data of one CpG (y) and exposure-covariates data (x)
    data <- merge(x = expo_cov, y = CpG, by = "id")
    
    # Create regression formula
    formula <-
      stats::as.formula(paste(
        "y ~",
        exposure,
        "+",
        paste(confounders, collapse = " + "),
        "+",
        paste(technical_confounders, collapse = " + ")
      ))
    
    # Fit the robust linear regression for one exposure at a time
    fit <- myRobustLinearRegression(formula = formula,
                                   data = data,
                                   maxit = maxit)
    
    # Obtain regression statistics
    sfit <- myObtainStatistics(fit = fit, data = data, exposure = exposure)
    
    return(sfit)
  }

#' Run linear regression for methylation at specific loci ~ 1 exposure at a time (ExWAS) adjusted for diferent set of confounders
#'
#' @param expo_cov A data.frame containing all variables needed for regression (exposures and confounders)
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
           expo_cov,
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
    
    if(!("id" %in% colnames(expo_cov)))  expo_cov <- data.frame(id = rownames(expo_cov), expo_cov)
    # Create a subset of data containing methylation data of one CpG (y) and exposure-covariates data (x)
    data <- merge(x = expo_cov, y = CpG, by = "id")
    
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
#' Run Mixed Linear regression for methylation at specific loci ~ 1 exposure at a time (ExWAS) adjusted 
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

myMLRegressionLoci <-
  function(CpG,
           exposure,
           expo_cov,
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
    
    if(!("id" %in% colnames(expo_cov)))  expo_cov <- data.frame(id = rownames(expo_cov), expo_cov)
    # Create a subset of data containing methylation data of one CpG (y) and exposure-covariates data (x)
    data <- merge(x = expo_cov, y = CpG, by = "id")
    data <- dplyr::filter(data, !is.na(y))
    
    # Create regression formula
    formula0 <-
      stats::as.formula(paste(
        exposure, " ~ ",
        paste(confounders, collapse = " + "),
        "+",
        paste(technical_confounders, collapse = " + "),
        " + (1 | id)"
      ))
    fit0 <- lme4::lmer(formula = formula0, data = data, REML = FALSE)
    
    formula1 <-
      stats::as.formula(paste(
        exposure, " ~ ",
        " y ",
        "+",
        paste(confounders, collapse = " + "),
        "+",
        paste(technical_confounders, collapse = " + "),
        " + (1 | id)"
      ))
    fit1 <- lme4::lmer(formula = formula1, data = data, REML = FALSE)
    
    # Calculate CIs
    suppressMessages(
      CIs <- stats::confint(object = fit1, parm = "y", level = 0.95) %>%
      as.data.frame() %>%
      dplyr::rename(conf_low = `2.5 %`, conf_high = `97.5 %`)
    )
    
    # Obtain estimate
    Estimate <- summary(fit1)$coefficients["y", "Estimate"]
    
    # Obtain standard error of the estimate
    SE <- summary(fit1)$coefficients["y", "Std. Error"]
    
    LRT <- anova(fit0,fit1)
    pval <- as.numeric(LRT$`Pr(>Chisq)`[2])
    
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

################################################################################
#' Run Mixed Linear regression for methylation at specific loci ~ 1 exposure at a time (ExWAS) adjusted 
#'
#' @param expo_cov A data.frame containing all variables needed for regression (exposures and confounders)
#' @param exposure A character string with the exposure name
#' @param technical_confounders A character vector defining technical confounders the regression will be adjusted to
#' @param CpG A numeric vector containing methylation data from 1 CpG
#' @param confounders A character vector containing names of the confounders
#'
#' @importFrom stats as.formula p.adjust sd
#' @importFrom lme4 lmer ranef
#' @importFrom survey regTermTest
#' @importFrom dplyr filter distinct
#'
#' @return A named vector of p values for regressions for each CpG

twoStageMLM <-
  function(CpG,
           exposure,
           expo_cov,
           confounders = 0,
           technical_confounders = 0,
           transformToMvalue = transformToMvalue,
           maxit = 25){
    
    #--------------------------------------------------------------------------#
    # Prepare data
    
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
    
    if(!("id" %in% colnames(expo_cov)))  expo_cov <- data.frame(id = rownames(expo_cov), expo_cov)

    #--------------------------------------------------------------------------#
    # Two-stage model
    #--------------------------------------------------------------------------#
    ## Stage one: mixed model 
    ## (assuming 'constant' exposure, repeatedly measured)
    
    ### ML regression formula: exposure ~ 1 + (1 | id)
    formula.repeat <- stats::as.formula(paste(exposure, " ~ (1 | id)"))
    ### Fit mixed linear model
    fit <- lme4::lmer(formula = formula.repeat, data = expo_cov)
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
    expo_cov <- expo_cov %>% dplyr::distinct(id, .keep_all = TRUE)
    expo_cov <- merge(x = expo_cov, y = CpG, by = "id")
    expo_cov <- dplyr::filter(expo_cov, !is.na(y))
    expo_cov <- merge(expo_cov, exposure.estim, by ='id')
    
    ### Create regression formula
    formula.ewas <-
      stats::as.formula(paste(
        "y", " ~ exposure + ",
        paste(confounders, collapse = " + "),
        "+",
        paste(technical_confounders, collapse = " + ")
        ))
    
    ### Fit linear regression   
    fit <- lm(formula = formula.ewas, data = expo_cov)
    
    #--------------------------------------------------------------------------#
    ## Compute estimates
    
    ### Calculate 95% CIs
    suppressMessages(
      CIs <- stats::confint(object = fit, parm = "exposure", level = 0.95) %>%
        as.data.frame() %>%
        dplyr::rename(conf_low = `2.5 %`, conf_high = `97.5 %`)
    )
    
    ### Obtain beta estimates, standard errors, p-values
    
    Intercept <- summary(fit)$coefficients["(Intercept)", "Estimate"]

    Estimate <- summary(fit)$coefficients["exposure", "Estimate"]
    
    MeanBeta <- abs(m2beta(Intercept) - m2beta(Intercept+Estimate))
    
    SE <- summary(fit)$coefficients["exposure", "Std. Error"]
    
    pval <- summary(fit)$coefficients["exposure","Pr(>|t|)"]
    
    #### Calculate p-value using Wald test and LRT
    wald <- survey::regTermTest(model = fit,
                          test.terms = "exposure",
                          null = NULL,
                          df = Inf,
                          method = "Wald")$p
    
    fit.rob <- MASS::rlm(formula = formula.ewas, data = expo_cov, maxit = maxit)
      
    robust_p_value <- survey::regTermTest(model = fit.rob,
                                          test.terms = "exposure",
                                          null = NULL,
                                          df = Inf,
                                          method = "Wald")$p
      
      sfit <- cbind(Estimate, SE, CIs, 
                    raw_p_value = wald,
                    robust_p_value = robust_p_value,
                    MeanBeta)

    return(sfit)
  }

################################################################################
#' Run Mixed Linear regression for methylation at specific loci ~ 1 exposure at a time (ExWAS) adjusted 
#'
#' @param expo_cov A data.frame containing all variables needed for regression (exposures and confounders)
#' @param exposure A character string with the exposure name
#' @param technical_confounders A character vector defining technical confounders the regression will be adjusted to
#' @param CpG A numeric vector containing methylation data from 1 CpG
#' @param confounders A character vector containing names of the confounders
#'
#' @import ggplot2
#' @importFrom stats as.formula p.adjust sd
#' @importFrom jtools center 
#' @importFrom interactions interact_plot
#' @importFrom sjPlot plot_model
#' @importFrom cowplot ggdraw draw_plot draw_plot_label
#' @importFrom lme4 lmer
#' @importFrom dplyr filter
#'
#' @return A named vector of p values for regressions for each CpG

TimeMLReg <- function(CpG,
           covariates.df,
           mainEffect,
           time,
           testTimeInteraction = TRUE,
           threshold.interaction = 0.05,
           clinical_confounders = 0,
           technical_confounders = 0,
           transformToMvalue = transformToMvalue,
           plot = FALSE, linearity.check = FALSE, plot.points = FALSE){
    
    if( transformToMvalue & max(CpG, na.rm = T)<1 & min(CpG, na.rm = T)>0){
      CpG <- logit2(CpG)
    }
  
    CpG <- jtools::center(CpG)
    CpG <- as.matrix(CpG)
    colnames(CpG) <- "y"
    CpG <- data.frame(id = rownames(CpG), CpG)
    
    if(!("id" %in% colnames(covariates.df)))  covariates.df <- data.frame(id = rownames(covariates.df), 
                                                                          covariates.df)
    # Create a subset of data containing methylation data of one CpG (y) and exposure-covariates data (x)
    data <- merge(x = covariates.df, y = CpG, by = "id")
    data <- dplyr::filter(data, !is.na(y))
    
    # Create regression formulae
    formula0 <-
      stats::as.formula(paste(
        mainEffect, " ~ ",
        time, "+",
        paste(clinical_confounders, collapse = " + "),
        "+",
        paste(technical_confounders, collapse = " + "),
        " + (1 | id)"
      ))
    fit0 <- myMLRegression(formula = formula0, data = data)
    
    #With mainEffect
    formula1 <-
      stats::as.formula(paste(
        mainEffect, " ~ ",
        " y ", "+", time, "+", 
        paste(clinical_confounders, collapse = " + "),
        "+",
        paste(technical_confounders, collapse = " + "),
        " + (1 | id)"
      ))
    fit1 <- myMLRegression(formula = formula1, data = data)
    
    #With interaction term mainEffect*Time
    formula2 <-
      stats::as.formula(paste(
        mainEffect, " ~ ",
        " y ", "*", time, "+", 
        paste(clinical_confounders, collapse = " + "),
        "+",
        paste(technical_confounders, collapse = " + "),
        " + (1 | id)"
      ))
    fit2 <- myMLRegression(formula = formula2, data = data)
    
    LRT1 <- anova(fit1,fit2)
    pval1 <- as.numeric(LRT1$`Pr(>Chisq)`[2])
    
    LRT0 <- anova(fit0,fit1)
    pval0 <- as.numeric(LRT0$`Pr(>Chisq)`[2])
    
    sfit <- data.frame(Estimate.interaction = summary(fit2)$coefficients[paste0("y:",time), "Estimate"],
                  Estimate.y2 = summary(fit2)$coefficients["y", "Estimate"],
                  Pval.interaction = pval1,
                  Estimate.y =  summary(fit1)$coefficients["y", "Estimate"],
                  Pval.y = pval0)
    
    if( plot ){
      if(pval1<threshold.interaction){
        pMod <- sjPlot::plot_model(fit2, terms = c("y", time, clinical_confounders, paste0("y:",time))) +
          theme_minimal()
        p1 <- interactions::interact_plot(fit2, pred = !!time, modx = "y",
                                         linearity.check = TRUE, plot.points = TRUE) 
        p2 <- interactions::interact_plot(fit2, pred = !!time, modx = "y") 
     
        p <- ggdraw() +
          draw_plot(pMod, 0, 0.5, .5, .5) +
          draw_plot(p2, .5, 0.5, .5, .5) +
          draw_plot(p1, 0, 0, 1, 0.5) +
          draw_plot_label(c("A", "B", "C"), c(0, 0.5, 0), c(1, 1, 0.5), size = 15) 
        
        print( p )
        
      }else{
        pMod <- sjPlot::plot_model(fit1, c("y", time, clinical_confounders))
        print( pMod )
      }
    }

    
    return( sfit )
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
#' @return A 2-entry list with synchronized data (ie: same samples, matching order)
#'

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