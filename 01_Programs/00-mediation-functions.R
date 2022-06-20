################################################################################
### Author: Max Rubinstein (modified from Alec McLean)
### Purpose: Cleaned estimation functions for double robust estimators
###
################################################################################

library(tidyverse)
library(tidymodels)
library(stacks)
library(tune)

#' Estimate nuisance functions for double robust estimator with missing treatment
#' and outcome data.  See McClean et al. (2021?) for more details on algorithm
#' and theory.
#' 
#' @param outcome_tx_dat Dataset with outcome and treatment variables.
#' It is expected that this dataframe contains an ID column ("id"),
#' an outcome and outcome missingness column ("Y", "R_Y") and a
#' treatment and treatment missingness column ("A", R_A"), and a 
#' complete cases column "cc"
#' 
#' @param X_dat Covariate data.
#' It is expected that this dataframe contains all the covariates 
#' 
#' @param outer_folds Number of folds for cross-validation and splitting
#' 
#' @param sl_lib List of SuperLearner libraries for SuperLearner
#' 
#' @param sl_folds Number of folds for the SuperLearner
#' 
#' @return Dataset with estimated nuisance functions
EstimateNuisanceFunctions <- function(outcome_tx_dat, 
                                      X_dat, 
                                      outer_folds, 
                                      sl_lib,
                                      sl_folds,
                                      mediation) {
  
  library(tidyverse)
  library(SuperLearner)
  library(testit)

  #######################
  ### Data Cleaning
  #######################
  
  ### Create folds
  folds <- sample(outer_folds, nrow(X_dat), replace = T)
  
  ### Turn covariate data into model matrix 
  X_dat <- X_dat %>% 
    stats::model.matrix(~ ., data = .) %>%
    as.data.frame(.)
  
  ### Column naming hack to make models work with SuperLearner
  colnames(X_dat) <- paste0("V", seq(1, ncol(X_dat)))
  
  ### Initialize nuisance function columns
  
  if (mediation == FALSE) {
    dat <- dplyr::select(outcome_tx_dat, id)
    dat$mu_1 <- NA; dat$mu_0 <- NA
    dat$eta <- NA; dat$pi <- NA; 
  }
  
#  if (mediation == FALSE & pscore == TRUE) {
#    dat <- dplyr::select(outcome_tx_dat, id)
#    dat$mu_1 <- NA; dat$mu_0 <- NA
#    dat$eta <- NA; dat$pi <- NA; 
#  }
  
  
  a_vals <- sort(unique(outcome_tx_dat$A[outcome_tx_dat$cc == 1]))
  
  if (mediation == TRUE) {
    
    dat <- dplyr::select(outcome_tx_dat, id)
    
    m_vals <- sort(unique(outcome_tx_dat$M[outcome_tx_dat$cc == 1]))
    
    for (m in m_vals) {
      dat[[paste0("gamma.m", m, ".a1")]] <- NA; dat[[paste0("gamma.m", m, ".a0")]] <- NA
      dat[[paste0("mu_a1_m", m)]] <- NA; dat[[paste0("mu_a0_m", m)]] <- NA
    }
  }
  
  ###################################
  ### Nuisance Function Estimation
  ###################################
  
  for (FOLD in seq(1, outer_folds, 1)) {
    
    cat("\n------------------------------------------------",
        "\n-- Estimating fold ", FOLD, " out of ", outer_folds,
        "\n------------------------------------------------")
    
    TEST_ROWS <- folds == FOLD
    TRAIN_ROWS <- !TEST_ROWS
    
    if (outer_folds == 1) {
      TRAIN_ROWS <- TEST_ROWS
    }
    
    if (mediation == FALSE) {

      #############################
      ### Eta = P(CC = 1 | X)
      #############################
      
      cat("\nEstimating Eta")
      
      # If we observe all outcomes, set eta = 1
      if (sum(outcome_tx_dat$cc[TRAIN_ROWS] != 1) == 0) {
        
        dat$eta[TEST_ROWS] <- 1
        
      } else {
        
        outcome <- outcome_tx_dat$cc[TRAIN_ROWS]
        
          mod_eta <-
            SuperLearner::SuperLearner(
              Y = outcome,   # Outcome
              X = X_dat[TRAIN_ROWS, -ncol(X_dat)], # X value training data (remove id)
              SL.library = sl_lib,                 # Algorithms to use
              family = binomial,                   # Binomial outcome
              cvControl = SuperLearner.CV.control(V = sl_folds, stratifyCV = TRUE)
            )
          
          dat$eta[TEST_ROWS] <-
            predict(mod_eta, newdata = X_dat[TEST_ROWS, -ncol(X_dat)])$pred[,1]
      }
        ########################
        ### Pi: P(A = 1 | X) ###
        ########################
        
        cat("\nEstimating Pi\n")
        
        outcome <- outcome_tx_dat$A[TRAIN_ROWS & outcome_tx_dat$cc == 1]

          mod_pi <- 
            SuperLearner::SuperLearner(
              Y = outcome,   # Outcome
              X = X_dat[TRAIN_ROWS & outcome_tx_dat$cc == 1, -ncol(X_dat)], # X value training data (remove id)
              SL.library = sl_lib,                 # Algorithms to use
              family = binomial,                   # Binomial outcome
              cvControl = SuperLearner.CV.control(V = sl_folds, stratifyCV = TRUE)
            )

          dat[["pi"]][TEST_ROWS] <- predict(mod_pi, newdata = as.matrix(X_dat[TEST_ROWS, -ncol(X_dat)]))$pred[,1]
          
    }
    
    if (mediation == FALSE) {

      ############################
      ### Mu: P(Y | X, A = a) ####
      ############################
      
      cat("\nEstimating Mu")
      
      ### Indicator for each treatment value
      for (a in a_vals) {
        outcome <- outcome_tx_dat$Y[outcome_tx_dat$cc == 1 & outcome_tx_dat$A == a & TRAIN_ROWS]
        
          mod_mu <-
            SuperLearner::SuperLearner(
              Y = outcome,
              X = X_dat[outcome_tx_dat$cc == 1 & outcome_tx_dat$A == a & TRAIN_ROWS, -ncol(X_dat)],
              SL.library = sl_lib,
              family = binomial,
              cvControl = SuperLearner.CV.control(V = sl_folds, stratifyCV = TRUE)
            )
          
          dat[[paste0("mu_", a)]][TEST_ROWS] <- predict(mod_mu, newdata = X_dat[TEST_ROWS, -ncol(X_dat)])$pred[,1]
      }
    }
      
    if (mediation == TRUE) {
      
      ####################################
      ### Mu-NFX: P(Y | X, A = a, M = m) #
      ####################################
      
      cat("\nEstimating Mu-NFX")
      
      ### Indicator for each treatment value
      for (a in a_vals) {
        for(m in m_vals) {
          
          outcome <- outcome_tx_dat$Y[outcome_tx_dat$cc == 1 & outcome_tx_dat$A == a & 
                                        outcome_tx_dat$M == m & TRAIN_ROWS]
          
          mod_mu_med <-
            SuperLearner::SuperLearner(
              Y = outcome,
              X = X_dat[outcome_tx_dat$cc == 1 & outcome_tx_dat$A == a & 
                          outcome_tx_dat$M == m & TRAIN_ROWS, -ncol(X_dat)],
              SL.library = sl_lib,
              family = binomial,
              cvControl = SuperLearner.CV.control(V = sl_folds, stratifyCV = TRUE)
            )
          
          dat[[paste0("mu_a", a, "_m", m)]][TEST_ROWS] <- predict(mod_mu_med, newdata = X_dat[TEST_ROWS, -ncol(X_dat)])$pred[,1]
        }
      }
      
      ################################
      ### Gamma: P(M = m | X, A = a) #
      ################################
      
      cat("\nEstimating Gamma")
      
      if (length(m_vals) <= 2) {
        
        for (a in a_vals) {
          
          outcome <- outcome_tx_dat$M[TRAIN_ROWS & outcome_tx_dat$cc == 1 & outcome_tx_dat$A == a]
          
          mod_gamma <-
            SuperLearner::SuperLearner(
              Y = outcome,
              X = X_dat[TRAIN_ROWS & outcome_tx_dat$cc == 1 & outcome_tx_dat$A == a, -ncol(X_dat)],
              SL.library = sl_lib,
              family = binomial,
              cvControl = SuperLearner.CV.control(V = sl_folds, stratifyCV = TRUE)
            )
          
          pm1hat <- predict(mod_gamma, newdata = X_dat[TEST_ROWS, -ncol(X_dat)])$pred[,1]
          
          dat[[paste0("gamma.m1.a", a)]][TEST_ROWS] <- pm1hat
          dat[[paste0("gamma.m0.a", a)]][TEST_ROWS] <- 1 - pm1hat
        }
      }
      
      if (length(m_vals) > 2) {
        for (a in a_vals) {
          
          outcome <- outcome_tx_dat$M[TRAIN_ROWS & outcome_tx_dat$cc == 1 & outcome_tx_dat$A == a]
          
          mdat <- cbind(
            tibble(M = outcome),
            X_dat[TRAIN_ROWS & outcome_tx_dat$cc == 1 & outcome_tx_dat$A == a, -ncol(X_dat)]
          ) 
          
          mod_gamma <- nnet::multinom(M ~ ., data = mdat, MaxNWts = 3000, maxit = 1000)
          
          preds <- predict(mod_gamma, newdata = X_dat[TEST_ROWS, -ncol(X_dat)], type = "prob")
          
          for (i in 1:length(m_vals)) {
            dat[[paste0("gamma.m", m_vals[i], ".a", a)]][TEST_ROWS] <- preds[,i]
          }
        }
      }
    }
  }
  outcome_tx_dat <- inner_join(outcome_tx_dat, dat, by = "id")
  assert("We lost people!?", nrow(outcome_tx_dat) == nrow(X_dat))
  return(outcome_tx_dat)
}

GenerateTrimWeights <- function(pi, eta, trim, epsilon = 1) {
  plogis(pi - trim, scale = epsilon) * plogis(1 - pi - trim, scale = epsilon) * plogis(eta - trim, scale = epsilon)
}

TrimWeightPrime <- function(R, A, pi, eta, trim, epsilon = 1) {
  a2 <- 1 - trim
  dlogis(pi - trim, scale = epsilon) * plogis(1 - pi - trim, scale = epsilon) * plogis(eta - trim, scale = epsilon) *
    ((R / eta) * (A - pi))
  - plogis(pi - trim, scale = epsilon) * dlogis(1 - pi - trim, scale = epsilon) * plogis(eta - trim, scale = epsilon) *
    (((R / eta) * (A - pi)))
  + plogis(pi - trim, scale = epsilon) * plogis(1 - pi - trim, scale = epsilon) * dlogis(eta - trim, scale = epsilon) *
    (R - eta)
}

#' Calculates the rest of the influence function values that were left 
#' uncalculated in Estimate Nuisance Functions.  Can return a clean or 
#' verbose version of the data
#' 
#' @param data Output from EstimateNuisanceFunctions.
#' @param clean Whether user wants just phi_0 and phi_1 or all the nuisance
#' function values. Assumed FALSE
#' 
#' @return data with full influence function values

CalculateNuisFuncs <- function(data) {
  
    # Bias-correcting weight
    data$weight_y1 <- (data$cc / data$eta) * (data$A / data$pi) 
    data$weight_y0 <- (data$cc / data$eta) * ((1 - data$A) / (1 - data$pi)) 

    # Residuals
    data$ry1 <- (data$Y - data$mu_1)
    data$ry0 <- (data$Y - data$mu_0)
    
    #replace NAs
    data <- data %>%
      replace_na(list(weight_y1 = 0, weight_y0 = 0,
                      ry1 = 0, ry0 = 0))
    
    # EIFs
    data$eif_y1 <- (data$mu_1 + data$weight_y1 * data$ry1) * data$iweights + data$mu_1 * data$iweights.p
    data$eif_y0 <- (data$mu_0 + data$weight_y0 * data$ry0) * data$iweights + data$mu_0 * data$iweights.p
    data$eif_d  <- data$iweights + data$iweights.p

    # contrasts
    data$te <- (data$eif_y1 - data$eif_y0) 
    
  return(data)
}

#' Attaches grouping variables to the CalcNuisFuncs output 
#' 
#' @param models Output from EstimateNuisFuncs
#' 
#' @param analytic_data list of dataframes containing analytic dataset and
#' any alternative subsets of the original data (for sensitivity analyses)
#' 
#' @return data frame containing EIF values and relevant subgroup variables
#' for heterogeneous treatment effect analysis

JoinSubgroups <- function(models, analytic_data, subgroup_effects, mediation) {
  subgroup_data <- map(analytic_data, ~mutate(.x, all = 1) %>%
                         dplyr::select(id, any_of(subgroup_effects))) 

  if (mediation == "FALSE") {
    models <- map(models, CalculateNuisFuncs)
  }
  if (mediation == "TRUE") {
    models <- map(models, CalculateIDEEIFs)
  }
  
  all_data <- map2(models, subgroup_data, ~left_join(.x, .y, by = "id")) 

  data_labels <- c("Analytic data", "Full data", "Sensitivity data", "Analytic alternate")
  
  names(all_data) <- data_labels[1:length(analytic_data)]
  
  return(all_data)
}

delta_ci <- function(proportion, numerator, denominator) {
  res <- proportion*sqrt((var(numerator) / mean(numerator)^2) + 
                           var(denominator) / mean(denominator)^2 - 
                           2*cov(numerator, denominator) / (mean(numerator)*mean(denominator)))
  res <- ifelse(res < 0, -1*res, res)
  return(res)
}

#' Calculates the TE along with the standard errors and 95 percent 
#' condidence intervals for outcomes analysis.
#'  
#' @param data Output from CalculateNuisFuncs
#' 
#' @return data frame with all the estimands and confidence intervals. 

CalculateYAtxfx <- function(data) {
    mu_1 <- mean(data$eif_y1) / mean(data$eif_d); mu_1_tmle <- mean(data$mu_1_tmle)
    mu_0 <- mean(data$eif_y0) / mean(data$eif_d); mu_0_tmle <- mean(data$mu_0_tmle)
    te <- mu_1 - mu_0; te_tmle <- mu_1_tmle - mu_0_tmle
    rr.01 <- mean(data$eif_y0) / mean(data$eif_y1)
    rr.01_tmle <- mean(data$mu_0_tmle) / mean(data$mu_1_tmle)
    rr.10 <- mean(data$eif_y1) / mean(data$eif_y0)
    rr.10_tmle <- mean(data$mu_1_tmle) / mean(data$mu_0_tmle)
    
    se.te   <- delta_ci(te, data$te, data$eif_d) / sqrt(nrow(data))
    se.mu_1 <- delta_ci(mu_1, data$eif_y1, data$eif_d) / sqrt(nrow(data))
    se.mu_0 <- delta_ci(mu_0, data$eif_y0, data$eif_d) / sqrt(nrow(data))

    se.rr.01      <- delta_ci(rr.01, data$eif_y0, data$eif_y1) / sqrt(nrow(data))
    se.rr.01_tmle <- delta_ci(rr.01_tmle, data$mu_0_tmle, data$mu_1_tmle) / sqrt(nrow(data))
    se.rr.10      <- delta_ci(rr.10, data$eif_y1, data$eif_y0) / sqrt(nrow(data))
    se.rr.10_tmle <- delta_ci(rr.10_tmle, data$mu_1_tmle, data$mu_0_tmle) / sqrt(nrow(data))
    
    est.list <- c(te, te_tmle, 
                  rr.01, rr.01_tmle,
                  rr.10, rr.10_tmle,
                  mu_1, mu_1_tmle, 
                  mu_0, mu_0_tmle)
    
    se.list <- c(se.te, sqrt(var(data$eif_tmle)/nrow(data)),
                 se.rr.01, se.rr.01_tmle,
                 se.rr.10, se.rr.10_tmle,
                 se.mu_1, 0,
                 se.mu_0, 0)
    
    res <- tibble(est = est.list, 
                  se = se.list,
                  lci = est - 1.96*se, uci = est + 1.96*se,
                  estimand = rep(c("rd", "rr.01", "rr.10", "mu1", "mu0"), each = 2),
                  estimator = rep(c("OSE", "TMLE"), 5),
                  ess    = sum(data$iweights),
                  ess.dr = sum(data$eif_d),
                  nfull  = nrow(data))
  return(res)
}

EstimateTMLE <- function(dat) {
  dat$H_1 <- with(dat, (cc * A) / (pi * eta))
  dat$H_0 <- with(dat, - (cc * (1 - A)) / (eta * (1 - pi)))
  dat$H_1[is.na(dat$H_1)] <- 0
  dat$H_0[is.na(dat$H_0)] <- 0
  dat$H_a <- dat$H_1 + dat$H_0
  dat$pred <- with(dat, case_when(A == 1 ~ mu_1, A == 0 ~ mu_0, is.na(A) ~ 1))
  dat$pred.logit <- with(dat, log(pred / (1 - pred)))
  dat$pred.logit.1 <- with(dat, log(mu_1 / (1 - mu_1)))
  dat$pred.logit.0 <- with(dat, log(mu_0 / (1 - mu_0)))

  IterFolds <- function(fold, model) {
    fold$mu_1_tmle <- plogis(fold$pred.logit.1 + coef(model)*fold$H_1)
    fold$mu_0_tmle <- plogis(fold$pred.logit.0 + coef(model)*fold$H_0)
    fold$mu_tmle <- with(fold, case_when(A == 1 ~ mu_1_tmle, A == 0 ~ mu_0_tmle, is.na(A) ~ 1)) 
    fold$te_tmle <- fold$mu_1_tmle - fold$mu_0_tmle
    fold$eif_y1_tmle <- (fold$Y - fold$mu_1_tmle) * fold$H_1 + fold$mu_1_tmle
    fold$eif_y0_tmle <- (fold$Y - fold$mu_0_tmle) * fold$H_0 + fold$mu_0_tmle
    fold$eif_y1_tmle[is.na(fold$eif_y1_tmle)] <- fold$mu_1_tmle[is.na(fold$eif_y1_tmle)]
    fold$eif_y0_tmle[is.na(fold$eif_y0_tmle)] <- fold$mu_0_tmle[is.na(fold$eif_y0_tmle)]
    fold$eif_tmle <- fold$eif_y1_tmle - fold$eif_y0_tmle
    return(fold)
  }
  
  tmle <- glm(Y ~ -1 + H_a, offset = pred.logit, family = binomial(), data = dat[dat$cc == 1,])
  res <- IterFolds(dat, tmle)
  
  return(res)
}

#' Calculates the TE, NIE, NDE, and the TR-version of the TE (which is a less
#' efficient estimate) along with the standard errors and 95 percent 
#' confidence intervals. Note that these confidence intervals assume
#' that all models are correctly specified in the parametric case, or that
#' they converge at sufficient rates in the non-parametric case
#' 
#' @param data output from JoinSubgroups
#' 
#' @param group_list vector of grouping variables to estimate effect by
#' 
#' @return data frame with all the estimands and confidence intervals by sub-group 
CalculateTxfx <- function(analytic_data, models, group_list, mediation, cc,
                          trim = 0, s = 0.01, xg = FALSE) {
  
  if (cc == FALSE) {
    models <- map(models, ~mutate(.x, iweights    = GenerateTrimWeights(pi, eta, trim, epsilon = s),
                                      iweights.p  = TrimWeightPrime(cc, A, pi, eta, trim, epsilon = s)))
  }
  if (cc == TRUE) {
    models <- map(models, ~filter(.x, cc == 1) %>%
                    mutate(eta = 1, 
                           iweights   = GenerateTrimWeights(pi, eta, trim, s),
                           iweights.p = TrimWeightPrime(cc, A, pi, eta, trim, s)))
  }
  
  if (all(group_list %in% names(models[[1]]))) {
    data <- map(models, CalculateIDEEIFs)
  }
  
  if (any(!group_list %in% names(models[[1]]))) {
    glist <- if (any(grepl("^cc$", group_list))) group_list[-grep("^cc$", group_list)] else group_list
    models <- map(models, ~.x %>% select(-any_of(glist)))
    data <- JoinSubgroups(models, analytic_data, group_list, mediation = mediation) 
  }

  ResultsFun <- switch(mediation,
                       "FALSE" = CalculateYAtxfx,
                       "TRUE"  = CalculateIDEtxfx)

  subfx <- function(data, group){
    map(data, ~.x %>%
                   nest(-sym(!!group)) %>%
                   mutate(data = map(data, ~ResultsFun(.x))) %>%
                   unnest(cols = data)) %>%
    map2(names(data), ~mutate(.x, Dataset = .y)) %>%
    invoke(rbind, .)
  }
  
  effects <- map(group_list, ~subfx(data, .x))
  
  res <- map2(effects, group_list, ~mutate(.x, variable_name = .y)) %>%
      map(~set_names(.x, c("value", names(.x)[-1]))) %>%
      invoke(rbind, .) %>%
      mutate_at("value", as.numeric) 
  
  if (any(group_list %in% names(analytic_data))) {
    res <- res %>%
      left_join(codebook, by = "variable_name") %>%
      left_join(labels %>% mutate_at("value", as.numeric)) %>%
      replace_na(list(variable = "all", description = "All"))
  }
  
  if (mediation == "TRUE") {
    final <- list(res = res, data = data)
  }
  
  if (mediation == "FALSE") {
    final <- res
  }
  return(final)
}

JoinMedTable <- function(results1, results2) {
  map(results1, ~.x$res) %>%
    map2(rep(c("Binomial", "Multinomial"), 3), ~mutate(.x, Specification = .y)) %>%
    map2(c("M2", "M2", "M1", "M1", "M1M2", "M1M2"), ~mutate(.x, mediator = .y)) %>%
    invoke(rbind, .) %>%
    bind_rows(results2)  
}


#' Calculates model diagnostics for binary classification models used to calculate effects
#' where the predicted values are among the observed factual outcomes
#' 
#' @param data output from CalcNuisFuncs
#' #' 
#' @return dataframe containing AUC, TPR, TNR, PPV, NPV, Gini, F1, and
#' Brier scores for all of the associated nuisance functions 
CalcDiagnostics <- function(data) {
  # Note: need to edit this to handle multinomial mediator values
  library(ModelMetrics)
  calc_metric <- function(metric, data, outcome, pred) {
    metric(data[[outcome]], data[[pred]])
  }
  
  yd.a1 <- filter(data, A == 1, !is.na(Y)) 
  yd.a0 <- filter(data, A == 0, !is.na(Y)) 
  md.a1 <- filter(data, A == 1, !is.na(M)) 
  md.a0 <- filter(data, A == 0, !is.na(M)) 
  yd.a1m1 <- filter(data, A == 1, M == 1, !is.na(Y), !is.na(M)) 
  yd.a0m0 <- filter(data, A == 0, M == 0, !is.na(Y), !is.na(M)) 
  yd.a1m0 <- filter(data, A == 1, M == 0, !is.na(Y), !is.na(M)) 
  yd.a0m1 <- filter(data, A == 0, M == 1, !is.na(Y), !is.na(M)) 
  ad <- filter(data, !is.na(A)) 
  
  eval_models <- tibble(
    model = c("mu_1", "mu_0", "plug_in_y1m1", "plug_in_y0m0", 
              "mu_a1_m1", "mu_a1_m0", "mu_a0_m1", "mu_a0_m0", 
              "gamma_1", "gamma_0", "pi", "eta"),
    data = list(yd.a1, yd.a0, yd.a1, yd.a0, yd.a1m1, yd.a1m0,
                yd.a0m1, yd.a0m0, md.a1, md.a0, ad, data)
  ) %>%
    mutate(outcome = c(rep("Y", 8), "M", "M", "A", "cc"),
           pred = c("mu_1", "mu_0", "plug_in_y1m1", "plug_in_y0m0",
                    "mu_a1_m1", "mu_a1_m0", "mu_a0_m1", "mu_a0_m0",
                    "gamma1", "gamma0", "pi", "eta")) %>%
    mutate(auc = pmap(list(data = data, outcome = outcome, pred = pred), calc_metric, auc),
           tpr = pmap(list(data = data, outcome = outcome, pred = pred), calc_metric, sensitivity),
           tnr = pmap(list(data = data, outcome = outcome, pred = pred), calc_metric, specificity),
           ppv = pmap(list(data = data, outcome = outcome, pred = pred), calc_metric, ppv),
           npv = pmap(list(data = data, outcome = outcome, pred = pred), calc_metric, npv),
           gini = pmap(list(data = data, outcome = outcome, pred = pred), calc_metric, gini),
           f1 = pmap(list(data = data, outcome = outcome, pred = pred), calc_metric, f1Score),
           brier = pmap(list(data = data, outcome = outcome, pred = pred), calc_metric, brier))
  
  eval_models %>%
    select(-data) %>%
    unnest()
}

########################################
### Difference in means analysis
########################################

#' Calculates the differences in means for a treatment A and
#' outcome Y
#' 
#' @param data Output from EstimateNuisanceFunctions.
#' @param weights returns a weighted mean if a weight
#' vector is supplied. Assumed NULL
#' 
#' @return dataframe with difference in means and relative risk reduction

diff_means <- function(data, weights = NULL) {
  if (is.null(weights)) {
    weights <- rep(1, nrow(data))
  }
  
  model <- lm(Y ~ A, data = data, weights = weights)
  
  model_vcov <- sandwich::vcovHC(model, type = "HC1")
  y0 <- coef(model)[1]
  y1 <- coef(model)[2] + y0
  var.y1 <- (c(1, 1) %*% model_vcov %*% c(1, 1))[1,1]
  var.y0 <- model_vcov[1, 1]
  cov.y0y1 <- var.y0 + model_vcov[1, 2]
  
  rr.01 <- y0 / y1
  rr.10 <- y1 / y0
  
  se.rr.01 <- rr.01 * sqrt(var.y0/(y0^2) + var.y1/(y1^2) - 2 * cov.y0y1/(y0 * y1))
  se.rr.10 <- rr.10 * sqrt(var.y0/(y0^2) + var.y1/(y1^2) - 2 * cov.y0y1/(y0 * y1))
  
  lmtest::coeftest(model, vcov = model_vcov) %>%
    broom::tidy() %>%
    select(estimate, SE = `std.error`) %>%
    mutate(estimand = c("y0", "rd")) %>%
    bind_rows(tibble(estimate = c(rr.01, rr.10, y1), 
                     SE = c(se.rr.01, se.rr.10, sqrt(var.y1)), 
                     estimand = c("rr.01", "rr.10", "y1")))
}

#' Calculates the difference in means for subgroups 
#' contained in a dataframe
#' 
#' @param data dataset containing outcome variable Y,
#' treatment A, and subgroups specified by the category
#' argument
#' 
#' @param weights returns a weighted mean using a column
#' named eta if TRUE. Assumed FALSE
#' 
#' @return list of dataframes containing dim estimates by subgroup

dim_estimators <- function(data, category, weights = FALSE) {
  if (weights == TRUE) {
    wts <- map(data, ~1/.x$eta)
  }
  if (weights == FALSE) {
    wts <- map(data, ~rep(1, nrow(.x)))
  }
  map2(data, wts, ~.x %>%
         dplyr::select(!!category, Y, A) %>%
         mutate(weights = .y) %>%
         nest(data = any_of(c("Y", "A", "weights"))) %>%
         mutate(model = map(data, ~diff_means(.x, .x$weights))) %>%
         select(!!category, model) %>%
         unnest(cols = c(model)) %>%
         mutate(l95ci = estimate - 1.96*SE, u95ci = estimate + 1.96*SE))
}

iter_dim <- function(data_list, subgroups) {
  dim_prim <- map(subgroups, ~dim_estimators(data_list, .x, weights = FALSE))
  dim_prim_w <- map(subgroups, ~dim_estimators(data_list, .x, weights = TRUE))
  list(dim_prim, dim_prim_w)
}


#' Calculates the difference in means for subgroups 
#' contained in a dataframe
#' 
#' @param results output from the dim_estimators function
#' 
#' @param subgroups a vector containing names of variables that 
#' specify the relevant subgroup effects
#' 
#' @return dataframe containing dim estimates by subgroup

subgroup_results <- function(results, subgroups) {
  results %>%
    map2(subgroups, ~mutate(.x, variable_name = .y)) %>%
    map(~set_names(.x, c("value", names(.x)[-1]))) %>%
    invoke(rbind, .) %>%
    mutate_at("value", ~as.numeric(as.character(.)))
}


# estimate parametric models
EstimateGLMmodels <- function(outcome_data, out_file, mediation) {
  #assertthat::assert_that(grepl("\\.rds", out_file))
  glm_models <- map2(outcome_data$outcome_data, outcome_data$X_data, 
                     ~EstimateNuisanceFunctions(.x, .y, outer_folds = 1,
                                                sl_lib = "SL.glm", sl_folds = 2,
                                                mediation = mediation))

  saveRDS(glm_models, paste0("03_ModelOutput/", out_file))
  glm_models
}

# subgroup differences
gen_figure_data <- function(sub_results, subgroups) {
  data <- map(transpose(sub_results), ~subgroup_results(.x, subgroups)) %>%
    map2(data_labels, ~mutate(.x, Dataset = .y)) %>%
    map(~left_join(.x, labels, by = c("variable_name", "value"))) %>%
    map(~filter(.x, coef == 2)) 
  
  orders <- labels %>%
    right_join(data[[1]] %>% dplyr::select(variable_name, label_full), by = c("variable_name", "label_full")) %>%
    filter(!label_full %in% c("No response", "Missing")) %>%
    .$label_full
  
  data %>%
    map(~dplyr::filter(.x, !label_full %in% c("Missing", "No response")) %>%  
          mutate(label_full = factor(label_full, levels = orders))) %>%
    invoke(rbind, .) 
}

format_cells <- function(num, lci, uci, digits, mult = 1) {
  paste0(format(round(mult*num, digits), nsmall = digits), " (", 
         format(round(mult*lci, digits), nsmall = digits), ", ", 
         format(round(mult*uci, digits), nsmall = digits), ")")
}

# calculate classification diagnostics for a model
DiagnosticTable <- function(data, mediate = TRUE) {
  a_vals <- sort(unique(data$A[data$cc == 1]))
  
  # eta diagnostics
  eta_pred  <- factor(ifelse(data$eta > 0.5, 1, 0), levels = c(0, 1))
  #  eta_auc <- auc(roc(data$cc, eta_pred))[[1]]
  eta_diagnostics <- ml_test(eta_pred, data$cc, output.as.table = TRUE) %>%
    rownames_to_column() %>%
    filter(rowname != 0) %>%
    mutate(outcome = "CC", subset = "All", N = length(eta_pred), class_pct = mean(data$cc))
  
  # pi diagnostics
  pi_pred  <- factor(ifelse(data$pi[!is.na(data$A)] > 0.5, 1, 0), levels = c(0, 1))
  #  pi_auc <- auc(roc(data$A[!is.na(data$A)], pi_pred))[[1]]
  pi_diagnostics <- ml_test(pi_pred, data$A[!is.na(data$A)], output.as.table = TRUE) %>%
    rownames_to_column() %>%
    filter(rowname != 0) %>%
    mutate(outcome = "A", subset = "All", N = length(pi_pred), class_pct = mean(data$pi))
  
  # mu_a diagnostics
  outcome_diagnostics <- list()
  for (a in a_vals) {
    sub_data <- data %>%
      filter(A == a, !is.na(Y))
    class.pred <- factor(ifelse(sub_data[[paste0("mu_", a)]] > 0.5, 1, 0), levels = c(0, 1))
    #    mu_auc <- auc(roc(sub_data$Y, class.pred))[[1]]
    out_diag <- ml_test(class.pred, sub_data$Y, output.as.table = TRUE) %>%
      rownames_to_column() %>%
      filter(rowname != 0) %>%
      mutate(outcome = "Y_A", subset = paste0("A = ", a), N = length(class.pred), class_pct = mean(sub_data$Y))
    outcome_diagnostics <- append(outcome_diagnostics,
                                  list(out_diag))
  }
  
  if (mediate == TRUE) {
    m_vals <- sort(unique(data$M[data$cc == 1]))
    # mu_{a, m} diagnostics
    nde_diagnostics <- list()
    for (a in a_vals) {
      for (m in m_vals) {
        am_string <- paste0("mu_a", a, "_m", m)
        sub_data <- data %>%
          filter(A == a, M == m, !is.na(Y)) 
        class.pred <- factor(ifelse(sub_data[[paste0("mu_a", a, "_m", m)]] > 0.5, 1, 0), levels = c(0, 1))
        #        mu_auc <- auc(roc(sub_data$Y, class.pred))[[1]]
        out_diag <- ml_test(class.pred, sub_data$Y, output.as.table = TRUE) %>%
          rownames_to_column() %>%
          filter(rowname != 0) %>%
          mutate(outcome = "Y_AM", subset = paste0("A = ", a, ", M = ", m), 
                 N = length(class.pred), class_pct = mean(sub_data$Y))
        
        nde_diagnostics <- append(nde_diagnostics, list(out_diag))
      }
    }
    
    # gamma diagnostics
    mediator_diagnostics <- list()
    
    for (a in a_vals) {
      a_string <- paste0("a", a)
      sub_data <- data %>%
        filter(A == a, !is.na(M)) 
      outcome <- sub_data$M
      pcts <- as.numeric(prop.table(table(sub_data$M)))
      Ns <- as.numeric(table(sub_data$M))
      sub_data <- sub_data %>%
        select(starts_with("gamma")) %>%
        select(ends_with(a_string)) %>%
        select(sort(names(.))) 
      class.pred <- factor(max.col(sub_data, ties.method = "first"), 
                           levels = sort(unique(outcome)))
      bin.tru <- map(sort(unique(outcome)), ~ifelse(outcome == .x, 1, 0))
      bin.prd <- map(sort(unique(outcome)), ~ifelse(class.pred == .x, 1, 0))
      #      auc.list <- unlist(map2(bin.tru, bin.prd, ~auc(roc(.x, .y))[[1]]))
      out_diag <- ml_test(class.pred, outcome, output.as.table = TRUE) %>%
        mutate(outcome = "M", subset = paste0("A = ", a),
               N = Ns, class_pct = pcts)
      mediator_diagnostics <- append(mediator_diagnostics, list(out_diag))
    }
    res <- bind_rows(
      eta_diagnostics,
      pi_diagnostics,
      invoke(rbind, mediator_diagnostics),
      invoke(rbind, nde_diagnostics),
      invoke(rbind, outcome_diagnostics)
    )
  }
  
  if (mediate == FALSE) {
    res <- bind_rows(
      eta_diagnostics,
      pi_diagnostics,
      invoke(rbind, outcome_diagnostics))
  }
  return(res)
}

# Summary Table for Mediation analysis
DiagnosticsSummary <- function(diagnostic_table_list, mediator_specification) {
  map(diagnostic_table_list, ~.x %>%
        select(BA = balanced.accuracy, F1, FDR, FNR, FPR, precision, recall, specificity, outcome, N, class_pct) %>%
        group_by(outcome) %>%
        summarize_if(is.numeric, ~sum(.*N)/sum(N)) %>%
        select(-N)) %>%
    map2(rep(c("Worried", "Isolated", "Both"), each = 2), ~mutate(.x, Mediator = .y)) %>%
    map2(rep(c("Binomial", "Multinomial"), 3), ~mutate(.x, Specification = .y)) %>%
    invoke(rbind, .) %>%
    mutate(outcome = case_when(
      outcome == "CC" ~ "R",
      outcome == "M" ~ paste0("M: ", Mediator),
      outcome == "Y_AM" ~ paste0("Y_AM", ": ", Mediator),
      TRUE ~ outcome
    )) %>%
    filter(Specification == mediator_specification) %>%
    distinct(Outcome = outcome, `Balanced Accuracy` = BA, Precision = precision, Sensitivity = recall,
             Specificity = specificity) %>%
    mutate_at("Outcome", ~factor(., levels = c("R", "A", "Y_A", "Y_AM: Worried",
                                               "Y_AM: Isolated", "Y_AM: Both",
                                               "M: Worried", "M: Isolated", "M: Both"))) %>%
    arrange(Outcome)
}

# Summary Table for Outcomes Analysis
OutcomeDiagTable <- function(outcome_diagnostics, outcome_order) {
  outcome_diagnostics %>%
    map2(outcome_order, ~mutate(.x, Outcome = .y)) %>%
    invoke(rbind, .) %>%
    filter(outcome == "Y_A") %>%
    select(Outcome, `Balanced Accuracy` = balanced.accuracy, Precision = precision, Sensitivity = recall,
           Specificity = specificity, N) %>%
    group_by(Outcome) %>%
    summarize_if(is.numeric, ~sum(.*N)/sum(N)) %>%
    select(-N)
}

################################################################################
########################### INTERVENTIONAL EFFECTS #############################
################################################################################
MarginalizeJointDensity <- function(data, marginalize) {
  GetCols <- function(col_string) {
    unlist(map(col_string, ~grep(.x, names(data), value = TRUE)))
  }
  
  ValueList <- function(col_values) {
    Reduce(`+`, map(col_values, ~data[[.x]]))
  }
  
  GenXWalk <- function(data, marginal_mediator_column) {
    data %>%
      filter(cc == 1) %>%
      select(all_of(c("M", marginal_mediator_column))) %>%
      distinct() %>%
      arrange(marginal_mediator_column) %>%
      split(f = .[[marginal_mediator_column]]) %>%
      map(~.x$M)
  }
  
  if (marginalize == TRUE) {
    M1_values <- GenXWalk(data, "M1")
    M2_values <- GenXWalk(data, "M2")
    
    M1_v <- sort(unique(data$M1))
    M2_v <- sort(unique(data$M2))
    
    for(a in c(0, 1)) {
      m1.chars <- map(M1_values, ~paste0("m", .x, ".a", a))
      m2.chars <- map(M2_values, ~paste0("m", .x, ".a", a))
      m1_col_vals <- map(m1.chars, GetCols)
      m2_col_vals <- map(m2.chars, GetCols)
      
      m1_list <- map(m1_col_vals, ValueList)
      m1_labels <- unlist(map(M1_v, ~paste0("gamma.M1", .x, ".a", a)))
      m2_list <- map(m2_col_vals, ValueList)
      m2_labels <- unlist(map(M2_v, ~paste0("gamma.M2", .x, ".a", a)))
      
      for (i in 1:length(M1_values)) {
        data[[m1_labels[i]]] <- m1_list[[i]]
        data[[m2_labels[i]]] <- m2_list[[i]]
      }
    }
  }
  return(data)
}

MarginalizeOutcomesJointDensityInt <- function(data, a_val, m1_val, m2_val) {
  xwalk <- data %>%
    filter(cc == 1) %>%
    select(M, M1, M2) %>%
    distinct() %>%
    arrange(M) %>%
    mutate(M1_label = paste0("gamma.M1", M1, ".a", m1_val),
           M2_label = paste0("gamma.M2", M2, ".a", m2_val),
           Y_label = paste0("mu_a", a_val, "_m", M))
  
  sum_list <- list()
  
  for(i in 1:nrow(xwalk)) {
    sum_list[[i]] <- data[[xwalk$Y_label[i]]]*data[[xwalk$M1_label[i]]]*data[[xwalk$M2_label[i]]]
  }
  
  data[[paste0("mu_a", a_val, "_m1.", m1_val, "xm2.", m2_val)]] <- Reduce(`+`, sum_list)
  
  return(data)
}

TauXAA <- function(data, mediator, a_val, m_val) {
  
  medstr <- ifelse(mediator == "M", "m", mediator)
  
  xwalk <- data %>%
    filter(cc == 1) %>%
    select(mediator) %>%
    distinct() %>%
    arrange(mediator) %>%
    mutate(M_label = paste0("gamma.", medstr, !!sym(mediator), ".a", m_val),
           Y_label  = case_when(
             mediator == "M"  ~ paste0("mu_a", a_val, "_m", !!sym(mediator)),
             mediator == "M1" ~ paste0("mu_a", a_val, "_marg_m2a", a_val, "_m1", !!sym(mediator)),
             mediator == "M2" ~ paste0("mu_a", a_val, "_marg_m1a", a_val, "_m2", !!sym(mediator)))
             )

  sum_list <- list()
  
  for(i in 1:nrow(xwalk)) {
    sum_list[[i]] <- data[[xwalk$Y_label[i]]]*data[[xwalk$M_label[i]]]
  }
  
  data[[paste0("tau.", mediator, ".", a_val, m_val)]] <- Reduce(`+`, sum_list)
  
  return(data)
}

MarginalizeOutcomesJointDensityNat <- function(data, a_val, ma_val) {
  xwalk <- data %>%
    filter(cc == 1) %>%
    select(M) %>%
    distinct() %>%
    arrange(M) %>%
    mutate(M_label = paste0("gamma.m", M, ".a", ma_val),
           Y_label = paste0("mu_a", a_val, "_m", M))
  
  sum_list <- list()
  
  for(i in 1:nrow(xwalk)) {
    sum_list[[i]] <- data[[xwalk$Y_label[i]]]*data[[xwalk$M_label[i]]]
  }
  
  data[[paste0("mu_a", a_val, "_m1m2.a", ma_val)]] <- Reduce(`+`, sum_list)
  
  return(data)
}

MarginalizeOutcomesMarginalDensity <- function(data, a_val_y, a_val_m, mediator) {
  med.obs <- ifelse(mediator == "M1", "M2", "M1")
  
  med.obs.vals <- sort(unique(data[[med.obs]][data$cc == 1]))
  med.vals <- sort(unique(data[[mediator]][data$cc == 1]))
  
  for (i in 1:length(med.obs.vals)) {
    med_vars_m1 <- map(med.vals, ~paste0("jmu_a", a_val_y, "_m1", .x, "_m2", med.obs.vals[i]))
    med_vars_m2 <- map(med.vals, ~paste0("jmu_a", a_val_y, "_m1", med.obs.vals[i], "_m2", .x))
    med_vars    <- switch(mediator, "M1" = med_vars_m1, "M2" = med_vars_m2)
    med_density_vars <- map(med.vals, 
                            ~grep(paste0("gamma\\.", mediator, .x, "\\.a", a_val_m), 
                                  names(data), value = TRUE))
    sum_list <- list()
    
    for(j in 1:length(med_vars)) {
      sum_list[[j]] <- data[[med_vars[[j]]]]*data[[med_density_vars[[j]]]]
    }
    
    data[[paste0("mu_a", a_val_y, "_marg_", tolower(mediator), "a", 
                 a_val_m, "_", tolower(med.obs), med.obs.vals[i])]] <- 
      Reduce(`+`, sum_list)
  }
  return(data)
}

ReturnObservedM <- function(data, mediator) {
  m_vals <- sort(unique(data[[mediator]][data$cc == 1]))
  
  med_val <- mediator
  if (mediator == "M") med_val <- "m"
  
  gamma_a0_val <- unlist(map(paste0("gamma.", med_val, m_vals, ".a0"), ~grep(.x, names(data), value = TRUE)))
  gamma_a1_val <- unlist(map(paste0("gamma.", med_val, m_vals, ".a1"), ~grep(.x, names(data), value = TRUE)))
  
  for (i in 1:length(m_vals)) {
    data[[paste0("gamma.a0.", med_val, ".obs")]][data[[mediator]] == m_vals[i] | is.na(data[[mediator]])] <- data[[gamma_a0_val[i]]][data[[mediator]] == m_vals[i]]
    data[[paste0("gamma.a1.", med_val, ".obs")]][data[[mediator]] == m_vals[i] | is.na(data[[mediator]])] <- data[[gamma_a1_val[i]]][data[[mediator]] == m_vals[i]]
  }
  return(data)
}

ReturnObservedY <- function(data) {
  m_vals <- sort(unique(data$M[data$cc == 1]))
  m1_vals <- sort(unique(data$M1[data$cc == 1]))
  m2_vals <- sort(unique(data$M2[data$cc == 1]))
  
  m1_range <- paste0("m1[", m1_vals[1], "-", m1_vals[length(m1_vals)], "]$")
  m2_range <- paste0("m2[", m2_vals[1], "-", m1_vals[length(m2_vals)], "]$")
  
  mu_0M_val <- unlist(map(paste0("mu_a0_m", m_vals, "$"), ~grep(.x, names(data), value = TRUE)))
  mu_1M_val <- unlist(map(paste0("mu_a1_m", m_vals, "$"), ~grep(.x, names(data), value = TRUE)))
  mu_0M1a0_val <- unlist(map(paste0("mu_a0_marg_m2a0_", m1_range), ~grep(.x, names(data), value = TRUE)))
  mu_1M1a0_val <- unlist(map(paste0("mu_a1_marg_m2a0_", m1_range), ~grep(.x, names(data), value = TRUE)))
  mu_0M2a0_val <- unlist(map(paste0("mu_a0_marg_m1a0_", m2_range), ~grep(.x, names(data), value = TRUE)))
  mu_1M2a0_val <- unlist(map(paste0("mu_a1_marg_m1a0_", m2_range), ~grep(.x, names(data), value = TRUE)))
  mu_0M1a1_val <- unlist(map(paste0("mu_a0_marg_m2a1_", m1_range), ~grep(.x, names(data), value = TRUE)))
  mu_1M1a1_val <- unlist(map(paste0("mu_a1_marg_m2a1_", m1_range), ~grep(.x, names(data), value = TRUE)))
  mu_0M2a1_val <- unlist(map(paste0("mu_a0_marg_m1a1_", m2_range), ~grep(.x, names(data), value = TRUE)))
  mu_1M2a1_val <- unlist(map(paste0("mu_a1_marg_m1a1_", m2_range), ~grep(.x, names(data), value = TRUE)))
  
  for (i in 1:length(m_vals)) {
    data[[paste0("plug_in_mu1m.obs")]][data$M == m_vals[i] | is.na(data$M)] <- data[[mu_1M_val[i]]][data$M == m_vals[i]]
    data[[paste0("plug_in_mu0m.obs")]][data$M == m_vals[i] | is.na(data$M)] <- data[[mu_0M_val[i]]][data$M == m_vals[i]]
  }
  
  for(i in 1:length(m1_vals)) {
    data[[paste0("plug_in_mu1m1.marg.m2a0.obs")]][data$M1 == m1_vals[i] | is.na(data$M1)] <- data[[mu_1M1a0_val[i]]][data$M1 == m1_vals[i]]
    data[[paste0("plug_in_mu0m1.marg.m2a0.obs")]][data$M1 == m1_vals[i] | is.na(data$M1)] <- data[[mu_0M1a0_val[i]]][data$M1 == m1_vals[i]]
    data[[paste0("plug_in_mu1m1.marg.m2a1.obs")]][data$M1 == m1_vals[i] | is.na(data$M1)] <- data[[mu_1M1a1_val[i]]][data$M1 == m1_vals[i]]
    data[[paste0("plug_in_mu0m1.marg.m2a1.obs")]][data$M1 == m1_vals[i] | is.na(data$M1)] <- data[[mu_0M1a1_val[i]]][data$M1 == m1_vals[i]]
  }
  
  for(i in 1:length(m2_vals)) {
    data[[paste0("plug_in_mu1m2.marg.m1a0.obs")]][data$M2 == m2_vals[i] | is.na(data$M2)] <- data[[mu_1M2a0_val[i]]][data$M2 == m2_vals[i]]
    data[[paste0("plug_in_mu0m2.marg.m1a0.obs")]][data$M2 == m2_vals[i] | is.na(data$M2)] <- data[[mu_0M2a0_val[i]]][data$M2 == m2_vals[i]]
    data[[paste0("plug_in_mu1m2.marg.m1a1.obs")]][data$M2 == m2_vals[i] | is.na(data$M2)] <- data[[mu_1M2a1_val[i]]][data$M2 == m2_vals[i]]
    data[[paste0("plug_in_mu0m2.marg.m1a1.obs")]][data$M2 == m2_vals[i] | is.na(data$M2)] <- data[[mu_0M2a1_val[i]]][data$M2 == m2_vals[i]]
  }
  return(data)
}

CalculateIDEEIFs <- function(data, stabilize = FALSE) {
  
  hybrd.replace <- function(x) { mutate_all(x, ~replace(., is.na(.), 0)) }
  data <- hybrd.replace(data)
  
  # EIF: Denominator
  data$eif_d <- with(data, iweights + iweights.p)
  
  # EIF: Total Effect
  data$weight_a1 <- with(data, ((A * cc) / (pi * eta)))
  data$weight_a0 <- with(data, (((1 - A) * cc) / ((1 - pi) * eta)))
  
  if (stabilize == TRUE) {
    data$weight_a1 <- data$weight_a1 / mean(data$weight_a1)
    data$weight_a0 <- data$weight_a0 / mean(data$weight_a0)
  }
  
  data$mu_A <- with(data, ifelse(A == 1, mu_1, mu_0))
  data$eif_te_1 <- with(data, (Y - mu_A) * (weight_a1 - weight_a0))
  data$eif_te_2 <- with(data, mu_1 - mu_0)
  data$eif_te  <- with(data, (eif_te_1 + eif_te_2) * iweights + (mu_1 - mu_0) * iweights.p)

  # EIF: Interventional Direct Effect (note: equal to NDE of M1-M2)
  
  # Ref: A = 1 (pure)
  data$ratio.1 <- with(data, ifelse(gamma.a1.m.obs != 0, gamma.a0.m.obs / gamma.a1.m.obs, 0))
  data$eif_de_t1_a1 <- with(data, weight_a1 * ratio.1 * (Y - plug_in_mu1m.obs))
  data$eif_de_t2_a1 <- with(data, weight_a0 * ((plug_in_mu1m.obs - mu_a1_m1m2.a0) - (Y - mu_0)))
  data$eif_de_t3_a1 <- with(data, mu_a1_m1m2.a0 - mu_0)
  
  # Ref: A = 0 (total)
  data$ratio.0 <- with(data, ifelse(gamma.a0.m.obs != 0, gamma.a1.m.obs / gamma.a0.m.obs, 0))
  data$eif_de_t1_a0 <- with(data, weight_a0 * ratio.0 * (Y - plug_in_mu0m.obs))
  data$eif_de_t2_a0 <- with(data, weight_a1 * ((plug_in_mu0m.obs - mu_a0_m1m2.a1) - (Y - mu_1)))
  data$eif_de_t3_a0 <- with(data, mu_a0_m1m2.a1 - mu_1)

  if (stabilize == TRUE) {
    data$eif_de_t1_a1 <- data$eif_de_t1_a1 / with(data, mean(weight_a1 * ratio.1))
    data$eif_de_t2_a1 <- data$eif_de_t2_a1 / mean(data$weight_a0)
    data$eif_de_t1_a0 <- data$eif_de_t1_a0 / with(data, mean(weight_a0 * ratio.0))
    data$eif_de_t2_a0 <- data$eif_de_t2_a0 / mean(data$weight_a1)
  }
  
  data$eif_de_a1 <- with(data, (eif_de_t1_a1 + eif_de_t2_a1 + eif_de_t3_a1) * iweights 
                         + (mu_a1_m1m2.a0 - mu_0) * iweights.p) 
  data$eif_de_a0 <- with(data, -1 * ((eif_de_t1_a0 + eif_de_t2_a0 + eif_de_t3_a0) * iweights 
                                     + (mu_a0_m1m2.a1 - mu_1) * iweights.p)) 

  # EIF: Interventional Indirect Effect: M1
  
  # Ref: A = 1
  data$ratio.m1.1 <- with(data, ifelse(gamma.a1.m.obs != 0, gamma.a0.M2.obs / gamma.a1.m.obs, 0))
  data$eif_m1_t1_a1 <- with(data, weight_a1 * (gamma.a1.M1.obs - gamma.a0.M1.obs) * ratio.m1.1
                            * (Y - plug_in_mu1m.obs))
  data$eif_m1_t2_a1 <- with(data,  weight_a1  * (plug_in_mu1m1.marg.m2a0.obs - mu_a1_m1.1xm2.0))
  data$eif_m1_t3_a1 <- with(data,  -weight_a0 * (plug_in_mu1m1.marg.m2a0.obs - mu_a1_m1.0xm2.0))
  data$eif_m1_t4_a1 <- with(data,  weight_a0  * ((plug_in_mu1m2.marg.m1a1.obs - plug_in_mu1m2.marg.m1a0.obs) -
                                                   (mu_a1_m1.1xm2.0 - mu_a1_m1.0xm2.0)))
  data$eif_m1_t5_a1 <- with(data, mu_a1_m1.1xm2.0 - mu_a1_m1.0xm2.0)

  # Breaking up into two component estimands:
  data$eif_m1_t1_a1.1 <- with(data, weight_a1 * gamma.a1.M1.obs * ratio.m1.1 * (Y - plug_in_mu1m.obs))
  data$eif_m1_t2_a1.1 <- with(data, weight_a1 * (plug_in_mu1m1.marg.m2a0.obs  - mu_a1_m1.1xm2.0))
  data$eif_m1_t3_a1.1 <- with(data, weight_a0 * ((plug_in_mu1m2.marg.m1a1.obs - mu_a1_m1.1xm2.0)))
  data$eif_m1_t4_a1.1 <- with(data, mu_a1_m1.1xm2.0)

  data$eif_m1_t1_a1.0 <- with(data, weight_a1 * gamma.a0.M1.obs * ratio.m1.1 * (Y - plug_in_mu1m.obs))
  data$eif_m1_t2_a1.0 <- with(data, weight_a0 * (plug_in_mu1m1.marg.m2a0.obs  - mu_a1_m1.0xm2.0))
  data$eif_m1_t3_a1.0 <- with(data, weight_a0 * ((plug_in_mu1m2.marg.m1a0.obs - mu_a1_m1.0xm2.0)))
  data$eif_m1_t4_a1.0 <- with(data, mu_a1_m1.0xm2.0)
    
  # Ref: A = 0
  data$ratio.m1.0 <- with(data, ifelse(gamma.a0.m.obs!= 0, gamma.a1.M2.obs / gamma.a0.m.obs, 0))
  data$eif_m1_t1_a0 <- with(data, weight_a0 * (gamma.a0.M1.obs - gamma.a1.M1.obs) * ratio.m1.0
                            * (Y  - plug_in_mu0m.obs))
  data$eif_m1_t2_a0 <- with(data,  weight_a0  * (plug_in_mu0m1.marg.m2a1.obs - mu_a0_m1.0xm2.1))
  data$eif_m1_t3_a0 <- with(data,  -weight_a1 * (plug_in_mu0m1.marg.m2a1.obs - mu_a0_m1.1xm2.1))
  data$eif_m1_t4_a0 <- with(data,  weight_a1  * ((plug_in_mu0m2.marg.m1a0.obs - plug_in_mu0m2.marg.m1a1.obs) -
                                                   (mu_a0_m1.0xm2.1 - mu_a0_m1.1xm2.1)))
  data$eif_m1_t5_a0 <- with(data, mu_a0_m1.0xm2.1 - mu_a0_m1.1xm2.1)
  
  # EIF: Interventional Indirect Effect: M2
  # Ref: A = 1
  data$ratio.m2.1 <- with(data, ifelse(gamma.a1.m.obs != 0, gamma.a1.M1.obs / gamma.a1.m.obs, 0))
  data$eif_m2_t1_a1 <- with(data, weight_a1 * (gamma.a1.M2.obs - gamma.a0.M2.obs) * ratio.m2.1
                            * (Y  - plug_in_mu1m.obs))
  data$eif_m2_t2_a1 <- with(data,  weight_a1 * (plug_in_mu1m2.marg.m1a1.obs - mu_a1_m1.1xm2.1))
  data$eif_m2_t3_a1 <- with(data,  -weight_a0 * (plug_in_mu1m2.marg.m1a1.obs - mu_a1_m1.1xm2.0))
  data$eif_m2_t4_a1 <- with(data,  weight_a1 * ((plug_in_mu1m1.marg.m2a1.obs - plug_in_mu1m1.marg.m2a0.obs) -
                                                   (mu_a1_m1.1xm2.1 - mu_a1_m1.1xm2.0)))
  data$eif_m2_t5_a1 <- with(data,  mu_a1_m1.1xm2.1 - mu_a1_m1.1xm2.0)
  
  # Ref: A = 0
  data$ratio.m2.0 <- with(data, ifelse(gamma.a0.m.obs != 0, gamma.a0.M1.obs / gamma.a0.m.obs, 0))
  data$eif_m2_t1_a0 <- with(data, weight_a0 * (gamma.a0.M2.obs - gamma.a1.M2.obs) * ratio.m2.0
                            * (Y  - plug_in_mu0m.obs))
  data$eif_m2_t2_a0 <- with(data,  weight_a0  * (plug_in_mu0m2.marg.m1a0.obs - mu_a0_m1.0xm2.0))
  data$eif_m2_t3_a0 <- with(data,  -weight_a1 * (plug_in_mu0m2.marg.m1a0.obs - mu_a0_m1.0xm2.1))
  data$eif_m2_t4_a0 <- with(data,  weight_a0  * ((plug_in_mu0m1.marg.m2a0.obs - plug_in_mu0m1.marg.m2a1.obs) -
                                                   (mu_a0_m1.0xm2.0 - mu_a0_m1.0xm2.1)))
  data$eif_m2_t5_a0 <- with(data,  mu_a0_m1.0xm2.0 - mu_a0_m1.0xm2.1)
  
  data <- data %>%
    replace_na(list(eif_m1_t1_a0 = 0, eif_m1_t2_a0 = 0, eif_m1_t3_a0 = 0,
                    eif_m1_t4_a0 = 0, eif_m1_t5_a0 = 0,
                    eif_m2_t1_a0 = 0, eif_m2_t2_a0 = 0, eif_m2_t3_a0 = 0,
                    eif_m2_t4_a0 = 0, eif_m2_t5_a0 = 0,
                    eif_m1_t1_a1 = 0, eif_m1_t2_a1 = 0, eif_m1_t3_a1 = 0,
                    eif_m1_t4_a1 = 0, eif_m1_t5_a1 = 0,
                    eif_m2_t1_a1 = 0, eif_m2_t2_a1 = 0, eif_m2_t3_a1 = 0,
                    eif_m2_t4_a1 = 0, eif_m2_t5_a1 = 0,
                    eif_m1_t1_a1.1 = 0, eif_m1_t2_a1.1 = 0, eif_m1_t3_a1.1 = 0,
                    eif_m1_t4_a1.1 = 0,
                    eif_m1_t1_a1.0 = 0, eif_m1_t2_a1.0 = 0, eif_m1_t3_a1.0 = 0,
                    eif_m1_t4_a1.0 = 0))
  
  if (stabilize == TRUE) {
    data$eif_m1_t1_a1 <- data$eif_m1_t1_a1 / with(data, mean(weight_a1 * ratio.m1.1))
    data$eif_m1_t2_a1 <- data$eif_m1_t2_a1 / mean(data$weight_a1)
    data$eif_m1_t3_a1 <- data$eif_m1_t3_a1 / mean(data$weight_a0)
    data$eif_m1_t4_a1 <- data$eif_m1_t4_a1 / mean(data$weight_a0)
    data$eif_m1_t1_a0 <- data$eif_m1_t1_a0 / with(data, mean(weight_a0 * ratio.m1.0))
    data$eif_m1_t2_a0 <- data$eif_m1_t2_a0 / mean(data$weight_a0)
    data$eif_m1_t3_a0 <- data$eif_m1_t3_a0 / mean(data$weight_a1)
    data$eif_m1_t4_a0 <- data$eif_m1_t4_a0 / mean(data$weight_a1)
    data$eif_m2_t1_a1 <- data$eif_m2_t1_a1 / with(data, mean(weight_a1 * ratio.m2.1))
    data$eif_m2_t2_a1 <- data$eif_m2_t2_a1 / mean(data$weight_a1)
    data$eif_m2_t3_a1 <- data$eif_m2_t3_a1 / mean(data$weight_a0)
    data$eif_m2_t4_a1 <- data$eif_m2_t4_a1 / mean(data$weight_a1)
    data$eif_m2_t1_a0 <- data$eif_m2_t1_a0 / with(data, mean(weight_a0 * ratio.m2.0))
    data$eif_m2_t2_a0 <- data$eif_m2_t2_a0 / mean(data$weight_a0)
    data$eif_m2_t3_a0 <- data$eif_m2_t3_a0 / mean(data$weight_a1)
    data$eif_m2_t4_a0 <- data$eif_m2_t4_a0 / mean(data$weight_a0)
  }
  
  data$eif_m1_a0 <- with(data, -1 * ((eif_m1_t1_a0 + eif_m1_t2_a0 + 
                                       eif_m1_t3_a0 + eif_m1_t4_a0 + eif_m1_t5_a0) * iweights +
                                       eif_m1_t5_a0 * iweights.p)) 
  
  data$eif_m2_a0 <- with(data, -1 * ((eif_m2_t1_a0 + eif_m2_t2_a0 + 
                                       eif_m2_t3_a0 + eif_m2_t4_a0 + eif_m2_t5_a0) * iweights +
                                       eif_m2_t5_a0 * iweights.p)) 
  
  data$eif_m1_a1 <- with(data, (eif_m1_t1_a1 + eif_m1_t2_a1 + 
                           eif_m1_t3_a1 + eif_m1_t4_a1 + eif_m1_t5_a1) * iweights + 
                           eif_m1_t5_a1 * iweights.p) 
  
  data$eif_m1_a1.1 <- with(data, (eif_m1_t1_a1.1 + eif_m1_t2_a1.1 + 
                                  eif_m1_t3_a1.1 + eif_m1_t4_a1.1) * iweights + 
                           eif_m1_t4_a1.1 * iweights.p)

  data$eif_m1_a1.0 <- with(data, (eif_m1_t1_a1.0 + eif_m1_t2_a1.0 + 
                                  eif_m1_t3_a1.0 + eif_m1_t4_a1.0) * iweights + 
                             eif_m1_t4_a1.0 * iweights.p)
  
  data$eif_m2_a1 <- with(data, (eif_m2_t1_a1 + eif_m2_t2_a1 + 
                           eif_m2_t3_a1 + eif_m2_t4_a1 + eif_m2_t5_a1) * iweights +
                           eif_m2_t5_a1 * iweights.p) 
  

  # EIF: Covariant effect
  data$eif_cov_a0 <- with(data, eif_te - eif_de_a0 - eif_m1_a0 - eif_m2_a0)
  data$eif_cov_a1 <- with(data, eif_te - eif_de_a1 - eif_m1_a1 - eif_m2_a1)

  # EIF: Joint Indirect effect
  data$eif_m1m2_a1 <- with(data, eif_cov_a1 + eif_m1_a1 + eif_m2_a1)
  data$eif_m1m2_a0 <- with(data, eif_cov_a0 + eif_m1_a0 + eif_m2_a0)

  return(data)
}

ProcessAllData <- function(data, marg_jd) {
  RelabOutcomes <- function(data) {
    xwalk <- data %>%
      select(M, M1, M2) %>%
      distinct() %>%
      arrange(M)
    
    a_vals <- sort(unique(data$A[data$cc == 1]))
    for (a in a_vals) {
      for(i in 1:nrow(xwalk)) {
        M <- xwalk$M[i]
        M1 <- xwalk$M1[i]
        M2 <- xwalk$M2[i]
        data[[paste0("jmu_a", a, "_m1", M1, "_m2", M2)]] <- data[[paste0("mu_a", a, "_m", M)]]
      }
    }
    return(data)
  }

  data %>%
    RelabOutcomes() %>%
    MarginalizeJointDensity(marginalize = marg_jd) %>%
    MarginalizeOutcomesJointDensityNat(0, 0) %>%
    MarginalizeOutcomesJointDensityNat(1, 1) %>%
    MarginalizeOutcomesJointDensityNat(0, 1) %>%
    MarginalizeOutcomesJointDensityNat(1, 0) %>%
    MarginalizeOutcomesJointDensityInt(1, 1, 1) %>%
    MarginalizeOutcomesJointDensityInt(0, 0, 0) %>%
    MarginalizeOutcomesJointDensityInt(1, 1, 0) %>%
    MarginalizeOutcomesJointDensityInt(0, 0, 1) %>%
    MarginalizeOutcomesJointDensityInt(1, 0, 0) %>%
    MarginalizeOutcomesJointDensityInt(0, 1, 1) %>%
    MarginalizeOutcomesJointDensityInt(0, 1, 0) %>%
    MarginalizeOutcomesJointDensityInt(1, 0, 1) %>%
    MarginalizeOutcomesMarginalDensity(0, 0, "M1") %>%
    MarginalizeOutcomesMarginalDensity(1, 1, "M1") %>%
    MarginalizeOutcomesMarginalDensity(0, 0, "M2") %>%
    MarginalizeOutcomesMarginalDensity(1, 1, "M2") %>%
    MarginalizeOutcomesMarginalDensity(0, 1, "M1") %>%
    MarginalizeOutcomesMarginalDensity(1, 0, "M1") %>%
    MarginalizeOutcomesMarginalDensity(0, 1, "M2") %>%
    MarginalizeOutcomesMarginalDensity(1, 0, "M2") %>%
    ReturnObservedM("M") %>%
    ReturnObservedM("M1") %>%
    ReturnObservedM("M2") %>%
    ReturnObservedY()
}

CalculateIDEtxfx <- function(data) {
  te <- mean(data$eif_te) / mean(data$eif_d) 
  se.te <- ProportionSE(data, te, "eif_te", "eif_d")
  
  de.a1 <- mean(data$eif_de_a1) / mean(data$eif_d)
  de.a0 <- mean(data$eif_de_a0) / mean(data$eif_d)
  se.de.a1 <- ProportionSE(data, de.a1, "eif_de_a1", "eif_d") 
  se.de.a0 <- ProportionSE(data, de.a0, "eif_de_a0", "eif_d")
  
  prop.de.a1 <- mean(data$eif_de_a1) / mean(data$eif_te)
  prop.de.a0 <- mean(data$eif_de_a0) / mean(data$eif_te)
  se.prop.de.a1 <- ProportionSE(data, prop.de.a1, "eif_de_a1", "eif_te")
  se.prop.de.a0 <- ProportionSE(data, prop.de.a0, "eif_de_a0", "eif_te")
  
  m1.a1 <- mean(data$eif_m1_a1) / mean(data$eif_d)
  m1.a0 <- mean(data$eif_m1_a0) / mean(data$eif_d)
  se.m1.a1 <- ProportionSE(data, m1.a1, "eif_m1_a1", "eif_d") 
  se.m1.a0 <- ProportionSE(data, m1.a0, "eif_m1_a0", "eif_d") 
  
  prop.m1.a1 <- mean(data$eif_m1_a1) / mean(data$eif_te)
  prop.m1.a0 <- mean(data$eif_m1_a0) / mean(data$eif_te)
  se.prop.m1.a1 <- ProportionSE(data, prop.m1.a1, "eif_m1_a1", "eif_te")
  se.prop.m1.a0 <- ProportionSE(data, prop.m1.a0, "eif_m1_a0", "eif_te")
  
  m2.a1 <- mean(data$eif_m2_a1) / mean(data$eif_d)
  m2.a0 <- mean(data$eif_m2_a0) / mean(data$eif_d)
  se.m2.a1 <- ProportionSE(data, m2.a1, "eif_m2_a1", "eif_d") 
  se.m2.a0 <- ProportionSE(data, m2.a0, "eif_m2_a0", "eif_d")
  
  prop.m2.a1 <- mean(data$eif_m2_a1) / mean(data$eif_te)
  prop.m2.a0 <- mean(data$eif_m2_a0) / mean(data$eif_te)
  se.prop.m2.a1 <- ProportionSE(data, prop.m2.a1, "eif_m2_a1", "eif_te")
  se.prop.m2.a0 <- ProportionSE(data, prop.m2.a0, "eif_m2_a0", "eif_te")
  
  cov.a1 <- mean(data$eif_cov_a1) / mean(data$eif_d)
  cov.a0 <- mean(data$eif_cov_a0) / mean(data$eif_d)
  se.cov.a1 <- ProportionSE(data, cov.a1, "eif_cov_a1", "eif_d")
  se.cov.a0 <- ProportionSE(data, cov.a0, "eif_cov_a0", "eif_d")
  
  prop.c.a1 <- mean(data$eif_cov_a1) / mean(data$eif_te)
  prop.c.a0 <- mean(data$eif_cov_a0) / mean(data$eif_te)
  se.prop.c.a1 <- ProportionSE(data, prop.c.a1, "eif_cov_a1", "eif_te")
  se.prop.c.a0 <- ProportionSE(data, prop.c.a0, "eif_cov_a0", "eif_te")
  
  ide.a1 <- mean(data$eif_m1m2_a1) / mean(data$eif_d)
  ide.a0 <- mean(data$eif_m1m2_a0) / mean(data$eif_d)
  se.ide.a1 <- ProportionSE(data, ide.a1, "eif_m1m2_a1", "eif_d")
  se.ide.a0 <- ProportionSE(data, ide.a0, "eif_m1m2_a0", "eif_d")
  
  prop.ide.a1 <- mean(data$eif_m1m2_a1) / mean(data$eif_te)
  prop.ide.a0 <- mean(data$eif_m1m2_a0) / mean(data$eif_te)
  se.prop.ide.a1 <- ProportionSE(data, prop.ide.a1, "eif_m1m2_a1", "eif_te")
  se.prop.ide.a0 <- ProportionSE(data, prop.ide.a0, "eif_m1m2_a0", "eif_te")
  
  tibble(
    ests = c(te, 
             de.a1, m1.a1, m2.a1, cov.a1, ide.a1,
             prop.de.a1, prop.m1.a1, prop.m2.a1, prop.c.a1, prop.ide.a1,
             de.a0, m1.a0, m2.a0, cov.a0, ide.a0,
             prop.de.a0, prop.m1.a0, prop.m2.a0, prop.c.a0, prop.ide.a0),
    se = c(se.te, 
           se.de.a1, se.m1.a1, se.m2.a1, se.ide.a1, se.cov.a1,
           se.prop.de.a1, se.prop.m1.a1, se.prop.m2.a1, se.prop.c.a1, se.prop.ide.a1,
           se.de.a0, se.m1.a0, se.m2.a0, se.cov.a0, se.ide.a0,
           se.prop.de.a0, se.prop.m1.a0, se.prop.m2.a0, se.prop.c.a0, se.prop.ide.a0),
    n = nrow(data),
    reference = c(NA, 0, rep(1, 10), rep(0, 9)),
    effect = c("TE", rep(c("IDE/NDE", "IIE - Isolated", "IIE - Health", "IIE - Covariant", "IIE/NIE"), 4)),
    type = c(rep("diff", 6), rep("prop", 5), rep("diff", 5), rep("prop", 5))
  ) %>%
    mutate(lci = ests - 1.96*se, uci = ests + 1.96*se)
}

ProportionSE <- function(data, proportion, numerator, denominator) {
  se <- abs(proportion) * sqrt((var(data[[numerator]])/mean(data[[numerator]])^2 + 
                           var(data[[denominator]])/mean(data[[denominator]])^2 - 
                           2*cov(data[[numerator]], data[[denominator]])/
                           (mean(data[[numerator]]*mean(data[[denominator]]))))/nrow(data))
  return(se)
}

CalculateXGNuisFuncs <- function(data, cc) {
  if (sum(grepl("^M1$", names(data))) == 1) {
    data <- data
  }
  if (sum(grepl("^M1$", names(data))) == 0) {
    data <- data %>%
      left_join(mxwalk, by = "M") 
  }
  
  mxwalk <- data %>%
    filter(cc == 1) %>%
    select(M, M1, M2) %>%
    distinct() %>%
    mutate_all(factor)
  
  data <- data %>%
    MarginalizeJointDensity() %>%
    MarginalizeOutcomesMarginalDensity(1, 1, "M1") %>%
    MarginalizeOutcomesMarginalDensity(0, 0, "M1") %>%
    MarginalizeOutcomesMarginalDensity(1, 1, "M2") %>%
    MarginalizeOutcomesMarginalDensity(0, 0, "M2") %>%
    ReturnObservedM("M") %>%
    ReturnObservedM("M1") %>%
    ReturnObservedM("M2") %>%
    ReturnObservedY() %>%
    TauXAA("M1", 1, 0) %>%
    TauXAA("M1", 0, 1) %>%
    TauXAA("M2", 1, 0) %>%
    TauXAA("M2", 0, 1) %>%
    TauXAA("M",  1, 0) %>%
    TauXAA("M",  0, 1) %>%
    CalculateNaturalEIFs(cc = cc) %>%
    CalculateEIFPropMed()
  return(data)
}

#' CreateJointData: processes mediator density models for interventional effects. this function
#' is essentially a byproduct of having previously calculated natural effects sequentially and
#' exists to repurpose this output to calculate the interventional effects, which we prefer for 
#' our final analysis
#'
#' @param glm_models list of glm models for pscore, outcome, mediator density
#' @param dataset specify which dataset to run join on: full (1), non-hesitant (2); SA (3); (4) NH + NA-Hesitant
#' @param mediator_specification multinomial or binomial specification of mediator
#'
#' @return dataset with all relevant functions
CreateJointData <- function(glm_models, dataset, mediator_specification = "Multinomial") {
  spec_vector <- c(2, 4, 6)
  
  if (mediator_specification == "Binomial") {
    spec_vector <- c(1, 3, 5)
  }
  
  joint_data <- glm_models[[spec_vector[3]]][[dataset]] # mediator is M1 x M2
  m1_data <- glm_models[[spec_vector[1]]][[dataset]] # mediator is isolation (M1)
  m2_data <- glm_models[[spec_vector[2]]][[dataset]] # mediator is health (M2)
  
  joint_data <- joint_data %>%
    mutate(M1 = m1_data$M, M2 = m2_data$M) %>%
    select(id, Y, A, M, M1, M2, cc, everything()) 
  
  return(joint_data)
}

#' CalculateAllIntFx: calculate interventional treatment effect estimates across
#' different data specifications with and without propensity score trimming/complete
#' case estimates
#'
#' @param joint_data output of CreateJointData, ProcessAllData, CalculateIDEEIFs
#' @param subgroup_effects discrete subgroups to calcuate effects within
#'
#' @return list of treatment effect estimates across all data: results include four
#' specifications: untrimmed, trimmed (0.001), complete-case, and complete-case trimmed 
#' estimates
CalculateAllIntFx <- function(joint_data, subgroup_effects, trim, s) {
  intfx_all <- CalculateTxfx(analytic_data, joint_data, subgroup_effects, 
                             mediation = "TRUE", cc = FALSE, trim = 0, s = 1e-10)

  intfx_all_trim <- CalculateTxfx(analytic_data, joint_data, subgroup_effects, 
                                  mediation = "TRUE", cc = FALSE, trim, s)
  
  intfx_all_cc <- CalculateTxfx(analytic_data, joint_data, subgroup_effects, 
                                mediation = "TRUE", cc = TRUE, 
                                trim = 0, s = 1e-10)
  
  intfx_all_trim_cc <- CalculateTxfx(analytic_data, joint_data, subgroup_effects, 
                                     mediation = "TRUE", cc = TRUE, trim, s)
  
  list(Full = intfx_all, Trimmed = intfx_all_trim, CC = intfx_all_cc, `CC-Trimmed` = intfx_all_trim_cc)
}


