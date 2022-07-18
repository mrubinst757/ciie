# program: simuation-funs.R
# purpose: functions for simulations of projection, DRL, and plug-in approaches to
# estimation of heterogeneous interventional indirect effects

library(tidyverse); expit <- function(x){ exp(x)/(1+exp(x)) }; logit <- function(x){ log(x/(1-x)) }

################################################################################
######################## SIMULATE POPULATION DATA ##############################
################################################################################

#' GeneratePopulation: creates population of data and associated nuisance functions
#'
#' @param N population size
#' @param outcome if "binary" makes all outcomes fall within 0 and 1
#'
#' @return dataframe containing population information
GeneratePopulation <- function(N, outcome = NULL) {
  
  pop    <- data.frame(cc = rep(1, N), R = rep(1, N), eta = 1, 
                       iweights = 1, iweights.p = 0)
  pop$X  <- rnorm(N, 1, sqrt(0.5))
  pop$X  <- ifelse(pop$X < -2, -2, pop$X)
  pop$X  <- ifelse(pop$X > 4, 4, pop$X)
  
  pop$pi1 <- with(pop, 0.2 * (X < -1))
  pop$pi2 <- with(pop, (0.2 + 0.55 * abs(X + 1)) * (X < 0 & X > -1))
  pop$pi3 <- with(pop, (0.75 - 0.25*X) * (X > 0 & X < 1))
  pop$pi4 <- with(pop, (0.5 - 0.25*(X - 1)^2) * (X > 1 & X < 2))
  pop$pi5 <- with(pop, (0.25 + 0.5 * (X - 2)) * (X > 2 & X < 3))
  pop$pi6 <- with(pop, 0.75 * (X > 3))
  pop$pi  <- with(pop, pi1 + pi2 + pi3 + pi4 + pi5 + pi6)
  
  pop$A  <- rbinom(N, 1, pop$pi)
  pop$U  <- rbinom(N, 1, 0.5)
  
  # Counterfactual mediator densities --------------------------------------------
  
  # density values setting U = u
  gamma.M1.a0.u0 <- with(pop, 0.15 + (X + 1) * 0.1)
  gamma.M1.a1.u0 <- with(pop, 0.55 + (X + 1) * 0.05)
  gamma.M1.a0.u1 <- with(pop, 0.1)
  gamma.M1.a1.u1 <- with(pop, 0.8)
  
  gamma.M2.a0.u1 <- with(pop, 0.15 + (X + 1) * 0.125)
  gamma.M2.a1.u1 <- with(pop, 0.4 + (X + 0.5) * 0.1)
  gamma.M2.a0.u0 <- with(pop, 0.1)
  gamma.M2.a1.u0 <- with(pop, 0.8)
  
  gamma.m1.a0.u0 <- (1 - gamma.M1.a0.u0) * (1 - gamma.M2.a0.u0)
  gamma.m2.a0.u0 <- gamma.M1.a0.u0 * (1 - gamma.M2.a0.u0)
  gamma.m3.a0.u0 <- (1 - gamma.M1.a0.u0) * gamma.M2.a0.u0
  gamma.m4.a0.u0 <- gamma.M1.a0.u0 * gamma.M2.a0.u0
  gamma.m1.a1.u0 <- (1 - gamma.M1.a1.u0) * (1 - gamma.M2.a1.u0)
  gamma.m2.a1.u0 <- gamma.M1.a1.u0 * (1 - gamma.M2.a1.u0)
  gamma.m3.a1.u0 <- (1 - gamma.M1.a1.u0) * gamma.M2.a1.u0
  gamma.m4.a1.u0 <- gamma.M1.a1.u0 * gamma.M2.a1.u0
  gamma.m1.a0.u1 <- (1 - gamma.M1.a0.u1) * (1 - gamma.M2.a0.u1)
  gamma.m2.a0.u1 <- gamma.M1.a0.u1 * (1 - gamma.M2.a0.u1)
  gamma.m3.a0.u1 <- (1 - gamma.M1.a0.u1) * gamma.M2.a0.u1
  gamma.m4.a0.u1 <- gamma.M1.a0.u1 * gamma.M2.a0.u1
  gamma.m1.a1.u1 <- (1 - gamma.M1.a1.u1) * (1 - gamma.M2.a1.u1)
  gamma.m2.a1.u1 <- gamma.M1.a1.u1 * (1 - gamma.M2.a1.u1)
  gamma.m3.a1.u1 <- (1 - gamma.M1.a1.u1) * gamma.M2.a1.u1
  gamma.m4.a1.u1 <- gamma.M1.a1.u1 * gamma.M2.a1.u1
  
  # Observed mediator values -----------------------------------------------------
  pop$M1.1 <- rbinom(N, 1, gamma.M1.a1.u1 * pop$U + gamma.M1.a1.u0 * (1 - pop$U))
  pop$M2.1 <- rbinom(N, 1, gamma.M2.a1.u1 * pop$U + gamma.M2.a1.u0 * (1 - pop$U))
  pop$M1.0 <- rbinom(N, 1, gamma.M1.a0.u1 * pop$U + gamma.M1.a0.u0 * (1 - pop$U))
  pop$M2.0 <- rbinom(N, 1, gamma.M2.a0.u1 * pop$U + gamma.M2.a0.u0 * (1 - pop$U))
  pop$M1 <- ifelse(pop$A == 1, pop$M1.1, pop$M1.0)
  pop$M2 <- ifelse(pop$A == 1, pop$M2.1, pop$M2.0)
  pop$M  <- with(pop, as.numeric(interaction(M1, M2)))
  
  # counterfactual mediator densities (marginalized over U)
  pop$gamma.M1.a0 <- gamma.M1.a0.u0 * 0.5 + gamma.M1.a0.u1 * 0.5
  pop$gamma.M1.a1 <- gamma.M1.a1.u0 * 0.5 + gamma.M1.a1.u1 * 0.5
  pop$gamma.M2.a0 <- gamma.M2.a0.u0 * 0.5 + gamma.M2.a0.u1 * 0.5
  pop$gamma.M2.a1 <- gamma.M2.a1.u0 * 0.5 + gamma.M2.a1.u1 * 0.5
  
  pop$gamma.M11.a1 <- pop$gamma.M1.a1
  pop$gamma.M11.a0 <- pop$gamma.M1.a0
  pop$gamma.M10.a1 <- 1 - pop$gamma.M1.a1
  pop$gamma.M10.a0 <- 1 - pop$gamma.M1.a0
  
  pop$gamma.M21.a1 <- pop$gamma.M2.a1
  pop$gamma.M21.a0 <- pop$gamma.M2.a0
  pop$gamma.M20.a1 <- 1 - pop$gamma.M2.a1
  pop$gamma.M20.a0 <- 1 - pop$gamma.M2.a0
  
  pop$gamma.m1.a1 <- gamma.m1.a1.u0 * 0.5 + gamma.m1.a1.u1 * 0.5
  pop$gamma.m2.a1 <- gamma.m2.a1.u0 * 0.5 + gamma.m2.a1.u1 * 0.5
  pop$gamma.m3.a1 <- gamma.m3.a1.u0 * 0.5 + gamma.m3.a1.u1 * 0.5
  pop$gamma.m4.a1 <- gamma.m4.a1.u0 * 0.5 + gamma.m4.a1.u1 * 0.5
  pop$gamma.m1.a0 <- gamma.m1.a0.u0 * 0.5 + gamma.m1.a0.u1 * 0.5
  pop$gamma.m2.a0 <- gamma.m2.a0.u0 * 0.5 + gamma.m2.a0.u1 * 0.5
  pop$gamma.m3.a0 <- gamma.m3.a0.u0 * 0.5 + gamma.m3.a0.u1 * 0.5
  pop$gamma.m4.a0 <- gamma.m4.a0.u0 * 0.5 + gamma.m4.a0.u1 * 0.5
  
  # Counterfactual outcome model -------------------------------------------------
  pop$mu1 <- with(pop, (X - X^2) * (X < - 0.5))
  pop$mu2 <- with(pop, (-2 +  X) * (X > -0.5) * (X < 0))
  pop$mu3 <- with(pop, (- 12 + 10*sin(X^2) + 10*cos(X^2)) * (X > 0) * (X < 1))
  pop$mu4 <- with(pop, (-12 + 10*sin(1) + 10*cos(1) - 5*(X - 1) - 5*(X - 1)^2) * (X > 1 & X < 1.5))
  pop$mu5 <- with(pop, (- 12 + 10*sin(1) + 10*cos(1) - 5*(0.5) - 5*(0.5)^2 + 
                          0.5 * (X - 1.5) - (X - 1.5)^2 + 3*(X - 1.5)^3) * (X > 1.5 & X < 2.5))
  pop$mu6 <- with(pop, (-4 + 2*X) * (X > 2.5))
  pop$mu <- with(pop, mu1 + mu2 + mu3 + mu4 + mu5 + mu6) 
  
  #pop$mu_a1_m4 <- 10 + pop$mu + 2 * (pop$X + pop$X^2)
  #pop$mu_a1_m3 <- 5 + pop$mu + 1 * (pop$X + pop$X^2)
  #pop$mu_a1_m2 <- 5 + pop$mu + 0.5 * (pop$X + pop$X^2)
  #pop$mu_a1_m1 <- pop$mu + 0 * (pop$X + pop$X^2)
  
  pop$mu_a1_m4 <- 10 + pop$mu + 2*pop$X + 0.5*pop$X^2 #0 * (pop$X + pop$X^2)
  pop$mu_a1_m3 <- 4 + pop$mu #0 * (pop$X + pop$X^2)
  pop$mu_a1_m2 <- 8 + pop$mu + 2*pop$X + 0.5*pop$X^2#0 * (pop$X + pop$X^2)
  pop$mu_a1_m1 <- pop$mu #0 * (pop$X + pop$X^2)
  pop$mu_a0_m4 <- 10 + pop$mu + 2*pop$X + 0.5*pop$X^2#0 * (pop$X + pop$X^2)
  pop$mu_a0_m3 <- 4 + pop$mu + 0*pop$X #0 * (pop$X + pop$X^2)
  pop$mu_a0_m2 <- 8 + pop$mu + 2*pop$X + 0.5*pop$X^2#0 * (pop$X + pop$X^2)
  pop$mu_a0_m1 <- pop$mu + 0*pop$X# 0 * (pop$X + pop$X^2)
  
  pop$mu_1 <- case_when(
    pop$M == 1 ~ pop$mu_a1_m1,
    pop$M == 2 ~ pop$mu_a1_m2,
    pop$M == 3 ~ pop$mu_a1_m3,
    pop$M == 4 ~ pop$mu_a1_m4
  )
  
  pop$mu_0 <- case_when(
    pop$M == 1 ~ pop$mu_a0_m1,
    pop$M == 2 ~ pop$mu_a0_m2,
    pop$M == 3 ~ pop$mu_a0_m3,
    pop$M == 4 ~ pop$mu_a0_m4
  )
  
  pop$mu_1 <- with(pop, mu_a1_m1 * gamma.m1.a1 + mu_a1_m2 * gamma.m2.a1 + 
                     mu_a1_m3 * gamma.m3.a1 + mu_a1_m4 * gamma.m4.a1)
  
  pop$mu_0 <- with(pop, mu_a0_m1 * gamma.m1.a0 + mu_a0_m2 * gamma.m2.a0 + 
                     mu_a0_m3 * gamma.m3.a0 + mu_a0_m4 * gamma.m4.a0)
  
  if (outcome == "binary") {
    max_val_a1 <- max(pop$mu_a1_m1, pop$mu_a1_m2, pop$mu_a1_m3, pop$mu_a1_m4)

    max_val_a0 <- max(pop$mu_a0_m1, pop$mu_a0_m2, pop$mu_a0_m3, pop$mu_a0_m4)
    
    min_val_a1 <- min(pop$mu_a1_m1, pop$mu_a1_m2, pop$mu_a1_m3, pop$mu_a1_m4)
    
    min_val_a0 <- min(pop$mu_a0_m1, pop$mu_a0_m2, pop$mu_a0_m3, pop$mu_a0_m4)

    pop <- mutate_at(pop, vars(matches("mu_a0_m[1-4]")), ~(1 / (max_val_a0 + 10 - min_val_a0 + 10)) * (. - min_val_a0 + 10))
    pop <- mutate_at(pop, vars(matches("mu_a1_m[1-4]")), ~(1 / (max_val_a1 + 10 - min_val_a1 + 10)) * (. - min_val_a1 + 10))

    pop$mu_1 <- with(pop, mu_a1_m1 * gamma.m1.a1 + mu_a1_m2 * gamma.m2.a1 + 
                       mu_a1_m3 * gamma.m3.a1 + mu_a1_m4 * gamma.m4.a1)
    
    pop$mu_0 <- with(pop, mu_a0_m1 * gamma.m1.a0 + mu_a0_m2 * gamma.m2.a0 + 
                       mu_a0_m3 * gamma.m3.a0 + mu_a0_m4 * gamma.m4.a0)
    
  }
  
  pop$Y111 <- rbinom(nrow(pop), 1, pop$mu_a1_m4)
  pop$Y101 <- rbinom(nrow(pop), 1, pop$mu_a1_m3)
  pop$Y110 <- rbinom(nrow(pop), 1, pop$mu_a1_m2)
  pop$Y100 <- rbinom(nrow(pop), 1, pop$mu_a1_m1)
  pop$Y011 <- rbinom(nrow(pop), 1, pop$mu_a0_m4)
  pop$Y001 <- rbinom(nrow(pop), 1, pop$mu_a0_m3)
  pop$Y010 <- rbinom(nrow(pop), 1, pop$mu_a0_m2)
  pop$Y000 <- rbinom(nrow(pop), 1, pop$mu_a0_m1)

  pop$Y14 <- pop$Y111; pop$Y04 <- pop$Y011
  pop$Y13 <- pop$Y101; pop$Y03 <- pop$Y001
  pop$Y12 <- pop$Y110; pop$Y02 <- pop$Y010
  pop$Y11 <- pop$Y100; pop$Y01 <- pop$Y000
  
  pop <- mutate(pop, Y1 = case_when(
    M1.1 == 0 & M2.1 == 0 ~ Y11,
    M1.1 == 1 & M2.1 == 0 ~ Y12,
    M1.1 == 0 & M2.1 == 1 ~ Y13,
    M1.1 == 1 & M2.1 == 1 ~ Y14
  ))
  
  pop <- mutate(pop, Y0 = case_when(
    M1.0 == 0 & M2.0 == 0 ~ Y01,
    M1.0 == 1 & M2.0 == 0 ~ Y02,
    M1.0 == 0 & M2.0 == 1 ~ Y03,
    M1.0 == 1 & M2.0 == 1 ~ Y04
  ))
  
  #pop$Y1 <- rbinom(nrow(pop), 1, pop$mu_1)
  #pop$Y0 <- rbinom(nrow(pop), 1, pop$mu_0)
  
  pop <- mutate(pop, Y = ifelse(A == 1, Y1, Y0))
  
  return(pop)
}

################################################################################
################## METHOD ONE: ESTIMATION VIA SUPERLEARNER #####################
################################################################################

#' EstimateCATEFuns: nuisance estimation associated with estimating interventional effects via M1
#'
#' @param folds training folds for nuisance estimation; 6 folds if independent nuisance estimation learning and
#' only one fold otherwise (the former is to train the DRL with theoretic guarantees on convergence)
#' @param sl.libs libraries for superlearner
#' @param sample_split option whether to learn each nuisance function on a separate sample or to learn
#' all on one sample split. the former is useful for the DR learner; the latter for a projection estimator
#' @param return_data option of whether to return the test data associated with nuisance estimation
#'
#' @return a list containing projection estimators; DRL models; plug0in models, (and the test data if specified)

EstimateCATEFuns <- function(folds, sl.libs, sample_split, return_data = FALSE) {
  
  if (sample_split == FALSE) {
    test_dat <- folds[[2]]
    folds[[2]] <- folds[[1]]
    folds[[3]] <- folds[[1]]
    folds[[4]] <- folds[[1]]
    folds[[5]] <- folds[[1]]
    folds[[6]] <- test_dat
  }
  
  # Learn E(Y | X, M1, M2, A = 1)
  X <- folds[[1]][,c("X", "A", "M1", "M2")]
  Y <- folds[[1]]$Y
  y1_mod <- SuperLearner(Y[X$A==1], 
                         X[X$A==1, -grep("A", names(X))], 
                         SL.library = sl.libs,
                         family = binomial())

  # Learn E(Y | X, A = 1)
  X <- folds[[1]][,c("X", "A")]
  y1_mod.te <- SuperLearner(Y[X$A==1], 
                         data.frame(X = X[X$A == 1, -grep("A", names(X))]), 
                         SL.library = sl.libs,
                         family = binomial())

  # Learn E(Y | X, A = 0)
  y0_mod.te <- SuperLearner(Y[X$A==0], 
                         data.frame(X = X[X$A == 0, -grep("A", names(X))]), 
                         SL.library = sl.libs,
                         family = binomial())
  
  # Learn P(A = 1 | X)
  A <- folds[[2]]$A
  X <- folds[[2]]$X
  pi_mod <- SuperLearner(A, 
                         tibble(X = X), 
                         SL.library = sl.libs,
                         family = binomial())
  
  
  # Learn P(M1 = 1 | X, A = 1); P(M1 = 1 | X, A = 0)
  M1 <- folds[[3]]$M1
  X  <- folds[[3]][,c("X", "A")]
  
  m1.1_mod <- SuperLearner(M1[X$A == 1], 
                           X = data.frame(X = X[X$A == 1, -grep("A", names(X))]), 
                           SL.library = sl.libs,
                           family = binomial())
  
  m1.0_mod <- SuperLearner(M1[X$A == 0], 
                           X = data.frame(X = X[X$A == 0, -grep("A", names(X))]), 
                           SL.library = sl.libs,
                           family = binomial())
  
  # Learn P(M2 = 1 | X, A = 0)
  M2 <- folds[[4]]$M2
  X  <- folds[[4]][,c("X", "A")]
  
  m2.0_mod <- SuperLearner(M2[X$A == 0], 
                           X = data.frame(X = X[X$A == 0, -grep("A", names(X))]), 
                           SL.library = sl.libs,
                           family = binomial())
  
  # Learn P(M1 = m1, M2 = m2 | X, A = 1) 
  # Factorize as: P(M1 = m1 | X, M2 = m2, A = 1) * P(M2 = m2 | X, A = 1)
  M1  <- folds[[5]]$M1
  M2  <- folds[[5]]$M2
  X1  <- data.frame(folds[[5]][,c("X", "A", "M2")])
  X2  <- data.frame(folds[[5]][,c("X", "A")])
  
  m1_1m2.mod <- SuperLearner(M1[X1$A == 1], 
                             X = X1[X1$A == 1, -grep("A", names(X1))], 
                             SL.library = sl.libs,
                             family = binomial())
  
  m2_1.mod <- SuperLearner(M2[X2$A == 1], 
                           X = data.frame(X = X2[X2$A == 1, -grep("A", names(X2))]), 
                           SL.library = sl.libs,
                           family = binomial())
  
  test_dat <- folds[[6]]
  test_XM1M2 <- select(test_dat, X, M1, M2)
  test_XM2 <- select(test_dat, X, M2)
  test_X <- select(test_dat, X)
  
  test_dat$pi <- predict(pi_mod, test_dat)$pred
  test_dat$mu_a1_m1 <- predict(y1_mod, mutate(test_XM1M2, M1 = 0, M2 = 0))$pred
  test_dat$mu_a1_m2 <- predict(y1_mod, mutate(test_XM1M2, M1 = 1, M2 = 0))$pred
  test_dat$mu_a1_m3 <- predict(y1_mod, mutate(test_XM1M2, M1 = 0, M2 = 1))$pred
  test_dat$mu_a1_m4 <- predict(y1_mod, mutate(test_XM1M2, M1 = 1, M2 = 1))$pred

  test_dat$mu_1 <- predict(y1_mod.te, test_X)$pred
  test_dat$mu_0 <- predict(y0_mod.te, test_X)$pred

  test_dat$gamma.m1.a1 <- (1 - predict(m1_1m2.mod, mutate(test_XM2, M2 = 0))$pred) * (1 - predict(m2_1.mod, test_X)$pred)
  test_dat$gamma.m2.a1 <- predict(m1_1m2.mod, mutate(test_XM2, M2 = 0))$pred * (1 - predict(m2_1.mod, test_X)$pred)
  test_dat$gamma.m3.a1 <- (1 - predict(m1_1m2.mod, mutate(test_XM2, M2 = 1))$pred) * predict(m2_1.mod, test_X)$pred
  test_dat$gamma.m4.a1 <- predict(m1_1m2.mod, mutate(test_XM2, M2 = 1))$pred * predict(m2_1.mod, test_X)$pred
  
  test_dat$gamma.M11.a1 <- predict(m1.1_mod, test_X)$pred
  test_dat$gamma.M10.a1 <- 1 - predict(m1.1_mod, test_X)$pred
  test_dat$gamma.M11.a0 <- predict(m1.0_mod, test_X)$pred
  test_dat$gamma.M10.a0 <- 1 - predict(m1.0_mod, test_X)$pred
  
  test_dat$gamma.M21.a0 <- predict(m2.0_mod, test_X)$pred
  test_dat$gamma.M20.a0 <- 1 - predict(m2.0_mod, test_X)$pred
  test_dat$gamma.M21.a1 <- predict(m2_1.mod, test_X)$pred
  test_dat$gamma.M20.a1 <- 1 - predict(m2_1.mod, test_X)$pred
  
  mjd <- ifelse(sample_split == TRUE, FALSE, TRUE)
  
  test_dat <- CalculateIDEEIFs(ProcessAllData(test_dat, marg_jd = mjd))
  test_dat$eif_prop <- with(test_dat, eif_m1_a1 / eif_te_2 - 
                              eif_te * (eif_m1_t5_a1 / eif_te_2^2) + 
                              eif_m1_t5_a1 / eif_te_2)

  test_dat$prop <- with(test_dat, eif_m1_t5_a1 / eif_te_2)
  
  Xmat <- cbind(1, test_dat$X)
  Xmat2 <- cbind(1, test_dat$X, test_dat$X^2)
  
  XXi  <- solve(t(Xmat) %*% Xmat)
  XXi2  <- solve(t(Xmat2) %*% Xmat2)
  
  # Run Projection estimator
  model1 <- lm(eif_m1_a1 ~ X, test_dat)
  model2 <- lm(eif_m1_a1 ~ poly(X, 2, raw = TRUE), test_dat)
  model3 <- lm(eif_m1_t5_a1 ~ X, test_dat)
  model4 <- lm(eif_m1_t5_a1 ~ poly(X, 2, raw = TRUE), test_dat)
  
  model5 <- lm(eif_te ~ X, test_dat)
  model6 <- lm(eif_te ~ poly(X, 2, raw = TRUE), test_dat)
  model7 <- lm(eif_te_2 ~ X, test_dat)
  model8 <- lm(eif_te_2 ~ poly(X, 2, raw = TRUE), test_dat)
  
  model9  <- lm(eif_prop ~ X, test_dat)
  model10 <- lm(eif_prop ~ poly(X, 2, raw = TRUE), test_dat)
  model11 <- lm(prop ~ X, test_dat)
  model12 <- lm(prop ~ poly(X, 2, raw = TRUE), test_dat)
    
  varc1  <- sandwich::vcovHC(model1, type = "HC1")
  varc2  <- sandwich::vcovHC(model2, type = "HC1")
  varc3  <- sandwich::vcovHC(model3, type = "HC1")
  varc4  <- sandwich::vcovHC(model4, type = "HC1")

  varc5  <- sandwich::vcovHC(model5, type = "HC1")
  varc6  <- sandwich::vcovHC(model6, type = "HC1")
  varc7  <- sandwich::vcovHC(model7, type = "HC1")
  varc8  <- sandwich::vcovHC(model8, type = "HC1")

  varc9  <- sandwich::vcovHC(model9, type = "HC1")
  varc10 <- sandwich::vcovHC(model10, type = "HC1")
  varc11 <- sandwich::vcovHC(model11, type = "HC1")
  varc12  <- sandwich::vcovHC(model12, type = "HC1")
  
  avgest1 <- coef(model1); avgest2 <- coef(model2); avgest3 <- coef(model3); avgest4 <- coef(model4)
  avgest5 <- coef(model5); avgest6 <- coef(model6); avgest7 <- coef(model7); avgest8 <- coef(model8)
  avgest9 <- coef(model9); avgest10 <- coef(model10); avgest11 <- coef(model11); avgest12 <- coef(model12)
  
  resid1 <- resid(model1); resid2 <- resid(model2); resid3 <- resid(model3); resid4 <- resid(model4)
  resid5 <- resid(model5); resid6 <- resid(model6); resid7 <- resid(model7); resid8 <- resid(model8)

  res.projection <- list(model = list(avgest1, avgest2, avgest3, avgest4, 
                                      avgest5, avgest6, avgest7, avgest8,
                                      avgest9, avgest10, avgest11, avgest12), 
                         vcov  = list(varc1, varc2, varc3, varc4,
                                      varc5, varc6, varc7, varc8,
                                      varc9, varc10, varc11, varc12),
                         resid = list(resid1, resid2, resid3, resid4,
                                      resid5, resid6, resid7, resid8),
                         XXi  = list(XXi, XXi2),
                         Xmat = list(Xmat, Xmat2))
  
  # Run DR Learner
  drl_model    <- smooth.spline(test_dat$X, test_dat$eif_m1_a1)
  drl_model_te <- smooth.spline(test_dat$X, test_dat$eif_te)
  drl_model_prop <- smooth.spline(test_dat$X, test_dat$eif_prop)
  
  plugin_model <- list(mu.mod = y1_mod, 
                       m1.1.mod = m1.1_mod,
                       m1.0.mod = m1.0_mod, 
                       m2.0.mod = m2.0_mod,
                       mu1.mod.te = y1_mod.te,
                       mu0.mod.te = y0_mod.te)
  
  statistics <- tibble(
    max_pi1_inv = max(test_dat$A / test_dat$pi),
    max_pi0_inv = max((1 - test_dat$A) / (1-test_dat$pi)),
    max_gamma_pi_inv = max(test_dat$A / (test_dat$pi * test_dat$gamma.a1.m.obs)),
    max_gamma_pi_prop_inv = max(test_dat$A / (test_dat$pi * test_dat$gamma.a1.m.obs * test_dat$eif_te_2)),
    max_pi1_prop_inv = max(test_dat$A / (test_dat$pi * (test_dat$eif_te_2^2))),
    max_pi0_prop_inv = max((1 - test_dat$A) / ((1 - test_dat$pi) * (test_dat$eif_te_2^2))),
    max_te_inv  = max(1/test_dat$eif_te_2),
    max_te_inv2 = max(1/(test_dat$eif_te_2^2))
  )
  
  if (return_data == TRUE) {
    res <- list(proj = res.projection, 
                drl_model = list(drl_model = drl_model, 
                                 drl_model_te = drl_model_te,
                                 drl_model_prop = drl_model_prop), 
                plugin = plugin_model, 
                data = test_dat)
  }
  if (return_data == FALSE) {
    res <- list(proj = res.projection, 
                drl_model = list(drl_model = drl_model, 
                                 drl_model_te = drl_model_te,
                                 drl_model_prop = drl_model_prop), 
                plugin = plugin_model,
                stats = statistics)
  }
  return(res)
}

#' ProcessResults: takes output from EstimateCATEFuns and returns point estimates and
#' variance estimates. this function is used in the final step of RunExperiments
#'
#' @param result_list output from EstimateCATEFuns
#' @param point point to predict X = x
#' @param test_data test set independent of nuisance estimation data
#' @param truth_vector vector of true population values
#'
#' @return dataframe containing results
ProcessResults <- function(result_list, point, test_data, truth_vector) {
  #res1 <- result_list[[1]]
  res1 <- result_list[[1]] #non-parametric nuisance estimation
  res2 <- result_list[[2]] #parametric nuisance estimation
  stats_np  <- result_list[[1]]$stats
  stats_glm <- result_list[[2]]$stats
  
  # IIE via M1
  #drl.pred.ss  <- predict(res1$drl_model$drl_model, point)$y
  drl.pred.nss <- predict(res1$drl_model$drl_model, point)$y
  #plugin.pred.ss  <- PredictPlugIn(res1$plugin, point) # add plug-in total effect?
  plugin.pred.nss <- PredictPlugIn(res1$plugin, point)
  
  proj.pred.lin.glm <- c(1, point$X) %*% res2$proj$model[[1]] 
  proj.pred.quad.glm <- c(1, point$X, point$X^2) %*% res2$proj$model[[2]]
  proj.pred.lin.np <- c(1, point$X) %*% res1$proj$model[[1]] 
  proj.pred.quad.np <- c(1, point$X, point$X^2) %*% res1$proj$model[[2]]
  
  proj.pred.lin.glm.plugin <- c(1, point$X) %*% res2$proj$model[[3]]
  proj.pred.quad.glm.plugin <- c(1, point$X, point$X^2) %*% res2$proj$model[[4]]
  proj.pred.lin.np.plugin <- c(1, point$X) %*% res1$proj$model[[3]] 
  proj.pred.quad.np.plugin <- c(1, point$X, point$X^2) %*% res1$proj$model[[4]]
  
  #drl.pred.ss.var <- CalcSSVar(res1$drl_model$drl_model, test_data$X, point$X)
  drl.pred.nss.var <- CalcSSVar(res1$drl_model$drl_model, test_data$X, point$X)
  #plugin.pred.ss.var <- NA
  plugin.pred.nss.var <- NA
  proj.lin.glm.var <- c(1, point$X) %*% res2$proj$vcov[[1]] %*% c(1, point$X)
  proj.quad.glm.var <- c(1, point$X, point$X^2) %*% res2$proj$vcov[[2]] %*% c(1, point$X, point$X^2)
  proj.lin.np.var <- c(1, point$X) %*% res1$proj$vcov[[1]] %*% c(1, point$X)
  proj.quad.np.var <- c(1, point$X, point$X^2) %*% res1$proj$vcov[[2]] %*% c(1, point$X, point$X^2)
  
  proj.lin.glm.var.pi <- c(1, point$X) %*% res2$proj$vcov[[3]] %*% c(1, point$X)
  proj.quad.glm.var.pi <- c(1, point$X, point$X^2) %*% res2$proj$vcov[[4]] %*% c(1, point$X, point$X^2)
  proj.lin.np.var.pi <- c(1, point$X) %*% res1$proj$vcov[[3]] %*% c(1, point$X)
  proj.quad.np.var.pi <- c(1, point$X, point$X^2) %*% res1$proj$vcov[[4]] %*% c(1, point$X, point$X^2)
  
  res_med <- tibble(
    ests = c(proj.pred.lin.glm, proj.pred.quad.glm, 
             proj.pred.lin.np, proj.pred.quad.np,
             #drl.pred.ss, 
             drl.pred.nss, 
             #plugin.pred.ss, 
             plugin.pred.nss,
             proj.pred.lin.glm.plugin, proj.pred.quad.glm.plugin, 
             proj.pred.lin.np.plugin, proj.pred.quad.np.plugin),
    var = c(proj.lin.glm.var, proj.quad.glm.var, 
            proj.lin.np.var, proj.quad.np.var, 
            #drl.pred.ss.var, 
            drl.pred.nss.var,
            #plugin.pred.ss.var, 
            plugin.pred.nss.var,
            proj.lin.glm.var.pi, proj.quad.glm.var.pi, 
            proj.lin.np.var.pi, proj.quad.np.var.pi)) %>%
    unnest() %>%
    mutate(
      estimand = "IIE-M1",
      lci = ests - 1.96*sqrt(var), uci = ests + 1.96*sqrt(var),
      type = c("Projection-Linear", "Projection-Quad", 
               "Projection-Linear", "Projection-Quad",
               #"DRLearner", "DRLearner", 
               #"Plugin-DRL", "Plugin-DRL", 
               "DRLearner", "Plugin-DRL",
               "Proj-Lin-PI", "Proj-Quad-PI", 
               "Proj-Lin-PI", "Proj-Quad-PI"),
      nuisance_est = c("GLM", "GLM", "NP", "NP", "NP", 
                       "NP", "GLM", "GLM", "NP", "NP"),
      #sample_split = c(NA, NA, NA, NA, "Yes", "No", "Yes", "No", NA, NA, NA, NA),
      truth = truth_vector$mediation)
  
  # Total effect
  #drl.pred.ss.te  <- predict(res1$drl_model$drl_model_te, point)$y
  drl.pred.nss.te <- predict(res1$drl_model$drl_model_te, point)$y
  #plugin.pred.ss.te  <- PredictPlugIn(res1$plugin, point, "cate") # add plug-in total effect?
  plugin.pred.nss.te <- PredictPlugIn(res1$plugin, point, "cate")
  
  proj.pred.lin.glm.te <- c(1, point$X) %*% res2$proj$model[[5]] 
  proj.pred.quad.glm.te <- c(1, point$X, point$X^2) %*% res2$proj$model[[6]]
  proj.pred.lin.np.te <- c(1, point$X) %*% res1$proj$model[[5]] 
  proj.pred.quad.np.te <- c(1, point$X, point$X^2) %*% res1$proj$model[[6]]
  
  proj.pred.lin.glm.plugin.te <- c(1, point$X) %*% res2$proj$model[[7]]
  proj.pred.quad.glm.plugin.te <- c(1, point$X, point$X^2) %*% res2$proj$model[[8]]
  proj.pred.lin.np.plugin.te <- c(1, point$X) %*% res1$proj$model[[7]] 
  proj.pred.quad.np.plugin.te <- c(1, point$X, point$X^2) %*% res1$proj$model[[8]]
  
  #drl.pred.ss.var.te <- CalcSSVar(res1$drl_model$drl_model_te, test_data$X, point$X)
  drl.pred.nss.var.te <- CalcSSVar(res1$drl_model$drl_model_te, test_data$X, point$X)
  #plugin.pred.ss.var.te <- NA
  plugin.pred.nss.var.te <- NA
  
  proj.lin.glm.var.te <- c(1, point$X) %*% res2$proj$vcov[[5]] %*% c(1, point$X)
  proj.quad.glm.var.te <- c(1, point$X, point$X^2) %*% res2$proj$vcov[[6]] %*% c(1, point$X, point$X^2)
  proj.lin.np.var.te <- c(1, point$X) %*% res1$proj$vcov[[5]] %*% c(1, point$X)
  proj.quad.np.var.te <- c(1, point$X, point$X^2) %*% res1$proj$vcov[[6]] %*% c(1, point$X, point$X^2)
  
  proj.lin.glm.var.pi.te <- c(1, point$X) %*% res2$proj$vcov[[7]] %*% c(1, point$X)
  proj.quad.glm.var.pi.te <- c(1, point$X, point$X^2) %*% res2$proj$vcov[[8]] %*% c(1, point$X, point$X^2)
  proj.lin.np.var.pi.te <- c(1, point$X) %*% res1$proj$vcov[[7]] %*% c(1, point$X)
  proj.quad.np.var.pi.te <- c(1, point$X, point$X^2) %*% res1$proj$vcov[[8]] %*% c(1, point$X, point$X^2)
  
  res_te <- tibble(
    ests = c(proj.pred.lin.glm.te, proj.pred.quad.glm.te, 
             proj.pred.lin.np.te, proj.pred.quad.np.te,
             #drl.pred.ss.te, 
             drl.pred.nss.te,
             #plugin.pred.ss.te, 
             plugin.pred.nss.te,
             proj.pred.lin.glm.plugin.te, proj.pred.quad.glm.plugin.te, 
             proj.pred.lin.np.plugin.te, proj.pred.quad.np.plugin.te),
    var = c(proj.lin.glm.var.te, proj.quad.glm.var.te, 
            proj.lin.np.var.te, proj.quad.np.var.te, 
            #drl.pred.ss.var.te, 
            drl.pred.nss.var.te,
            #plugin.pred.ss.var.te, 
            plugin.pred.nss.var.te,
            proj.lin.glm.var.pi.te, proj.quad.glm.var.pi.te, 
            proj.lin.np.var.pi.te, proj.quad.np.var.pi.te)) %>%
    unnest() %>%
    mutate(
      lci = ests - 1.96*sqrt(var), uci = ests + 1.96*sqrt(var),
      type = c("Projection-Linear", "Projection-Quad", 
               "Projection-Linear", "Projection-Quad",
               #"DRLearner", "DRLearner", 
               #"Plugin-DRL", "Plugin-DRL",
               "DRLearner", "Plugin-DRL",
               "Proj-Lin-PI", "Proj-Quad-PI", 
               "Proj-Lin-PI", "Proj-Quad-PI"),
      estimand = "Total effect",
      nuisance_est = c("GLM", "GLM", "NP", "NP", "NP",
                       "NP", "GLM", "GLM", "NP", "NP"),
      #sample_split = c(NA, NA, NA, NA, "Yes", "No", "Yes", "No", NA, NA, NA, NA),
      truth = truth_vector$total_effect)

  # Proportion mediated
  
  #drl.pred.ss.prop  <- predict(res1$drl_model$drl_model_prop, point)$y
  drl.pred.nss.prop <- predict(res1$drl_model$drl_model_prop, point)$y
  #plugin.pred.ss.prop  <- PredictPlugIn(res1$plugin, point, "prop.med") 
  plugin.pred.nss.prop <- PredictPlugIn(res1$plugin, point, "prop.med")
  
  proj.pred.lin.glm.prop <- c(1, point$X) %*% res2$proj$model[[9]] 
  proj.pred.quad.glm.prop <- c(1, point$X, point$X^2) %*% res2$proj$model[[10]]
  proj.pred.lin.np.prop <- c(1, point$X) %*% res1$proj$model[[9]] 
  proj.pred.quad.np.prop <- c(1, point$X, point$X^2) %*% res1$proj$model[[10]]
  
  proj.pred.lin.glm.plugin.prop <- c(1, point$X) %*% res2$proj$model[[11]]
  proj.pred.quad.glm.plugin.prop <- c(1, point$X, point$X^2) %*% res2$proj$model[[12]]
  proj.pred.lin.np.plugin.prop <- c(1, point$X) %*% res1$proj$model[[11]] 
  proj.pred.quad.np.plugin.prop <- c(1, point$X, point$X^2) %*% res1$proj$model[[12]]
  
  #drl.pred.ss.var.prop <- CalcSSVar(res1$drl_model$drl_model_prop, test_data$X, point$X)
  drl.pred.nss.var.prop <- CalcSSVar(res1$drl_model$drl_model_prop, test_data$X, point$X)
  #plugin.pred.ss.var.prop <- NA
  plugin.pred.nss.var.prop <- NA
  
  proj.lin.glm.var.prop <- c(1, point$X) %*% res2$proj$vcov[[9]] %*% c(1, point$X)
  proj.quad.glm.var.prop <- c(1, point$X, point$X^2) %*% res2$proj$vcov[[10]] %*% c(1, point$X, point$X^2)
  proj.lin.np.var.prop <- c(1, point$X) %*% res1$proj$vcov[[9]] %*% c(1, point$X)
  proj.quad.np.var.prop <- c(1, point$X, point$X^2) %*% res1$proj$vcov[[10]] %*% c(1, point$X, point$X^2)
  
  proj.lin.glm.var.pi.prop <- c(1, point$X) %*% res2$proj$vcov[[11]] %*% c(1, point$X)
  proj.quad.glm.var.pi.prop <- c(1, point$X, point$X^2) %*% res2$proj$vcov[[12]] %*% c(1, point$X, point$X^2)
  proj.lin.np.var.pi.prop <- c(1, point$X) %*% res1$proj$vcov[[11]] %*% c(1, point$X)
  proj.quad.np.var.pi.prop <- c(1, point$X, point$X^2) %*% res1$proj$vcov[[12]] %*% c(1, point$X, point$X^2)
  
  res_prop <- tibble(
    ests = c(proj.pred.lin.glm.prop, proj.pred.quad.glm.prop, 
             proj.pred.lin.np.prop, proj.pred.quad.np.prop,
             #drl.pred.ss.prop, 
             drl.pred.nss.prop,
             #plugin.pred.ss.prop, 
             plugin.pred.nss.prop,
             proj.pred.lin.glm.plugin.prop, proj.pred.quad.glm.plugin.prop, 
             proj.pred.lin.np.plugin.prop, proj.pred.quad.np.plugin.prop),
    var = c(proj.lin.glm.var.prop, proj.quad.glm.var.prop, 
            proj.lin.np.var.prop, proj.quad.np.var.prop, 
            #drl.pred.ss.var.prop, 
            drl.pred.nss.var.prop,
            #plugin.pred.ss.var.prop, 
            plugin.pred.nss.var.prop,
            proj.lin.glm.var.pi.prop, proj.quad.glm.var.pi.prop, 
            proj.lin.np.var.pi.prop, proj.quad.np.var.pi.prop)) %>%
    unnest() %>%
    mutate(
      lci = ests - 1.96*sqrt(var), uci = ests + 1.96*sqrt(var),
      type = c("Projection-Linear", "Projection-Quad", 
               "Projection-Linear", "Projection-Quad",
               #"DRLearner", "DRLearner", 
               #"Plugin-DRL", "Plugin-DRL",
               "DRLearner", "Plugin-DRL",
               "Proj-Lin-PI", "Proj-Quad-PI", 
               "Proj-Lin-PI", "Proj-Quad-PI"),
      estimand = "Proportion mediated",
      nuisance_est = c("GLM", "GLM", "NP", "NP", "NP",
                       "NP", "GLM", "GLM", "NP", "NP"),
      #sample_split = c(NA, NA, NA, NA, "Yes", "No", "Yes", "No", NA, NA, NA, NA),
      truth = truth_vector$proportion)
  
  # Proportion mediated -- alternate
  X <- res1$proj$Xmat[[1]]
  XXi <- res1$proj$XXi[[1]]

  X2 <- res1$proj$Xmat[[2]]
  XXi2 <- res1$proj$XXi[[2]]
  
  covmat.np.dr.lin <- XXi %*% t(X) %*% diag(res1$proj$resid[[1]]*res1$proj$resid[[5]]) %*% X %*% XXi
  covmat.np.dr.qad <- XXi2 %*% t(X2) %*% diag(res1$proj$resid[[2]]*res1$proj$resid[[6]]) %*% X2 %*% XXi2
  covmat.np.pi.lin <- XXi %*% t(X) %*% diag(res1$proj$resid[[3]]*res1$proj$resid[[7]]) %*% X %*% XXi
  covmat.np.pi.qad <- XXi2 %*% t(X2) %*% diag(res1$proj$resid[[4]]*res1$proj$resid[[8]]) %*% X2 %*% XXi2

  covmat.glm.dr.lin <- XXi %*% t(X) %*% diag(res2$proj$resid[[1]]*res2$proj$resid[[5]]) %*% X %*% XXi
  covmat.glm.dr.qad <- XXi2 %*% t(X2) %*% diag(res2$proj$resid[[2]]*res2$proj$resid[[6]]) %*% X2 %*% XXi2
  covmat.glm.pi.lin <- XXi %*% t(X) %*% diag(res2$proj$resid[[3]]*res2$proj$resid[[7]]) %*% X %*% XXi
  covmat.glm.pi.qad <- XXi2 %*% t(X2) %*% diag(res2$proj$resid[[4]]*res2$proj$resid[[8]]) %*% X2 %*% XXi2
  
  proj.lin.np.dr.cov  <- c(1, point$X) %*% covmat.np.dr.lin %*% c(1, point$X)
  proj.lin.np.pi.cov  <- c(1, point$X) %*% covmat.np.pi.lin %*% c(1, point$X)
  proj.lin.glm.dr.cov <- c(1, point$X) %*% covmat.glm.dr.lin %*% c(1, point$X)
  proj.lin.glm.pi.cov <- c(1, point$X) %*% covmat.glm.pi.lin %*% c(1, point$X)

  proj.qad.np.dr.cov <- c(1, point$X, point$X^2) %*% covmat.np.dr.qad %*% c(1, point$X, point$X^2)
  proj.qad.np.pi.cov <- c(1, point$X, point$X^2) %*% covmat.np.pi.qad %*% c(1, point$X, point$X^2)
  proj.qad.glm.dr.cov <- c(1, point$X, point$X^2) %*% covmat.glm.dr.qad %*% c(1, point$X, point$X^2)
  proj.qad.glm.pi.cov <- c(1, point$X, point$X^2) %*% covmat.glm.pi.qad %*% c(1, point$X, point$X^2)
  
  #drl.pred.ss.prop2  <- drl.pred.ss / drl.pred.ss.te  
  drl.pred.nss.prop2 <- drl.pred.nss / drl.pred.nss.te
  #plugin.pred.ss.prop2  <- plugin.pred.ss / plugin.pred.ss.te
  plugin.pred.nss.prop2 <- plugin.pred.nss / plugin.pred.nss.te
  
  proj.pred.lin.glm.prop2  <- proj.pred.lin.glm / proj.pred.lin.glm.te
  proj.pred.quad.glm.prop2 <- proj.pred.quad.glm / proj.pred.quad.glm.te 
  proj.pred.lin.np.prop2   <-  proj.pred.lin.np / proj.pred.lin.np.te
  proj.pred.quad.np.prop2  <- proj.pred.quad.np / proj.pred.quad.np.te
  
  proj.pred.lin.glm.plugin.prop2  <- proj.pred.lin.glm.plugin / proj.pred.lin.glm.plugin.te
  proj.pred.quad.glm.plugin.prop2 <- proj.pred.quad.glm.plugin / proj.pred.quad.glm.plugin.te
  proj.pred.lin.np.plugin.prop2   <- proj.pred.lin.np.plugin / proj.pred.lin.np.plugin.te  
  proj.pred.quad.np.plugin.prop2  <- proj.pred.quad.np.plugin / proj.pred.quad.np.plugin.te
  
  #drl.pred.ss.var.prop2 <- drl.pred.ss.prop2^2 * (drl.pred.ss.var / drl.pred.ss^2 + drl.pred.ss.var.te / drl.pred.ss.te^2)
  drl.pred.nss.var.prop2 <- drl.pred.nss.prop2^2 * (drl.pred.nss.var / drl.pred.nss^2 + drl.pred.nss.var.te / drl.pred.nss.te^2)
  #plugin.pred.ss.var.prop2 <- NA
  plugin.pred.nss.var.prop2 <- NA
  
  proj.lin.glm.var.prop2  <- proj.pred.lin.glm.prop2^2 * (proj.lin.glm.var / proj.pred.lin.glm^2 + 
                                                            proj.lin.glm.var.te / proj.pred.lin.glm.te^2 -
                                                            2 * proj.lin.glm.dr.cov / (proj.pred.lin.glm * proj.pred.lin.glm.te))
  proj.quad.glm.var.prop2 <- proj.pred.quad.glm.prop2^2 * (proj.quad.glm.var / proj.pred.quad.glm^2 + 
                                                             proj.quad.glm.var.te / proj.pred.quad.glm.te^2 -
                                                             2 * proj.qad.glm.dr.cov / (proj.pred.quad.glm*proj.pred.quad.glm.te))
  proj.lin.np.var.prop2   <- proj.pred.lin.np.prop2^2 * (proj.lin.np.var / proj.pred.lin.np^2 + 
                                                           proj.lin.np.var.te / proj.pred.lin.np.te^2 -
                                                           2 * proj.lin.np.dr.cov / (proj.pred.lin.np * proj.pred.lin.np.te))
  proj.quad.np.var.prop2  <- proj.pred.quad.np.prop2^2 * (proj.quad.np.var / proj.pred.quad.np^2 + 
                                                            proj.quad.np.var.te / proj.pred.quad.np.te^2 -
                                                            2 * proj.qad.np.dr.cov / (proj.pred.quad.np * proj.pred.quad.np.te))
  
  proj.lin.glm.var.pi.prop2  <- proj.pred.lin.glm.plugin.prop2^2 * (proj.lin.glm.var.pi / proj.pred.lin.glm.plugin^2 + 
                                                                      proj.lin.glm.var.pi.te / proj.pred.lin.glm.plugin.te^2 -
                                                                      2 * proj.lin.glm.pi.cov / (proj.pred.lin.glm.plugin * proj.pred.lin.glm.plugin.te))
  proj.quad.glm.var.pi.prop2 <- proj.pred.quad.glm.plugin.prop2^2 * (proj.quad.glm.var.pi / proj.pred.quad.glm.plugin^2 + 
                                                                       proj.quad.glm.var.pi.te / proj.pred.quad.glm.plugin.te^2 -
                                                                       2 * proj.qad.glm.pi.cov / (proj.pred.quad.glm.plugin * proj.pred.quad.glm.plugin.te)) 
  proj.lin.np.var.pi.prop2   <- proj.pred.lin.np.plugin.prop2^2 * (proj.lin.np.var.pi / proj.pred.lin.np.plugin^2 + 
                                                                     proj.lin.np.var.pi.te / proj.pred.lin.np.plugin.te^2 -
                                                                     2 * proj.lin.np.pi.cov / (proj.pred.lin.np.plugin * proj.pred.lin.np.plugin.te))
  proj.quad.np.var.pi.prop2  <- proj.pred.quad.np.plugin.prop2^2 * (proj.quad.np.var.pi / proj.pred.quad.np.plugin^2 + 
                                                                      proj.quad.np.var.pi.te / proj.pred.quad.np.plugin.te^2 -
                                                                      2 * proj.qad.np.pi.cov / (proj.pred.quad.np.plugin * proj.pred.quad.np.plugin.te)) 
  
  res_prop2 <- tibble(
    ests = c(proj.pred.lin.glm.prop2, proj.pred.quad.glm.prop2, 
             proj.pred.lin.np.prop2, proj.pred.quad.np.prop2,
             #drl.pred.ss.prop2, 
             drl.pred.nss.prop2,
             #plugin.pred.ss.prop2, 
             plugin.pred.nss.prop2,
             proj.pred.lin.glm.plugin.prop2, proj.pred.quad.glm.plugin.prop2, 
             proj.pred.lin.np.plugin.prop2, proj.pred.quad.np.plugin.prop2),
    var = c(proj.lin.glm.var.prop2, proj.quad.glm.var.prop2, 
            proj.lin.np.var.prop2, proj.quad.np.var.prop2, 
            #drl.pred.ss.var.prop2, 
            drl.pred.nss.var.prop2,
            #plugin.pred.ss.var.prop2, 
            plugin.pred.nss.var.prop2,
            proj.lin.glm.var.pi.prop2, proj.quad.glm.var.pi.prop2, 
            proj.lin.np.var.pi.prop2, proj.quad.np.var.pi.prop2)) %>%
    unnest() %>%
    mutate(
      lci = ests - 1.96*sqrt(var), uci = ests + 1.96*sqrt(var),
      type = c("Projection-Linear", "Projection-Quad", 
               "Projection-Linear", "Projection-Quad",
               #"DRLearner", "DRLearner", 
               #"Plugin-DRL", "Plugin-DRL",
               "DRLearner", "Plugin-DRL",
               "Proj-Lin-PI", "Proj-Quad-PI", 
               "Proj-Lin-PI", "Proj-Quad-PI"),
      estimand = "Proportion mediated alt",
      nuisance_est = c("GLM", "GLM", "NP", "NP", "NP",
                       "NP", "GLM", "GLM", "NP", "NP"),
      #sample_split = c(NA, NA, NA, NA, "Yes", "No", "Yes", "No", NA, NA, NA, NA),
      truth = truth_vector$proportion)
  
  
  res_main <- bind_rows(res_med, res_te, res_prop, res_prop2)
  res_stat <- list(stats_glm = stats_glm, stats_np = stats_np)
  res <- list(res_main = res_main, res_stat = res_stat)
    
  return(res)
}

#' RunExperiments: iterate EstimateCATEFuns a specified number of times over repeated sampling
#' from a population; returns simulation results
#'
#' @param nsims number of simulations to run
#' @param population dataframe containing a full population of observations to sample from
#' @param sample_size sample size for test set
#' @param point_value predict at X = x
#' @param sl.libs SuperLearner libraries to use
#'
#' @return dataframe containing bias, MSE, and variance associated with estimates
#' averaged over all nsims
RunExperiments <- function(nsims, population, sample_size, point_value, sl.libs) {
  truth_vector <- CalcTruth(population, point_value)
  res <- list()
  for (i in 1:nsims) {
    data_sample <- sample_n(population, sample_size * 6)
    
    folds <- ntile(runif(sample_size * 2, 0, 1), 6)
    folds <- map(1:2, ~folds == .x)
    folds <- map(folds, ~data_sample[.x, ]) 
    
    res1 <- EstimateCATEFuns(folds, sl.libs, sample_split = FALSE)
    res2 <- EstimateCATEFuns(folds, sl.libs = "SL.glm", sample_split = FALSE)
    
    res_list <- list(res1, res2)
    res[[i]] <- ProcessResults(res_list, point_value, folds[[2]], truth_vector)
  }
  return(res)
}

#' CalcSSVar: calculate estimated variance at X = x from smoothing spline model
#'
#' @param model smoothing spline model
#' @param x complete vector of X values in teh data
#' @param point_value point to predict variance at
#'
#' @return variance estimate of smoothing spline at X = x
CalcSSVar <- function(model, x, point_value) {
  yi <- c(1, rep(0, length(x)))
  w <- predict(smooth.spline(c(point_value, x), yi, lambda = model$lambda), x)$y
  w <- w / sum(w)
  varcomp <- (w * resid(model))^2
  varcomp <- varcomp[!is.na(varcomp)]
  varest <- sum(varcomp)
  return(varest)
}

#' PredictPlugIn: generates estimates for non-parametric plug-in estimator for 
#' mediation functionals (via M1)
#'
#' @param models list of nuisance functions
#' @param point point to predict
#'
#' @return plug-in point estimate at point X = x
PredictPlugIn <- function(models, point, estimand = "ciie") {
  
  if (estimand == "ciie") {
    mu_a1_m1 <- predict(models$mu.mod, mutate(point, M1 = 0, M2 = 0))$pred
    mu_a1_m2 <- predict(models$mu.mod, mutate(point, M1 = 1, M2 = 0))$pred
    mu_a1_m3 <- predict(models$mu.mod, mutate(point, M1 = 0, M2 = 1))$pred
    mu_a1_m4 <- predict(models$mu.mod, mutate(point, M1 = 1, M2 = 1))$pred
    
    pm.1.1 <- (1 - predict(models$m1.1.mod, point)$pred) * (1 - predict(models$m2.0.mod, point)$pred)
    pm.2.1 <- predict(models$m1.1.mod, point)$pred * (1 - predict(models$m2.0.mod, point)$pred)
    pm.3.1 <- (1 - predict(models$m1.1.mod, point)$pred) * predict(models$m2.0.mod, point)$pred
    pm.4.1 <- predict(models$m1.1.mod, point)$pred * predict(models$m2.0.mod, point)$pred
    
    pm.1.0 <- (1 - predict(models$m1.0.mod, point)$pred) * (1 - predict(models$m2.0.mod, point)$pred)
    pm.2.0 <- predict(models$m1.0.mod, point)$pred * (1 - predict(models$m2.0.mod, point)$pred)
    pm.3.0 <- (1 - predict(models$m1.0.mod, point)$pred) * predict(models$m2.0.mod, point)$pred
    pm.4.0 <- predict(models$m1.0.mod, point)$pred * predict(models$m2.0.mod, point)$pred
    
    res.1 <- mu_a1_m1 * pm.1.1 + mu_a1_m2 * pm.2.1 + mu_a1_m3 * pm.3.1 + mu_a1_m4 * pm.4.1
    res.0 <- mu_a1_m1 * pm.1.0 + mu_a1_m2 * pm.2.0 + mu_a1_m3 * pm.3.0 + mu_a1_m4 * pm.4.0
    res <- res.1 - res.0
  }
  
  if (estimand == "cate") {
    res.1 <- predict(models$mu1.mod.te, mutate(point, A = 1))$pred
    res.0 <- predict(models$mu0.mod.te, mutate(point, A = 0))$pred
    res <- res.1 - res.0
  }
  
  if (estimand == "prop.med") {
    mu_a1_m1 <- predict(models$mu.mod, mutate(point, M1 = 0, M2 = 0))$pred
    mu_a1_m2 <- predict(models$mu.mod, mutate(point, M1 = 1, M2 = 0))$pred
    mu_a1_m3 <- predict(models$mu.mod, mutate(point, M1 = 0, M2 = 1))$pred
    mu_a1_m4 <- predict(models$mu.mod, mutate(point, M1 = 1, M2 = 1))$pred
    
    pm.1.1 <- (1 - predict(models$m1.1.mod, point)$pred) * (1 - predict(models$m2.0.mod, point)$pred)
    pm.2.1 <- predict(models$m1.1.mod, point)$pred * (1 - predict(models$m2.0.mod, point)$pred)
    pm.3.1 <- (1 - predict(models$m1.1.mod, point)$pred) * predict(models$m2.0.mod, point)$pred
    pm.4.1 <- predict(models$m1.1.mod, point)$pred * predict(models$m2.0.mod, point)$pred
    
    pm.1.0 <- (1 - predict(models$m1.0.mod, point)$pred) * (1 - predict(models$m2.0.mod, point)$pred)
    pm.2.0 <- predict(models$m1.0.mod, point)$pred * (1 - predict(models$m2.0.mod, point)$pred)
    pm.3.0 <- (1 - predict(models$m1.0.mod, point)$pred) * predict(models$m2.0.mod, point)$pred
    pm.4.0 <- predict(models$m1.0.mod, point)$pred * predict(models$m2.0.mod, point)$pred
    
    res.1.m <- mu_a1_m1 * pm.1.1 + mu_a1_m2 * pm.2.1 + mu_a1_m3 * pm.3.1 + mu_a1_m4 * pm.4.1
    res.0.m <- mu_a1_m1 * pm.1.0 + mu_a1_m2 * pm.2.0 + mu_a1_m3 * pm.3.0 + mu_a1_m4 * pm.4.0
  
    res.1.t <- predict(models$mu1.mod.te, mutate(point, A = 1))$pred
    res.0.t <- predict(models$mu0.mod.te, mutate(point, A = 0))$pred
    res <- (res.1.m - res.0.m) / (res.1.t - res.0.t)
  }
  return(res)
}

#' CalcTruth: calculates true population parameters
#'
#' @param population population of data
#' @param point point X = x
#'
#' @return returns true population values of mediation function at X = x and
#' the true best linear / quadratic projections
CalcTruth <- function(population, point) {
  population$plugin.1 <- with(population, mu_a1_m4 * gamma.M1.a1 * gamma.M2.a0 +
                                mu_a1_m3 * (1 - gamma.M1.a1) * gamma.M2.a0 +
                                mu_a1_m2 * gamma.M1.a1 * (1 - gamma.M2.a0) +
                                mu_a1_m1 * (1 - gamma.M1.a1) * (1 - gamma.M2.a0))
  
  population$plugin.0 <- with(population, mu_a1_m4 * gamma.M1.a0 * gamma.M2.a0 +
                                mu_a1_m3 * (1 - gamma.M1.a0) * gamma.M2.a0 +
                                mu_a1_m2 * gamma.M1.a0 * (1 - gamma.M2.a0) +
                                mu_a1_m1 * (1 - gamma.M1.a0) * (1 - gamma.M2.a0))
  
  population$plugin <- with(population, plugin.1 - plugin.0)
  population$plugin.te <- with(population, mu_1 - mu_0)
  
  truth <- mean(population$plugin[population$X > point$X - 0.0001 & 
                                  population$X < point$X + 0.0001])
  
  truth.te <- mean(population$plugin.te[population$X > point$X - 0.0001 & 
                                        population$X < point$X + 0.0001])
  
  truth.prop <- truth / truth.te
  
  lin.proj.truth  <- lm(plugin ~ X, population)
  quad.proj.truth <- lm(plugin ~ poly(X, 2, raw = TRUE), population)
  lin.proj.truth.te  <- lm(plugin.te ~ X, population)
  quad.proj.truth.te <- lm(plugin.te ~ poly(X, 2, raw = TRUE), population)
  
  point.truth.lin <- c(1, point$X) %*% coef(lin.proj.truth)
  point.truth.qad <- c(1, point$X, point$X^2) %*% coef(quad.proj.truth)
  point.truth.lin.te <- c(1, point$X) %*% coef(lin.proj.truth.te)
  point.truth.qad.te <- c(1, point$X, point$X^2) %*% coef(quad.proj.truth.te)
  point.truth.lin.prop <- point.truth.lin / point.truth.lin.te
  point.truth.qad.prop <- point.truth.qad / point.truth.qad.te
  
  mediation <- c(point.truth.lin, point.truth.qad, 
           point.truth.lin, point.truth.qad, 
           truth, truth, #truth, truth, 
           point.truth.lin, point.truth.qad, 
           point.truth.lin, point.truth.qad)

  total_effect <- c(point.truth.lin.te, point.truth.qad.te, 
                 point.truth.lin.te, point.truth.qad.te, 
                 truth.te, truth.te, #truth.te, truth.te, 
                 point.truth.lin.te, point.truth.qad.te, 
                 point.truth.lin.te, point.truth.qad.te)
  
  proportion <- c(point.truth.lin.prop, point.truth.qad.prop, 
                  point.truth.lin.prop, point.truth.qad.prop, 
                  truth.prop, truth.prop, #truth.prop, truth.prop, 
                  point.truth.lin.prop, point.truth.qad.prop, 
                  point.truth.lin.prop, point.truth.qad.prop)
  
  res <- list(mediation = mediation, total_effect = total_effect, proportion = proportion)
  
  return(res)
}

################################################################################
############# METHOD TWO: SIMULATE NUISANCE ESTIMATION ERROR ###################
################################################################################


#' EstimateDRL: simulate nuisance estimation of DRL by specifying rate of convergence
#' as guassian noise for each nuisance function
#'
#' @param population population
#' @param sample_size sample size
#' @param mu_rate RMSE of outcome model
#' @param pi_rate RMSE of p-score model
#' @param m_rate RMSE of joint-mediator probability
#' @param m1_rate RMSE of M1-probability
#' @param m2_rate RMSE of M2-probability
#' @param bin_constant constant to multiply estimation by
#' @param return_data option whether to return the data or not
#'
#' @return list of results
EstimateDRL <- function(population, sample_size, mu_rate, pi_rate, m_rate, m1_rate, m2_rate,
                         bin_constant = 1, mu_constant = 8, return_data = FALSE) {
  n <- sample_size
  sample_rows <- sample(1:nrow(population), 2*n, replace = FALSE)
  samp   <- population[sample_rows,]
  psuedo <- CalculateIDEEIFs(ProcessAllData(samp, marg_jd = TRUE))[, c("eif_m1_a1", "eif_m1_t1_a1", "eif_m1_t2_a1", "eif_m1_t3_a1", 
                                                                        "eif_m1_t4_a1", "eif_m1_t5_a1")]
  
  mu.sigma <- matrix(rep(0.15/n^(2*mu_rate), 16), 4, 4)
  diag(mu.sigma) <- 1/n^(2*mu_rate)
  mu.sigma <- mu_constant * mu.sigma
  
  m.sigma <- bin_constant * matrix(rep(0.15/n^(2*m_rate), 16), 4, 4)
  diag(m.sigma) <- bin_constant/n^(2*m_rate)
  
  m1.sigma <- bin_constant * matrix(rep(0.15/n^(2*m1_rate), 4), 2, 2)
  diag(m1.sigma) <- bin_constant/n^(2*m1_rate)
  
  dat <- tibble(ones = rep(1, 2*n))
  dat$pi <- with(samp, expit(logit(pi) + rnorm(2*n, 
                                               mean = sqrt(bin_constant)/n^pi_rate, 
                                               sd   = sqrt(bin_constant)/n^pi_rate)))
  
  mu.ests <- expit(logit(samp[, c("mu_a1_m1", "mu_a1_m2", "mu_a1_m3", "mu_a1_m4")]) + MASS::mvrnorm(2*n, rep(sqrt(bin_constant)/n^mu_rate, 4), mu.sigma))

  m.ests <- expit(logit(samp[, c("gamma.m1.a1", "gamma.m2.a1", "gamma.m3.a1", "gamma.m4.a1")]) + MASS::mvrnorm(2*n, rep(sqrt(bin_constant)/n^m_rate, 4), m.sigma))
  m.ests <- m.ests / rowSums(m.ests)
  
  m1.ests <- expit(logit(samp[, c("gamma.M1.a1", "gamma.M1.a0")]) + MASS::mvrnorm(2*n, rep(sqrt(bin_constant)/n^m1_rate, 2), m1.sigma))
  
  dat$gamma.M2.a0 <- with(samp, expit(logit(gamma.M2.a0) + rnorm(2*n, 
                                                                 mean = sqrt(bin_constant) / n^m2_rate, 
                                                                 sd   = sqrt(bin_constant) / n^m2_rate)))
  
  dat <- cbind(dat, mu.ests, m.ests, m1.ests, 
               tibble(X = samp$X), tibble(psuedo = psuedo$eif_m1_a1))
  
  dat <- cbind(samp[, !(names(samp) %in% names(dat))], dat)
  
  dat$gamma.M11.a1 <- dat$gamma.M1.a1
  dat$gamma.M11.a0 <- dat$gamma.M1.a0
  dat$gamma.M10.a1 <- 1 - dat$gamma.M1.a1
  dat$gamma.M10.a0 <- 1 - dat$gamma.M1.a0
  
  dat$gamma.M21.a1 <- dat$gamma.M2.a1
  dat$gamma.M21.a0 <- dat$gamma.M2.a0
  dat$gamma.M20.a1 <- 1 - dat$gamma.M2.a1
  dat$gamma.M20.a0 <- 1 - dat$gamma.M2.a0
  
  dat$plugin.1 <- with(dat, mu_a1_m4 * gamma.M1.a1 * gamma.M2.a0 +
                         mu_a1_m3 * (1 - gamma.M1.a1) * gamma.M2.a0 +
                         mu_a1_m2 * gamma.M1.a1 * (1 - gamma.M2.a0) +
                         mu_a1_m1 * (1 - gamma.M1.a1) * (1 - gamma.M2.a0))
  
  dat$plugin.0 <- with(dat, mu_a1_m4 * gamma.M1.a0 * gamma.M2.a0 +
                         mu_a1_m3 * (1 - gamma.M1.a0) * gamma.M2.a0 +
                         mu_a1_m2 * gamma.M1.a0 * (1 - gamma.M2.a0) +
                         mu_a1_m1 * (1 - gamma.M1.a0) * (1 - gamma.M2.a0)) 
  
  samp$plugin.1 <- with(samp, mu_a1_m4 * gamma.M1.a1 * gamma.M2.a0 +
                          mu_a1_m3 * (1 - gamma.M1.a1) * gamma.M2.a0 +
                          mu_a1_m2 * gamma.M1.a1 * (1 - gamma.M2.a0) +
                          mu_a1_m1 * (1 - gamma.M1.a1) * (1 - gamma.M2.a0)) + rnorm(2*n, 1/n^0.5, 1/n^0.5)
  
  samp$plugin.0 <- with(samp, mu_a1_m4 * gamma.M1.a0 * gamma.M2.a0 +
                          mu_a1_m3 * (1 - gamma.M1.a0) * gamma.M2.a0 +
                          mu_a1_m2 * gamma.M1.a0 * (1 - gamma.M2.a0) +
                          mu_a1_m1 * (1 - gamma.M1.a0) * (1 - gamma.M2.a0)) + rnorm(2*n, 1/n^0.5, 1/n^0.5)
  
  dat <- CalculateIDEEIFs(ProcessAllData(dat, marg_jd = FALSE))
  S <- sample(c(rep(1, n), rep(2, n)), 2*n, replace = FALSE)
  plugin.1 <- dat$plugin.1[S==2] 
  plugin.0 <- dat$plugin.0[S==2] 
  plugin.1.s <- samp$plugin.1[S==2] 
  plugin.0.s <- samp$plugin.0[S==2] 
  
  plugin  <- plugin.1 - plugin.0
  plugin.s  <- plugin.1.s - plugin.0.s
  drl     <- predict(smooth.spline(dat$X[S == 1], dat$eif_m1_a1[S == 1]), dat$X)$y[S == 2]
  oracle  <- predict(smooth.spline(dat$X[S == 1], dat$psuedo[S == 1]), dat$X)$y[S == 2]
  tau <- psuedo$eif_m1_t5_a1[S == 2]
  
  res.df <- c("drl"     = mean((drl - tau)^2), 
              "oracle"  = mean((oracle - tau)^2),
              "plugin"  = mean((plugin - tau)^2),
              "plugin.s" = mean((plugin.s - tau)^2))
  
  if (return_data == TRUE) {
    res <- list(res = res.df, cdat = tibble(plugin = plugin, 
                                            drl = drl, 
                                            oracle = oracle, 
                                            tau = tau,
                                            X = dat$X[S==2]))
  }
  if (return_data == FALSE) {
    res <- list(res = res.df)
  }
  return(res)
}


#' DRLSimulation: iterates EstimateDRL for a sequence of fixed rates for all nuisance estimation
#'
#' @param population population of data
#' @param sample_size sample_size to draw from population
#' @param rate_seq sequence of convergence rates
#'
#' @return output of EstimateDRL iterated over a rate sequence
DRLSimulation <- function(population, sample_size, rate_seq) {
  res <- map(rate_seq, ~EstimateDRL(population, sample_size, .x, .x, .x, .x, .x)) %>%
    invoke(rbind, .) %>%
    as_tibble() %>%
    mutate(rate = rate_seq)
  return(res)
}

#' DRLConvergenceTest: compares DRL to a plug-in approach across different sample sizes
#'
#' @param population population of data to sample from
#' @param nsims number of simulations
#' @param mu_rate RMSE of outcome model
#' @param pi_rate RMSE of propensity score model
#' @param m_rate RMSE of joint mediator probability model
#' @param m1_rate RMSE of M1 model
#' @param m2_rate RMSE of M2 model
#'
#' @return list of output from EstimateDRL
DRLConvergenceTest <- function(population, nsims, mu_rate, pi_rate, m_rate, m1_rate, m2_rate) {
  test1 <- map(1:nsims, ~EstimateDRL(population, 500, mu_rate = mu_rate, pi_rate = pi_rate, m_rate = m_rate, 
                                      m1_rate = m1_rate, m2_rate = m2_rate,
                                      bin_constant = 15, return_data = FALSE))
  test2 <- map(1:nsims, ~EstimateDRL(population, 1000, mu_rate = mu_rate, pi_rate = pi_rate, m_rate = m_rate, 
                                      m1_rate = m1_rate, m2_rate = m2_rate, 
                                      bin_constant = 15, return_data = FALSE))
  test3 <- map(1:nsims, ~EstimateDRL(population, 5000, mu_rate = mu_rate, pi_rate = pi_rate, m_rate = m_rate, 
                                      m1_rate = m1_rate, m2_rate = m2_rate, 
                                      bin_constant = 15, return_data = FALSE))
  test4 <- map(1:nsims, ~EstimateDRL(population, 10000, mu_rate = mu_rate, pi_rate = pi_rate, m_rate = m_rate, 
                                      m1_rate = m1_rate, m2_rate = m2_rate, 
                                      bin_constant = 15, return_data = FALSE))
  test5 <- map(1:nsims, ~EstimateDRL(population, 50000, mu_rate = mu_rate, pi_rate = pi_rate, m_rate = m_rate, 
                                      m1_rate = m1_rate, m2_rate = m2_rate, 
                                      bin_constant = 15, return_data = FALSE))
  test6 <- map(1:nsims, ~EstimateDRL(population, 100000, mu_rate = mu_rate, pi_rate = pi_rate, m_rate = m_rate, 
                                      m1_rate = m1_rate, m2_rate = m2_rate, 
                                      bin_constant = 15, return_data = FALSE))
  res <- list(test1, test2, test3, test4, test5, test6)
  return(res)
}

