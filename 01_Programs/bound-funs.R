
#' AddConfounding: adds confounding on risk-ratio scale to true regression function
#'
#' @param data dataset containing true parameter values
#' @param scale controls the amount of confounding: scale = 1 roughly implies the minimum ratio of unobserved
#' counterfactuals to biased estimate is 0.97. To be precise, we have that
#' min E[Y(m) | A, M != m]/E[Y(m) | A, M = m] = 0.97 and
#' min (1-E[Y(m) | A, X, M != m])/(1-E[Y(m) | A, X, M = m]) =  0.97. 
#' Scale = 2 would, for example, change this number equal to 0.94. We also define taustar as 1 minus this number
#' so that for scale = 1 it is 0.03.
#'
#' @return list containing dataset with confounded vector and the 
AddConfounding <- function(data, scale) {
  
  sfun <- scale * ((data$X < -1) * 0.01 + (data$X > -1 & data$X < 0) * 0.02 +
                     (data$X > 0 & data$X < 1) * 0.03 + (data$X > 1 & data$X < 2) * 0.02 +
                     (data$X > 2 & data$X < 3) * 0.01 + (data$X > 3) * 0.03)
  
  data$mu_a1_m1.cu <- data$mu_a1_m1 / (1 - (1 - data$gamma.m1.a1) * sfun)
  data$mu_a1_m2.cu <- data$mu_a1_m2 / (1 - (1 - data$gamma.m2.a1) * sfun)
  data$mu_a1_m3.cu <- data$mu_a1_m3 / (1 - (1 - data$gamma.m3.a1) * sfun)
  data$mu_a1_m4.cu <- data$mu_a1_m4 / (1 - (1 - data$gamma.m4.a1) * sfun)
  
  data$mu_a1_m1.cl <- (data$mu_a1_m1 - sfun * (1 - data$gamma.m1.a1)) / (1 - (1 - data$gamma.m1.a1) * sfun)
  data$mu_a1_m2.cl <- (data$mu_a1_m2 - sfun * (1 - data$gamma.m2.a1)) / (1 - (1 - data$gamma.m2.a1) * sfun)
  data$mu_a1_m3.cl <- (data$mu_a1_m3 - sfun * (1 - data$gamma.m3.a1)) / (1 - (1 - data$gamma.m3.a1) * sfun)
  data$mu_a1_m4.cl <- (data$mu_a1_m4 - sfun * (1 - data$gamma.m4.a1)) / (1 - (1 - data$gamma.m4.a1) * sfun)
  
  data$mu_a1_m1.c <- ifelse(data$X < 1, data$mu_a1_m1.cl, data$mu_a1_m1.cu)
  data$mu_a1_m2.c <- ifelse(data$X < 1, data$mu_a1_m2.cl, data$mu_a1_m2.cu)
  data$mu_a1_m3.c <- ifelse(data$X < 1, data$mu_a1_m3.cl, data$mu_a1_m3.cu)
  data$mu_a1_m4.c <- ifelse(data$X < 1, data$mu_a1_m4.cl, data$mu_a1_m4.cu)
  
  data$mu_a1_m1.mpl <- data$mu_a1_m1.c * (1 - sfun)
  data$mu_a1_m2.mpl <- data$mu_a1_m2.c * (1 - sfun)
  data$mu_a1_m3.mpl <- data$mu_a1_m3.c * (1 - sfun)
  data$mu_a1_m4.mpl <- data$mu_a1_m4.c * (1 - sfun)
  
  data$mu_a1_m1.mpu <- data$mu_a1_m1.c * (1 - sfun) + sfun
  data$mu_a1_m2.mpu <- data$mu_a1_m2.c * (1 - sfun) + sfun
  data$mu_a1_m3.mpu <- data$mu_a1_m3.c * (1 - sfun) + sfun
  data$mu_a1_m4.mpu <- data$mu_a1_m4.c * (1 - sfun) + sfun
  
  data$mu_a1_m1.mp <- ifelse(data$X < 1, data$mu_a1_m1.mpl, data$mu_a1_m1.mpu)
  data$mu_a1_m2.mp <- ifelse(data$X < 1, data$mu_a1_m2.mpl, data$mu_a1_m1.mpu)
  data$mu_a1_m3.mp <- ifelse(data$X < 1, data$mu_a1_m3.mpl, data$mu_a1_m1.mpu)
  data$mu_a1_m4.mp <- ifelse(data$X < 1, data$mu_a1_m4.mpl, data$mu_a1_m1.mpu)
  
  tau.star <- max(sfun)
  
  return(list(data = data, tau.star = tau.star))
}

#' CalculateEIFBound: calculates EIF associated with bound on risk ratio scale. Requires use of
#' AddBoundVector to calculate actual bound for a given sensitivity parameter
#'
#' @param data dataset containing biased estimate
#'
#' @return dataset with EIFs of bounds -- caution: needs to be further processed for a given sensitivity parameter
CalculateEIFBound <- function(data) {
  data$t1 <- with(data, mu_a1_m1 * ((1 - gamma.M1.a1) - (1 - gamma.M1.a0)) * (1 - gamma.M2.a0) * gamma.m1.a1 +
                    mu_a1_m2 * (gamma.M1.a1 - gamma.M1.a0) * (1 - gamma.M2.a0) * gamma.m2.a1 +
                    mu_a1_m3 * ((1 - gamma.M1.a1) - (1 - gamma.M1.a0)) * gamma.M2.a0 * gamma.m3.a1 +
                    mu_a1_m4 * (gamma.M1.a1 - gamma.M1.a0) * gamma.M2.a0 * gamma.m4.a1)
  
  data$t1.1 <- with(data, mu_a1_m1 * (1 - gamma.M1.a1) * (1 - gamma.M2.a0) * gamma.m1.a1 +
                      mu_a1_m2 * gamma.M1.a1 * (1 - gamma.M2.a0) * gamma.m2.a1 +
                      mu_a1_m3 * (1 - gamma.M1.a1) * gamma.M2.a0 * gamma.m3.a1 +
                      mu_a1_m4 * gamma.M1.a1 * gamma.M2.a0 * gamma.m4.a1)
  
  data$t1.0 <- with(data, (mu_a1_m1 * (1 - gamma.M1.a0) * (1 - gamma.M2.a0) * gamma.m1.a1 +
                             mu_a1_m2 * gamma.M1.a0 * (1 - gamma.M2.a0) * gamma.m2.a1 +
                             mu_a1_m3 * (1 - gamma.M1.a0) * gamma.M2.a0 * gamma.m3.a1 +
                             mu_a1_m4 * gamma.M1.a0 * gamma.M2.a0 * gamma.m4.a1))
  
  data$t2 <- with(data, (1 - gamma.M1.a1) * (1 - gamma.M2.a0) * gamma.m1.a1 +
                    gamma.M1.a1 * (1 - gamma.M2.a0) * gamma.m2.a1 +
                    (1 - gamma.M1.a1) * gamma.M2.a0 * gamma.m3.a1 +
                    gamma.M1.a1 * gamma.M2.a0 * gamma.m4.a1)
  
  data$t3 <- with(data, (1 - gamma.M1.a0) * (1 - gamma.M2.a0) * gamma.m1.a1 +
                    gamma.M1.a0 * (1 - gamma.M2.a0) * gamma.m2.a1 +
                    (1 - gamma.M1.a0) * gamma.M2.a0 * gamma.m3.a1 +
                    gamma.M1.a0 * gamma.M2.a0 * gamma.m4.a1)
  
  
  data$gamma.M1.a1.obs <- with(data, ifelse(M1 == 1, gamma.M1.a1, 1 - gamma.M1.a1))
  data$gamma.M1.a0.obs <- with(data, ifelse(M1 == 1, gamma.M1.a0, 1 - gamma.M1.a0))
  data$gamma.M2.a1.obs <- with(data, ifelse(M2 == 1, gamma.M2.a1, 1 - gamma.M2.a1))
  data$gamma.M2.a0.obs <- with(data, ifelse(M2 == 1, gamma.M2.a0, 1 - gamma.M2.a0))
  
  data$gamma.M.M1.1.M2.obs <- with(data, ifelse(M2 == 1, gamma.m4.a1, gamma.m2.a1))
  data$gamma.M.M1.0.M2.obs <- with(data, ifelse(M2 == 1, gamma.m3.a1, gamma.m1.a1))
  
  data$gamma.M.M2.1.M1.obs <- with(data, ifelse(M1 == 1, gamma.m4.a1, gamma.m3.a1))
  data$gamma.M.M2.0.M1.obs <- with(data, ifelse(M1 == 1, gamma.m2.a1, gamma.m1.a1))
  
  data$gamma.m.a1.marg.M1.a1 <- with(data, gamma.M1.a1 * gamma.M.M1.1.M2.obs + (1 - gamma.M1.a1) * gamma.M.M1.0.M2.obs)
  data$gamma.m.a1.marg.M1.a0 <- with(data, gamma.M1.a0 * gamma.M.M1.1.M2.obs + (1 - gamma.M1.a0) * gamma.M.M1.0.M2.obs)
  data$gamma.m.a1.marg.M2.a0 <- with(data, gamma.M2.a0 * gamma.M.M2.1.M1.obs + (1 - gamma.M2.a0) * gamma.M.M2.0.M1.obs)

  data$eif.t1.1.1 <- with(data, (A / pi) * (Y * gamma.M1.a1.obs  * gamma.M2.a0.obs - t1.1))
  data$eif.t1.1.2 <- with(data, (A / pi) * (mu_a1.marg.M2.a0X.M.a1.M1.obs - t1.1))
  data$eif.t1.1.3 <- with(data, (A / pi) * (mu_a1.marg.M1.a1X.M.a1.M2.obs - t1.1))
  data$eif.t1.1   <- with(data, eif.t1.1.1 + eif.t1.1.2 + eif.t1.1.3 + t1.1)
  
  data$eif.t1.0.1 <- with(data, (A / pi) * (Y * gamma.M1.a0.obs * gamma.M2.a0.obs - t1.0))
  data$eif.t1.0.2 <- with(data, ((1 - A) / (1 - pi)) * (mu_a1.marg.M2.a0X.M.a1.M1.obs - t1.0))
  data$eif.t1.0.3 <- with(data, (A / pi) * (mu_a1.marg.M1.a0X.M.a1.M2.obs - t1.0))
  data$eif.t1.0   <- with(data, eif.t1.0.1 + eif.t1.0.2 + eif.t1.0.3 + t1.0)
  
  data$eif.t1.1.a <- with(data, (A / pi) * (Y * (gamma.M1.a1.obs - gamma.M1.a0.obs) * gamma.M2.a0.obs - t1))
  data$eif.t1.2.a <- with(data, ((A / pi) * (mu_a1.marg.M2.a0X.M.a1.M1.obs - t1.1) - ((1 - A) / (1 - pi)) * (mu_a1.marg.M2.a0X.M.a1.M1.obs - t1.0)))
  data$eif.t1.3.a <- with(data, (A / pi) * ((mu_a1.marg.M1.a1X.M.a1.M2.obs - mu_a1.marg.M1.a0X.M.a1.M2.obs) - t1))
  data$eif.t1   <- with(data, eif.t1.1.a + eif.t1.2.a + eif.t1.3.a + t1)
  
  data$eif.t2.1 <- with(data, (A / pi) * (gamma.M1.a1.obs * gamma.M2.a0.obs - t2))
  data$eif.t2.2 <- with(data, (A / pi) * (gamma.m.a1.marg.M2.a0 - t2))
  data$eif.t2.3 <- with(data, ((1 - A) / (1 - pi)) * (gamma.m.a1.marg.M1.a1 - t2))
  data$eif.t2   <- with(data, eif.t2.1 + eif.t2.2 + eif.t2.3 + t2)
  
  data$eif.t3.1 <- with(data, (A / pi) * (gamma.M1.a0.obs * gamma.M2.a0.obs - t3)) 
  data$eif.t3.2 <- with(data, (1 - A) / (1 - pi) * (gamma.m.a1.marg.M2.a0 - t3)) 
  data$eif.t3.3 <- with(data, (1 - A) / (1 - pi) * (gamma.m.a1.marg.M1.a0 - t3))
  data$eif.t3   <- with(data, eif.t3.1 + eif.t3.2 + eif.t3.3 + t3)

  return(data)
}

#' ProcessBoundData: applies functions above to generate EIFs for bounds
#'
#' @param data dataset
#' @param confounded TRUE/FALSE indicator of whether the data is confounded. If not then
#' this function simply calculates EIFs associated with the true quantities rather than
#' the biased targets
#'
#' @return processed dataset containing EIFs of the target estimand and the associated bounds
ProcessBoundData <- function(data, confounded = FALSE) {
  if (confounded == FALSE) {
    data <- data 
  }
  if (confounded == TRUE) {
    data <- data %>%
      select(-matches("mu_a1_m[1-4]$")) %>% 
      set_names(gsub("\\.c", "", names(.)))
  }
  data <- data %>%
    ProcessAllData(marg_jd = FALSE) %>%
    CalculateIDEEIFs() %>%
    MarginalizeOutcomesMargDensityBounds() %>%
    CalculateEIFBound()
  
  return(data)
}


#' AddBoundVector: returns EIF of bound for given sensitivity parameters
#'
#' @param data output ofr CalculateEIFBound
#' @param tau0 min E[Y(a, m) | A = a, M ! = m, X] / E[Y(a, m) | A = a, M = m, X]
#' @param tau1 min (1-E[Y(a, m) | A = a, M ! = m, X]) / (1-E[Y(a, m) | A = a, M = m, X])
#'
#' @return dataset with EIFs for upper and lower bounds for given sensitivity parameters
AddBoundVector <- function(data, tau0, tau1) {
  data$upper_true <- data$eif_m1_t5_a1 + # biased estimand
    tau1 * (1 - data$t2) + # tau scaled by sum of products of mediator probabilities
    tau0 * (data$eif_m1_t4_a1.0 - (tau1 / tau0) * data$eif_m1_t4_a1.1) + # scaled effect estimate
    tau0 * ((tau1 / tau0) * data$t1.1 - data$t1.0) # estimand times joint probability
  
  data$lower_true <- data$eif_m1_t5_a1 + # biased estimand
    tau1 * (data$t3 - 1) + # tau scaled by sum of products of mediator probabilities
    tau1 * (data$eif_m1_t4_a1.0 - (tau0 / tau1) * data$eif_m1_t4_a1.1) + # scaled effect estimate
    tau1 * ((tau0 / tau1) * data$t1.1 - data$t1.0) # estimand times joint probability
  
  return(data)
}

#' MarginalizeOutcomesMargDensityBounds: marginalizes quantities for CalculateEIFBound
#'
#' @param data dataset
#'
#' @return data with columns used in CalculateEIFBound
MarginalizeOutcomesMargDensityBounds <- function(data) {
  xwalk.m1.m20 <- data %>%
    filter(cc == 1) %>%
    select(M, M1, M2) %>%
    distinct() %>%
    arrange(M) %>%
    filter(M2 == 0) %>%
    mutate(M1_label.a1 = paste0("gamma.M1", M1, ".a1"),
           M1_label.a0 = paste0("gamma.M1", M1, ".a0"),
           M_label.a1  = paste0("gamma.m", M, ".a1"),
           Y_label = paste0("mu_a1_m", M))
  
  xwalk.m1.m21 <- data %>%
    filter(cc == 1) %>%
    select(M, M1, M2) %>%
    distinct() %>%
    arrange(M) %>%
    filter(M2 == 1) %>%
    mutate(M1_label.a1 = paste0("gamma.M1", M1, ".a1"),
           M1_label.a0 = paste0("gamma.M1", M1, ".a0"),
           M_label.a1  = paste0("gamma.m", M, ".a1"),
           Y_label = paste0("mu_a1_m", M))
  
  xwalk.m2.m10 <- data %>%
    filter(cc == 1) %>%
    select(M, M1, M2) %>%
    distinct() %>%
    filter(M1 == 0) %>%
    arrange(M) %>%
    mutate(M2_label.a0 = paste0("gamma.M2", M2, ".a0"),
           M_label.a1  = paste0("gamma.m", M, ".a1"),
           Y_label = paste0("mu_a1_m", M))
  
  xwalk.m2.m11 <- data %>%
    filter(cc == 1) %>%
    select(M, M1, M2) %>%
    distinct() %>%
    filter(M1 == 1) %>%
    arrange(M) %>%
    mutate(M2_label.a0 = paste0("gamma.M2", M2, ".a0"),
           M_label.a1  = paste0("gamma.m", M, ".a1"),
           Y_label = paste0("mu_a1_m", M))
  
  
  sum_list1 <- list(); sum_list5 <- list()
  sum_list2 <- list(); sum_list6 <- list()
  sum_list3 <- list() 
  sum_list4 <- list() 
  
  for(i in 1:nrow(xwalk.m1.m20)) {
    # M2 = 0
    sum_list1[[i]] <- data[[xwalk.m1.m20$Y_label[i]]]*data[[xwalk.m1.m20$M1_label.a1[i]]]*data[[xwalk.m1.m20$M_label.a1[i]]]
    sum_list2[[i]] <- data[[xwalk.m1.m20$Y_label[i]]]*data[[xwalk.m1.m20$M1_label.a0[i]]]*data[[xwalk.m1.m20$M_label.a1[i]]]
    # M2 = 1
    sum_list3[[i]] <- data[[xwalk.m1.m21$Y_label[i]]]*data[[xwalk.m1.m21$M1_label.a1[i]]]*data[[xwalk.m1.m21$M_label.a1[i]]]
    sum_list4[[i]] <- data[[xwalk.m1.m21$Y_label[i]]]*data[[xwalk.m1.m21$M1_label.a0[i]]]*data[[xwalk.m1.m21$M_label.a1[i]]]
    # M1 = 0
    sum_list5[[i]] <- data[[xwalk.m2.m10$Y_label[i]]]*data[[xwalk.m2.m10$M2_label.a0[i]]]*data[[xwalk.m2.m10$M_label.a1[i]]]
    # M1 = 1
    sum_list6[[i]] <- data[[xwalk.m2.m11$Y_label[i]]]*data[[xwalk.m2.m11$M2_label.a0[i]]]*data[[xwalk.m2.m11$M_label.a1[i]]]
  }
  
  data[["mu_a1.marg.M1.a1xM.a1.M2e0"]] <- Reduce(`+`, sum_list1)
  data[["mu_a1.marg.M1.a0xM.a1.M2e0"]] <- Reduce(`+`, sum_list2)
  data[["mu_a1.marg.M1.a1xM.a1.M2e1"]] <- Reduce(`+`, sum_list3)
  data[["mu_a1.marg.M1.a0xM.a1.M2e1"]] <- Reduce(`+`, sum_list4)
  
  data[["mu_a1.marg.M2.a0xM.a1.M1e0"]] <- Reduce(`+`, sum_list5)
  data[["mu_a1.marg.M2.a0xM.a1.M1e1"]] <- Reduce(`+`, sum_list6)
  
  data[["mu_a1.marg.M1.a1X.M.a1.M2.obs"]] <- ifelse(data$M2 == 1, data$mu_a1.marg.M1.a1xM.a1.M2e1, data$mu_a1.marg.M1.a1xM.a1.M2e0)
  data[["mu_a1.marg.M1.a0X.M.a1.M2.obs"]] <- ifelse(data$M2 == 1, data$mu_a1.marg.M1.a0xM.a1.M2e1, data$mu_a1.marg.M1.a0xM.a1.M2e0)
  data[["mu_a1.marg.M2.a0X.M.a1.M1.obs"]] <- ifelse(data$M1 == 1, data$mu_a1.marg.M2.a0xM.a1.M1e1, data$mu_a1.marg.M2.a0xM.a1.M1e0)
  
  return(data)
}

#' CalculateBound
#'
#' @param data dataset
#' @param tau0 min E[Y(a, m) | A = a, M ! = m, X] / E[Y(a, m) | A = a, M = m, X]
#' @param tau1 min (1-E[Y(a, m) | A = a, M ! = m, X]) / (1-E[Y(a, m) | A = a, M = m, X])
#'
#' @return dataframe with estimates of the bound using the EIF and associated CIs as
#' well as the truth
CalculateBound <- function(data, tau0, tau1) {
  
  upper_true <- data$eif_m1_t5_a1 + 
    tau1 * ((tau0 / tau1) * data$eif_m1_t4_a1.0 - data$eif_m1_t4_a1.1) + 
    tau1 * (1 + data$t1.1 - (tau0 / tau1) * data$t1.0 - data$t2)

  lower_true <- data$eif_m1_t5_a1 +
    tau0 * ((tau1 / tau0) * data$eif_m1_t4_a1.0 - data$eif_m1_t4_a1.1) +
    tau0 * (-1 + data$t1.1 - (tau1 / tau0) * data$t1.0 + (tau1 / tau0) * data$t3)

  eif_upper <- data$eif_m1_a1 + 
    tau1 * ((tau0 / tau1) * data$eif_m1_a1.0 - data$eif_m1_a1.1) + 
    tau1 * (1 + data$eif.t1.1 - (tau0 / tau1) * data$eif.t1.0 - data$eif.t2)
  
  eif_lower <- data$eif_m1_a1 +
    tau0 * ((tau1 / tau0) * data$eif_m1_t4_a1.0 - data$eif_m1_t4_a1.1) +
    tau0 * (-1 + data$eif.t1.1 - (tau1 / tau0) * data$eif.t1.0 + (tau1 / tau0) * data$eif.t3)

  res <- tibble(
    bound  = c(mean(eif_upper), mean(eif_lower)),
    se.est = sqrt(c(var(eif_upper)/nrow(data), var(eif_lower)/nrow(data))),
    lci = bound - 1.96 * se.est,
    uci = bound + 1.96 * se.est,
    truth = c(mean(upper_true), mean(lower_true)),
    type = c("upper", "lower"),
    tau = c(tau1, tau0)
  )
  return(res)
}

MediationBoundDRL <- function(n, mu_rate, pi_rate, m_rate, m1_rate, m2_rate, population, 
                         bin_constant = 1, mu_constant = 8, return_data = FALSE) {
  sample_rows <- sample(1:nrow(population), 2*n, replace = FALSE)
  samp   <- population[sample_rows,]
  psuedo <- BoundDeltaEIF(CalculateEIFbound(CalculateIDEEIFs(ProcessAllData(samp, 
                                                                            marg_jd = TRUE))), delta)[, c("eif_m1_a1", "eif_m1_t1_a1", "eif_m1_t2_a1", "eif_m1_t3_a1", 
                                                                                                                 "eif_m1_t4_a1", "eif_m1_t5_a1",
                                                                                                                 "eif_bound_u", "eif_bound_l")]
  
  mu.sigma <- matrix(rep(0.15/n^(2*mu_rate), 16), 4, 4)
  diag(mu.sigma) <- 1/n^(2*mu_rate)
  mu.sigma <- mu_constant * mu.sigma
  
  m.sigma <- bin_constant * matrix(rep(0.15/n^(2*m_rate), 9), 3, 3)
  diag(m.sigma) <- bin_constant/n^(2*m_rate)
  
  m1.sigma <- bin_constant * matrix(rep(0.15/n^(2*m1_rate), 4), 2, 2)
  diag(m1.sigma) <- bin_constant/n^(2*m1_rate)
  
  dat <- tibble(ones = rep(1, 2*n))
  dat$pi <- with(samp, expit(logit(pi) + rnorm(2*n, 
                                               mean = sqrt(bin_constant)/n^pi_rate, 
                                               sd   = sqrt(bin_constant)/n^pi_rate)))
  
  mu.ests <- samp[, c("mu_a1_m1", "mu_a1_m2", "mu_a1_m3", "mu_a1_m4")] + MASS::mvrnorm(2*n, rep(sqrt(mu_constant)/n^mu_rate, 4), mu.sigma)
  
  m.ests <- expit(logit(samp[, c("gamma.m1.a1", "gamma.m2.a1", "gamma.m3.a1")]) + MASS::mvrnorm(2*n, rep(sqrt(bin_constant)/n^m_rate, 3), m.sigma))
  
  m.ests$gamma.m4.a1 <- 1 - rowSums(m.ests)
  
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
  
  dat <- CalculateIDEEIFs(ProcessAllData(dat, marg_jd = FALSE)) %>%
    CalculateEIFbound() %>%
    BoundDeltaEIF(delta)
  
  S <- sample(c(rep(1, n), rep(2, n)), 2*n, replace = FALSE)
  
  # average effects bounds
  lb.avg.true <- mean(psuedo$eif_bound_u[S==2])
  ub.avg.true <- mean(psuedo$eif_bound_l[S==2])
  lb.avg.est  <- mean(dat$eif_bound_u[S==2])
  ub.avg.est  <- mean(dat$eif_bound_l[S==2])
  
  lb.se.true <- sd(psuedo$eif_bound_u[S==2]) / sqrt(length(S[S==2]))
  ub.se.true <- sd(psuedo$eif_bound_l[S==2]) / sqrt(length(S[S==2]))
  lb.se.est  <- sd(dat$eif_bound_u[S==2]) / sqrt(length(S[S==2]))
  ub.se.est  <- sd(dat$eif_bound_l[S==2]) / sqrt(length(S[S==2]))
  
  tau.avg <- mean(tau)
  
  res.avg <- tibble(
    truth = c(tau.avg, tau.avg),
    estimator = c("true eif", "est eif"),
    lower_bound = c(lb.avg.true, lb.avg.est),
    upper_bound = c(ub.avg.true, ub.avg.true),
    lb.se = c(lb.se.true, lb.se.est),
    ub.se = c(ub.se.true, ub.se.est)
  ) %>%
    mutate(
      lci.lb = lower_bound - 1.96 * lb.se,
      uci.ub = upper_bound + 1.96 * ub.se
    )
  
  # conditional effects bounds
  
  return(res.avg)
}

