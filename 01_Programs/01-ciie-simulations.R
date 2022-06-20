source("../../Depression/Programs/01_Analysis/02_AnalyzeData/00-mediation-functions.R")
source("01_Programs/simulation-funs.R")
library(SuperLearner)

TProc <- function(result, rate) {
  t1 <- map(result, ~map(.x, ~.x$res)) %>%
    map(~invoke(rbind, .x) %>% as_tibble()) %>%
    map2(c(500, 1000, 5000, 10000, 50000, 100000), ~mutate(.x, n = .y)) %>%
    invoke(rbind, .) %>%
    group_by(n) %>%
    summarize_at(c("drl", "oracle", "plugin"), ~sqrt(mean(.)))
  
  t1 %>%
    mutate_at(c("drl", "oracle", "plugin"), ~.*n^(rate))
}

ProcRe <- function(result) {
  map(result, ~.x$res) %>%
    invoke(rbind, .) %>%
    as_tibble() %>%
    summarize_all(mean)
}

set.seed(100)
population  <- GeneratePopulation(1000000, "binary")

# Test 1: Projection and DRL learner at points X = {0, 2} using SuperLearner
test.p2 <- RunExperiments(1000, population, 1000, tibble(X = 2), c("SL.glm", "SL.ranger"))
test.p0 <- RunExperiments(1000, population, 1000, tibble(X = 0), c("SL.glm", "SL.ranger"))
test.p2.n2k <- RunExperiments(1000, population, 2000, tibble(X = 2), c("SL.glm", "SL.ranger"))
test.p0.n2k <- RunExperiments(1000, population, 2000, tibble(X = 0), c("SL.glm", "SL.ranger"))

saveRDS(list(test.p0 = test.p0, test.p2 = test.p2), "02_Output/simulations-1k-p2.rds")
saveRDS(list(test.p0.n2k = test.p0.n2k, test.p2.n2k = test.p2.n2k), "02_Output/simulations-2k-p1.rds")

# Test 2: Compare DRL and Oracle across range of convergence rates
set.seed(NULL)
test1 <- DRLConvergenceTest(population, 900, mu_rate = 0.4, pi_rate = 0.4, m_rate = 0.4, m1_rate = 0.4, m2_rate = 0.4) # All 0.4
test2 <- DRLConvergenceTest(population, 900, mu_rate = 0.1, pi_rate = 0.4, m_rate = 0.4, m1_rate = 0.4, m2_rate = 0.4) # mu slow
test3 <- DRLConvergenceTest(population, 900, mu_rate = 0.4, pi_rate = 0.4, m_rate = 0.4, m1_rate = 0.1, m2_rate = 0.4) # m1 slow
test4 <- DRLConvergenceTest(population, 900, mu_rate = 0.4, pi_rate = 0.4, m_rate = 0.4, m1_rate = 0.4, m2_rate = 0.1) # m2 slow
test5 <- DRLConvergenceTest(population, 900, mu_rate = 0.1, pi_rate = 0.1, m_rate = 0.1, m1_rate = 0.1, m2_rate = 0.1) # All 0.1
test6 <- DRLConvergenceTest(population, 900, mu_rate = 0.4, pi_rate = 0.1, m_rate = 0.4, m1_rate = 0.4, m2_rate = 0.4) # pi slow
test7 <- DRLConvergenceTest(population, 900, mu_rate = 0.4, pi_rate = 0.4, m_rate = 0.1, m1_rate = 0.4, m2_rate = 0.4) # m slow

# --

test <- list(test1, test2, test3, test4, test5, test6, test7)

saveRDS(test, "02_Output/sim-test-1-p2.rds")

test <- readRDS("02_Output/sim-test-1.rds")

TProc(test[[1]], 0.5)

TProc(test2, 0.25)
TProc(test3, 0.3)
TProc(test4, 0.4)
TProc(test5, 0.5)
TProc(test5, 0.55)

TProc(test4, 0.5)


#-----------
SimSeqPanel <- function(rate_seq, fixed_rate, n, population) {
  res1 <- map(rate_seq, ~DRLMediation(n, pi_rate = fixed_rate, med_rate = fixed_rate, 
                                      mu_rate = .x, population)) %>%
    invoke(rbind, .) %>%
    as_tibble() %>%
    mutate(mu_rate = rate_seq, pi_rate = fixed_rate, med_rate = fixed_rate,
           fixed = "Pi/Med")
  
  res2 <- map(rate_seq, ~DRLMediation(n, pi_rate = fixed_rate, med_rate = .x, 
                                      mu_rate = fixed_rate, population)) %>%
    invoke(rbind, .) %>%
    as_tibble() %>%
    mutate(med_rate = rate_seq, pi_rate = fixed_rate, mu_rate = fixed_rate,
           fixed = "Pi/Mu")
  
  res3 <- map(rate_seq, ~DRLMediation(n, pi_rate = .x, med_rate = fixed_rate, 
                                      mu_rate = fixed_rate, population)) %>%
    invoke(rbind, .) %>%
    as_tibble() %>%
    mutate(pi_rate = rate_seq, med_rate = fixed_rate, mu_rate = fixed_rate,
           fixed = "Mu/Med")
  
  bind_rows(res1, res2, res3)
}

test.fr30 <- map(1:250, ~SimSeqPanel(rateseq, 0.3,  n, population = population))
test.fr25 <- map(1:250, ~SimSeqPanel(rateseq, 0.25, n, population = population))
test.fr20 <- map(1:250, ~SimSeqPanel(rateseq, 0.2,  n, population = population))

test.fr25 %>%
  invoke(rbind, .) %>%
  mutate(convergence = rep(rateseq, 3*250)) %>%
  group_by(convergence, fixed) %>%
  summarize(drl = mean(drl), oracle = mean(oracle)) %>%
  gather(key, value, drl, oracle) %>%
  filter(value < 1) %>%
  ggplot(aes(x = convergence, y = value, color = key)) +
  geom_point() +
  geom_line() +
  facet_wrap(~fixed) 
