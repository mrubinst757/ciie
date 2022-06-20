# verified that my proposed EIF estimates the quantity it is supposed to estimate
truth.c1 <- mean(population.c1$eif_m1_t5_a1)
truth.c2 <- mean(population.c2$eif_m1_t5_a1)
truth.u1 <- mean(population.u1$eif_m1_t5_a1)

population1.proc <- CalculateEIFBound(population1.proc) 
population2.proc <- CalculateEIFBound(population2.proc) 

samp_test <- function(population, sample_size, tau1, tau0) {
  truth <- CalculateBound(population, tau1, tau0)$truth
  samp <- sample_n(population, sample_size)
  test <- CalculateBound(samp, tau1, tau0) %>%
    mutate(truth = truth, 
           captured = ifelse(truth > lci & truth < uci, 1, 0))
  return(test)
}

test.c1 %>%
  group_by(type) %>%
  summarize_all(mean)

test.c1 <- map(1:250, ~samp_test(population1.proc, 1000, 
                                 tau1 = 1-tau.star.c1$upper, 
                                 tau0 = 1-tau.star.c1$lower)) %>%
  invoke(rbind, .) 

test.c2 <- map(1:250, ~samp_test(population2.proc, 1000, 
                                 tau1 = 1-tau.star.c2$upper, 
                                 tau0 = 1-tau.star.c2$lower)) %>%
  invoke(rbind, .) 

test.c2a <- map(1:250, ~samp_test(population2.proc, 1000, 
                                  tau1 = 0.05, 
                                  tau0 = 0.05)) %>%
  invoke(rbind, .) 


test.c2a %>%
  select(-captured, -se.est, - truth) %>%
  nest(data = c(bound, lci, uci, tau)) %>%
  spread(type, data) %>%
  unnest() %>%
  mutate(truth = truth.c2,
         captured = ifelse(truth > lci & truth < uci1, 1, 0)) %>%
  summarize_all(mean)
