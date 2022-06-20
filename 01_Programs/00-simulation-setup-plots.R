source("../../Depression/Programs/01_Analysis/02_AnalyzeData/00-mediation-functions.R")
source("01_Programs/simulation-funs.R")
source("01_Programs/bound-funs.R")
library(RColorBrewer)

# Nuisance function and estimand plots ------------------------------------------------
population.u1 <- GeneratePopulation(200000, "binary")

population.u1 <- population.u1 %>%
  ProcessBoundData(confounded = FALSE)

set.seed(12)
samp <- sample_n(population.u1, 2000) %>%
  mutate(te = mu_1 - mu_0)

# Nuisance function plots
nuis.fun.plot <- samp %>%
  select(X, matches("mu_a1_m[1-4]$"), pi, matches("gamma.m[1-4].a1$", ignore.case = FALSE),
         matches("gamma.M1.a1$"), gamma.M1.a0, gamma.M2.a1, gamma.M2.a0) %>%
  gather(Function, value, -X) %>%
  mutate(type = case_when(
    grepl("mu", Function) ~ "Outcome model",
    grepl("pi", Function) ~ "Propensity score",
    grepl("M1", Function) ~ "Probability of M1", 
    grepl("M2", Function) ~ "Probability of M2", 
    grepl("m[1-4]", Function) ~ "Joint mediator probability" 
  )) %>%
  mutate_at("type", ~factor(., levels = c("Outcome model", "Joint mediator probability",
                                          "Propensity score", "Probability of M1",
                                          "Probability of M2"))) %>%
  mutate_at("Function", ~stringr::str_replace_all(.,
    c("gamma.M1.a0" = "P(M1(0) = 1 | X)",
      "gamma.M1.a1" = "P(M1(1) = 1 | X)",
      "gamma.M2.a0" = "P(M2(0) = 1 | X)",
      "gamma.M2.a1" = "P(M2(1) = 1 | X)",
      "mu_a1_m1" = "E(Y(a, M1 = 0, M2 = 0) | X)",
      "mu_a1_m2" = "E(Y(a, M1 = 0, M2 = 1) | X)",
      "mu_a1_m3" = "E(Y(a, M1 = 1, M2 = 0) | X)",
      "mu_a1_m4" = "E(Y(a, M1 = 1, M2 = 1) | X)",
      "gamma.m1.a1" = "P(M1(1) = 0, M2(1) = 0 | X)",
      "gamma.m2.a1" = "P(M1(1) = 0, M2(1) = 1 | X)",
      "gamma.m3.a1" = "P(M1(1) = 1, M2(1) = 0 | X)",
      "gamma.m4.a1" = "P(M1(1) = 1, M2(1) = 1 | X)",
      "pi" = "P(A = 1 | X)")
  )) %>%
  ggplot(aes(x = X, y = value, color = Function)) +
  geom_point() +
  facet_wrap(~type) +
  theme_minimal() +
  ylab("")

# Estimands plot
estimands.plot <- samp %>%
  mutate(te = mu_1 - mu_0,
         eif_m1_t5_a0 = -eif_m1_t5_a0, eif_m2_t5_a0 = -eif_m2_t5_a0,
         eif_de_t3_a0 = -eif_de_t3_a0,
         cov_a1 = te - eif_de_t3_a0 - eif_m1_t5_a1 - eif_m2_t5_a1,
         cov_a0 = te - eif_de_t3_a1 - eif_m1_t5_a0 - eif_m2_t5_a0) %>%
  mutate(prop_m1_a0 = 100*eif_m1_t5_a0/te,
         prop_m2_a0 = 100*eif_m2_t5_a0/te,
         prop_m1_a1 = 100*eif_m1_t5_a1/te,
         prop_m2_a1 = 100*eif_m2_t5_a1/te) %>%
  select(X, te, eif_m1_t5_a1, eif_m2_t5_a1, prop_m1_a1, prop_m2_a1) %>%
  gather(`Effect type`, value, -X) %>%
  mutate(type = ifelse(grepl("prop", `Effect type`), "Proportion of total effect", "Effect")) %>%
  mutate_at("Effect type", ~stringr::str_replace_all(.,
                  c("eif_m1_t5_a1" = "Via M1",
                  "eif_m2_t5_a1" = "Via M2",
                  "te" = "Total effect",
                  "prop_m1_a1" = "Via M1",
                  "prop_m2_a1" = "Via M2"))) %>%
  ggplot(aes(x = X, y = value, fill = `Effect type`, color = `Effect type`)) +
  geom_point() +
  ylab("") +
  facet_wrap(~type, scales = "free") +
  theme_minimal() +
  scale_color_manual(values = c("#E65955", "#015799", "#3E9DE6"))

# Inverse probability weight plot
inv.weight.plot <- samp %>%
  rowwise() %>%
  mutate(min.gamma  = min(gamma.m1.a1, gamma.m2.a1, gamma.m3.a1, gamma.m4.a1),
         max.weight = max(1/(pi * min.gamma), 1/(1-pi), 1/pi)) %>%
  ggplot(aes(x = X, y = max.weight)) +
  geom_point() +
  theme_minimal() +
  ylab("Maximum inverse probability weight")

ggsave("../paper/plots/sim-nuis-funs.png", nuis.fun.plot)
ggsave("../paper/plots/sim-estimands.png", estimands.plot)
ggsave("../paper/plots/sim-inv-weights.png", inv.weight.plot)


# Confounded data and sensitivity analysis plots ---------------------------------------
population.u1.c1 <- population.u1 %>%
  AddConfounding(10/3) 

# bias of regression function plot
bias.plot <- population.u1.c1$data %>%
  sample_n(1000) %>%
  select(X, matches("mu_a1_m[1-4]$"), matches("mu_a1_m[1-4].c$")) %>%
  gather(key, value, -X) %>%
  mutate(Mclass = as.character(stringr::str_extract_all(key, "m[1-4]")),
         `Target function` = case_when(
           grepl("c", key) ~ "Biased",
           TRUE ~ "True"
         )) %>%
  mutate_at("Mclass", ~stringr::str_replace_all(.,
                    c("m1" = "E[Y(1, M1 = 0, M2 = 0) | X]",
                      "m2" = "E[Y(1, M1 = 0, M2 = 1) | X]",
                      "m3" = "E[Y(1, M1 = 1, M2 = 0) | X]",
                      "m4" = "E[Y(1, M1 = 1, M2 = 1) | X]"))) %>%
  ggplot(aes(x = X, y = value, fill = `Target function`, color = `Target function`)) +
  geom_point() +
  facet_wrap(~Mclass, scales = "free") +
  theme_minimal() +
  ylab("")

ggsave("../paper/plots/sim-bias-tau-10.png", bias.plot)

population.u1.c1$data <- population.u1.c1$data %>%
  ProcessBoundData(confounded = TRUE)

# bound as function of tau plot
set.seed(12)
sample_rows <- sample(seq(1:nrow(population.u1)), 5000, replace = FALSE)
samp.c1.list <- map(seq(0.01, 0.1, 0.01), ~AddBoundVector(population.u1.c1$data, .x, .x)) %>%
  map(~.x[sample_rows,])
samp <- population.u1[sample_rows,]

dat <- samp.c1.list %>%
  map(~select(.x, X, upper_true, lower_true)) %>%
  map2(seq(0.01, 0.1, 0.01), ~mutate(.x, Tau = .y)) %>%
  map(~mutate(.x, Truth = samp$eif_m1_t5_a1)) %>%
  invoke(rbind, .) %>%
  gather(key, value, upper_true, lower_true) %>%
  mutate_at("Tau", factor)

blues.alt <- colorRampPalette(brewer.pal(9, "Blues")[3:7])(10)
greens.alt <- colorRampPalette(brewer.pal(9, "Greens")[3:7])(10)

library(ggnewscale)

bound.plot <- dat %>%
  ggplot() +
  geom_point(data = samp, aes(x = X, y = eif_m1_t5_a1, alpha = 0.5)) +
  geom_point(data = filter(dat, key == "upper_true"), aes(x = X, y = value, color = Tau, alpha = 0.5)) +
  scale_color_manual(values = blues.alt) +
  new_scale_color() +
  geom_point(data = filter(dat, key == "lower_true"), aes(x = X, y = value, color = Tau, alpha = 0.5)) +
  scale_color_manual(values = greens.alt) +
  theme_minimal() +
  ylab("") +
  guides(alpha = "none")

ggsave("../Mediation-Project/plots/sim-bound.png", bound.plot)


