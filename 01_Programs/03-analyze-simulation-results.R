library(tidyverse); library(gt)

# 1. Prediction at X = \{0, 2\}

# process results --------------------------------------------------------------
results1 <- list.files("02_Output/", "simulations-1k-", full.names = TRUE) %>%
  map(readRDS) 

results1 <- map2(results1[[1]], results1[[2]], append)

results2 <- list.files("02_Output/", "simulations-2k-", full.names = TRUE) %>%
  map(readRDS) 

results2 <- map2(results2[[1]], results2[[2]], append)

#results2 <- readRDS("02_Output/simulations-1k-p2.rds")
results <- append(results1, results2)

all_results <- map(results, ~invoke(rbind, .x)) %>%
  map2(c(0, 2, 0, 2), ~mutate(.x, point = .y)) %>% # add when ready
  map2(c(1000, 1000, 2000, 2000), ~mutate(.x, n = .y)) %>%
  invoke(rbind, .) %>%
  group_by(estimand, nuisance_est, type, n, point) %>%
  mutate(covered = ifelse(truth > lci & truth < uci, 1, 0), 
         bias = ests - truth,
         var  = (ests - mean(ests))^2,
         mse = (ests - truth)^2) %>%
  summarize_all(mean) %>%
  mutate_at(c("point", "n"), as.integer)

prop_med_comp <- all_results %>%
  filter(grepl("Proportion", estimand), nuisance_est == "NP") %>%
  filter(!grepl("PI|Plugin", type), n == 1000) %>%
  nest(data = c(ests, var, lci, uci, truth, covered, bias, mse)) %>%
  spread(estimand, data) %>%
  unnest() %>%
  ungroup() %>%
  select(Estimator = type, Point = point, `MSE: EIF Ratio` = mse, `Coverage: EIF Ratio` = covered, 
         `MSE: EIF Separate` = mse1, `Coverage: EIF Separate` = covered1)

print(xtable::xtable(prop_med_comp, caption = "Proportion mediated: comparison of two approaches",
                     label = "tab:propmedcomp"), include.rownames = FALSE,
      caption.placement = "top", file = "../Mediation-Project/tables/prop-med-comp-table.tex")

# projection results
projection_table1 <- all_results %>%
  ungroup() %>%
  #filter(grepl("Proj", type), estimand != "Proportion mediated") %>%
  filter(grepl("Proj", type)) %>%
  filter(n == 1000, nuisance_est == "NP") %>%
  mutate(rmse = sqrt(mse), Std = sqrt(var)) %>%
  select(Truth = truth, Estimand = estimand, Point = point, type, 
         Bias = bias, Std, RMSE = rmse, Coverage = covered) %>%
  nest(data = c(Point, Truth, type, Bias, Std, RMSE, Coverage)) %>%
  spread(Estimand, data) %>%
  select(`IIE-M1`, `Total effect`, `Proportion mediated` = `Proportion mediated alt`) %>%
  #select(`IIE-M1`, `Total effect`) %>%
  unnest() %>%
  #select(-type1, -Point1, -type2, -Point2) %>%
  #select(-type1, -Point1) %>%
  mutate(Projection = ifelse(grepl("Lin", type), "Linear", "Quadratic")) %>%
  mutate(Strategy = ifelse(grepl("PI", type), "Plugin", "DR")) %>%
  select(Point, Strategy, Projection, Truth, Bias, Std, RMSE, Coverage, 
         Truth1, Bias1, Std1, RMSE1, Coverage1,#) %>%
         Truth2, Bias2, Std2, RMSE2, Coverage2) %>%
  mutate_at(vars(matches("coverage")), ~.*100)

gt_projection_table1 <- gt::gt(projection_table1) %>%
  tab_spanner(label = "CIIE-M1", columns = c(Truth, Bias, Std, RMSE, Coverage)) %>%
  tab_spanner(label = "CATE", columns = c(Truth1, Bias1, Std1, RMSE1, Coverage1)) %>%
  tab_spanner(label = "Proportion mediated", columns = c(Truth2, Bias2, Std2, RMSE2, Coverage2)) %>%
  cols_label(Truth1 = "Truth", Bias1 = "Bias", Std1 = "Std", RMSE1 = "RMSE", Coverage1 = "Coverage") %>%
  cols_label(Truth2 = "Truth", Bias2 = "Bias", Std2 = "Std", RMSE2 = "RMSE", Coverage2 = "Coverage") %>%
  fmt_number(columns = matches("Truth|Bias|Std|RMSE"), decimals = 3) %>%
  fmt_number(columns = matches("Cover"), decimals = 1) %>%
  tab_header(title = "Projection estimators: simulation performance (N = 1000)")

projection_table1a <- projection_table1 %>%
  select(-matches("[1-2]$")) %>%
  mutate_at(vars(matches("Truth|Bias|Std|RMSE")), ~round(., 3)) %>%
  mutate_at(vars(matches("Cover")), ~round(., 2))
  
gtsave(gt_projection_table1, "projection-table-gt.tex", path = "../Mediation-Project/tables/")

print(xtable::xtable(projection_table1a, caption = "Projection estimators: simulation performance (N = 1000)",
                     label = "tab:sim1"), include.rownames = FALSE,
      caption.placement = "top", file = "../Mediation-Project/tables/projection-table.tex")

drl_table2 <- all_results %>%
  ungroup() %>%
  filter(grepl("DRL", type)) %>%
  filter(nuisance_est == "NP") %>%
  mutate(rmse = sqrt(mse), Std = sqrt(var)) %>%
  select(Estimand = estimand, type, Point = point, Truth = truth,
         Bias = bias, Std, RMSE = rmse, 
         Coverage = covered, N = n) %>%
  nest(data = c(N, Point, Truth, type, Bias, Std, RMSE, Coverage)) %>%
  spread(Estimand, data) %>%
  select(`IIE-M1`, `Total effect`) %>%
  unnest() %>%
  select(-type1, -Point1, -N1) %>%
  mutate(Strategy = ifelse(grepl("Plugin", type), "Plug-in", "DRLearner")) %>%
  select(Point, Strategy, N, Truth, Bias, Std, RMSE, Coverage, Truth1, Bias1, Std1, RMSE1, Coverage1) %>%
  mutate_at(vars(matches("coverage")), ~.*100) %>%
  arrange(Point, Strategy)

gt_drl_table2 <- drl_table2 %>%
  mutate_at(vars(matches("Cover")), ~as.character(round(., 2))) %>%
  replace_na(list(Coverage = "-", Coverage1 = "-")) %>%
  gt() %>%
  tab_spanner(label = "CIIE-M1", columns = c(Truth, Bias, Std, RMSE, Coverage)) %>%
  tab_spanner(label = "CATE", columns = c(Truth1, Bias1, Std1, RMSE1, Coverage1)) %>%
  cols_label(Truth1 = "Truth", Bias1 = "Bias", Std1 = "Std", RMSE1 = "RMSE", Coverage1 = "Coverage") %>%
  fmt_number(columns = matches("Truth|Bias|Std|RMSE"), decimals = 3) %>%
  #fmt_number(columns = matches("Coverage"), decimals = 1) %>%
  tab_header(title = "Non-parametric estimators: simulation performance")

drl_table2a <- drl_table2 %>%
  select(-ends_with("1")) %>%
  mutate_at(vars(matches("Truth|Bias|Std|RMSE")), ~round(., 3)) %>%
  mutate_at(vars(matches("Cover")), ~as.character(round(., 2))) %>%
  replace_na(list(Coverage = "-")) %>%
  arrange(Point)

gtsave(gt_drl_table2, "drl-table-gt.tex", path = "../Mediation-Project/tables/")

print(xtable::xtable(drl_table2a, caption = "Non-parametric estimators: simulation performance",
                     label = "tab:sim2"), include.rownames = FALSE,
      caption.placement = "top", file = "../Mediation-Project/tables/drl-table.tex")

#2. integrated MSE
#sims_imse <- 

readRDS("02_Output/sim-test-1.rds") %>%
  map(flatten) %>%
  map(~map(.x, ~as_tibble(.x$res)))
  map(~invoke(rbind, .x)) %>%
  map(as_tibble) %>%
  map2(c("All 0.5", "Slow Mu", "Slow M1", "Slow M2", "All 0.1", "Slow Pi", "Slow M"),
       ~mutate(.x, Specification = .y)) %>%
  map(~mutate(.x, n = rep(c(500, 1000, 5000, 10000, 50000, 100000), 1000))) %>%
  map(unnest) %>%
  map(~mutate(.x, Estimator = rep(c("drl", "oracle", "plugin", "plugin.s"), nrow(.) / 4))) %>%
  invoke(rbind, .) %>%
  group_by(n, Specification, Estimator) %>%
  summarize_all(~mean(.)*sqrt(n))

imse_paper_plot <- sims_imse %>%
  distinct() %>%
  spread(Estimator, res) %>%
  filter(n > 1000, Specification != "Slow M", Specification != "Slow Pi", Specification != "All 0.1") %>%
  rename(`DRLearner` = drl, `Oracle` = oracle, Plugin = plugin) %>%
  gather(Estimator, value, DRLearner, Oracle, Plugin) %>%
  ggplot(aes(y = value, x = n, fill = Estimator, color = Estimator)) +
  geom_point() +
  geom_line() +
  facet_wrap(~Specification, scales = "free_y") +
  theme_minimal() +
  ylab("n^(0.5) * MSE") +
  xlab("Sample size") +
  scale_color_brewer(palette = "Set1")

ggsave("../Mediation-Project/plots/sim-imse.png", imse_paper_plot)

imse_data %>%
  distinct() %>%
  filter(n > 1000, Specification %in% c("Slow M", "Slow Pi")) %>%
  gather(key, value, drl, oracle, plugin) %>%
  ggplot(aes(y = value, x = n, fill = key, color = key)) +
  geom_point() +
  geom_line() +
  facet_wrap(~Specification, scales = "free_y") +
  theme_minimal() +
  ylab("n^(0.5) * RMSE") +
  xlab("Sample size")

