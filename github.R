library(dplyr)
library(glmmTMB)
library(readr)
library(sandwich)
library(clubSandwich)
library(lmtest)
library(parameters)
library(nnet)
library(ggplot2)
library(dplyr)
library(future)
library(furrr)
library(performance)

set.seed(123)

# Table 2

df <- read_csv("data_without_price.csv")

df <- df %>%
  mutate(
    CPV = gsub("'", "", CPV),
    CPV = trimws(CPV)
  ) %>%
  filter(substr(CPV, 3, 3) != "0") %>%
  mutate(CPV = substr(CPV, 1, 3))

df <- df %>%
  group_by(CPV) %>%
  filter(n() >= 30) %>%
  ungroup()

df <- df %>%
  mutate(
    G_WEIGHT_SQUARED = G_WEIGHT^2
  )

df <- df %>%
  mutate(
    CPV = as.factor(CPV),
    REGION = as.factor(REGION),
    STATUS = as.factor(STATUS)
  )

df <- df %>%
  mutate(
    CPV = gsub("'", "", CPV),
    CPV = trimws(CPV)
  ) %>%
  filter(substr(CPV, 3, 3) != "0") %>%
  mutate(CPV = substr(CPV, 1, 3)) %>%
  mutate(CPV2 = substr(CPV, 1, 2))

df <- df %>%
  group_by(CAE_SIREN) %>%
  mutate(n_siren = n()) %>%
  ungroup() %>%
  mutate(w_inv = 1 / n_siren)
df <- df %>%
  mutate(w_samp = w_inv * (nrow(df) / sum(w_inv)))


fixed_trunc_model <- glmmTMB(
  OFFERS ~ G_CLAUSE + G_WEIGHT + ALLOTMENT + FRAMEWORK_AGREEMENT + 
    P_CRITERION_WEIGHT_ENV + REGION + STATUS + factor(CPV2),
  data = df,
  weights = w_samp,
  family = truncated_nbinom2(),
  ziformula = ~0
)
robust_fixed_trunc_model <- model_parameters(fixed_trunc_model, robust = TRUE)
print(robust_fixed_trunc_model)
r2(fixed_trunc_model)

random_trunc_model <- glmmTMB(
  OFFERS ~  G_CLAUSE + G_WEIGHT + ALLOTMENT + FRAMEWORK_AGREEMENT + 
    P_CRITERION_WEIGHT_ENV + REGION + STATUS + (1 | CPV),
  data = df,
  weights = w_samp,
  family = truncated_nbinom2(),
  ziformula = ~0  
)
robust_random_trunc_model <- model_parameters(random_trunc_model, robust = TRUE)
r2(random_trunc_model)
print(robust_random_trunc_model)

random_trunc_model_quad <- glmmTMB(
  OFFERS ~  G_CLAUSE + G_WEIGHT + G_WEIGHT_SQUARED + ALLOTMENT + FRAMEWORK_AGREEMENT + 
    P_CRITERION_WEIGHT_ENV + REGION + STATUS + (1 | CPV),
  data = df,
  weights = w_samp,
  family = truncated_nbinom2(),
  ziformula = ~0  
)
robust_random_trunc_model_quad <- model_parameters(random_trunc_model_quad, robust = TRUE)
print(robust_random_trunc_model_quad)
r2(random_trunc_model_quad)




# Appendix D

set.seed(123)

df <- read_csv("data_without_price.csv")

df <- df %>%
  mutate(
    CPV = gsub("'", "", CPV),
    CPV = trimws(CPV)
  ) %>%
  filter(substr(CPV, 3, 3) != "0") %>%
  mutate(CPV = substr(CPV, 1, 3))

df <- df %>%
  group_by(CPV) %>%
  filter(n() >= 30) %>%
  ungroup()

df <- df %>%
  mutate(
    G_WEIGHT_SQUARED = G_WEIGHT^2
  )

df <- df %>%
  mutate(
    CPV = as.factor(CPV),
    REGION = as.factor(REGION),
    STATUS = as.factor(STATUS)
  )

df <- df %>%
  mutate(
    CPV = gsub("'", "", CPV),
    CPV = trimws(CPV)
  ) %>%
  filter(substr(CPV, 3, 3) != "0") %>%
  mutate(CPV = substr(CPV, 1, 3)) %>%
  mutate(CPV2 = substr(CPV, 1, 2))

df <- df %>%
  group_by(CAE_SIREN) %>%
  mutate(n_siren = n()) %>%
  ungroup() %>%
  mutate(w_inv = 1 / n_siren)
df <- df %>%
  mutate(w_samp = w_inv * (nrow(df) / sum(w_inv)))

plan(multisession, workers = 30) # adapt to your computer

sample_one_per_contract <- function(data, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  data %>%
    group_by(ID_CONTRACT) %>%
    slice_sample(n = 1) %>%
    ungroup()
}

estimate_model <- function(sampled_data) {
  tryCatch({
    model <- glmmTMB(
      OFFERS ~ G_CLAUSE + 
        G_WEIGHT + G_WEIGHT_SQUARED +
        ALLOTMENT + FRAMEWORK_AGREEMENT + P_CRITERION_WEIGHT_ENV + 
        REGION + STATUS + (1 | CPV),
      data = sampled_data,
      family = truncated_nbinom2(),
      weights = w_samp,
      ziformula = ~0
    )
    params <- model_parameters(model, robust = TRUE)
    vars_interet <- c("G_CLAUSE",
                      "G_WEIGHT", 
                      "G_WEIGHT_SQUARED")
    params %>%
      filter(Parameter %in% vars_interet) %>%
      select(Parameter, Coefficient, SE, p) %>%
      mutate(converged = TRUE, n_obs = nrow(sampled_data))
  }, error = function(e) {
    tibble(
      Parameter = c("ENVIRONMENTAL_CLAUSE",
                    "ENVIRONMENTAL_WEIGHT", "ENV_WEIGHT_SQUARED"),
      Coefficient = NA_real_, SE = NA_real_, p = NA_real_,
      converged = FALSE, n_obs = NA_integer_
    )
  })
}

bootstrap_results <- future_map_dfr(
  1:1000,
  ~{
    sampled <- sample_one_per_contract(df, seed = .x)
    estimate_model(sampled)
  },
  .options = furrr_options(seed = TRUE),
  .progress = TRUE
)

results_clean <- bootstrap_results %>%
  filter(converged == TRUE, !is.na(Coefficient))

summary_bootstrap <- results_clean %>%
  group_by(Parameter) %>%
  summarise(
    mean    = mean(Coefficient),
    sd      = sd(Coefficient),
    ci95_lo = quantile(Coefficient, 0.025),
    ci95_hi = quantile(Coefficient, 0.975),
    ci99_lo = quantile(Coefficient, 0.005),
    ci99_hi = quantile(Coefficient, 0.995),
    sig_rate = mean(p < 0.05, na.rm = TRUE) * 100,
    n_iter  = n()
  )

print(summary_bootstrap)

param_labels <- c(
  "G_CLAUSE"    = "G_clause",
  "G_WEIGHT"    = "G_weight",
  "G_WEIGHT_SQUARED"      = "G_weight²"
)

results_clean_plot <- results_clean %>%
  mutate(Parameter = factor(Parameter, levels = names(param_labels)))

summary_bootstrap_plot <- summary_bootstrap %>%
  mutate(Parameter = factor(Parameter, levels = names(param_labels)))

x_limits <- results_clean_plot %>%
  group_by(Parameter) %>%
  summarise(
    xmin = min(min(Coefficient), 0),
    xmax = max(max(Coefficient), 0)
  )

ggplot(results_clean_plot, aes(x = Coefficient)) +
  geom_histogram(bins = 50, fill = "grey60", color = "white", alpha = 0.8) +
  geom_vline(
    data = summary_bootstrap_plot,
    aes(xintercept = mean),
    color = "black", linetype = "solid", linewidth = 0.8
  ) +
  geom_vline(
    data = summary_bootstrap_plot,
    aes(xintercept = ci95_lo), color = "black", linetype = "dashed"
  ) +
  geom_vline(
    data = summary_bootstrap_plot,
    aes(xintercept = ci95_hi), color = "black", linetype = "dashed"
  ) +
  geom_blank(data = x_limits, aes(x = xmin)) +
  geom_blank(data = x_limits, aes(x = xmax)) +
  facet_wrap(~Parameter, scales = "free", labeller = as_labeller(param_labels)) +
  labs(
    x = "Estimated coefficient", y = "Frequency"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(color = "black"),
    plot.subtitle = element_text(color = "black"),
    panel.grid.minor = element_blank()
  )


# Appendices E & F

set.seed(123)

df <- read_csv("data_without_price.csv")

df <- df %>%
  mutate(
    CPV = gsub("'", "", CPV),
    CPV = trimws(CPV)
  ) %>%
  filter(substr(CPV, 3, 3) != "0") %>%
  mutate(CPV = substr(CPV, 1, 3))

df <- df %>%
  group_by(CPV) %>%
  filter(n() >= 30) %>%
  ungroup()

df <- df %>%
  mutate(
    G_WEIGHT_SQUARED = G_WEIGHT^2
  )

df <- df %>%
  mutate(
    CPV = as.factor(CPV),
    REGION = as.factor(REGION),
    STATUS = as.factor(STATUS)
  )

df <- df %>%
  mutate(
    CPV = gsub("'", "", CPV),
    CPV = trimws(CPV)
  ) %>%
  filter(substr(CPV, 3, 3) != "0") %>%
  mutate(CPV = substr(CPV, 1, 3)) %>%
  mutate(CPV2 = substr(CPV, 1, 2))

df <- df %>%
  group_by(CAE_SIREN) %>%
  mutate(n_siren = n()) %>%
  ungroup() %>%
  mutate(w_inv = 1 / n_siren)
df <- df %>%
  mutate(w_samp = w_inv * (nrow(df) / sum(w_inv)))

df <- df %>%
  mutate(
    treat_gc = case_when(
      G_CLAUSE == 0 & G_CRITERION == 0 ~ 0,
      G_CLAUSE == 1 & G_CRITERION == 0 ~ 1,
      G_CLAUSE == 0 & G_CRITERION == 1 ~ 2,
      G_CLAUSE == 1 & G_CRITERION == 1 ~ 3
    ),
    treat_gc = factor(treat_gc, 
                      levels = 0:3, 
                      labels = c("none", "clause only", "criterion only", "both"))
  )

mlogit_fit <- multinom(
  treat_gc ~ ALLOTMENT  + FRAMEWORK_AGREEMENT + P_CRITERION_WEIGHT_ENV + 
    factor(CPV) + factor(REGION) + factor(STATUS),
  data = df
)

pr <- predict(mlogit_fit, type = "probs")
df <- df %>%
  mutate(
    pr0 = pr[, "none"],
    pr1 = pr[, "clause only"],
    pr2 = pr[, "criterion only"],
    pr3 = pr[, "both"]
  )

p_treat <- prop.table(table(df$treat_gc))
df <- df %>%
  mutate(
    ipw_gc = case_when(
      treat_gc == "none"          ~ p_treat["none"] / pr0,
      treat_gc == "clause only"   ~ p_treat["clause only"] / pr1,
      treat_gc == "criterion only"~ p_treat["criterion only"] / pr2,
      treat_gc == "both"          ~ p_treat["both"] / pr3
    ),
    w_final_gc = ipw_gc * w_samp,
    w_final_norm_gc = w_final_gc * (n() / sum(w_final_gc, na.rm = TRUE))
  )

random_intercept_dr_gc <- glmmTMB(
  OFFERS ~ G_CLAUSE + G_CRITERION + ALLOTMENT +
    FRAMEWORK_AGREEMENT + P_CRITERION_WEIGHT_ENV +
    factor(REGION) + factor(STATUS) +
    (1 | CPV),
  data = df,
  family = nbinom2(),
  weights = w_final_norm_gc,
  ziformula = ~0
)
robust_random_intercept_dr_gc <- model_parameters(random_intercept_dr_gc, robust = TRUE)
robust_random_intercept_dr_gc

df <- df %>%
  mutate(
    weight_bins = case_when(
      G_WEIGHT == 0 ~ "no_weight",
      G_WEIGHT > 0 & G_WEIGHT <= 10 ~ "low",
      G_WEIGHT > 10 & G_WEIGHT <= 30 ~ "medium",
      G_WEIGHT > 30 ~ "high"
    ),
    weight_bins = factor(weight_bins, levels = c("no_weight", "low", "medium", "high"))
  )

df <- df %>%
  mutate(
    treat_clause_cont = case_when(
      G_CLAUSE == 0 & weight_bins == "no_weight" ~ 0,
      G_CLAUSE == 1 & weight_bins == "no_weight" ~ 1,
      G_CLAUSE == 0 & weight_bins == "low" ~ 2,
      G_CLAUSE == 1 & weight_bins == "low" ~ 3,
      G_CLAUSE == 0 & weight_bins == "medium" ~ 4,
      G_CLAUSE == 1 & weight_bins == "medium" ~ 5,
      G_CLAUSE == 0 & weight_bins == "high" ~ 6,
      G_CLAUSE == 1 & weight_bins == "high" ~ 7
    ),
    treat_clause_cont = factor(treat_clause_cont,
                               levels = 0:7,
                               labels = c("none", "clause_only", "weight_low", "clause_weight_low",
                                          "weight_medium", "clause_weight_medium", "weight_high", 
                                          "clause_weight_high"))
  )

table(df$treat_clause_cont)

mlogit_fit_cont <- multinom(
  treat_clause_cont ~ ALLOTMENT + FRAMEWORK_AGREEMENT + P_CRITERION_WEIGHT_ENV + 
    factor(CPV2) + factor(REGION) + factor(STATUS),
  data = df,
  maxit = 1000
)

pr_cont <- predict(mlogit_fit_cont, type = "probs")

if (is.vector(pr_cont)) {
  pr_cont <- matrix(pr_cont, ncol = 1)
  colnames(pr_cont) <- levels(df$treat_clause_cont)[2]
}

pr_df_cont <- as.data.frame(pr_cont)

all_levels_cont <- levels(df$treat_clause_cont)

for (level in all_levels_cont) {
  if (!level %in% colnames(pr_df_cont)) {
    pr_df_cont[[level]] <- 0.001
  }
}

pr_df_cont <- pr_df_cont[, all_levels_cont, drop = FALSE]

colnames(pr_df_cont) <- paste0(colnames(pr_df_cont), "_cont")

df <- df %>%
  bind_cols(pr_df_cont)

p_treat_cont <- prop.table(table(df$treat_clause_cont))

new_cols <- colnames(df)[grepl("_cont$", colnames(df))]

df <- df %>%
  mutate(
    ipw_clause_cont = case_when(
      treat_clause_cont == "none" ~ p_treat_cont["none"] / pmax(none_cont, 0.001),
      treat_clause_cont == "clause_only" ~ p_treat_cont["clause_only"] / pmax(clause_only_cont, 0.001),
      treat_clause_cont == "weight_low" ~ p_treat_cont["weight_low"] / pmax(weight_low_cont, 0.001),
      treat_clause_cont == "clause_weight_low" ~ p_treat_cont["clause_weight_low"] / pmax(clause_weight_low_cont, 0.001),
      treat_clause_cont == "weight_medium" ~ p_treat_cont["weight_medium"] / pmax(weight_medium_cont, 0.001),
      treat_clause_cont == "clause_weight_medium" ~ p_treat_cont["clause_weight_medium"] / pmax(clause_weight_medium_cont, 0.001),
      treat_clause_cont == "weight_high" ~ p_treat_cont["weight_high"] / pmax(weight_high_cont, 0.001),
      treat_clause_cont == "clause_weight_high" ~ p_treat_cont["clause_weight_high"] / pmax(clause_weight_high_cont, 0.001),
      TRUE ~ 1
    ),
    ipw_clause_cont = pmin(pmax(ipw_clause_cont, 0.1), 10),
    w_final_clause_cont = ipw_clause_cont * w_samp,
    w_final_norm_clause_cont = w_final_clause_cont * (n() / sum(w_final_clause_cont, na.rm = TRUE))
  )

random_intercept_dr_cont <- glmmTMB(
  OFFERS ~ G_CLAUSE + G_WEIGHT + 
    ALLOTMENT + FRAMEWORK_AGREEMENT + P_CRITERION_WEIGHT_ENV +
    factor(REGION) + factor(STATUS) + (1 | CPV),
  data = df,
  family = truncated_nbinom2(),
  weights = w_final_norm_clause_cont,
  ziformula = ~0
)

robust_random_intercept_dr_cont <- model_parameters(random_intercept_dr_cont, robust = TRUE)
robust_random_intercept_dr_cont

random_intercept_dr_cont_quad <- glmmTMB(
  OFFERS ~ G_CLAUSE + G_WEIGHT + G_WEIGHT_SQUARED +
    ALLOTMENT + FRAMEWORK_AGREEMENT + P_CRITERION_WEIGHT_ENV +
    factor(REGION) + factor(STATUS) + (1 | CPV),
  data = df,
  family = truncated_nbinom2(),
  weights = w_final_norm_clause_cont,
  ziformula = ~0
)

robust_random_intercept_dr_cont_quad <- model_parameters(random_intercept_dr_cont_quad, robust = TRUE)
robust_random_intercept_dr_cont_quad



# Appendix G

set.seed(123)

df <- read_csv("data_with_price.csv")

df <- df %>%
  mutate(
    CPV = gsub("'", "", CPV),
    CPV = trimws(CPV)
  ) %>%
  filter(substr(CPV, 3, 3) != "0") %>%
  mutate(CPV = substr(CPV, 1, 3))

df <- df %>%
  group_by(CPV) %>%
  filter(n() >= 30) %>%
  ungroup()

df <- df %>%
  mutate(
    G_WEIGHT_SQUARED = G_WEIGHT^2
  )

df <- df %>%
  mutate(
    CPV = as.factor(CPV),
    REGION = as.factor(REGION),
    STATUS = as.factor(STATUS)
  )

df <- df %>%
  mutate(
    CPV = gsub("'", "", CPV),
    CPV = trimws(CPV)
  ) %>%
  filter(substr(CPV, 3, 3) != "0") %>%
  mutate(CPV = substr(CPV, 1, 3)) %>%
  mutate(CPV2 = substr(CPV, 1, 2))

df <- df %>%
  group_by(CAE_SIREN) %>%
  mutate(n_siren = n()) %>%
  ungroup() %>%
  mutate(w_inv = 1 / n_siren)
df <- df %>%
  mutate(w_samp = w_inv * (nrow(df) / sum(w_inv)))


fixed_trunc_model_price <- glmmTMB(
  OFFERS ~ G_CLAUSE + G_WEIGHT + AWARD_PRICE + ALLOTMENT + FRAMEWORK_AGREEMENT + 
    P_CRITERION_WEIGHT_ENV + REGION + STATUS + factor(CPV2),
  data = df,
  weights = w_samp,
  family = truncated_nbinom2(),
  ziformula = ~0
)
robust_fixed_trunc_model_price <- model_parameters(fixed_trunc_model_price, robust = TRUE)
print(robust_fixed_trunc_model_price)
r2(fixed_trunc_model_price)

random_trunc_model_price <- glmmTMB(
  OFFERS ~  G_CLAUSE + G_WEIGHT + AWARD_PRICE + ALLOTMENT + FRAMEWORK_AGREEMENT + 
    P_CRITERION_WEIGHT_ENV + REGION + STATUS + (1 | CPV),
  data = df,
  weights = w_samp,
  family = truncated_nbinom2(),
  ziformula = ~0  
)
robust_random_trunc_model_price <- model_parameters(random_trunc_model_price, robust = TRUE)
r2(random_trunc_model_price)
print(robust_random_trunc_model_price)

