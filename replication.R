# Table-by-table replication
# Place this script and the three CSV files in the same folder.
# Run the blocks in R/RStudio, or execute:
# source("replication.R", print.eval = TRUE)
# No automatic exports: results and fitted models remain in the R session.

library(dplyr)
library(glmmTMB)
library(readr)
library(parameters)
library(nnet)
library(ggplot2)
library(future)
library(furrr)
library(performance)

set.seed(123)

# TABLE 2 — Main results

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


# TABLE 1 — Descriptive statistics for the baseline sample

variables_descriptives <- df %>% select(
  OFFERS, G_CLAUSE, G_CRITERION, G_WEIGHT, P_CRITERION_WEIGHT_ENV,
  ALLOTMENT, FRAMEWORK_AGREEMENT
)

statistiques_descriptives <- data.frame(
  Variable = names(variables_descriptives),
  Moyenne = sapply(variables_descriptives, mean),
  Ecart_type = sapply(variables_descriptives, sd),
  Observations = sapply(variables_descriptives, function(x) sum(!is.na(x))),
  row.names = NULL
)
print(statistiques_descriptives)
mean(df$G_WEIGHT[df$G_CRITERION == 1])
table(df$STATUS)
table(df$REGION)

# TABLE 2 — Main estimates

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





# TABLE 3 — Sector heterogeneity
# Baseline sample; three-digit CPV fixed effects.

# CPV 45
df_45 <- df %>% filter(CPV2 == "45")
df_45 <- df_45 %>%
  group_by(CAE_SIREN) %>% mutate(n_siren = n()) %>% ungroup() %>%
  mutate(w_inv = 1 / n_siren,
         w_samp = w_inv * (nrow(df_45) / sum(w_inv)))

sector_45 <- glmmTMB(
  OFFERS ~ G_CLAUSE + G_WEIGHT + ALLOTMENT + FRAMEWORK_AGREEMENT +
    P_CRITERION_WEIGHT_ENV + REGION + STATUS + factor(CPV),
  data = df_45, weights = w_samp,
  family = truncated_nbinom2(), ziformula = ~0
)
robust_sector_45 <- model_parameters(sector_45, robust = TRUE)
print(robust_sector_45)
nobs(sector_45)

# CPV 90
df_90 <- df %>% filter(CPV2 == "90")
df_90 <- df_90 %>%
  group_by(CAE_SIREN) %>% mutate(n_siren = n()) %>% ungroup() %>%
  mutate(w_inv = 1 / n_siren,
         w_samp = w_inv * (nrow(df_90) / sum(w_inv)))

sector_90 <- glmmTMB(
  OFFERS ~ G_CLAUSE + G_WEIGHT + ALLOTMENT + FRAMEWORK_AGREEMENT +
    P_CRITERION_WEIGHT_ENV + REGION + STATUS + factor(CPV),
  data = df_90, weights = w_samp,
  family = truncated_nbinom2(), ziformula = ~0
)
robust_sector_90 <- model_parameters(sector_90, robust = TRUE)
print(robust_sector_90)
nobs(sector_90)

# CPV 33
df_33 <- df %>% filter(CPV2 == "33")
df_33 <- df_33 %>%
  group_by(CAE_SIREN) %>% mutate(n_siren = n()) %>% ungroup() %>%
  mutate(w_inv = 1 / n_siren,
         w_samp = w_inv * (nrow(df_33) / sum(w_inv)))

sector_33 <- glmmTMB(
  OFFERS ~ G_CLAUSE + G_WEIGHT + ALLOTMENT + FRAMEWORK_AGREEMENT +
    P_CRITERION_WEIGHT_ENV + REGION + STATUS + factor(CPV),
  data = df_33, weights = w_samp,
  family = truncated_nbinom2(), ziformula = ~0
)
robust_sector_33 <- model_parameters(sector_33, robust = TRUE)
print(robust_sector_33)
nobs(sector_33)

# APPENDIX B — Correlation matrix
correlation_matrix <- cor(df %>% select(
  OFFERS, G_CLAUSE, G_WEIGHT, P_CRITERION_WEIGHT_ENV,
  ALLOTMENT, FRAMEWORK_AGREEMENT
))
print(round(correlation_matrix, 2))

# APPENDIX C — Heterogeneity by buyer type
# Weights recalculated within each buyer type; CPV random intercept.
# STATUS is omitted from the formula because it is constant within the subsample.

# Department
df_departement <- df %>% filter(STATUS == "Departement") %>%
  mutate(CPV = droplevels(factor(CPV)), REGION = droplevels(REGION)) %>%
  group_by(CAE_SIREN) %>% mutate(n_siren = n()) %>% ungroup() %>%
  mutate(w_inv = 1 / n_siren,
         w_samp = w_inv * (n() / sum(w_inv)))

buyer_departement <- glmmTMB(
  OFFERS ~ G_CLAUSE + G_WEIGHT + ALLOTMENT + FRAMEWORK_AGREEMENT +
    P_CRITERION_WEIGHT_ENV + REGION + (1 | CPV),
  data = df_departement, weights = w_samp,
  family = truncated_nbinom2(), ziformula = ~0
)
# Same optimizer fallback as in the supplied heterogeneity script.
if (!isTRUE(buyer_departement$sdr$pdHess)) {
  buyer_departement_retry <- tryCatch(
    update(buyer_departement, control = glmmTMBControl(
      optimizer = optim, optArgs = list(method = "BFGS"))),
    error = function(e) NULL
  )
  if (!is.null(buyer_departement_retry) && isTRUE(buyer_departement_retry$sdr$pdHess)) {
    buyer_departement <- buyer_departement_retry
  }
}
robust_buyer_departement <- model_parameters(buyer_departement, robust = TRUE)
print(robust_buyer_departement)
nobs(buyer_departement)

# Local agency
df_local_agency <- df %>% filter(STATUS == "Local agency") %>%
  mutate(CPV = droplevels(factor(CPV)), REGION = droplevels(REGION)) %>%
  group_by(CAE_SIREN) %>% mutate(n_siren = n()) %>% ungroup() %>%
  mutate(w_inv = 1 / n_siren,
         w_samp = w_inv * (n() / sum(w_inv)))

buyer_local_agency <- glmmTMB(
  OFFERS ~ G_CLAUSE + G_WEIGHT + ALLOTMENT + FRAMEWORK_AGREEMENT +
    P_CRITERION_WEIGHT_ENV + REGION + (1 | CPV),
  data = df_local_agency, weights = w_samp,
  family = truncated_nbinom2(), ziformula = ~0
)
# Same optimizer fallback as in the supplied heterogeneity script.
if (!isTRUE(buyer_local_agency$sdr$pdHess)) {
  buyer_local_agency_retry <- tryCatch(
    update(buyer_local_agency, control = glmmTMBControl(
      optimizer = optim, optArgs = list(method = "BFGS"))),
    error = function(e) NULL
  )
  if (!is.null(buyer_local_agency_retry) && isTRUE(buyer_local_agency_retry$sdr$pdHess)) {
    buyer_local_agency <- buyer_local_agency_retry
  }
}
robust_buyer_local_agency <- model_parameters(buyer_local_agency, robust = TRUE)
print(robust_buyer_local_agency)
nobs(buyer_local_agency)

# Municipal federation
df_municipal_federation <- df %>% filter(STATUS == "Municipal federation") %>%
  mutate(CPV = droplevels(factor(CPV)), REGION = droplevels(REGION)) %>%
  group_by(CAE_SIREN) %>% mutate(n_siren = n()) %>% ungroup() %>%
  mutate(w_inv = 1 / n_siren,
         w_samp = w_inv * (n() / sum(w_inv)))

buyer_municipal_federation <- glmmTMB(
  OFFERS ~ G_CLAUSE + G_WEIGHT + ALLOTMENT + FRAMEWORK_AGREEMENT +
    P_CRITERION_WEIGHT_ENV + REGION + (1 | CPV),
  data = df_municipal_federation, weights = w_samp,
  family = truncated_nbinom2(), ziformula = ~0
)
# Same optimizer fallback as in the supplied heterogeneity script.
if (!isTRUE(buyer_municipal_federation$sdr$pdHess)) {
  buyer_municipal_federation_retry <- tryCatch(
    update(buyer_municipal_federation, control = glmmTMBControl(
      optimizer = optim, optArgs = list(method = "BFGS"))),
    error = function(e) NULL
  )
  if (!is.null(buyer_municipal_federation_retry) && isTRUE(buyer_municipal_federation_retry$sdr$pdHess)) {
    buyer_municipal_federation <- buyer_municipal_federation_retry
  }
}
robust_buyer_municipal_federation <- model_parameters(buyer_municipal_federation, robust = TRUE)
print(robust_buyer_municipal_federation)
nobs(buyer_municipal_federation)

# Municipality
df_municipality <- df %>% filter(STATUS == "Municipality") %>%
  mutate(CPV = droplevels(factor(CPV)), REGION = droplevels(REGION)) %>%
  group_by(CAE_SIREN) %>% mutate(n_siren = n()) %>% ungroup() %>%
  mutate(w_inv = 1 / n_siren,
         w_samp = w_inv * (n() / sum(w_inv)))

buyer_municipality <- glmmTMB(
  OFFERS ~ G_CLAUSE + G_WEIGHT + ALLOTMENT + FRAMEWORK_AGREEMENT +
    P_CRITERION_WEIGHT_ENV + REGION + (1 | CPV),
  data = df_municipality, weights = w_samp,
  family = truncated_nbinom2(), ziformula = ~0
)
# Same optimizer fallback as in the supplied heterogeneity script.
if (!isTRUE(buyer_municipality$sdr$pdHess)) {
  buyer_municipality_retry <- tryCatch(
    update(buyer_municipality, control = glmmTMBControl(
      optimizer = optim, optArgs = list(method = "BFGS"))),
    error = function(e) NULL
  )
  if (!is.null(buyer_municipality_retry) && isTRUE(buyer_municipality_retry$sdr$pdHess)) {
    buyer_municipality <- buyer_municipality_retry
  }
}
robust_buyer_municipality <- model_parameters(buyer_municipality, robust = TRUE)
print(robust_buyer_municipality)
nobs(buyer_municipality)

# National agency
df_national_agency <- df %>% filter(STATUS == "National agency") %>%
  mutate(CPV = droplevels(factor(CPV)), REGION = droplevels(REGION)) %>%
  group_by(CAE_SIREN) %>% mutate(n_siren = n()) %>% ungroup() %>%
  mutate(w_inv = 1 / n_siren,
         w_samp = w_inv * (n() / sum(w_inv)))

buyer_national_agency <- glmmTMB(
  OFFERS ~ G_CLAUSE + G_WEIGHT + ALLOTMENT + FRAMEWORK_AGREEMENT +
    P_CRITERION_WEIGHT_ENV + REGION + (1 | CPV),
  data = df_national_agency, weights = w_samp,
  family = truncated_nbinom2(), ziformula = ~0
)
# Same optimizer fallback as in the supplied heterogeneity script.
if (!isTRUE(buyer_national_agency$sdr$pdHess)) {
  buyer_national_agency_retry <- tryCatch(
    update(buyer_national_agency, control = glmmTMBControl(
      optimizer = optim, optArgs = list(method = "BFGS"))),
    error = function(e) NULL
  )
  if (!is.null(buyer_national_agency_retry) && isTRUE(buyer_national_agency_retry$sdr$pdHess)) {
    buyer_national_agency <- buyer_national_agency_retry
  }
}
robust_buyer_national_agency <- model_parameters(buyer_national_agency, robust = TRUE)
print(robust_buyer_national_agency)
nobs(buyer_national_agency)

# Region
df_region <- df %>% filter(STATUS == "Region") %>%
  mutate(CPV = droplevels(factor(CPV)), REGION = droplevels(REGION)) %>%
  group_by(CAE_SIREN) %>% mutate(n_siren = n()) %>% ungroup() %>%
  mutate(w_inv = 1 / n_siren,
         w_samp = w_inv * (n() / sum(w_inv)))

buyer_region <- glmmTMB(
  OFFERS ~ G_CLAUSE + G_WEIGHT + ALLOTMENT + FRAMEWORK_AGREEMENT +
    P_CRITERION_WEIGHT_ENV + REGION + (1 | CPV),
  data = df_region, weights = w_samp,
  family = truncated_nbinom2(), ziformula = ~0
)
# Same optimizer fallback as in the supplied heterogeneity script.
if (!isTRUE(buyer_region$sdr$pdHess)) {
  buyer_region_retry <- tryCatch(
    update(buyer_region, control = glmmTMBControl(
      optimizer = optim, optArgs = list(method = "BFGS"))),
    error = function(e) NULL
  )
  if (!is.null(buyer_region_retry) && isTRUE(buyer_region_retry$sdr$pdHess)) {
    buyer_region <- buyer_region_retry
  }
}
robust_buyer_region <- model_parameters(buyer_region, robust = TRUE)
print(robust_buyer_region)
nobs(buyer_region)

# State
df_state <- df %>% filter(STATUS == "State") %>%
  mutate(CPV = droplevels(factor(CPV)), REGION = droplevels(REGION)) %>%
  group_by(CAE_SIREN) %>% mutate(n_siren = n()) %>% ungroup() %>%
  mutate(w_inv = 1 / n_siren,
         w_samp = w_inv * (n() / sum(w_inv)))

buyer_state <- glmmTMB(
  OFFERS ~ G_CLAUSE + G_WEIGHT + ALLOTMENT + FRAMEWORK_AGREEMENT +
    P_CRITERION_WEIGHT_ENV + REGION + (1 | CPV),
  data = df_state, weights = w_samp,
  family = truncated_nbinom2(), ziformula = ~0
)
# Same optimizer fallback as in the supplied heterogeneity script.
if (!isTRUE(buyer_state$sdr$pdHess)) {
  buyer_state_retry <- tryCatch(
    update(buyer_state, control = glmmTMBControl(
      optimizer = optim, optArgs = list(method = "BFGS"))),
    error = function(e) NULL
  )
  if (!is.null(buyer_state_retry) && isTRUE(buyer_state_retry$sdr$pdHess)) {
    buyer_state <- buyer_state_retry
  }
}
robust_buyer_state <- model_parameters(buyer_state, robust = TRUE)
print(robust_buyer_state)
nobs(buyer_state)

# APPENDIX C — Coefficient heterogeneity Q tests
# Tests based on separate estimates, assumed to be independent.
buyer_coefficients <- bind_rows(
  as.data.frame(robust_buyer_departement) %>% mutate(buyer_type = "departement",
    convergence_code = buyer_departement$fit$convergence,
    positive_definite_hessian = isTRUE(buyer_departement$sdr$pdHess)),
  as.data.frame(robust_buyer_local_agency) %>% mutate(buyer_type = "local_agency",
    convergence_code = buyer_local_agency$fit$convergence,
    positive_definite_hessian = isTRUE(buyer_local_agency$sdr$pdHess)),
  as.data.frame(robust_buyer_municipal_federation) %>% mutate(buyer_type = "municipal_federation",
    convergence_code = buyer_municipal_federation$fit$convergence,
    positive_definite_hessian = isTRUE(buyer_municipal_federation$sdr$pdHess)),
  as.data.frame(robust_buyer_municipality) %>% mutate(buyer_type = "municipality",
    convergence_code = buyer_municipality$fit$convergence,
    positive_definite_hessian = isTRUE(buyer_municipality$sdr$pdHess)),
  as.data.frame(robust_buyer_national_agency) %>% mutate(buyer_type = "national_agency",
    convergence_code = buyer_national_agency$fit$convergence,
    positive_definite_hessian = isTRUE(buyer_national_agency$sdr$pdHess)),
  as.data.frame(robust_buyer_region) %>% mutate(buyer_type = "region",
    convergence_code = buyer_region$fit$convergence,
    positive_definite_hessian = isTRUE(buyer_region$sdr$pdHess)),
  as.data.frame(robust_buyer_state) %>% mutate(buyer_type = "state",
    convergence_code = buyer_state$fit$convergence,
    positive_definite_hessian = isTRUE(buyer_state$sdr$pdHess))
)

q_tests <- list()
for (parameter in c("G_CLAUSE", "G_WEIGHT")) {
  estimates <- buyer_coefficients %>% filter(
    Parameter == parameter, convergence_code == 0,
    positive_definite_hessian, is.finite(Coefficient), is.finite(SE), SE > 0
  )
  k <- nrow(estimates)
  if (k >= 2) {
    inverse_variance <- 1 / estimates$SE^2
    pooled <- sum(inverse_variance * estimates$Coefficient) / sum(inverse_variance)
    Q <- sum(inverse_variance * (estimates$Coefficient - pooled)^2)
    q_tests[[parameter]] <- data.frame(
      Parameter = parameter, Groups = k, Pooled_coefficient = pooled,
      Q = Q, df = k - 1, p = pchisq(Q, df = k - 1, lower.tail = FALSE)
    )
  } else {
    q_tests[[parameter]] <- data.frame(
      Parameter = parameter, Groups = k, Pooled_coefficient = NA_real_,
      Q = NA_real_, df = NA_integer_, p = NA_real_
    )
  }
}
print(bind_rows(q_tests))
# APPENDIX D — One-lot-per-contract draws (1,000 replications)

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

plan(multisession, workers = 4) # Adjust to the available number of cores and memory

# Weights calculated on the full sample are retained after each draw.
# The two functions below come from the original resampling block.
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
      Parameter = c("G_CLAUSE", "G_WEIGHT", "G_WEIGHT_SQUARED"),
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

# The solid line represents the mean, as in the original code.
# The manuscript caption referring to the median needs to be aligned separately.
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



plan(sequential)

# APPENDICES E AND F — Propensity models and reweighted estimates

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

# APPENDIX E — Multinomial regression: binary instruments
summary(mlogit_fit)

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

# APPENDIX F — Reweighted estimates: binary instruments
# The NONTRUNCATED family is retained, matching the original script.
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
print(robust_random_intercept_dr_gc)

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

# APPENDIX E — Multinomial regression: environmental-weight categories
summary(mlogit_fit_cont)

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

# APPENDIX F — Reweighted estimates: linear and quadratic environmental weight
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
print(robust_random_intercept_dr_cont)

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
print(robust_random_intercept_dr_cont_quad)



# APPENDIX G — Controlling for lot price
# AWARD_PRICE is already logged in the CSV: do not transform it again.

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

set.seed(123)

# APPENDIX H — Sensitivity to the bid-count restriction

df <- read_csv("data_without_price_extended.csv")

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


# Use the ENTIRE extended dataset, not only observations above ten offers.
# The upstream contract-level filter is not reapplied to this prepared CSV.

extended_fixed_trunc_model <- glmmTMB(
  OFFERS ~ G_CLAUSE + G_WEIGHT + ALLOTMENT + FRAMEWORK_AGREEMENT + 
    P_CRITERION_WEIGHT_ENV + REGION + STATUS + factor(CPV2),
  data = df,
  weights = w_samp,
  family = truncated_nbinom2(),
  ziformula = ~0
)
robust_extended_fixed_trunc_model <- model_parameters(extended_fixed_trunc_model, robust = TRUE)
print(robust_extended_fixed_trunc_model)
r2(extended_fixed_trunc_model)

extended_random_trunc_model <- glmmTMB(
  OFFERS ~  G_CLAUSE + G_WEIGHT + ALLOTMENT + FRAMEWORK_AGREEMENT + 
    P_CRITERION_WEIGHT_ENV + REGION + STATUS + (1 | CPV),
  data = df,
  weights = w_samp,
  family = truncated_nbinom2(),
  ziformula = ~0  
)
robust_extended_random_trunc_model <- model_parameters(extended_random_trunc_model, robust = TRUE)
r2(extended_random_trunc_model)
print(robust_extended_random_trunc_model)

extended_random_trunc_model_quad <- glmmTMB(
  OFFERS ~  G_CLAUSE + G_WEIGHT + G_WEIGHT_SQUARED + ALLOTMENT + FRAMEWORK_AGREEMENT + 
    P_CRITERION_WEIGHT_ENV + REGION + STATUS + (1 | CPV),
  data = df,
  weights = w_samp,
  family = truncated_nbinom2(),
  ziformula = ~0  
)
robust_extended_random_trunc_model_quad <- model_parameters(extended_random_trunc_model_quad, robust = TRUE)
print(robust_extended_random_trunc_model_quad)
r2(extended_random_trunc_model_quad)
