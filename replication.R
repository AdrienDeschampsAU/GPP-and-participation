# GPP-and-participation: unified replication script (initial analyses + revision).
# Run from the repository root: Rscript replication.R
# Or specify the folder: Rscript replication.R /path/to/repository
# In R/RStudio, set the working directory and run source("replication.R").
# Inputs: data_without_price.csv, data_with_price.csv,
#         data_without_price_extended.csv. Outputs: results_replication/.
# Requires R packages listed below. Original package used R 4.5.1.
# Numerical results have not yet been verified by executing this unified script.
#
# Specifications are preserved, including:
# - untruncated nbinom2 for the binary reweighted outcome model;
# - truncated_nbinom2 for the other count models;
# - baseline sampling weights retained after drawing one lot per contract;
# - AWARD_PRICE already logged in its CSV: do NOT log it a second time.
#
# Optional environment settings (defaults execute ALL analyses):
# Sys.setenv(GPP_WORKERS = "2", GPP_BOOTSTRAP_REPS = "1000")
# Sys.setenv(GPP_RUN_BOOTSTRAP = "false") # optional partial run only
# A complete replication uses GPP_RUN_BOOTSTRAP=true and 1000 replications.


packages <- c("dplyr", "readr", "glmmTMB", "parameters", "performance",
              "nnet", "ggplot2", "future", "furrr")
missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Install the required R packages first: ", paste(missing_packages, collapse = ", "))
}
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(glmmTMB)
})
set.seed(123)
args <- commandArgs(trailingOnly = TRUE)
repo_dir <- normalizePath(if (length(args)) args[[1]] else ".", mustWork = TRUE)
output_dir <- file.path(repo_dir, "results_replication")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

required <- c("OFFERS", "G_CLAUSE", "G_WEIGHT", "ALLOTMENT",
              "FRAMEWORK_AGREEMENT", "P_CRITERION_WEIGHT_ENV",
              "REGION", "STATUS", "CPV", "CAE_SIREN")
sample_audit <- list()
audit_sample <- function(data, dataset, stage) {
  data.frame(dataset = dataset, stage = stage, observations = nrow(data),
             observations_above_10 = sum(data$OFFERS > 10, na.rm = TRUE),
             cpv_groups = length(unique(data$CPV)), stringsAsFactors = FALSE)
}

read_sample <- function(filename, dataset, extra_required = character()) {
  required <- unique(c(required, extra_required))
  input <- file.path(repo_dir, filename)
  if (!file.exists(input)) stop("Input file not found: ", input)
  data <- read_csv(input, show_col_types = FALSE,
                   col_types = cols(.default = col_guess(), CPV = col_character(),
                                    CAE_SIREN = col_character()))
  if (nrow(problems(data))) stop("CSV parsing problems in ", filename)
  absent <- setdiff(required, names(data))
  if (length(absent)) stop("Missing columns in ", filename, ": ", paste(absent, collapse = ", "))
  audit <- list(audit_sample(data, dataset, "Imported"))
  # The CPV frequency threshold is applied to the full sample, before splitting
  # by sector or buyer, as in the supplied scripts and the existing replication.
  data <- data %>%
    mutate(CPV = trimws(gsub("'", "", as.character(CPV)))) %>%
    filter(substr(CPV, 3, 3) != "0") %>%
    mutate(CPV = substr(CPV, 1, 3))
  audit[[2]] <- audit_sample(data, dataset, "After CPV third-digit filter")
  data <- data %>% group_by(CPV) %>% filter(n() >= 30) %>% ungroup()
  audit[[3]] <- audit_sample(data, dataset, "After CPV minimum of 30 observations")
  data <- data %>% select(all_of(required)) %>%
    filter(if_all(everything(), ~ !is.na(.x))) %>% filter(OFFERS > 0)
  if (any(!is.finite(data$OFFERS)) || any(data$OFFERS != floor(data$OFFERS))) {
    stop("The outcome must contain finite positive integer counts: ", filename)
  }
  data <- data %>% mutate(CPV2 = substr(CPV, 1, 2), G_WEIGHT_SQUARED = G_WEIGHT^2,
                          CPV = factor(CPV), REGION = factor(REGION), STATUS = factor(STATUS))
  audit[[4]] <- audit_sample(data, dataset, "Final estimation sample")
  sample_audit[[dataset]] <<- bind_rows(audit)
  data
}

weight_sample <- function(data) {
  data %>%
    mutate(across(c(CPV, REGION, STATUS), droplevels)) %>%
    group_by(CAE_SIREN) %>% mutate(n_siren = n()) %>% ungroup() %>%
    mutate(w_inv = 1 / n_siren, w_samp = w_inv * (n() / sum(w_inv)))
}

all_estimates <- list()
all_diagnostics <- list()
models <- list()
model_logs <- list()

fit_model <- function(id, analysis, subgroup, data, formula, retry_optimizer = FALSE,
                      family = truncated_nbinom2(), weight_column = "w_samp") {
  diagnostic <- data.frame(
    model = id, analysis = analysis, subgroup = subgroup,
    observations = nrow(data), buyers = n_distinct(data$CAE_SIREN),
    cpv_groups = n_distinct(data$CPV), regions = n_distinct(data$REGION),
    formula = paste(deparse(formula), collapse = " "),
    family = family$family, weight_column = weight_column,
    convergence_code = NA_integer_, positive_definite_hessian = FALSE,
    valid_fit = FALSE, optimizer_retry_used = FALSE,
    logLik = NA_real_, AIC = NA_real_, error = NA_character_,
    stringsAsFactors = FALSE
  )
  tryCatch({
    if (!weight_column %in% names(data)) stop("Missing weight column: ", weight_column)
    data$.model_weight <- data[[weight_column]]
    if (any(!is.finite(data$.model_weight)) || any(data$.model_weight <= 0)) {
      stop("Model weights must be finite and positive.")
    }
    fit <- glmmTMB(formula, data = data, weights = .model_weight,
                   family = family, ziformula = ~ 0)
    # Retain the original buyer script's BFGS fallback for a non-positive Hessian.
    if (retry_optimizer && !isTRUE(fit$sdr$pdHess)) {
      retry <- tryCatch(update(fit, control = glmmTMBControl(
        optimizer = optim, optArgs = list(method = "BFGS"))), error = function(e) NULL)
      if (!is.null(retry) && isTRUE(retry$sdr$pdHess) && retry$fit$convergence == 0) {
        fit <- retry
        diagnostic$optimizer_retry_used <- TRUE
      }
    }
    models[[id]] <<- fit
    diagnostic$convergence_code <- fit$fit$convergence
    diagnostic$positive_definite_hessian <- isTRUE(fit$sdr$pdHess)
    diagnostic$valid_fit <- fit$fit$convergence == 0 && isTRUE(fit$sdr$pdHess)
    diagnostic$logLik <- as.numeric(logLik(fit))
    diagnostic$AIC <- AIC(fit)
    # Preserve the robust-SE call in the supplied replication scripts.
    estimates <- as.data.frame(parameters::model_parameters(fit, robust = TRUE, ci = 0.95))
    required_output <- c("Parameter", "Coefficient", "SE", "CI_low", "CI_high", "p")
    if (!all(required_output %in% names(estimates))) {
      stop("Unexpected model_parameters output; check your parameters package version.")
    }
    estimates$model <- id
    estimates$analysis <- analysis
    estimates$subgroup <- subgroup
    estimates$valid_fit <- diagnostic$valid_fit
    estimates$IRR <- exp(estimates$Coefficient)
    estimates$IRR_CI_low <- exp(estimates$CI_low)
    estimates$IRR_CI_high <- exp(estimates$CI_high)
    all_estimates[[id]] <<- estimates
    model_logs[[id]] <<- c(paste("MODEL:", id), capture.output(print(estimates)),
      "R-squared (performance package; interpretation depends on model type):",
      tryCatch(capture.output(print(performance::r2(fit))),
               error = function(e) paste("R-squared unavailable:", conditionMessage(e))))
  }, error = function(e) {
    diagnostic$error <<- conditionMessage(e)
    diagnostic$valid_fit <<- FALSE
    warning("Model ", id, " failed: ", conditionMessage(e), call. = FALSE)
  })
  all_diagnostics[[id]] <<- diagnostic
  invisible(NULL)
}

# Export partial results after each stage, so a later failure does not erase them.
stage_errors <- list()
checkpoint <- function() {
  write_csv(bind_rows(sample_audit), file.path(output_dir, "sample_audit.csv"))
  write_csv(bind_rows(all_diagnostics), file.path(output_dir, "model_diagnostics.csv"))
  if (length(all_estimates)) {
    coefficients <- bind_rows(all_estimates)
    write_csv(coefficients, file.path(output_dir, "all_coefficients.csv"))
    write_csv(coefficients %>% filter(Parameter %in%
      c("G_CLAUSE", "G_CRITERION", "G_WEIGHT", "G_WEIGHT_SQUARED")),
      file.path(output_dir, "green_coefficients.csv"))
  }
  saveRDS(models, file.path(output_dir, "models.rds"))
  writeLines(unlist(model_logs, use.names = FALSE), file.path(output_dir, "results_summary.txt"))
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
}
run_stage <- function(name, expression) {
  message("Running: ", name)
  tryCatch(force(expression), error = function(e) {
    stage_errors[[name]] <<- conditionMessage(e)
    warning("Stage ", name, " failed: ", conditionMessage(e), call. = FALSE)
  })
  checkpoint()
  invisible(NULL)
}

# 1. Input datasets. CPV filters are applied independently to each full dataset.
baseline <- weight_sample(read_sample("data_without_price.csv", "baseline",
                                     c("G_CRITERION", "ID_CONTRACT")))
price <- weight_sample(read_sample("data_with_price.csv", "price", "AWARD_PRICE"))
extended <- weight_sample(read_sample("data_without_price_extended.csv", "extended"))
if (any(baseline$OFFERS > 10)) stop("Baseline contains observations above ten offers.")
if (!any(extended$OFFERS > 10)) stop("Extended sample contains no observations above ten offers.")
if (any(!baseline$G_CLAUSE %in% c(0, 1)) || any(!baseline$G_CRITERION %in% c(0, 1))) {
  stop("G_CLAUSE and G_CRITERION must be binary.")
}
if (any(!is.finite(baseline$G_WEIGHT)) || any(baseline$G_WEIGHT < 0 | baseline$G_WEIGHT > 100)) {
  stop("G_WEIGHT must be finite and between 0 and 100.")
}
checkpoint()

linear_fixed <- OFFERS ~ G_CLAUSE + G_WEIGHT + ALLOTMENT + FRAMEWORK_AGREEMENT +
  P_CRITERION_WEIGHT_ENV + REGION + STATUS + factor(CPV2)
linear_random <- OFFERS ~ G_CLAUSE + G_WEIGHT + ALLOTMENT + FRAMEWORK_AGREEMENT +
  P_CRITERION_WEIGHT_ENV + REGION + STATUS + (1 | CPV)
quadratic_random <- OFFERS ~ G_CLAUSE + G_WEIGHT + G_WEIGHT_SQUARED + ALLOTMENT +
  FRAMEWORK_AGREEMENT + P_CRITERION_WEIGHT_ENV + REGION + STATUS + (1 | CPV)

# 2. Main estimates (Table 2).
run_stage("main", {
  fit_model("main_fixed_linear", "main", "baseline", baseline, linear_fixed)
  fit_model("main_random_linear", "main", "baseline", baseline, linear_random)
  fit_model("main_random_quadratic", "main", "baseline", baseline, quadratic_random)
})

# 3. Sector heterogeneity (Table 3): three-digit CPV FIXED effects.
run_stage("sector_heterogeneity", {
  for (sector in c("45", "90", "33")) {
    sample <- weight_sample(baseline %>% filter(CPV2 == sector))
    fit_model(paste0("sector_", sector), "sector_heterogeneity", sector, sample,
      OFFERS ~ G_CLAUSE + G_WEIGHT + ALLOTMENT + FRAMEWORK_AGREEMENT +
        P_CRITERION_WEIGHT_ENV + REGION + STATUS + factor(CPV))
  }
})

# 4. Buyer heterogeneity (Appendix C): three-digit CPV RANDOM intercepts.
run_stage("buyer_heterogeneity", {
  for (buyer in levels(baseline$STATUS)) {
    sample <- weight_sample(baseline %>% filter(STATUS == buyer))
    id <- paste0("buyer_", gsub("[^a-z0-9]+", "_", tolower(buyer)))
    fit_model(id, "buyer_heterogeneity", buyer, sample,
      OFFERS ~ G_CLAUSE + G_WEIGHT + ALLOTMENT + FRAMEWORK_AGREEMENT +
        P_CRITERION_WEIGHT_ENV + REGION + (1 | CPV), retry_optimizer = TRUE)
  }
  buyer_estimates <- bind_rows(all_estimates) %>% filter(analysis == "buyer_heterogeneity")
  cochran_q <- function(term) {
    x <- buyer_estimates %>% filter(Parameter == term, valid_fit,
      is.finite(Coefficient), is.finite(SE), SE > 0)
    k <- nrow(x)
    if (k < 2) return(data.frame(Parameter = term, groups = k,
      pooled_coefficient = NA_real_, Q = NA_real_, df = NA_integer_, p_value = NA_real_))
    weights <- 1 / x$SE^2
    pooled <- sum(weights * x$Coefficient) / sum(weights)
    q <- sum(weights * (x$Coefficient - pooled)^2)
    data.frame(Parameter = term, groups = k, pooled_coefficient = pooled,
      Q = q, df = k - 1, p_value = pchisq(q, k - 1, lower.tail = FALSE))
  }
  # Separate-fit inverse-variance tests, assuming independent subgroup estimates.
  q_tests <- bind_rows(lapply(c("G_CLAUSE", "G_WEIGHT"), cochran_q))
  write_csv(q_tests, file.path(output_dir, "buyer_heterogeneity_Q_tests.csv"))
})

# 5. Treatment reweighting (Appendices E-F).
# Preserve the original propensity models, stabilization and clipping choices.
fit_propensity <- function(id, data, treatment, formula, maxit = 100) {
  fit <- nnet::multinom(formula, data = data, maxit = maxit, trace = FALSE)
  models[[id]] <<- fit
  if (fit$convergence != 0) stop("Nonconverged propensity model: ", id)
  probabilities <- predict(fit, type = "probs")
  if (is.null(dim(probabilities))) stop("At least three treatment groups are required: ", id)
  levels_treatment <- levels(data[[treatment]])
  if (!all(levels_treatment %in% colnames(probabilities))) {
    stop("Missing treatment groups in propensity probabilities: ", id)
  }
  if (any(!is.finite(probabilities)) || any(probabilities < 0 | probabilities > 1)) {
    stop("Invalid propensity probabilities: ", id)
  }
  # Model-based multinomial standard errors, as in the original multinom output.
  summary_fit <- summary(fit)
  coefficients <- summary_fit$coefficients
  errors <- summary_fit$standard.errors
  coefficient_table <- bind_rows(lapply(seq_len(nrow(coefficients)), function(i) {
    z <- coefficients[i, ] / errors[i, ]
    data.frame(model = id, treatment_group = rownames(coefficients)[i],
      Parameter = colnames(coefficients), Coefficient = as.numeric(coefficients[i, ]),
      SE = as.numeric(errors[i, ]), z = as.numeric(z),
      p = 2 * pnorm(abs(z), lower.tail = FALSE), row.names = NULL)
  }))
  write_csv(coefficient_table, file.path(output_dir, paste0(id, "_coefficients.csv")))
  chosen <- probabilities[cbind(seq_len(nrow(data)),
                                match(as.character(data[[treatment]]), colnames(probabilities)))]
  shares <- prop.table(table(data[[treatment]]))
  numerator <- as.numeric(shares[as.character(data[[treatment]])])
  list(chosen = chosen, numerator = numerator, probabilities = probabilities)
}
export_propensity <- function(id, data, treatment, propensity, ipw) {
  write_csv(data.frame(observation = seq_len(nrow(data)),
    ID_CONTRACT = data$ID_CONTRACT, CAE_SIREN = data$CAE_SIREN,
    treatment = as.character(data[[treatment]]),
    propensity = propensity$chosen, treatment_share = propensity$numerator,
    stabilized_ipw = ipw, sampling_weight = data$w_samp,
    outcome_weight = data$w_outcome),
    file.path(output_dir, paste0(id, "_weights.csv")))
}

run_stage("binary_reweighting", {
  set.seed(123)
  binary <- baseline %>% mutate(
    treatment = factor(2 * G_CRITERION + G_CLAUSE, levels = 0:3,
      labels = c("none", "clause only", "criterion only", "both")))
  propensity <- fit_propensity("propensity_binary", binary, "treatment",
    treatment ~ ALLOTMENT + FRAMEWORK_AGREEMENT + P_CRITERION_WEIGHT_ENV +
      factor(CPV) + factor(REGION) + factor(STATUS))
  if (any(propensity$chosen <= 0)) stop("Zero estimated propensity in binary treatment model.")
  ipw <- propensity$numerator / propensity$chosen
  combined <- ipw * binary$w_samp
  binary$w_outcome <- combined * (nrow(binary) / sum(combined))
  export_propensity("propensity_binary", binary, "treatment", propensity, ipw)
  # Deliberately NONTRUNCATED: this matches replication.R's binary specification.
  fit_model("reweighted_binary", "treatment_reweighting", "binary", binary,
    OFFERS ~ G_CLAUSE + G_CRITERION + ALLOTMENT + FRAMEWORK_AGREEMENT +
      P_CRITERION_WEIGHT_ENV + factor(REGION) + factor(STATUS) + (1 | CPV),
    family = nbinom2(), weight_column = "w_outcome")
})

run_stage("weight_category_reweighting", {
  continuous <- baseline %>% mutate(
    weight_bin = case_when(G_WEIGHT == 0 ~ 0L, G_WEIGHT <= 10 ~ 1L,
                           G_WEIGHT <= 30 ~ 2L, TRUE ~ 3L),
    treatment = factor(2 * weight_bin + G_CLAUSE, levels = 0:7,
      labels = c("none", "clause_only", "weight_low", "clause_weight_low",
                 "weight_medium", "clause_weight_medium", "weight_high", "clause_weight_high")))
  propensity <- fit_propensity("propensity_weight_categories", continuous, "treatment",
    treatment ~ ALLOTMENT + FRAMEWORK_AGREEMENT + P_CRITERION_WEIGHT_ENV +
      factor(CPV2) + factor(REGION) + factor(STATUS), maxit = 1000)
  ipw <- propensity$numerator / pmax(propensity$chosen, 0.001)
  ipw <- pmin(pmax(ipw, 0.1), 10)
  combined <- ipw * continuous$w_samp
  continuous$w_outcome <- combined * (nrow(continuous) / sum(combined))
  export_propensity("propensity_weight_categories", continuous, "treatment", propensity, ipw)
  fit_model("reweighted_weight_linear", "treatment_reweighting", "weight_categories", continuous,
    OFFERS ~ G_CLAUSE + G_WEIGHT + ALLOTMENT + FRAMEWORK_AGREEMENT +
      P_CRITERION_WEIGHT_ENV + factor(REGION) + factor(STATUS) + (1 | CPV),
    weight_column = "w_outcome")
  fit_model("reweighted_weight_quadratic", "treatment_reweighting", "weight_categories", continuous,
    OFFERS ~ G_CLAUSE + G_WEIGHT + G_WEIGHT_SQUARED + ALLOTMENT + FRAMEWORK_AGREEMENT +
      P_CRITERION_WEIGHT_ENV + factor(REGION) + factor(STATUS) + (1 | CPV),
    weight_column = "w_outcome")
})

# 6. Award-price control (Appendix G). The CSV already contains log(price).
run_stage("price_control", {
  fit_model("price_fixed_linear", "price_control", "price_sample", price,
    OFFERS ~ G_CLAUSE + G_WEIGHT + AWARD_PRICE + ALLOTMENT + FRAMEWORK_AGREEMENT +
      P_CRITERION_WEIGHT_ENV + REGION + STATUS + factor(CPV2))
  fit_model("price_random_linear", "price_control", "price_sample", price,
    OFFERS ~ G_CLAUSE + G_WEIGHT + AWARD_PRICE + ALLOTMENT + FRAMEWORK_AGREEMENT +
      P_CRITERION_WEIGHT_ENV + REGION + STATUS + (1 | CPV))
})

# 7. Upper bid-count restriction (Appendix H): ALL observations in extended data.
# Input is the prepared sample. Do not reconstruct upstream contract exclusions
# using an export that may contain only part of each contract's original lots.
run_stage("upper_bid_count_restriction", {
  fit_model("extended_fixed_linear", "upper_bid_count_restriction", "extended", extended, linear_fixed)
  fit_model("extended_random_linear", "upper_bid_count_restriction", "extended", extended, linear_random)
  fit_model("extended_random_quadratic", "upper_bid_count_restriction", "extended", extended, quadratic_random)
})

# 8. Repeated one-lot-per-contract draws (Appendix D). Run last because expensive.
# Original sample weights are RETAINED after drawing, not recalculated.
positive_integer_setting <- function(name, default) {
  value <- suppressWarnings(as.numeric(Sys.getenv(name, unset = as.character(default))))
  if (length(value) != 1 || !is.finite(value) || value < 1 || value != floor(value)) {
    stop(name, " must be a positive integer.")
  }
  as.integer(value)
}
run_bootstrap <- function(data) {
  n_replications <- positive_integer_setting("GPP_BOOTSTRAP_REPS", 1000)
  n_workers <- positive_integer_setting("GPP_WORKERS", min(4L, future::availableCores()))
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  if (n_workers == 1L) future::plan(future::sequential)
  else future::plan(future::multisession, workers = n_workers)
  focal_terms <- c("G_CLAUSE", "G_WEIGHT", "G_WEIGHT_SQUARED")
  estimate_draw <- function(iteration) {
    tryCatch({
      set.seed(iteration)
      sampled <- data %>% group_by(ID_CONTRACT) %>% slice_sample(n = 1) %>% ungroup()
      fit <- glmmTMB(OFFERS ~ G_CLAUSE + G_WEIGHT + G_WEIGHT_SQUARED + ALLOTMENT +
        FRAMEWORK_AGREEMENT + P_CRITERION_WEIGHT_ENV + REGION + STATUS + (1 | CPV),
        data = sampled, weights = w_samp, family = truncated_nbinom2(), ziformula = ~ 0)
      valid <- fit$fit$convergence == 0 && isTRUE(fit$sdr$pdHess)
      estimates <- as.data.frame(parameters::model_parameters(fit, robust = TRUE))
      estimates <- estimates %>% filter(Parameter %in% focal_terms) %>%
        select(Parameter, Coefficient, SE, p)
      if (nrow(estimates) != length(focal_terms)) stop("Missing focal coefficients.")
      valid <- valid && all(is.finite(estimates$Coefficient)) &&
        all(is.finite(estimates$SE)) && all(is.finite(estimates$p))
      estimates %>% mutate(iteration = iteration, n_observations = nrow(sampled),
        valid_fit = valid, convergence_code = fit$fit$convergence,
        positive_definite_hessian = isTRUE(fit$sdr$pdHess), error = NA_character_)
    }, error = function(e) {
      data.frame(Parameter = focal_terms, Coefficient = NA_real_, SE = NA_real_, p = NA_real_,
        iteration = iteration, n_observations = NA_integer_, valid_fit = FALSE,
        convergence_code = NA_integer_, positive_definite_hessian = FALSE,
        error = conditionMessage(e), stringsAsFactors = FALSE)
    })
  }
  set.seed(123)
  results <- furrr::future_map_dfr(seq_len(n_replications), estimate_draw,
    .options = furrr::furrr_options(seed = TRUE, packages = c("dplyr", "glmmTMB", "parameters")),
    .progress = TRUE)
  write_csv(results, file.path(output_dir, "bootstrap_coefficients.csv"))
  clean <- results %>% filter(valid_fit, is.finite(Coefficient), is.finite(SE), is.finite(p))
  if (!nrow(clean)) stop("No valid subsampling estimates; inspect bootstrap_coefficients.csv.")
  summary_results <- clean %>% group_by(Parameter) %>% summarise(
    mean = mean(Coefficient), sd = sd(Coefficient),
    ci95_lo = quantile(Coefficient, 0.025), ci95_hi = quantile(Coefficient, 0.975),
    ci99_lo = quantile(Coefficient, 0.005), ci99_hi = quantile(Coefficient, 0.995),
    sig_rate = mean(p < 0.05) * 100, valid_replications = n(),
    requested_replications = n_replications, .groups = "drop")
  write_csv(summary_results, file.path(output_dir, "bootstrap_summary.csv"))
  # Keep the original plot's solid MEAN line and dashed empirical 95% interval.
  # The manuscript's caption says median; that caption is not implemented here.
  labels <- c(G_CLAUSE = "G_clause", G_WEIGHT = "G_weight", G_WEIGHT_SQUARED = "G_weight squared")
  clean$Parameter <- factor(clean$Parameter, levels = focal_terms)
  summary_results$Parameter <- factor(summary_results$Parameter, levels = focal_terms)
  limits <- clean %>% group_by(Parameter) %>% summarise(
    xmin = min(min(Coefficient), 0), xmax = max(max(Coefficient), 0), .groups = "drop")
  plot <- ggplot2::ggplot(clean, ggplot2::aes(x = Coefficient)) +
    ggplot2::geom_histogram(bins = 50, fill = "grey60", color = "white", alpha = 0.8) +
    ggplot2::geom_vline(data = summary_results, ggplot2::aes(xintercept = mean), linewidth = 0.8) +
    ggplot2::geom_vline(data = summary_results, ggplot2::aes(xintercept = ci95_lo), linetype = "dashed") +
    ggplot2::geom_vline(data = summary_results, ggplot2::aes(xintercept = ci95_hi), linetype = "dashed") +
    ggplot2::geom_blank(data = limits, ggplot2::aes(x = xmin)) +
    ggplot2::geom_blank(data = limits, ggplot2::aes(x = xmax)) +
    ggplot2::facet_wrap(~Parameter, scales = "free", labeller = ggplot2::as_labeller(labels)) +
    ggplot2::labs(x = "Estimated coefficient", y = "Frequency") +
    ggplot2::theme_minimal() + ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
  ggplot2::ggsave(file.path(output_dir, "bootstrap_distribution.png"), plot,
                  width = 10, height = 4, dpi = 300)
  ggplot2::ggsave(file.path(output_dir, "bootstrap_distribution.pdf"), plot,
                  width = 10, height = 4)
  print(summary_results)
  list(requested = n_replications,
       invalid = n_distinct(results$iteration[!results$valid_fit]), workers = n_workers)
}

bootstrap_flag <- tolower(Sys.getenv("GPP_RUN_BOOTSTRAP", unset = "true"))
if (!bootstrap_flag %in% c("true", "false")) stop("GPP_RUN_BOOTSTRAP must be true or false.")
bootstrap_status <- list(requested = 0L, invalid = 0L, workers = 0L)
if (bootstrap_flag == "true") {
  run_stage("within_contract_subsampling", { bootstrap_status <- run_bootstrap(baseline) })
}
checkpoint()
diagnostics <- bind_rows(all_diagnostics)
bad_models <- diagnostics %>% filter(!valid_fit | !is.na(error))
errors <- if (length(stage_errors)) data.frame(stage = names(stage_errors),
  error = unlist(stage_errors, use.names = FALSE)) else data.frame(stage = character(), error = character())
write_csv(errors, file.path(output_dir, "stage_errors.csv"))
status <- data.frame(expected_count_models = 21L, fitted_count_models = nrow(diagnostics),
  invalid_count_models = nrow(bad_models), failed_stages = nrow(errors),
  bootstrap_enabled = bootstrap_flag == "true", bootstrap_requested = bootstrap_status$requested,
  bootstrap_invalid = bootstrap_status$invalid, bootstrap_workers = bootstrap_status$workers,
  complete_design_run = bootstrap_flag == "true" && bootstrap_status$requested == 1000L &&
    nrow(bad_models) == 0L && nrow(errors) == 0L && nrow(diagnostics) == 21L)
write_csv(status, file.path(output_dir, "run_status.csv"))
print(status)
if (nrow(bad_models) || nrow(errors) || nrow(diagnostics) != 21L) {
  stop("Replication contains failed stages or invalid models. Outputs were saved; inspect diagnostics and stage_errors.csv.")
}
if (bootstrap_status$invalid > 0) {
  warning("Some subsampling fits were invalid and excluded; check bootstrap_coefficients.csv.", call. = FALSE)
}
message("Finished. Results saved in: ", output_dir)
