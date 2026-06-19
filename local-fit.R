# Local single-rep fit for the current LMM, BP, and GP workflow
library(here)
library(tidyverse)

# Initialise seed for reproducibility
global_seed <- 854098
source(here("Simulation-Workflow", "functions.R")) # Sources 'functions' folder

# Point to Stan files
stan_path_bp1 <- here("Stan", "Bernstein-Polynomials-JM-Hist.stan")
stan_path_bp2 <- here("Stan", "bernstein-polynomials.stan")
stan_path_gp <- here("Stan", "gaussian-process.stan")

# Use the 'fit_one_rep' function from 'functions.R'. This function represents the fitting of the models for a given simulated dataset, parameters are specified.
local_result <- fit_one_rep(
  rep_local_idx = 1, batch_id = 1, n_per_batch = 1, base_seed = global_seed,
  stan_path_bp1 = stan_path_bp1, stan_path_bp2 = stan_path_bp2, stan_path_gp = stan_path_gp,
  D = matrix(c(15^2, -0.10 * 15 * 0.20, -0.10 * 15 * 0.20, 0.20^2), 2, 2),
  beta_0 = 73, beta_1 = -0.04, beta_2 = 0.04, sigma_e = 12,
  log_AF = -0.90, alpha_AFT = 0.012,
  loglogistic_shape = 1.20, loglogistic_scale = 23,
  visit = c(0, 1, seq(3, 92, 3)), n_patients = 100, max_FU = 120,
  lambda_c = -1, aft_mode = "loglogistic",
  models_to_fit = c("LMM", "bp2", "gp"),
  n_stan_chains = 4, n_stan_warmup = 1000, n_stan_iter = 1000
)

# Summary formatting function for clean summary description
format_summary <- function(summary_tbl, y_mean = 0, model = c("bp", "gp")) {
  model <- match.arg(model)
  if (is.null(summary_tbl)) return(NULL)
  
  surv_var <- if (model == "gp") "gamma[1]" else "beta_surv[1]"
  
  summary_tbl %>%
    filter(variable %in% c(
      "beta_long[1]", "beta_long[2]", "beta_long[3]",
      surv_var, "alpha"
    )) %>%
    mutate(
      parameter = case_when(
        variable == "beta_long[1]" ~ "beta_0",
        variable == "beta_long[2]" ~ "beta_1",
        variable == "beta_long[3]" ~ "beta_2",
        variable == surv_var       ~ "gamma",
        variable == "alpha"        ~ "alpha",
        TRUE                       ~ variable
      ),
      mean = if_else(variable == "beta_long[1]", mean + y_mean, mean),
      q2.5 = if_else(variable == "beta_long[1]", q2.5 + y_mean, q2.5),
      q97.5 = if_else(variable == "beta_long[1]", q97.5 + y_mean, q97.5)
    ) %>%
    dplyr::select(parameter, mean, sd, q2.5, q97.5)
}

# Summary from the linear mixed model
lmm_summary <- local_result$LMM$summary %>%
  mutate(parameter = recode(variable,
                            `beta_long[1]` = "beta_0",
                            `beta_long[2]` = "beta_1",
                            `beta_long[3]` = "beta_2"),
         mean = if_else(variable == "beta_long[1]", mean + local_result$Y_long_mean, mean), # Adding back the longitudinal mean from centring
         q2.5 = if_else(variable == "beta_long[1]", q2.5 + local_result$Y_long_mean, q2.5),
         q97.5 = if_else(variable == "beta_long[1]", q97.5 + local_result$Y_long_mean, q97.5),
         sd = if ("sd" %in% names(.)) sd else (q97.5 - q2.5) / (2 * 1.96)) %>%
  dplyr::select(parameter, mean, sd, q2.5, q97.5)

# Direct comparison of results of models
result_comparison <- bind_rows(
  LMM = lmm_summary,
  BP = format_summary(local_result$bp2$summary, local_result$Y_long_mean, model = "bp"),
  GP = format_summary(local_result$gp$summary, local_result$Y_long_mean, model = "gp"),
  .id = "model"
)

result_comparison
