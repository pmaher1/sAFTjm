
# Pre-processing
##########

# Load relevant libraries
library(tidyverse) # for data manipulation and visualization
library(cmdstanr) # for Bayesian modeling
library(survival) # for survival analysis functions
library(nlme) # for linear mixed-effects models
library(MASS) # for mvrnorm function to simulate random effects
library(here) # for constructing file paths relative to project root

# Set seed for reproducibility
global_seed <- 412321
set.seed(global_seed)

# Source data generation function
# source(here("Programs", "Data-Gen-ENZAMET-01.R"))

source(here("Simulation-Workflow", "functions.R"))

##########



# Simulation scenarios
##########

### Study design
sample_size_scenarios <- c(n300 = 300, n500 = 500)
max_FU <- 120
visit <- c(0, 1, seq(3, 92, 3))

### Fixed longitudinal parameters (shared across scenarios)
beta_0 <- 73     # baseline HRQoL
beta_1 <- -0.04  # monthly rate of change (control arm)
sigma_e <- 12     # residual SD

### Random effects covariance (fixed)
sigma_0 <- 15
sigma_1 <- 0.20
D <- matrix(
  c(sigma_0^2, -0.10 * sigma_0 * sigma_1,
    -0.10 * sigma_0 * sigma_1, sigma_1^2),
  2, 2
)

alpha_AFT <- 0.012  # no longitudinal-TTE association (fixed across scenarios)

### Baseline hazard parameters
# Log-logistic
ll_shape <- 1.20
ll_scale <- 23.0

# Weibull (matched to give similar median survival)
wb_shape <- 0.90
wb_scale <- 38.0

### Scenario grid
# Treatment effect patterns x baseline hazard distributions x sample sizes.
# beta_2: treatment-by-time interaction on HRQoL
# log_AF: log acceleration factor for TTE (gamma in the manuscript table)
treatment_scenarios <- list(
  s1 = list(label = "S1", beta_2 =  0.00, log_AF =  0.00),
  s2 = list(label = "S2", beta_2 =  0.04, log_AF =  0.90),
  s3 = list(label = "S3", beta_2 = -0.04, log_AF =  0.90),
  s4 = list(label = "S4", beta_2 =  0.04, log_AF = -0.90),
  s5 = list(label = "S5", beta_2 = -0.04, log_AF = -0.90)
)

baseline_scenarios <- list(
  ll = list(label = "LL", description = "Log-Logistic",
            aft_mode = "loglogistic", baseline_shape = ll_shape,
            baseline_scale = ll_scale),
  wb = list(label = "WB", description = "Weibull",
            aft_mode = "weibull", baseline_shape = wb_shape,
            baseline_scale = wb_scale)
)

scenario_grid <- expand.grid(
  treatment_id = names(treatment_scenarios),
  baseline_id = names(baseline_scenarios),
  sample_size_id = names(sample_size_scenarios),
  stringsAsFactors = FALSE
)

make_scenario <- function(treatment_id, baseline_id, sample_size_id) {
  trt <- treatment_scenarios[[treatment_id]]
  base <- baseline_scenarios[[baseline_id]]
  n <- unname(sample_size_scenarios[[sample_size_id]])

  list(
    label = sprintf("%s (%s, n=%s)", trt$label, base$label, n),
    treatment_scenario = treatment_id,
    baseline_hazard = baseline_id,
    sample_size_id = sample_size_id,
    n_patients = n,
    beta_2 = trt$beta_2,
    log_AF = trt$log_AF,
    aft_mode = base$aft_mode,
    baseline_shape = base$baseline_shape,
    baseline_scale = base$baseline_scale
  )
}

scenario_ids <- paste(
  scenario_grid$treatment_id,
  scenario_grid$baseline_id,
  scenario_grid$sample_size_id,
  sep = "_"
)

scenarios <- setNames(
  lapply(seq_len(nrow(scenario_grid)), function(i) {
    make_scenario(scenario_grid$treatment_id[[i]],
                  scenario_grid$baseline_id[[i]],
                  scenario_grid$sample_size_id[[i]])
  }),
  scenario_ids
)

# Backward-compatible aliases for existing result folders and old CLI calls.
default_sample_size_id <- names(sample_size_scenarios)[[1]]
legacy_grid <- scenario_grid[scenario_grid$sample_size_id == default_sample_size_id, ]
legacy_ids <- paste(legacy_grid$treatment_id, legacy_grid$baseline_id, sep = "_")
legacy_scenarios <- setNames(
  lapply(seq_len(nrow(legacy_grid)), function(i) {
    make_scenario(legacy_grid$treatment_id[[i]],
                  legacy_grid$baseline_id[[i]],
                  legacy_grid$sample_size_id[[i]])
  }),
  legacy_ids
)
scenarios <- c(scenarios, legacy_scenarios)

### Select active scenario (default: scenario 5, log-logistic, n = 300)
# Overridden by hpc-simulate.R via the SCENARIO_ID environment variable / CLI arg
active_scenario_id <- Sys.getenv("SCENARIO_ID", unset = "s5_ll_n300")
if (!active_scenario_id %in% names(scenarios)) {
  stop(sprintf("Unknown SCENARIO_ID '%s'. Valid IDs: %s",
               active_scenario_id, paste(names(scenarios), collapse = ", ")))
}
sc <- scenarios[[active_scenario_id]]

# Expose scenario parameters as top-level variables (used by local-fit.R etc.)
n_patients <- sc$n_patients
beta_2 <- sc$beta_2
log_AF <- sc$log_AF
aft_mode <- sc$aft_mode
loglogistic_shape <- sc$baseline_shape   # name kept for backward-compat
loglogistic_scale <- sc$baseline_scale

### True parameter values (for performance metrics in result-formatting.R)
true_params <- c(
  n_patients = n_patients,
  beta_0 = beta_0,
  beta_1 = beta_1,
  beta_2 = beta_2,
  log_AF = log_AF,
  alpha_AFT = alpha_AFT
)

##########



