
# Pre-processing
##########

# Load relevant libraries
library(tidyverse) # for data manipulation and visualization
library(future) # for parallel execution
library(furrr) # for parallel map functions
library(here) # file path management
library(cmdstanr) # Stan interface

# Stan model paths (passed as strings; each worker compiles from disk)
stan_path_bp1 <- here("WPP","Stan", "JM", "Bernstein-Polynomials-JM-Hist.stan")
stan_path_bp2 <- here("WPP","Stan", "JM", "bernstein-polynomials.stan")
stan_path_gp <- here("WPP","Stan", "JM", "gaussian-process.stan")

##########



# Command-line arguments
##########

args <- commandArgs(trailingOnly = TRUE)

batch_id <- if (length(args) >= 1) as.integer(args[[1]]) else 1
n_per_batch <- if (length(args) >= 2) as.integer(args[[2]]) else NULL  # resolved after n_sim is sourced

# 4 parallel Stan chains per replicate
n_stan_chains <- 4

# CPUs available from SLURM
slurm_cpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
if (is.na(slurm_cpus)) {
  slurm_cpus <- parallelly::availableCores()
}

# Optional third command-line arg = requested max future workers
requested_workers <- if (length(args) >= 3) as.integer(args[[3]]) else Inf

# Model arguments begin at position four. Comma-separated and separate model
# arguments are both accepted; when a scenario is supplied it is always last.
# Examples: "LMM,jm2" s1_ll_n300, or LMM jm2 s1_ll_n300.
has_scenario_arg <- length(args) >= 4 &&
  grepl("^s[0-9]+_(ll|wb)(?:_n[0-9]+)?$", args[[length(args)]], perl = TRUE)

models_to_fit <- if (length(args) >= 4) {
  model_end <- if (has_scenario_arg) length(args) - 1 else length(args)
  model_args <- args[4:model_end]
  trimws(unlist(strsplit(model_args, ",", fixed = TRUE)))
} else {
  c("LMM", "bp2", "gp")
}

# Optional final command-line arg = scenario ID
# Example IDs: s1_ll_n300, s2_wb_n300. Legacy IDs such as s1_ll remain valid.
scenario_id <- if (has_scenario_arg) trimws(args[[length(args)]]) else "s1_ll_n300"
Sys.setenv(SCENARIO_ID = scenario_id)

# Source scenario parameters and simulation functions
# (datagen.R reads SCENARIO_ID from the environment and sets sc, beta_2, log_AF, etc.)
source(here("WPP", "Simulation-Workflow", "datagen.R"))

# Minimal simulation settings
n_sim <- 100  # replications (increase for production)
if (is.null(n_per_batch)) n_per_batch <- n_sim

base_seed <- 32134

# Output directory (scenario-specific subdirectory)
results_dir <- here("WPP","Simulation-Workflow", "Results", scenario_id)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

##########



# Worker count
##########

# Each future worker will itself run n_stan_chains parallel Stan chains.
# So n_workers must be capped by available CPUs / n_stan_chains.
max_workers_by_cpu <- max(1L, floor(slurm_cpus / n_stan_chains))

# JMbayes2 keeps considerably more per-chain state in memory than the CmdStan
# fits. Run one jm2 replicate at a time while retaining four parallel chains.
max_workers_by_model <- if ("jm2" %in% models_to_fit) 1L else Inf

n_workers <- min(
  n_per_batch,
  requested_workers,
  max_workers_by_cpu,
  max_workers_by_model
)

message(sprintf(
  "[Batch %03d | %s]  n_per_batch = %d  |  n_workers = %d  |  n_stan_chains = %d  |  cpus = %d",
  batch_id, scenario_id, n_per_batch, n_workers, n_stan_chains, slurm_cpus
))

##########



# Parallelised batch execution
##########

if (n_workers == 1L) {
  # Avoid an unnecessary outer R process around the four jm2 chain processes.
  future::plan(future::sequential)
} else {
  future::plan(future::multisession, workers = n_workers)
}

fit_one_rep_with_gc <- function(rep_local_idx, ...) {
  result <- fit_one_rep(rep_local_idx, ...)
  invisible(gc())
  result
}

results <- furrr::future_map(
  seq_len(n_per_batch),
  fit_one_rep_with_gc,
  batch_id = batch_id,
  n_per_batch = n_per_batch,
  base_seed = base_seed,
  stan_path_bp1 = stan_path_bp1,
  stan_path_bp2 = stan_path_bp2,
  stan_path_gp = stan_path_gp,
  D = D,
  beta_0 = beta_0,
  beta_1 = beta_1,
  beta_2 = beta_2,
  sigma_e = sigma_e,
  log_AF = log_AF,
  alpha_AFT = alpha_AFT,
  loglogistic_shape = loglogistic_shape,
  loglogistic_scale = loglogistic_scale,
  aft_mode = aft_mode,
  visit = visit,
  n_patients = n_patients,
  max_FU = max_FU,
  
  # important
  n_stan_chains = n_stan_chains,
  models_to_fit = models_to_fit,
  
  .options = furrr::furrr_options(seed = NULL)
)

future::plan(future::sequential)

n_ok <- sum(sapply(results, function(r) is.null(r$error)))
n_err <- n_per_batch - n_ok

message(sprintf(
  "[Batch %03d]  %d / %d succeeded%s",
  batch_id, n_ok, n_per_batch,
  if (n_err > 0) sprintf(" (%d failed)", n_err) else ""
))


##########



# Save batch output
##########

out_file <- file.path(results_dir, sprintf("WPP-%s-%s-%03d.rds", scenario_id, aft_mode, batch_id))
saveRDS(
  list(
    batch_id = batch_id,
    scenario_id = scenario_id,
    scenario = sc,
    n_sim = n_per_batch,
    true_params = true_params,
    results = results
  ),
  file = out_file
)
message(sprintf("[Batch %03d | %s]  Saved → %s", batch_id, scenario_id, out_file))

##########



