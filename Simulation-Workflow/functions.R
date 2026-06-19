
# Pre-processing
##########

# Load relevant libraries
library(tidyverse) # for data manipulation and visualization
library(cmdstanr) # for Bayesian modeling

# Set seed for reproducibility
# global_seed <- 4217312
set.seed(global_seed)

##########



# Data generation functions
##########

# Inverse CDF functions for survival time generation
ll_inv <- function(u, shape, scale) scale * ((u / (1 - u))^(1 / shape))
wb_inv <- function(u, shape, scale) scale * (-log(1 - u))^(1 / shape)

# Primary simulation function
simulate_joint_dataset <- function(D = matrix(c(15^2, -0.10*15*0.20, -0.10*15*0.20, 0.20^2), 2, 2), 
                                   beta_0 = 73, beta_1 = -0.04, beta_2 = 0.00, sigma_e = 12,
                                   log_AF = 0.00, alpha_AFT = 0.012,
                                   weibull_shape = 0.90, weibull_scale = 38,
                                   loglogistic_shape = 1.20, loglogistic_scale = 23,
                                   visit = c(0, 1, seq(3, 92, 3)), 
                                   seed, n_patients = 1100, max_FU = 120, lambda_c = -1,
                                   aft_mode = "loglogistic", link_type = "value",
                                   ...) {  # absorb unused PH/Weibull placeholders
  
  # Set seed for reproducibility
  set.seed(seed)
  
  # Initialise parameters
  n <- n_patients
  arm <- rep(0:1, each = n / 2)
  b <- mvrnorm(n, mu = c(0, 0), Sigma = D)
  b0 <- b[, 1]; b1 <- b[, 2]
  
  # Acceleration components (vectorized over participants)
  # ENZAMET-style current-value link uses the raw baseline value.
  Y_ref <- beta_0
  eta0 <- if (link_type == "value") beta_0 + b0 else rep(0, n)
  C1 <- log_AF * arm + alpha_AFT * eta0
  C2 <- if (link_type %in% c("value", "slope")) {
    alpha_AFT * (beta_1 + beta_2 * arm + b1)
  } else {
    rep(0, n)
  }
  
  # Inverse-sample base survival times
  U <- runif(n, 1e-6, 1 - 1e-6)
  kappa <- if (aft_mode == "loglogistic") {
    ll_inv(U, loglogistic_shape, loglogistic_scale)
  } else {
    wb_inv(U, weibull_shape, weibull_scale)  # reuse shape/scale args
  }
  
  # Survival times: T = -log(1 - C2*exp(C1)*kappa) / C2 when C2 != 0
  if (aft_mode == "weibull") {
    power_term <- (-log(1 - U))^(1 / weibull_shape)
    A <- pmin(C2 * weibull_scale * exp(C1) * power_term, 1 - 1e-8)
    T_i <- ifelse(abs(C2) < 1e-8,
                  weibull_scale * exp(C1) * power_term,
                  log(1 - A) / (-C2))
  } else {
    A <- pmin(C2 * exp(C1) * kappa, 1 - 1e-8)
    T_i <- ifelse(abs(C2) < 1e-8, exp(C1) * kappa, -log(1 - A) / C2)
  }
  
  # Censoring
  if (lambda_c <= 0) {
    T_obs <- pmin(T_i, max_FU)
    status <- as.integer(T_i <= max_FU)
  } else if (lambda_c > 0) {
    C_i <- rexp(n, rate = lambda_c)
    T_obs <- pmin(T_i, C_i, max_FU)
    status <- as.integer(T_i <= C_i & T_i <= max_FU)
  } 
  
  # Longitudinal data: expand each participant to their valid visit times
  long_data <- do.call(rbind, lapply(seq_len(n), function(i) {
    jit <- rnorm(length(visit), 0, 1)
    jit[visit == 0] <- 0
    vt <- pmax(visit + jit, 0)
    vt <- vt[vt < T_obs[i]]
    if (length(vt) == 0) return(NULL)
    y_true <- beta_0 + (beta_1 + b1[i]) * vt + beta_2 * arm[i] * vt + b0[i]
    data.frame(id = i, time = vt, Y_true = y_true,
               Y_obs = y_true + rnorm(length(vt), 0, sigma_e),
               arm = arm[i], T_obs = T_obs[i])
  }))
  
  surv_data <- data.frame(id = seq_len(n), arm = arm,
                          T_true = T_i, T_obs = T_obs, status = status,
                          b0 = b0, b1 = b1)
  
  trt_df <- data.frame(id = seq_len(n),randgrp = factor(arm, levels = c(0, 1), labels = c("Control", "Experimental")))
  long_data <- merge(long_data, trt_df, by = "id")
  surv_data <- merge(surv_data, trt_df, by = "id")
  
  list(longitudinal = long_data, survival = surv_data)
}

##########



# Analysis functions
##########

fit_one_rep <- function(rep_local_idx,
                        batch_id, n_per_batch, base_seed,
                        stan_path_bp1, stan_path_bp2, stan_path_gp = NULL,
                        D, beta_0, beta_1, beta_2, sigma_e,
                        log_AF, alpha_AFT,
                        loglogistic_shape, loglogistic_scale,
                        visit, n_patients, max_FU, lambda_c = -1,
                        aft_mode = "loglogistic",
                        models_to_fit = c("LMM", "bp1", "bp2"),
                        n_stan_chains = 4,
                        n_stan_warmup = 2000,
                        n_stan_iter   = 2000) {
  
  # Load packages and simulation functions
  library(tidyverse)
  library(cmdstanr)
  library(survival)
  library(nlme)
  library(MASS)
  library(here)
  library(JMbayes2)
  
  models_to_fit <- unique(models_to_fit)
  valid_models <- c("LMM", "bp1", "bp2", "gp", "jm2")
  invalid_models <- setdiff(models_to_fit, valid_models)
  if (length(invalid_models) > 0) {
    stop(sprintf("Invalid models_to_fit: %s", paste(invalid_models, collapse = ", ")))
  }
  
  # Compile selected Stan models in this worker (uses cached binary on disk)
  stan_mod_bp1 <- if ("bp1" %in% models_to_fit) {
    cmdstanr::cmdstan_model(stan_path_bp1, quiet = TRUE, force_recompile = FALSE)
  } else {
    NULL
  }
  stan_mod_bp2 <- if ("bp2" %in% models_to_fit) {
    cmdstanr::cmdstan_model(stan_path_bp2, quiet = TRUE, force_recompile = FALSE)
  } else {
    NULL
  }
  stan_mod_gp <- if ("gp" %in% models_to_fit) {
    if (is.null(stan_path_gp)) {
      stop("stan_path_gp must be provided when 'gp' is in models_to_fit")
    }
    cmdstanr::cmdstan_model(stan_path_gp, quiet = TRUE, force_recompile = FALSE)
  } else {
    NULL
  }
  
  rep_global_idx <- (batch_id - 1) * n_per_batch + rep_local_idx
  rep_seed       <- base_seed + rep_global_idx
  
  # ------------------------------------------------------------------
  # 1. Simulate dataset
  # ------------------------------------------------------------------
  sim <- tryCatch(
    simulate_joint_dataset(
      D = D, beta_0 = beta_0, beta_1 = beta_1, beta_2 = beta_2, sigma_e = sigma_e,
      log_AF = log_AF, alpha_AFT = alpha_AFT,
      loglogistic_shape = loglogistic_shape, loglogistic_scale = loglogistic_scale,
      visit = visit, seed = rep_seed, n_patients = n_patients,
      max_FU = max_FU, lambda_c = lambda_c,
      aft_mode = aft_mode, link_type = "value"
    ),
    error = function(e) list(error = conditionMessage(e))
  )
  if (!is.null(sim$error))
    return(list(rep_global_idx = rep_global_idx, rep_seed = rep_seed,
                error = paste("simulate:", sim$error)))
  
  sim$longitudinal <- sim$longitudinal %>% dplyr::mutate(time_by_arm = time * arm)
  sim$survival <- sim$survival     %>% dplyr::rename(time = T_obs)
  
  # ------------------------------------------------------------------
  # 2. LMM (provides fixed-effect estimates and prior anchors for Stan)
  # ------------------------------------------------------------------
  t_lmm <- system.time({
    fit_lmm <- tryCatch(
      nlme::lme(
        fixed     = Y_obs ~ time + time_by_arm,
        random    = ~ time | id,
        data      = sim$longitudinal,
        na.action = na.omit,
        control   = nlme::lmeControl(maxIter = 100, msMaxIter = 100, opt = "optim")
      ),
      error = function(e) list(error = conditionMessage(e))
    )
  })
  lmm_runtime <- unname(t_lmm["elapsed"])
  if (!is.null(fit_lmm$error))
    return(list(rep_global_idx = rep_global_idx, rep_seed = rep_seed,
                error = paste("lmm:", fit_lmm$error)))
  
  lmm_coef <- as.numeric(nlme::fixef(fit_lmm))
  lmm_data <- nlme::getData(fit_lmm)
  Y_long_mean <- mean(lmm_data$Y_obs, na.rm = TRUE)
  
  # Build an LMM summary table on the centred scale (matching Stan q2.5/q97.5 convention
  # of 2.5th–97.5th percentile, i.e. 95% CI) so extract_bp_plot_data can treat LMM
  # like the Bayesian models.  beta_long[1] is shifted by -Y_long_mean to match
  # the centred Y_long; adding Y_long_mean back in result-formatting recovers
  # the raw-scale intercept estimate.
  lmm_summ <- tryCatch({
    tt    <- summary(fit_lmm)$tTable          # rows: intercept, time, time_by_arm
    se    <- tt[1:3, "Std.Error"]
    df_v  <- tt[1:3, "DF"]
    t975  <- qt(0.975, df = df_v)             # 95% CI half-width multiplier
    means_c <- c(lmm_coef[1] - Y_long_mean, lmm_coef[2], lmm_coef[3])
    tibble::tibble(
      variable = c("beta_long[1]", "beta_long[2]", "beta_long[3]"),
      mean     = means_c,
      `q2.5`   = means_c - t975 * se,
      `q97.5`  = means_c + t975 * se
    )
  }, error = function(e) NULL)
  
  # ------------------------------------------------------------------
  # 3. Construct Stan data (mirroring local-fit.R)
  # ------------------------------------------------------------------
  surv_data <- sim$survival
  mf_long   <- lmm_data
  mf_long$time_by_arm <- mf_long$time * mf_long$arm
  
  X_long <- model.matrix(~ time + time_by_arm, data = mf_long)
  Y_long_raw <- mf_long$Y_obs
  Y_long <- Y_long_raw - Y_long_mean
  id_fac <- as.factor(mf_long$id)
  J_1_long <- as.integer(id_fac)
  N_1_long <- length(unique(J_1_long))
  
  X_surv <- model.matrix(~ 0 + arm, data = surv_data)
  J_1_unique <- match(unique(surv_data$id), levels(id_fac))
  surv_data$time_by_arm <- surv_data$time * surv_data$arm
  X_long_surv <- model.matrix(~ time + time_by_arm, data = surv_data)
  
  # Avoid coxph() on HPC because it is causing an illegal-instruction crash.
  # Use survreg only to provide rough survival prior anchors / initial values.
  
  fit_stage <- tryCatch(
    survival::survreg(
      survival::Surv(time, status) ~ arm,
      data = surv_data,
      dist = "exponential",
      x = TRUE
    ),
    error = function(e) list(error = conditionMessage(e))
  )
  
  if (!is.null(fit_stage$error)) {
    return(list(
      rep_global_idx = rep_global_idx,
      rep_seed = rep_seed,
      error = paste("survreg:", fit_stage$error)
    ))
  }
  
  b_fe <- as.numeric(nlme::fixef(fit_lmm))
  # Centre the intercept to match the centred Y_long; slope/treatment coefs
  # are unaffected by the intercept shift.
  b_fe_centred <- b_fe
  b_fe_centred[1] <- b_fe_centred[1] - Y_long_mean
  
  # survreg AFT coefficient for arm; used only as a weak prior anchor.
  b_sv <- as.numeric(stats::coef(fit_stage)[-1])
  
  # Initial values for Stan survival coefficients.
  # Use survreg estimate if finite; otherwise fall back to zero.
  surv_init <- rep(0, ncol(X_surv))
  
  if (length(b_sv) > 0 && all(is.finite(b_sv))) {
    surv_init[seq_len(min(length(b_sv), ncol(X_surv)))] <-
      b_sv[seq_len(min(length(b_sv), ncol(X_surv)))]
  }
  
  s_long <- sd(Y_long_raw, na.rm = TRUE)
  s_surv <- sd(surv_data$time)
  
  stan_data1 <- list(
    K_long = ncol(X_long), q = ncol(X_surv), m = 5,
    r = 2, tau_h = 0.01,
    beta_long_prior_mean  = b_fe_centred, beta_long_prior_scale = rep(10, length(b_fe_centred)),
    beta_surv_prior_mean  = b_sv, beta_surv_prior_scale = rep(10, length(b_sv)),
    s_long = s_long, s_surv = s_surv,
    alpha_tilde_sd = (log(2) / 1.96) * (s_long / s_surv),
    N_long = nrow(X_long), Y_long = Y_long, Y_long_mean = Y_long_mean, X_long = X_long,
    N_1_long = N_1_long, J_1_long = J_1_long,
    Z_1_1_long = rep(1, nrow(X_long)), Z_1_2_long = mf_long$time,
    n = nrow(X_surv), status = as.numeric(surv_data$status), time = surv_data$time,
    X_surv = X_surv, J_1_unique = J_1_unique, X_long_surv = X_long_surv
  )
  
  G <- 100
  phi_grid <- seq(1e-6, 1 - 1e-6, length.out = G)
  
  stan_data2 <- list(
    K_long = ncol(X_long), N_long = nrow(X_long),
    Y_long = Y_long, Y_long_mean = Y_long_mean, X_long = X_long,
    N_1_long = N_1_long, J_1_long = J_1_long,
    Z_1_1_long = rep(1, nrow(X_long)), Z_1_2_long = mf_long$time,
    N_surv = nrow(X_surv), K_surv = ncol(X_surv),
    y_surv = surv_data$time, X_surv = X_surv, delta = as.numeric(surv_data$status),
    m = 5, G = G, phi_grid = phi_grid,
    s_long = s_long, s_surv = s_surv,
    alpha_tilde_sd = (log(2) / 1.96) * (s_long / s_surv),
    J_1_unique = J_1_unique, X_long_surv = X_long_surv
  )
  
  stan_data_gp <- stan_data2
  stan_data_gp$knots <- seq(0, 1, length.out = stan_data_gp$m)
  
  init_bp1 <- replicate(n_stan_chains, list(
    beta_long = b_fe_centred,
    beta_surv = surv_init,
    sigma_long = max(1e-3, as.numeric(s_long)),
    sd_1_long = rep(1, 2),
    z_1_long = matrix(0, 2, N_1_long),
    L_1_long = diag(2),
    gamma = rep(0.1, stan_data1$m),
    alpha_tilde = 0
  ), simplify = FALSE)
  
  init_bp2 <- replicate(n_stan_chains, list(
    beta_long = b_fe_centred,
    beta_surv = surv_init,
    sigma_long = max(1e-3, as.numeric(s_long)),
    alpha_tilde = 0,
    sd_1_long = rep(1, 2),
    z_1_long = matrix(0, 2, N_1_long),
    L_1_long = diag(2),
    gamma = rep(0.1, stan_data2$m)
  ), simplify = FALSE)
  
  init_gp <- replicate(n_stan_chains, list(
    beta_long = b_fe_centred,
    gamma = surv_init,
    sigma_long = max(1e-3, as.numeric(s_long)),
    alpha_tilde = 0,
    sd_1_long = rep(1, 2),
    z_1_long = matrix(0, 2, N_1_long),
    L_1_long = diag(2),
    gp_mean = 0,
    gp_sigma = 1,
    gp_length_scale = max(1e-3, stats::median(diff(stan_data_gp$knots))),
    f_gp_raw = rep(0, stan_data_gp$m)
  ), simplify = FALSE)
  
  # ------------------------------------------------------------------
  # 4. Stan fits  (n_stan_chains = 1 avoids nested parallelism inside workers)
  # ------------------------------------------------------------------
  extract_summ <- function(fit, vars) {
    if (is.null(fit)) return(NULL)
    tryCatch(
      fit$summary(vars, mean, sd,
                  rhat = posterior::rhat,
                  ess_bulk = posterior::ess_bulk,
                  ess_tail = posterior::ess_tail,
                  ~quantile(.x, c(0.025, 0.975))) %>%
        dplyr::rename(`q2.5` = `2.5%`, `q97.5` = `97.5%`),
      error = function(e) NULL
    )
  }

  compute_ic <- function(fit) {
    if (is.null(fit)) return(NULL)
    tryCatch({
      loglik_mat <- posterior::as_draws_matrix(fit$draws("log_lik"))
      loo_obj <- loo::loo(loglik_mat)
      waic_obj <- loo::waic(loglik_mat)
      list(
        loo  = loo_obj,
        waic = waic_obj
      )
    }, error = function(e) NULL)
  }

  bp_fit_1 <- NULL; bp1_runtime <- NA_real_; bp1_summ <- NULL; bp1_ic <- NULL
  bp_fit_2 <- NULL; bp2_runtime <- NA_real_; bp2_summ <- NULL; bp2_ic <- NULL
  gp_fit <- NULL; gp_runtime <- NA_real_; gp_summ <- NULL; gp_ic <- NULL
  
  if ("bp1" %in% models_to_fit) {
    t_bp1 <- system.time({
      bp_fit_1 <- tryCatch(
        stan_mod_bp1$sample(
          data = stan_data1, seed = rep_seed,
          chains = n_stan_chains, parallel_chains = n_stan_chains, max_treedepth = 10, # or 12 maybe
          iter_warmup = n_stan_warmup, iter_sampling = n_stan_iter, refresh = 0, adapt_delta = .99,
          init = init_bp1
        ),
        error = function(e) NULL
      )
    })
    bp1_runtime <- unname(t_bp1["elapsed"])
    bp1_summ <- extract_summ(bp_fit_1,
                             c("beta_long[1]", "beta_long[2]", "beta_long[3]", "beta_surv[1]", "alpha"))
    bp1_ic <- compute_ic(bp_fit_1)
  }
  
  if ("bp2" %in% models_to_fit) {
    t_bp2 <- system.time({
      bp_fit_2 <- tryCatch(
        stan_mod_bp2$sample(
          data = stan_data2, seed = rep_seed,
          chains = n_stan_chains, parallel_chains = n_stan_chains, max_treedepth = 10, # or 12 maybe
          iter_warmup = n_stan_warmup, iter_sampling = n_stan_iter, refresh = 0, adapt_delta = .99,
          init = init_bp2
        ),
        error = function(e) NULL
      )
    })
    bp2_runtime <- unname(t_bp2["elapsed"])
    bp2_summ <- extract_summ(bp_fit_2,
                             c("beta_long[1]", "beta_long[2]", "beta_long[3]", "beta_surv[1]", "alpha",
                               "sigma_long", "sd_1_long[1]", "sd_1_long[2]",
                               paste0("gamma[", seq_len(stan_data2$m), "]")))
    bp2_ic <- compute_ic(bp_fit_2)
  }
  
  if ("gp" %in% models_to_fit) {
    t_gp <- system.time({
      gp_fit <- tryCatch(
        stan_mod_gp$sample(
          data = stan_data_gp, seed = rep_seed,
          chains = n_stan_chains, parallel_chains = n_stan_chains, max_treedepth = 10, # or 12 maybe
          iter_warmup = n_stan_warmup, iter_sampling = n_stan_iter, refresh = 0, adapt_delta = .99,
          init = init_gp
        ),
        error = function(e) NULL
      )
    })
    gp_runtime <- unname(t_gp["elapsed"])
    gp_summ <- extract_summ(gp_fit,
                            c("beta_long[1]", "beta_long[2]", "beta_long[3]", "gamma[1]", "alpha",
                              "sigma_long", "sd_1_long[1]", "sd_1_long[2]",
                              paste0("f_gp_raw[", seq_len(stan_data_gp$m), "]"),
                              "gp_length_scale", "gp_alpha", "gp_mean"))
    gp_ic <- compute_ic(gp_fit)
  }
  
  # ------------------------------------------------------------------
  # 5. JMBayes2 (Cox PH joint model, value association)
  # ------------------------------------------------------------------
  jm2_runtime <- NA_real_
  jm2_summ <- NULL
  jm2_fit_error <- NULL
  jm2_summary_error <- NULL
  if ("jm2" %in% models_to_fit) {
    # Forked chains avoid the full data copies made by PSOCK workers on Linux.
    # Retain the portable snow backend for local Windows runs.
    jm2_parallel <- if (.Platform$OS.type == "windows") "snow" else "multicore"

    fit_cox <- tryCatch(
      survival::coxph(Surv(time, status) ~ arm, data = surv_data, x = TRUE),
      error = function(e) list(error = conditionMessage(e))
    )
    
    t_jm2 <- system.time({
      fit_jm2 <- if (!is.null(fit_cox$error)) {
        jm2_fit_error <- paste("coxph:", fit_cox$error)
        NULL
      } else {
        tryCatch(
          JMbayes2::jm(fit_cox, fit_lmm, time_var = "time",
                       n_chains = n_stan_chains,
                       parallel = jm2_parallel, cores = n_stan_chains,
                       seed = rep_seed),
          error = function(e) {
            jm2_fit_error <<- conditionMessage(e)
            NULL
          }
        )
      }
    })
    jm2_runtime <- unname(t_jm2["elapsed"])
    
    if (!is.null(fit_jm2)) {
      jm2_summ <- tryCatch({
        s <- summary(fit_jm2)

        # JMbayes2 summary.jm stores these tables directly as Outcome1 and
        # Survival. Retain fallbacks for older summary layouts.
        long_summary <- if (!is.null(s$Outcome1)) s$Outcome1 else s$Longitudinal$Outcome1
        surv_summary <- if (!is.null(s$Survival)) s$Survival else s$Event

        if (is.null(long_summary) || is.null(surv_summary)) {
          stop("JMbayes2 summary did not contain Outcome1 and Survival tables")
        }

        long_tbl <- as.data.frame(long_summary) %>%
          tibble::rownames_to_column("variable") %>%
          dplyr::select(variable, mean = Mean, sd = StDev, q2.5 = `2.5%`, q97.5 = `97.5%`)
        surv_tbl <- as.data.frame(surv_summary) %>%
          tibble::rownames_to_column("variable") %>%
          dplyr::select(variable, mean = Mean, sd = StDev, q2.5 = `2.5%`, q97.5 = `97.5%`)
        dplyr::bind_rows(long_tbl, surv_tbl)
      }, error = function(e) {
        jm2_summary_error <<- conditionMessage(e)
        NULL
      })
    }
  }
  
  list(
    rep_global_idx = rep_global_idx,
    rep_seed = rep_seed,
    models_to_fit = models_to_fit,
    Y_long_mean = Y_long_mean,
    LMM = if ("LMM" %in% models_to_fit) list(coef = lmm_coef, summary = lmm_summ, runtime = lmm_runtime) else NULL,
    bp1 = if ("bp1" %in% models_to_fit) list(summary = bp1_summ, runtime = bp1_runtime, ic = bp1_ic) else NULL,
    bp2 = if ("bp2" %in% models_to_fit) list(summary = bp2_summ, runtime = bp2_runtime, ic = bp2_ic) else NULL,
    gp = if ("gp" %in% models_to_fit) list(summary = gp_summ, runtime = gp_runtime, ic = gp_ic) else NULL,
    jm2 = if ("jm2" %in% models_to_fit) {
      list(
        summary = jm2_summ,
        runtime = jm2_runtime,
        fit_error = jm2_fit_error,
        summary_error = jm2_summary_error
      )
    } else NULL
  )
}


##########





