// Functions
functions {

  // AFT joint survival log likelihood with GP
  vector survival_loglik_gp_jm(
    vector kappa_integral,
    vector delta,
    vector lin_pred,
    vector f_gp, // GP function values at knots (log-scale)
    vector knots, // knot locations in [0,1], length m
    real tau_aft
  )
  {
    int n = num_elements(delta);
    vector[n] log_h0;
    vector[n] H0;
    int m = num_elements(knots);

    real eps = 1e-6;

    // Scale integrated hazard to [0,1]
    vector[n] y_alt = fmin(fmax(kappa_integral / tau_aft, eps), 1 - eps);

    // Interpolate GP values (log hazard) using linear interpolation at each y_alt[i]
    for (i in 1:n) {
      real t = y_alt[i];

      int idx = 1;
      for (k in 1:(m - 1)) {
        if (t > knots[k]) {
          idx = k;
        }
      }

      if (idx < m) {
        real w = (t - knots[idx]) / (knots[idx + 1] - knots[idx]);
        log_h0[i] = (1-w)*f_gp[idx] + w*f_gp[idx + 1];
      } else {
        log_h0[i] = f_gp[m];
      }
    }

    // Cumulative hazard via trapezoidal rule
    for (i in 1:n) {
      real t = y_alt[i];
      real cum_haz = 0;

      for (k in 1:(m - 1)) {
        if (knots[k + 1] <= t) {
          // Full interval
          real h_left = exp(f_gp[k]);
          real h_right = exp(f_gp[k + 1]);
          cum_haz += 0.5 * (h_left + h_right) * (knots[k + 1] - knots[k]);
        } else if (knots[k] < t) {
          // Partial interval
          real w = (t - knots[k]) / (knots[k + 1] - knots[k]);
          real f_interp = (1 - w) * f_gp[k] + w * f_gp[k + 1];
          real h_left = exp(f_gp[k]);
          real h_right = exp(f_interp);
          cum_haz += 0.5 * (h_left + h_right) * (t - knots[k]);
        }
      }

      // H0[i] = cum_haz * tau_aft;
      H0[i] = cum_haz;
    }

    // log-lik: delta * (log_h0 - lin_pred) + (-H0)
    return delta .* (log_h0 - log(tau_aft) - lin_pred) - H0;
  }

}


// Data
data {

  //// longitudinal data
  int<lower=1> K_long; // number of population-level effects
  int<lower=1> N_long; // total number of observations
  vector[N_long] Y_long; // response variable
  matrix[N_long, K_long] X_long; // population-level design matrix
  int<lower=1> N_1_long; // number of individuals
  array[N_long] int<lower=1> J_1_long; // group indicator from ID
  
  // individual-level predictor values (longitudinal)
  vector[N_long] Z_1_1_long; // intercept
  vector[N_long] Z_1_2_long; // time

  //// survival data
  int<lower=1> N_surv; // number of observations
  int<lower=1> K_surv; // number of covariates
  vector[N_surv] y_surv; // observed times
  matrix[N_surv, K_surv] X_surv; // covariate matrix
  vector[N_surv] delta; // event indicators [1: event, 0: censored]

  //// GP-specific
  int<lower=2> m; // number of GP knots
  vector[m] knots; // knot locations in [0,1]

  //// association parameter prior calibration
  real s_long;
  real s_surv;
  real<lower=0> alpha_tilde_sd;

  //// plotting
  int<lower=1> G; // number of grid points
  vector[G] phi_grid; // grid for phi (scaled time)

  //// linking data
  array[N_surv] int<lower=1> J_1_unique;
  matrix[N_surv, K_long] X_long_surv;

}


// Parameters
parameters {

  //// longitudinal parameters
  vector[K_long] beta_long;
  real<lower=0> sigma_long;
  vector<lower=1e-3>[2] sd_1_long;
  matrix[2, N_1_long] z_1_long;
  cholesky_factor_corr[2] L_1_long;

  //// survival parameters
  vector[K_surv] gamma;

  //// GP parameters
  vector[m] f_gp_raw; // raw (non-centred) GP values at knots
  real<lower=0> gp_length_scale; // GP squared-exponential length-scale
  real<lower=0> gp_alpha; // GP marginal standard deviation
  real gp_mean; // GP mean (log baseline hazard intercept)

  //// linking parameter
  real alpha_tilde;

}


// Transformed parameters
transformed parameters {

  //// longitudinal random effects
  matrix[2, N_1_long] b_raw  = diag_pre_multiply(sd_1_long, L_1_long) * z_1_long;
  matrix[N_1_long, 2] b_long = b_raw';

  //// survival quantities
  real alpha = alpha_tilde * (s_surv / s_long);
  vector[N_surv] lin_pred;
  vector[N_surv] Y_long_surv;
  vector[N_surv] kappa_integral;

  {
    vector[N_surv] eta_surv = X_surv * gamma;

    for (i in 1:N_surv) {
      real C1 = eta_surv[i] + alpha * (beta_long[1] + b_long[J_1_unique[i], 1]);
      real C2 = alpha * (beta_long[2]
                         + X_surv[i, 1] * beta_long[3]
                         + b_long[J_1_unique[i], 2]);

      Y_long_surv[i] = beta_long[1]
                     + beta_long[2] * X_long_surv[i, 2]
                     + beta_long[3] * X_long_surv[i, 3]
                     + b_long[J_1_unique[i], 1]
                     + b_long[J_1_unique[i], 2] * y_surv[i];

      kappa_integral[i] =
        (abs(C2) > 1e-10) ? exp(-C1) * (-expm1(-C2 * y_surv[i])) / C2: exp(-C1) * y_surv[i];
    }

    lin_pred = eta_surv + alpha * Y_long_surv;
  }

  //// GP: squared-exponential kernel -> Cholesky -> non-centred parameterisation
  matrix[m, m] K_gp;
  for (i in 1:m) {
    for (j in 1:m) {
      K_gp[i, j] = square(gp_alpha) * exp(-0.5 * square((knots[i] - knots[j]) / gp_length_scale));
    }
    K_gp[i, i] += 1e-9; // jitter for numerical stability
  }
  matrix[m, m] L_gp = cholesky_decompose(K_gp);

  // GP values on log-hazard scale
  vector[m] f_gp = gp_mean + L_gp * f_gp_raw;

  real tau_aft = max(kappa_integral) + 1e-6;

}


// Model
model {

  //// longitudinal priors
  beta_long  ~ normal(0, 10);
  sd_1_long ~ cauchy(0, 5);
  sigma_long ~ cauchy(0, 5);
  L_1_long ~ lkj_corr_cholesky(2);
  to_vector(z_1_long) ~ std_normal();

  // longitudinal likelihood
  vector[N_long] mu_long = X_long * beta_long;
  for (i in 1:N_long) {
    mu_long[i] += b_long[J_1_long[i], 1] * Z_1_1_long[i]
                + b_long[J_1_long[i], 2] * Z_1_2_long[i];
  }
  Y_long ~ normal(mu_long, sigma_long);

  //// survival priors
  alpha_tilde ~ normal(0, alpha_tilde_sd);
  gamma ~ normal(0, 10);

  //// GP priors
  gp_length_scale ~ inv_gamma(5, 5); // weakly informative; centres around 1
  gp_alpha ~ cauchy(0, 0.5);
  gp_mean ~ normal(-1, 1);
  f_gp_raw ~ std_normal(); // non-centred parameterisation

  //// survival likelihood
  target += sum(survival_loglik_gp_jm(kappa_integral, delta, lin_pred, f_gp, knots, tau_aft));

}

// Generated quantities
generated quantities {
  // joint likelihood specification for WAIC, LOO and DIC statistics
  vector[N_surv] log_lik;
  {
    vector[N_surv] log_lik_surv = survival_loglik_gp_jm(kappa_integral, delta, lin_pred, f_gp, knots, tau_aft);
    vector[N_surv] log_lik_long = rep_vector(0, N_surv);
    vector[N_long] mu_long_ic = X_long * beta_long
                              + b_long[J_1_long, 1] .* Z_1_1_long
                              + b_long[J_1_long, 2] .* Z_1_2_long;

    for (j in 1:N_long) {
      // J_1_long maps each repeated longitudinal observation back to its patient.
      log_lik_long[J_1_long[j]] += normal_lpdf(Y_long[j] | mu_long_ic[j], sigma_long);
    }

    log_lik = log_lik_surv + log_lik_long;
  }
  
}

