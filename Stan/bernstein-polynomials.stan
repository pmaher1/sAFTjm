// Functions
functions {
  // Bernstein Polynomial joiny model
  vector survival_loglik_bp_jm(
    vector kappa, // accelerated time
    vector delta,
    vector lin_pred, // linear predictor
    vector gamma, // polynomial weights
    vector bp_pdf_coef, // bernstein basis coefficients (binomial coefficients)
    vector bp_pdf_to_cdf_coef,
    int m  // degree of bernstein polynomial
  ) {
    int N = num_elements(kappa);
    real eps = 1e-12;

    real tau_aft = max(kappa) + eps;
    vector[N] phi = kappa / tau_aft;
    vector[N] omphi = 1 - phi;

    vector[N] h0 = rep_vector(0.0, N);
    vector[N] H0 = rep_vector(0.0, N);

    vector[m] gamma_cum = cumulative_sum(gamma);
    vector[N] phi_pw = rep_vector(1.0, N); // phi^(k-1), starts at k=1
    vector[N] omphi_pw = pow(omphi, m - 1); // (1-phi)^(m-k), starts at k=1

    for (k in 1:m) {
      vector[N] bern = bp_pdf_coef[k] * phi_pw .* omphi_pw;
      h0 += gamma[k] * bern;
      H0 += gamma_cum[k] * bp_pdf_to_cdf_coef[k] * (phi .* bern);
      if (k < m) {
        phi_pw .*= phi;
        omphi_pw ./= omphi;
      }
    }
    h0 *= m / tau_aft;
    h0 = fmax(h0, eps);

    return delta .* (log(h0) - lin_pred) - H0;
  }
}


// Data 
data {
  
  //// longitudinal data
  int<lower=1> K_long; // number of population-level effects
  int<lower=1> N_long; // total number of observations
  vector[N_long] Y_long; // response variable
  matrix[N_long, K_long] X_long; // population-level design matrix
  int<lower=1> N_1_long; // number of individual levels (should equate to number of individuals)
  array[N_long] int<lower=1> J_1_long; // group indicator from ID (rank, ordered)
  
  // individual/group-level predictor values (longitudinal)
  vector[N_long] Z_1_1_long; // intercept
  vector[N_long] Z_1_2_long; // time 
  
  
  //// survival data
  int<lower=1> N_surv; // number of observations
  int<lower=1> K_surv; // number of covariates
  vector[N_surv] y_surv; // observed times
  matrix[N_surv,K_surv] X_surv; // covariate matrix
  vector[N_surv] delta; // censoring indicators [1: event, 0: censored]
  
  // bp specific
  int<lower=1> m; // order of bernstein polynomial
  
  // plotting purposes
  int<lower=1> G; // num of grid points
  vector[G] phi_grid; // grid for phi (scaled time) for plotting baseline hazards and cumulative hazards
  
  
  // association parameter prior calibration
  real s_long; 
  real s_surv;
  real<lower=0> alpha_tilde_sd;

  //// linking data
  array[N_surv] int<lower=1> J_1_unique;
  matrix[N_surv, K_long] X_long_surv;
} 

// Transformed Data
transformed data {
  // Bernstein basis constants precomputed once rather than at every HMC step.
  vector[m] bp_pdf_coef;
  vector[m] bp_pdf_to_cdf_coef;
  for (k in 1:m) {
    bp_pdf_coef[k] = choose(m - 1, k - 1);
    bp_pdf_to_cdf_coef[k] = m * 1.0 / k;
  }
}

// Parameters
parameters {
  //// longitudinal parameters
  vector[K_long] beta_long; // regression coefficients 
  real<lower=0> sigma_long; // dispersion parameter
  vector<lower=1e-3>[2] sd_1_long; // group/individual-level standard deviations
  matrix[2, N_1_long] z_1_long; // individual-level random effects
  cholesky_factor_corr[2] L_1_long; // cholesky factor of correlation matrix
  
  //// survival parameters
  vector[K_surv] beta_surv; // survival regression coefficients

  // bp params
  vector<lower=0>[m] gamma; // BP basis weights for baseline hazard

  //// linking parameter
  real alpha_tilde; // standardised association parameter
}

// Transformed parameters
transformed parameters{
  
  //// longitudinal transformed parameters
  matrix[2, N_1_long] b_raw = diag_pre_multiply(sd_1_long, L_1_long) * z_1_long;
  matrix[N_1_long, 2] b_long = b_raw';
  
  //// survival transformed parameters
  real alpha = alpha_tilde * (s_surv / s_long);

  vector[N_surv] lin_pred; // linear predictor at event/censoring time (current value)
  vector[N_surv] Y_long_surv; // predicted Y_i(t) at event/censoring time
  vector[N_surv] kappa_integral; // integral of exp(-lin_pred(t)) from 0 to event time

  {
    vector[N_surv] eta_surv = X_surv * beta_surv;
    for (i in 1:N_surv) {
      real C1 = eta_surv[i] + alpha * (beta_long[1] + b_long[J_1_unique[i], 1]);
      real C2 = alpha * (beta_long[2] + X_surv[i, 1] * beta_long[3] + b_long[J_1_unique[i], 2]);
      Y_long_surv[i] = beta_long[1]
                     + beta_long[2] * X_long_surv[i, 2]
                     + beta_long[3] * X_long_surv[i, 3]
                     + b_long[J_1_unique[i], 1]
                     + b_long[J_1_unique[i], 2] * y_surv[i];
      kappa_integral[i] = (abs(C2) > 1e-10)
                          ? exp(-C1) * (-expm1(-C2 * y_surv[i])) / C2
                          : exp(-C1) * y_surv[i];
    }
    lin_pred = eta_surv + alpha * Y_long_surv;
  }

}

// Model
model {
  
  //// longitudinal priors
  beta_long ~ normal(0, 10);
  sd_1_long ~ cauchy(0, 5);
  sigma_long ~ cauchy(0, 5);
  L_1_long ~ lkj_corr_cholesky(2);
  to_vector(z_1_long) ~ std_normal();
  
  // longitudinal likelihood
  vector[N_long] mu_long = X_long * beta_long
                         + b_long[J_1_long, 1] .* Z_1_1_long
                         + b_long[J_1_long, 2] .* Z_1_2_long;
  Y_long ~ normal(mu_long, sigma_long);
  
  //// time to event priors
  beta_surv ~ normal(0, 10);
  alpha_tilde ~ normal(0, alpha_tilde_sd);
  gamma ~ normal(0, 1) T[0, ];

  //// likelihood
  target += sum(survival_loglik_bp_jm(kappa_integral, delta, lin_pred, gamma, bp_pdf_coef, bp_pdf_to_cdf_coef, m));

}

// Generated quantities
generated quantities {
  // joint likelihood specification for WAIC, LOO and DIC statistics  
  vector[N_surv] log_lik;
  {
    vector[N_surv] log_lik_surv = survival_loglik_bp_jm(kappa_integral, delta, lin_pred, gamma, bp_pdf_coef, bp_pdf_to_cdf_coef, m);
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

