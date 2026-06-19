
# Pre-processing
##########

# Load relevant libraries
library(here) # file path management
library(tidyverse) # data manipulation and plotting
dir.create(here("Figures"), showWarnings = FALSE, recursive = TRUE)

##########


# Weibull baseline hazard
wb_h0 <- function(time, alpha, beta) {
  (alpha / beta) * (time / beta)^(alpha - 1)
}

# Weibull cumulative baseline hazard
wb_H0 <- function(time, alpha, beta) {
  (time / beta)^alpha
}

# Log-logistic baseline hazard
ll_h0 <- function(time, alpha, beta) {
  (alpha / beta) * (time / beta)^(alpha - 1) / (1 + (time / beta)^alpha)
}

# Log-logistic cumulative baseline hazard
ll_H0 <- function(time, alpha, beta) {
  log(1 + (time / beta)^alpha)
}

# Bernstein polynomial cumulative baseline hazard
bp_H0 <- function(phi, theta) {
  m <- length(theta) - 1
  H0 <- rep(0, length(phi))
  for (j in 0:m) {
    H0 <- H0 + theta[j + 1] * choose(m, j) * phi^j * (1 - phi)^(m - j)
  }
  H0
}

# Bernstein polynomial derivative on the scaled-time scale
bp_h0_phi <- function(phi, theta) {
  m <- length(theta) - 1
  h_phi <- rep(0, length(phi))
  for (j in 1:m) {
    h_phi <- h_phi + m * (theta[j + 1] - theta[j]) *
      choose(m - 1, j - 1) * phi^(j - 1) * (1 - phi)^(m - j)
  }
  h_phi
}

# Bernstein polynomial baseline hazard on the actual-time scale
bp_h0 <- function(phi, theta, tau = 1) {
  bp_h0_phi(phi, theta) / tau
}

# Bernstein polynomial hazard with respect to actual time under AFT
bp_h <- function(time, x, gamma, theta, tau = NULL) {
  eta <- as.numeric(gamma %*% x)
  kappa <- time * exp(-eta)
  if (is.null(tau)) {
    tau <- max(kappa)
  }
  phi <- kappa / tau
  phi <- pmin(pmax(phi, 1e-10), 1 - 1e-10)
  bp_h0_phi(phi, theta) * exp(-eta) / tau
}

# Pull posterior means for indexed parameters from a saved summary table
summary_vector <- function(summary_tbl, pattern) {
  if (is.null(summary_tbl) || !"variable" %in% names(summary_tbl) || !"mean" %in% names(summary_tbl)) {
    return(NULL)
  }
  rows <- grep(pattern, summary_tbl$variable)
  if (!length(rows)) {
    return(NULL)
  }
  vars <- summary_tbl$variable[rows]
  idx <- as.integer(sub(".*\\[([0-9]+)\\].*", "\\1", vars))
  summary_tbl$mean[rows][order(idx)]
}

summary_scalar <- function(summary_tbl, variable) {
  if (is.null(summary_tbl) || !"variable" %in% names(summary_tbl) || !"mean" %in% names(summary_tbl)) {
    return(NULL)
  }
  rows <- which(summary_tbl$variable == variable)
  if (!length(rows)) {
    return(NULL)
  }
  summary_tbl$mean[rows[1]]
}

cumtrapz <- function(x, y) {
  c(0, cumsum(diff(x) * (head(y, -1) + tail(y, -1)) / 2))
}

# Reconstruct the BP baseline from fitted basis weights if they were saved.
reconstruct_bp_baseline <- function(summary_tbl, time_grid, tau = max(time_grid)) {
  gamma <- summary_vector(summary_tbl, "^gamma\\[[0-9]+\\]$")
  if (!is.null(gamma)) {
    m <- length(gamma)
    phi <- pmin(pmax(time_grid / tau, 1e-6), 1 - 1e-6)
    h0 <- rep(0, length(phi))
    H0 <- rep(0, length(phi))
    gamma_cum <- cumsum(gamma)
    
    for (k in seq_len(m)) {
      bern <- choose(m - 1, k - 1) * phi^(k - 1) * (1 - phi)^(m - k)
      h0 <- h0 + gamma[k] * bern
      H0 <- H0 + gamma_cum[k] * (m / k) * phi * bern
    }
    
    return(data.frame(time = time_grid, H0 = H0, h0 = h0 * m / tau))
  }
  
  theta <- summary_vector(summary_tbl, "^theta\\[[0-9]+\\]$")
  if (!is.null(theta)) {
    phi <- pmin(pmax(time_grid / tau, 1e-6), 1 - 1e-6)
    return(data.frame(
      time = time_grid,
      H0 = bp_H0(phi, theta),
      h0 = bp_h0(phi, theta, tau = tau)
    ))
  }
  
  NULL
}

# Reconstruct the GP baseline from fitted log-hazard values if they were saved.
reconstruct_gp_baseline <- function(summary_tbl, time_grid, tau = max(time_grid)) {
  f_gp <- summary_vector(summary_tbl, "^f_gp\\[[0-9]+\\]$")
  
  if (is.null(f_gp)) {
    f_gp_raw <- summary_vector(summary_tbl, "^f_gp_raw\\[[0-9]+\\]$")
    gp_mean <- summary_scalar(summary_tbl, "gp_mean")
    gp_alpha <- summary_scalar(summary_tbl, "gp_alpha")
    gp_length_scale <- summary_scalar(summary_tbl, "gp_length_scale")
    
    if (is.null(f_gp_raw) || is.null(gp_mean) || is.null(gp_alpha) || is.null(gp_length_scale)) {
      return(NULL)
    }
    
    knots <- seq(0, 1, length.out = length(f_gp_raw))
    K_gp <- outer(knots, knots, function(a, b) {
      gp_alpha^2 * exp(-0.5 * ((a - b) / gp_length_scale)^2)
    })
    diag(K_gp) <- diag(K_gp) + 1e-6
    f_gp <- gp_mean + as.vector(t(chol(K_gp)) %*% f_gp_raw)
  }
  
  knots <- seq(0, 1, length.out = length(f_gp))
  phi <- pmin(pmax(time_grid / tau, 1e-6), 1 - 1e-6)
  phi_ext <- c(0, phi)
  h_scaled_ext <- exp(approx(knots, f_gp, xout = phi_ext, rule = 2)$y)
  
  data.frame(
    time = time_grid,
    H0 = cumtrapz(phi_ext, h_scaled_ext)[-1],
    h0 = h_scaled_ext[-1] / tau
  )
}

random_model_summary <- function(folder_pattern, model_name) {
  results_dir <- here("Simulation-Workflow", "Results")
  rds_files <- list.files(results_dir, pattern = "\\.rds$", recursive = TRUE, full.names = TRUE)
  rds_files <- rds_files[grepl(folder_pattern, rds_files)]
  if (!length(rds_files)) {
    return(NULL)
  }
  
  for (rds_file in sample(rds_files)) {
    rds_obj <- tryCatch(readRDS(rds_file), error = function(e) NULL)
    if (is.null(rds_obj) || is.null(rds_obj$results)) {
      next
    }
    
    for (result in rds_obj$results) {
      summary_tbl <- result[[model_name]]$summary
      if (!is.null(summary_tbl)) {
        return(list(file = rds_file, summary = summary_tbl))
      }
    }
  }
  
  NULL
}





# Sample baseline hazard plots with larger titles and labels
time_grid <- seq(0.01, 120, length.out = 100)

# Larger plot text theme
big_plot_theme <- theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 13)
  )

## Log-logistic baseline

ll_true_alpha <- 1.2 # shape
ll_true_beta  <- 23  # scale

ll_true_H0 <- ll_H0(time_grid, alpha = ll_true_alpha, beta = ll_true_beta)
ll_true_h0 <- ll_h0(time_grid, alpha = ll_true_alpha, beta = ll_true_beta)

ll_H0_plot <- ggplot(data.frame(time = time_grid, H0 = ll_true_H0), aes(x = time, y = H0)) +
  geom_line(color = "#3498DB", linewidth = 2) +
  labs(
    title = "True Cumulative Baseline Hazard\n(Log-Logistic)",
    x = "Time (months)",
    y = expression(Lambda[0](t))
  ) +
  big_plot_theme

ll_h0_plot <- ggplot(data.frame(time = time_grid, h0 = ll_true_h0), aes(x = time, y = h0)) +
  geom_line(color = "#3498DB", linewidth = 2) +
  labs(
    title = "True Baseline Hazard\n(Log-Logistic)",
    x = "Time (months)",
    y = expression(lambda[0](t))
  ) +
  big_plot_theme

ll_combined_plot <- cowplot::plot_grid(ll_H0_plot, ll_h0_plot, ncol = 2)
ll_combined_plot

ggsave(
  here("Figures", "log_logistic_baseline_hazard.pdf"),
  ll_combined_plot,
  width = 14,
  height = 6
)


## Weibull baseline

wb_true_alpha <- 0.9 # shape
wb_true_beta  <- 38  # scale

wb_true_H0 <- wb_H0(time_grid, alpha = wb_true_alpha, beta = wb_true_beta)
wb_true_h0 <- wb_h0(time_grid, alpha = wb_true_alpha, beta = wb_true_beta)

wb_H0_plot <- ggplot(data.frame(time = time_grid, H0 = wb_true_H0), aes(x = time, y = H0)) +
  geom_line(color = "#E74C3C", linewidth = 2) +
  labs(
    title = "True Cumulative Baseline Hazard\n(Weibull)",
    x = "Time (months)",
    y = expression(Lambda[0](t))
  ) +
  big_plot_theme

wb_h0_plot <- ggplot(data.frame(time = time_grid, h0 = wb_true_h0), aes(x = time, y = h0)) +
  geom_line(color = "#E74C3C", linewidth = 2) +
  labs(
    title = "True Baseline Hazard\n(Weibull)",
    x = "Time (months)",
    y = expression(lambda[0](t))
  ) +
  big_plot_theme

wb_combined_plot <- cowplot::plot_grid(wb_H0_plot, wb_h0_plot, ncol = 2)
wb_combined_plot

ggsave(
  here("Figures", "weibull_baseline_hazard.pdf"),
  wb_combined_plot,
  width = 14,
  height = 6
)


# Reconstruct one fitted BP and one fitted GP baseline hazard from random RDS outputs
bp_example <- random_model_summary("_bp[/\\\\]", "bp2")
gp_example <- random_model_summary("_gp[/\\\\]", "gp")

bp_fit_baseline <- if (!is.null(bp_example)) {
  reconstruct_bp_baseline(bp_example$summary, time_grid)
} else {
  NULL
}

gp_fit_baseline <- if (!is.null(gp_example)) {
  reconstruct_gp_baseline(gp_example$summary, time_grid)
} else {
  NULL
}

fit_baseline_df <- bind_rows(
  if (!is.null(bp_fit_baseline)) mutate(bp_fit_baseline, model = "BP"),
  if (!is.null(gp_fit_baseline)) mutate(gp_fit_baseline, model = "GP")
)

if (nrow(fit_baseline_df) > 0) {
  fit_H0_plot <- ggplot(fit_baseline_df, aes(x = time, y = H0, color = model)) +
    geom_line(linewidth = 2) +
    scale_color_manual(values = c(BP = "#2C7FB8", GP = "#31A354")) +
    labs(title = "Fitted Cumulative Baseline Hazard (BP and GP)", x = "time (months)", y = "H0(t)", color = "Model") +
    theme_bw()
  fit_H0_plot
  
  fit_h0_plot <- ggplot(fit_baseline_df, aes(x = time, y = h0, color = model)) +
    geom_line(linewidth = 2) +
    scale_color_manual(values = c(BP = "#2C7FB8", GP = "#31A354")) +
    labs(title = "Fitted Baseline Hazard (BP and GP)", x = "time (months)", y = "h0(t)", color = "Model") +
    theme_bw()
  fit_h0_plot
  
  fit_combined_plot <- cowplot::plot_grid(fit_H0_plot, fit_h0_plot, ncol = 2)
  fit_combined_plot
  
  ggsave(here("Figures", "bp_gp_reconstructed_baseline_hazard.pdf"), fit_combined_plot, width = 12, height = 6)
} else {
  warning(
    "Could not reconstruct BP/GP fitted baseline hazards from the sampled RDS outputs. ",
    "The saved summaries need BP gamma/theta terms and GP f_gp or f_gp_raw/gp_* terms."
  )
}








