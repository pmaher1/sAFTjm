

# Pre-processing
##########

# Load relevant libraries
library(here) # file path management
library(tidyverse) # data manipulation and plotting
library(purrr) # for functional programming (e.g. map)
library(tibble) # for tidy data frames
library(knitr) # for tables in writeup

# Creation of output directories for the writeup
dir.create(here("Writeup", "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(here("Writeup", "pictures"), showWarnings = FALSE, recursive = TRUE)
dir.create(here("Writeup", "tables"), showWarnings = FALSE, recursive = TRUE)

##########




# Scenario detection
##########

# Directory of results
base_results_dir <- here("Simulation-Workflow", "Results")

# Find all scenario folders, including optional sample-size and model suffixes:
#   sN_ll, sN_ll_n300, sN_ll_n300_bp, sN_wb_n300_gp, sN_ll_n300_ph, etc.
all_scenario_dirs <- list.dirs(base_results_dir, full.names = FALSE, recursive = FALSE)
scenario_dirs <- sort(grep("^s[0-9]+_(ll|wb)(_n[0-9]+)?(_bp|_gp)?$",
                           all_scenario_dirs, value = TRUE, perl = TRUE))
ph_scenario_dirs <- sort(grep("^s[0-9]+_(ll|wb)(_n[0-9]+)?_ph$",
                              all_scenario_dirs, value = TRUE, perl = TRUE))

# Also check for default (loose) files in Results/ root
default_files <- list.files(base_results_dir, pattern = "WPP-MWE-.*\\.rds",
                            full.names = TRUE, recursive = FALSE)

scenario_list <- list()

scenario_key <- function(x) sub("_(bp|gp|ph)$", "", x)

scenario_aft_mode <- function(base_sf) {
  if (grepl("_ll(_n[0-9]+)?$", base_sf, perl = TRUE)) {
    return("loglogistic")
  }
  "weibull"
}

scenario_label <- function(base_sf) {
  m <- regexec("^s([0-9]+)_(ll|wb)(?:_n([0-9]+))?$", base_sf, perl = TRUE)
  parts <- regmatches(base_sf, m)[[1]]
  snum <- parts[[2]]
  baseline <- if (parts[[3]] == "ll") "Log-Logistic" else "Weibull"
  n_label <- if (length(parts) >= 4 && nzchar(parts[[4]])) {
    paste0(" (n=", parts[[4]], ")")
  } else {
    ""
  }
  paste0("Scenario ", snum, ": ", baseline, n_label)
}

if (length(scenario_dirs) > 0) {
  base_scenarios <- unique(scenario_key(scenario_dirs))
  for (base_sf in base_scenarios) {
    variants <- scenario_dirs[scenario_key(scenario_dirs) == base_sf]
    aft_mode <- scenario_aft_mode(base_sf)
    files <- sort(unlist(lapply(variants, function(sf) {
      list.files(here("Simulation-Workflow", "Results", sf),
                 pattern = paste0("WPP-", base_sf, "-", aft_mode, "-.*\\.rds"),
                 full.names = TRUE)
    })))
    if (length(files) > 0) {
      scenario_list[[base_sf]] <- list(folder = base_sf,
                                       label = scenario_label(base_sf),
                                       files = files)
    }
  }
}

ph_scenario_list <- list()
if (length(ph_scenario_dirs) > 0) {
  for (sf in ph_scenario_dirs) {
    base_sf <- scenario_key(sf)
    aft_mode <- scenario_aft_mode(base_sf)
    files <- sort(list.files(
      here("Simulation-Workflow", "Results", sf),
      pattern = paste0("WPP-", base_sf, "-", aft_mode, "-.*\\.rds"),
      full.names = TRUE
    ))
    if (length(files) > 0) {
      ph_scenario_list[[base_sf]] <- list(
        folder = sf,
        label = scenario_label(base_sf),
        files = files
      )
    }
  }
}

if (length(scenario_list) == 0 && length(default_files) > 0) {
  scenario_list[["default"]] <- list(folder = "default", label = "Default", files = default_files)
}

if (length(scenario_list) == 0) {
  stop("No simulation result files found in any scenario folder under Simulation-Workflow/Results/")
}

message("Found ", length(scenario_list), " scenario(s): ",
        paste(sapply(scenario_list, `[[`, "label"), collapse = ", "))
message("Found ", length(ph_scenario_list), " PH/LMM scenario(s).")

# Set PH_ONLY=1 to refresh only the PH/LMM comparison artifacts without
# regenerating all BP/GP scenario outputs first.
run_standard_analysis <- Sys.getenv("PH_ONLY", unset = "0") != "1"

##########



# Scenario colour mapping (built after scenario_list is final)
##########

### Nice colour sets

## Colour set 1

colour_set1 <- c(
  "#C96480", # Blush Rose
  "#B47978", # Dusty Rose
  "#B1AE91", # Dry Sage
  "#95BF8F", # Muted Olive
  "#99D17B") # Willow Green

colour_set1_dark <- c(
  "#974B60", # Dark Blush Rose
  "#875B5A", # Dark Dusty Rose
  "#85836D", # Dark Dry Sage
  "#70936B", # Dark Muted Olive
  "#739D5C") # Dark Willow Green

colour_set1_light <- c(
  "#D78BA0", # Light Blush Rose
  "#C79B9A", # Light Dusty Rose
  "#C5C2AD", # Light Dry Sage
  "#B0CFAB", # Light Muted Olive
  "#B3DD9C") # Light Willow Green

## Colour set 2

colour_set2 <- c(
  "#E5F4E3", # Honeydew
  "#5DA9E9", # Cool Sky
  "#003F91", # French Blue
  "#FFFFFF", # White
  "#6D326D") # Velvet Purple

colour_set2_dark <- c(
  "#B8D9B3", # Dark Honeydew
  "#4A8CD1", # Dark Cool Sky
  "#002F6B", # Dark French Blue
  "#CCCCCC", # Light Grey (instead of dark white)
  "#4D244D") # Dark Velvet Purple

colour_set2_light <- c(
  "#E5F4E3", # Light Honeydew (same as base)
  "#A9D0F7", # Light Cool Sky
  "#6691E8", # Light French Blue
  "#FFFFFF", # White (same as base)
  "#9B5DA9") # Light Velvet Purple

# Fill-group palette for cross_scenario_boxplot_combined.
# Key format: "<ModelType>_<dgm>_S<N>"  (e.g. "BP_loglogistic_S2").
# Each model uses one colour; Weibull uses a lighter variant.
.make_pal <- function(col, model_type, dgm) {
  setNames(rep(col, 5), paste(model_type, dgm, paste0("S", 1:5), sep = "_"))
}
boxplot_cols <- c(
  .make_pal("#2E7D32", "LMM", "loglogistic"),
  .make_pal("#81C784", "LMM", "weibull"),
  .make_pal("#C62828", "BP", "loglogistic"),
  .make_pal("#EF9A9A", "BP", "weibull"),
  .make_pal("#1565C0", "GP", "loglogistic"),
  .make_pal("#90CAF9", "GP", "weibull")
)
rm(.make_pal)

primary_parameter_order <- c("beta_0", "beta_1", "beta_2", "alpha_AFT", "log_AF")

# Assign colour_set1 to ll scenarios and colour_set1_dark to wb scenarios,
# cycling within each group if there are more than 5 of either type.
# (Built here so it is available for the cross-scenario plots below.)
.ll_keys <- grep("_ll(_n[0-9]+)?$", names(scenario_list), value = TRUE, perl = TRUE)
.wb_keys <- grep("_wb(_n[0-9]+)?$", names(scenario_list), value = TRUE, perl = TRUE)

scenario_colour_map <- c(
  setNames(
    colour_set1[((seq_along(.ll_keys) - 1) %% length(colour_set1)) + 1],
    sapply(scenario_list[.ll_keys], `[[`, "label")
  ),
  setNames(
    colour_set1_dark[((seq_along(.wb_keys) - 1) %% length(colour_set1_dark)) + 1],
    sapply(scenario_list[.wb_keys], `[[`, "label")
  )
)
# Fallback colour for any scenario not matching _ll or _wb (e.g. "default")
.other_keys <- setdiff(names(scenario_list), c(.ll_keys, .wb_keys))
if (length(.other_keys) > 0) {
  scenario_colour_map <- c(
    scenario_colour_map,
    setNames(rep("grey40", length(.other_keys)),
             sapply(scenario_list[.other_keys], `[[`, "label"))
  )
}
rm(.ll_keys, .wb_keys, .other_keys)

##########



# Parameter specification
##########

include_models <- c("LMM", "bp2", "gp")
plot_models <- c("LMM","bp2", "gp")

# Model colour map: BP/LMM models → colour_set1, GP models → colour_set2
# (colour_set2 indexed from 2 to skip near-white Honeydew at position 1)
.bp_models <- plot_models[grepl("^(LMM|bp)", plot_models)]
.gp_models <- plot_models[grepl("^gp",        plot_models)]
model_colour_map <- c(
  setNames(colour_set1[seq_along(.bp_models)],          .bp_models),
  setNames(colour_set2[seq_along(.gp_models) + 1L],     .gp_models)
)
rm(.bp_models, .gp_models)

bp_param_map <- list(
  LMM = c(
    "beta_long[1]" = "beta_0",
    "beta_long[2]" = "beta_1",
    "beta_long[3]" = "beta_2"
  ),
  bp1 = c(
    "beta_long[1]" = "beta_0",
    "beta_long[2]" = "beta_1",
    "beta_long[3]" = "beta_2",
    "beta_surv[1]" = "log_AF",
    "alpha"  = "alpha_AFT"
  ),
  bp2 = c(
    "beta_long[1]" = "beta_0",
    "beta_long[2]" = "beta_1",
    "beta_long[3]" = "beta_2",
    "beta_surv[1]" = "log_AF",
    "alpha"  = "alpha_AFT"
  ),
  gp = c(
    "beta_long[1]" = "beta_0",
    "beta_long[2]" = "beta_1",
    "beta_long[3]" = "beta_2",
    "gamma[1]" = "log_AF",
    "alpha" = "alpha_AFT"
  )
)

bp_gp_sd_param_map <- list(
  bp2 = c(
    "sigma_long" = "sigma_long",
    "sd_1_long[1]" = "sd_intercept",
    "sd_1_long[2]" = "sd_slope"
  ),
  gp = c(
    "sigma_long" = "sigma_long",
    "sd_1_long[1]" = "sd_intercept",
    "sd_1_long[2]" = "sd_slope"
  )
)

##########


# Helper functions
##########

get_bp_summary <- function(r, model_name) {
  if (!is.null(r[[model_name]]) && !is.null(r[[model_name]]$summary)) {
    return(r[[model_name]]$summary)
  }
  if (identical(model_name, "bp1") && !is.null(r$bp) && !is.null(r$bp$summary)) {
    return(r$bp$summary)
  }
  NULL
}

get_true_value <- function(param_name, true_params) {
  if (param_name %in% names(true_params)) return(as.numeric(true_params[[param_name]]))
  if (param_name == "sigma_long") return(12)
  if (param_name == "sd_intercept") return(15)
  if (param_name == "sd_slope") return(0.20)
  NA_real_
}

# bp2 summaries store q5/q95; gp summaries store q2.5/q97.5.
# These helpers pick whichever column is present so both models work.
get_q_lower <- function(row) {
  if ("q2.5" %in% names(row)) return(as.numeric(row[["q2.5"]][[1]])[1L])
  if ("q5"   %in% names(row)) return(as.numeric(row[["q5"]][[1]])[1L])
  NA_real_
}
get_q_upper <- function(row) {
  if ("q97.5" %in% names(row)) return(as.numeric(row[["q97.5"]][[1]])[1L])
  if ("q95"   %in% names(row)) return(as.numeric(row[["q95"]][[1]])[1L])
  NA_real_
}

get_ic_value <- function(ic_obj, criterion, estimate) {
  if (is.null(ic_obj) || is.null(ic_obj[[criterion]])) return(NA_real_)
  x <- ic_obj[[criterion]]
  if (!is.null(x$estimates) && estimate %in% rownames(x$estimates)) {
    return(as.numeric(x$estimates[estimate, "Estimate"]))
  }
  NA_real_
}

extract_ic_data <- function(model_name, results) {
  map_dfr(seq_along(results), function(i) {
    model_result <- results[[i]][[model_name]]
    if (is.null(model_result) || is.null(model_result$ic)) return(NULL)
    ic <- model_result$ic
    tibble(
      model = model_name,
      replicate = i,
      elpd_loo = get_ic_value(ic, "loo", "elpd_loo"),
      p_loo = get_ic_value(ic, "loo", "p_loo"),
      looic = get_ic_value(ic, "loo", "looic"),
      elpd_waic = get_ic_value(ic, "waic", "elpd_waic"),
      p_waic = get_ic_value(ic, "waic", "p_waic"),
      waic = get_ic_value(ic, "waic", "waic")
    )
  })
}

make_ic_summary <- function(results, models) {
  mean_or_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  sd_or_na <- function(x) if (sum(!is.na(x)) <= 1) NA_real_ else sd(x, na.rm = TRUE)

  ic_long <- bind_rows(lapply(models, extract_ic_data, results = results))

  if (nrow(ic_long) == 0) {
    return(tibble(model = character(), n_ic = integer(),
                  mean_looic = numeric(), sd_looic = numeric(),
                  mean_waic = numeric(), sd_waic = numeric(),
                  mean_elpd_loo = numeric(), mean_elpd_waic = numeric(),
                  mean_p_loo = numeric(), mean_p_waic = numeric()))
  }

  ic_long <- ic_long %>%
    filter(if_any(c(elpd_loo, p_loo, looic, elpd_waic, p_waic, waic), ~ !is.na(.x)))

  if (nrow(ic_long) == 0) {
    return(tibble(model = character(), n_ic = integer(),
                  mean_looic = numeric(), sd_looic = numeric(),
                  mean_waic = numeric(), sd_waic = numeric(),
                  mean_elpd_loo = numeric(), mean_elpd_waic = numeric(),
                  mean_p_loo = numeric(), mean_p_waic = numeric()))
  }

  ic_long %>%
    group_by(model) %>%
    summarise(
      n_ic = n(),
      mean_looic = mean_or_na(looic),
      sd_looic = sd_or_na(looic),
      mean_waic = mean_or_na(waic),
      sd_waic = sd_or_na(waic),
      mean_elpd_loo = mean_or_na(elpd_loo),
      mean_elpd_waic = mean_or_na(elpd_waic),
      mean_p_loo = mean_or_na(p_loo),
      mean_p_waic = mean_or_na(p_waic),
      .groups = "drop"
    ) %>%
    arrange(mean_looic)
}

make_bp_metrics <- function(model_name, results, true_params) {
  param_map <- bp_param_map[[model_name]]

  get_y_long_mean <- function(r) {
    if (is.null(r$Y_long_mean)) return(NA_real_)
    as.numeric(r$Y_long_mean)
  }

  longitudinal_sample_mean <- mean(vapply(results, get_y_long_mean, numeric(1)), na.rm = TRUE)

  map_dfr(names(param_map), function(stan_var) {
    rows <- lapply(results, function(r) {
      summary_tbl <- get_bp_summary(r, model_name)
      if (is.null(summary_tbl)) {
        return(tibble(variable = character(), mean = numeric(), `q2.5` = numeric(), `q97.5` = numeric()))
      }
      summary_tbl %>% filter(variable == stan_var)
    })

    means <- sapply(seq_along(rows), function(i) {
      if (nrow(rows[[i]]) == 0) return(NA_real_)
      raw <- as.numeric(rows[[i]]$mean[[1]])[1L]
      if (stan_var == "beta_long[1]") raw + get_y_long_mean(results[[i]]) else raw
    })

    q2_5s <- sapply(seq_along(rows), function(i) {
      if (nrow(rows[[i]]) == 0) return(NA_real_)
      raw <- get_q_lower(rows[[i]])
      if (stan_var == "beta_long[1]") raw + get_y_long_mean(results[[i]]) else raw
    })

    q97_5s <- sapply(seq_along(rows), function(i) {
      if (nrow(rows[[i]]) == 0) return(NA_real_)
      raw <- get_q_upper(rows[[i]])
      if (stan_var == "beta_long[1]") raw + get_y_long_mean(results[[i]]) else raw
    })

    means  <- unlist(means)
    q2_5s  <- unlist(q2_5s)
    q97_5s <- unlist(q97_5s)

    tv <- get_true_value(param_map[[stan_var]], true_params)
    tibble(
      parameter   = param_map[[stan_var]],
      true        = tv,
      sample_mean = if (stan_var == "beta_long[1]") longitudinal_sample_mean else NA_real_,
      bias        = if (is.na(tv)) NA_real_ else mean(means, na.rm = TRUE) - tv,
      rmse        = if (is.na(tv)) NA_real_ else sqrt(mean((means - tv)^2, na.rm = TRUE)),
      coverage    = if (is.na(tv)) NA_real_ else mean(q2_5s <= tv & q97_5s >= tv, na.rm = TRUE)
    )
  }) %>% mutate(model = model_name)
}

extract_bp_plot_data <- function(model_name, results, true_params) {
  param_map <- bp_param_map[[model_name]]
  map_dfr(seq_along(results), function(i) {
    summary_tbl <- get_bp_summary(results[[i]], model_name)
    if (is.null(summary_tbl)) return(NULL)
    # Safe Y_long_mean for older result objects where the centring constant may
    # be missing. Fall back to 0 rather than silently dropping beta_0 rows.
    .ylm <- { v <- results[[i]]$Y_long_mean; if (is.null(v) || length(v) == 0) 0 else as.numeric(v) }
    map_dfr(names(param_map), function(stan_var) {
      row <- summary_tbl %>% filter(variable == stan_var)
      if (nrow(row) == 0) return(NULL)
      tibble(
        model      = model_name,
        replicate  = i,
        parameter  = param_map[[stan_var]],
        estimate   = if (stan_var == "beta_long[1]") row$mean[[1]] + .ylm else row$mean[[1]],
        sample_mean = if (stan_var == "beta_long[1]") .ylm else NA_real_,
        posterior_sd = if ("sd" %in% names(row)) {
          as.numeric(row$sd[[1]])
        } else {
          (get_q_upper(row) - get_q_lower(row)) / (2 * 1.96)
        },
        ci_lower   = if (stan_var == "beta_long[1]") get_q_lower(row) + .ylm else get_q_lower(row),
        ci_upper   = if (stan_var == "beta_long[1]") get_q_upper(row) + .ylm else get_q_upper(row),
        true       = get_true_value(param_map[[stan_var]], true_params)
      )
    })
  })
}

extract_bp_sd_plot_data <- function(model_name, results, true_params) {
  param_map <- bp_gp_sd_param_map[[model_name]]
  if (is.null(param_map)) return(tibble())

  map_dfr(seq_along(results), function(i) {
    summary_tbl <- get_bp_summary(results[[i]], model_name)
    if (is.null(summary_tbl)) return(NULL)

    map_dfr(names(param_map), function(stan_var) {
      row <- summary_tbl %>% filter(variable == stan_var)
      if (nrow(row) == 0) return(NULL)
      tibble(
        model = model_name,
        replicate = i,
        parameter = param_map[[stan_var]],
        estimate = as.numeric(row$mean[[1]])[1L],
        sample_mean = NA_real_,
        posterior_sd = if ("sd" %in% names(row)) {
          as.numeric(row$sd[[1]])[1L]
        } else {
          (get_q_upper(row) - get_q_lower(row)) / (2 * 1.96)
        },
        ci_lower = get_q_lower(row),
        ci_upper = get_q_upper(row),
        true = get_true_value(param_map[[stan_var]], true_params)
      )
    })
  })
}

extract_mcmc_diagnostics <- function(model_name, results) {
  map_dfr(seq_along(results), function(i) {
    summary_tbl <- get_bp_summary(results[[i]], model_name)
    if (is.null(summary_tbl)) return(NULL)

    diagnostic_cols <- intersect(c("rhat", "ess_bulk", "ess_tail"), names(summary_tbl))
    if (length(diagnostic_cols) == 0) return(NULL)

    summary_tbl %>%
      dplyr::select(variable, any_of(c("rhat", "ess_bulk", "ess_tail"))) %>%
      mutate(
        model = model_name,
        replicate = i,
        .before = 1
      )
  })
}

make_mcmc_diagnostic_table <- function(results, models) {
  diagnostic_long <- bind_rows(lapply(models, extract_mcmc_diagnostics,
                                      results = results))

  if (nrow(diagnostic_long) == 0) {
    return(tibble(
      model = character(), variable = character(), n = integer(),
      mean_rhat = numeric(), max_rhat = numeric(),
      min_bulk_ess = numeric(), median_bulk_ess = numeric(),
      min_tail_ess = numeric(), median_tail_ess = numeric()
    ))
  }

  diagnostic_long %>%
    group_by(model, variable) %>%
    summarise(
      n = n(),
      mean_rhat = mean(rhat, na.rm = TRUE),
      max_rhat = max(rhat, na.rm = TRUE),
      min_bulk_ess = min(ess_bulk, na.rm = TRUE),
      median_bulk_ess = median(ess_bulk, na.rm = TRUE),
      min_tail_ess = min(ess_tail, na.rm = TRUE),
      median_tail_ess = median(ess_tail, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(model, variable)
}

make_mcmc_ess_histogram <- function(diagnostic_long, scenario_label,
                                    ess_col = "ess_bulk",
                                    ess_label = "Bulk ESS") {
  common_parameter_label <- function(model, variable, model_label) {
    mapped <- dplyr::case_when(
      variable == "beta_long[1]" ~ "beta_0",
      variable == "beta_long[2]" ~ "beta_1",
      variable == "beta_long[3]" ~ "beta_2",
      variable == "alpha" ~ "alpha",
      model == "bp2" & variable == "beta_surv[1]" ~ "gamma",
      model == "gp" & variable == "gamma[1]" ~ "gamma",
      variable == "sigma_long" ~ "sigma_long",
      variable == "sd_1_long[1]" ~ "sd_intercept",
      variable == "sd_1_long[2]" ~ "sd_slope",
      TRUE ~ NA_character_
    )

    label <- ifelse(is.na(mapped), variable, mapped)
    list(
      label = ifelse(is.na(mapped), paste(model_label, label, sep = ": "), label),
      group = dplyr::case_when(
        !is.na(mapped) ~ "Shared",
        model == "gp" ~ "GP-specific",
        model == "bp2" ~ "BP-specific",
        TRUE ~ "Model-specific"
      )
    )
  }

  if (!ess_col %in% names(diagnostic_long)) return(NULL)

  ess_plot_data <- diagnostic_long %>%
    mutate(
      model_label = recode(model, bp2 = "BP", gp = "GP"),
      model_label = factor(model_label, levels = c("BP", "GP")),
      ess_value = .data[[ess_col]]
    ) %>%
    filter(is.finite(ess_value))

  if (nrow(ess_plot_data) == 0) return(NULL)

  facet_info <- common_parameter_label(
    ess_plot_data$model,
    ess_plot_data$variable,
    as.character(ess_plot_data$model_label)
  )

  ess_plot_data <- ess_plot_data %>%
    mutate(
      facet_label = facet_info$label,
      facet_group = factor(
        facet_info$group,
        levels = c("Shared", "GP-specific", "BP-specific", "Model-specific")
      )
    ) %>%
    arrange(facet_group, facet_label, model_label, variable) %>%
    mutate(facet_label = factor(facet_label, levels = unique(facet_label)))

  ggplot(ess_plot_data, aes(x = ess_value, fill = model_label, colour = model_label)) +
    geom_histogram(
      bins = 30,
      position = "identity",
      linewidth = 0.2,
      alpha = 0.45
    ) +
    facet_wrap(~ facet_label, scales = "free_y", ncol = 5) +
    scale_fill_manual(values = c(BP = "#D73027", GP = "#2166AC"), name = "Model") +
    scale_colour_manual(values = c(BP = "#D73027", GP = "#2166AC"), name = "Model") +
    labs(
      title = paste0(ess_label, " Distributions: ", scenario_label),
      x = ess_label,
      y = "Replicate count"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 9, colour = "black"),
      strip.text = element_text(size = 8),
      strip.background = element_rect(fill = "grey90", colour = "grey50"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 16)
    )
}

write_mcmc_diagnostics_parameter_tex <- function(all_scenarios_mcmc_diagnostics,
                                                output_file) {
  if (nrow(all_scenarios_mcmc_diagnostics) == 0) return(invisible(NULL))

  format_num <- function(x, digits) {
    ifelse(
      is.na(x) | !is.finite(x),
      "--",
      formatC(x, format = "f", digits = digits)
    )
  }

  format_parameter <- function(x) {
    escaped <- x
    escaped <- gsub("\\", "\\textbackslash{}", escaped, fixed = TRUE)
    escaped <- gsub("_", "\\_", escaped, fixed = TRUE)
    escaped <- gsub("%", "\\%", escaped, fixed = TRUE)
    escaped <- gsub("&", "\\&", escaped, fixed = TRUE)
    escaped <- gsub("#", "\\#", escaped, fixed = TRUE)
    escaped <- gsub("$", "\\$", escaped, fixed = TRUE)
    paste0("\\texttt{", escaped, "}")
  }

  diagnostics_wide <- all_scenarios_mcmc_diagnostics %>%
    mutate(
      scenario_number = as.integer(
        regmatches(scenario, regexpr("[0-9]+", scenario))
      ),
      distribution_key = ifelse(
        grepl("log.logistic|Log-Logistic", scenario, ignore.case = TRUE),
        "ll", "wb"
      ),
      model = recode(model, bp2 = "BP", gp = "GP"),
      model = factor(model, levels = c("BP", "GP")),
      parameter = variable,
      mean_rhat = round(mean_rhat, 3),
      max_rhat = round(max_rhat, 3),
      min_bulk_ess = round(min_bulk_ess, 0),
      median_bulk_ess = round(median_bulk_ess, 0),
      min_tail_ess = round(min_tail_ess, 0),
      median_tail_ess = round(median_tail_ess, 0)
    ) %>%
    dplyr::select(
      scenario_number, distribution_key, model, parameter,
      mean_rhat, max_rhat, min_bulk_ess, median_bulk_ess,
      min_tail_ess, median_tail_ess
    ) %>%
    tidyr::pivot_wider(
      names_from = distribution_key,
      values_from = c(
        mean_rhat, max_rhat, min_bulk_ess, median_bulk_ess,
        min_tail_ess, median_tail_ess
      ),
      names_glue = "{distribution_key}_{.value}"
    ) %>%
    arrange(scenario_number, model, parameter)

  make_split_table <- function(scenarios, caption_scenarios, label_suffix = "") {
    split_data <- diagnostics_wide %>%
      filter(scenario_number %in% scenarios)

    if (nrow(split_data) == 0) return(character())

    split_data <- split_data %>%
      group_by(scenario_number) %>%
      mutate(display_scenario = ifelse(row_number() == 1L,
                                       as.character(scenario_number), "")) %>%
      ungroup()

    row_lines <- pmap_chr(split_data, function(scenario_number, model, parameter,
                                               ll_mean_rhat, wb_mean_rhat,
                                               ll_max_rhat, wb_max_rhat,
                                               ll_min_bulk_ess, wb_min_bulk_ess,
                                               ll_median_bulk_ess, wb_median_bulk_ess,
                                               ll_min_tail_ess, wb_min_tail_ess,
                                               ll_median_tail_ess, wb_median_tail_ess,
                                               display_scenario) {
      row_prefix <- ifelse(as.character(model) == "GP", "\\gprow ", "")
      paste0(
        row_prefix,
        display_scenario, " & ",
        as.character(model), " & ",
        format_parameter(parameter), " & ",
        format_num(ll_mean_rhat, 3), " & ",
        format_num(ll_max_rhat, 3), " & ",
        format_num(ll_min_bulk_ess, 0), " & ",
        format_num(ll_median_bulk_ess, 0), " & ",
        format_num(ll_min_tail_ess, 0), " & ",
        format_num(ll_median_tail_ess, 0), " & ",
        format_num(wb_mean_rhat, 3), " & ",
        format_num(wb_max_rhat, 3), " & ",
        format_num(wb_min_bulk_ess, 0), " & ",
        format_num(wb_median_bulk_ess, 0), " & ",
        format_num(wb_min_tail_ess, 0), " & ",
        format_num(wb_median_tail_ess, 0), " \\\\"
      )
    })

    row_lines_by_scenario <- split(row_lines, split_data$scenario_number)
    row_lines <- unlist(
      Map(function(lines, idx) {
        if (idx < length(row_lines_by_scenario)) {
          c(lines, "\\addlinespace")
        } else {
          lines
        }
      }, row_lines_by_scenario, seq_along(row_lines_by_scenario)),
      use.names = FALSE
    )

    c(
      "\\begin{table}[!htbp]",
      "\\centering",
      "\\fontsize{7}{8}\\selectfont",
      paste0(
        "\\caption{\\large Cross-scenario MCMC diagnostics by parameter for the BP ",
        "and GP sAFT-JM models for ", caption_scenarios, ".}"
      ),
      paste0("\\label{tab:cross-scenario-mcmc-diagnostics-by-parameter",
             label_suffix, "}"),
      "\\definecolor{gpgray}{HTML}{F2F2F2}",
      "\\providecommand{\\gprow}{\\rowcolor{gpgray}}",
      "\\makebox[\\textwidth][c]{\\resizebox{1.018\\textwidth}{!}{%",
      "\\begin{tabular}{lllrrrrrrrrrrrr}",
      "\\toprule",
      paste0(
        "Scen. & Model & Par. & \\multicolumn{6}{c}{Log-logistic} & ",
        "\\multicolumn{6}{c}{Weibull} \\\\"
      ),
      "\\cmidrule(lr){4-9} \\cmidrule(lr){10-15}",
      paste0(
        " & & & Mean R-hat & Max R-hat & Min bulk ESS & Med. bulk ESS & ",
        "Min tail ESS & Med. tail ESS & Mean R-hat & Max R-hat & ",
        "Min bulk ESS & Med. bulk ESS & Min tail ESS & Med. tail ESS \\\\"
      ),
      "\\midrule",
      row_lines,
      "\\bottomrule",
      "\\end{tabular}%",
      "}}",
      "\\vspace{0.5ex}",
      "\\begin{minipage}{0.95\\linewidth}",
      paste0(
        "\\footnotesize Notes: BP and GP denote the Bernstein-polynomial and ",
        "Gaussian-process sAFT joint models. ESS denotes effective sample size. ",
        "Diagnostics are summarised by fitted parameter across simulation ",
        "replicates; the replicate count column is omitted here to keep the ",
        "table focused on R-hat and ESS quantities."
      ),
      "\\end{minipage}",
      "\\end{table}"
    )
  }

  tex_lines <- c(
    make_split_table(1:3, "Scenarios 1--3"),
    "",
    make_split_table(4:5, "Scenarios 4--5", "-4-5")
  )

  writeLines(tex_lines, output_file)
  invisible(output_file)
}

write_mcmc_diagnostic_summary_tex <- function(all_scenarios_mcmc_diagnostics,
                                              output_file) {
  if (nrow(all_scenarios_mcmc_diagnostics) == 0) {
    return(invisible(NULL))
  }

  format_num <- function(x, digits = 3) {
    ifelse(
      is.na(x) | !is.finite(x),
      "--",
      formatC(x, format = "f", digits = digits)
    )
  }

  format_parameter <- function(x) {
    escaped <- x
    escaped <- gsub("\\", "\\textbackslash{}", escaped, fixed = TRUE)
    escaped <- gsub("_", "\\_", escaped, fixed = TRUE)
    escaped <- gsub("%", "\\%", escaped, fixed = TRUE)
    escaped <- gsub("&", "\\&", escaped, fixed = TRUE)
    escaped <- gsub("#", "\\#", escaped, fixed = TRUE)
    escaped <- gsub("$", "\\$", escaped, fixed = TRUE)
    paste0("\\texttt{", escaped, "}")
  }

  summary_data <- all_scenarios_mcmc_diagnostics %>%
    mutate(
      scenario_number = as.integer(
        regmatches(scenario, regexpr("[0-9]+", scenario))
      ),
      distribution = ifelse(
        grepl("log.logistic|Log-Logistic", scenario, ignore.case = TRUE),
        "Log-logistic",
        "Weibull"
      ),
      model = recode(model, bp2 = "BP", gp = "GP"),
      parameter = variable
    ) %>%
    arrange(distribution, scenario_number, model, parameter)

  make_distribution_table <- function(distribution_name, label_suffix) {
    table_data <- summary_data %>%
      filter(distribution == distribution_name) %>%
      transmute(
        Scenario = scenario_number,
        Model = model,
        Parameter = parameter,
        `Mean R-hat` = mean_rhat,
        `Max R-hat` = max_rhat,
        `Min bulk ESS` = min_bulk_ess,
        `Med. bulk ESS` = median_bulk_ess,
        `Min tail ESS` = min_tail_ess,
        `Med. tail ESS` = median_tail_ess
      )

    if (nrow(table_data) == 0) return(character())

    row_lines <- pmap_chr(table_data, function(...) {
      values <- list(...)
      row_prefix <- ifelse(values$Model == "GP", "\\gprow ", "")
      paste0(
        row_prefix,
        values$Scenario, " & ",
        values$Model, " & ",
        format_parameter(values$Parameter), " & ",
        format_num(values[["Mean R-hat"]]), " & ",
        format_num(values[["Max R-hat"]]), " & ",
        format_num(values[["Min bulk ESS"]], 0), " & ",
        format_num(values[["Med. bulk ESS"]], 0), " & ",
        format_num(values[["Min tail ESS"]], 0), " & ",
        format_num(values[["Med. tail ESS"]], 0), " \\\\"
      )
    })

    c(
      "\\begingroup",
      "\\fontsize{7}{8}\\selectfont",
      "\\setlength{\\LTleft}{-0.04\\textwidth}",
      "\\setlength{\\LTright}{-0.04\\textwidth}",
      "\\setlength{\\tabcolsep}{4pt}",
      "\\renewcommand{\\arraystretch}{1.0}",
      "\\definecolor{gpgray}{HTML}{F2F2F2}",
      "\\providecommand{\\gprow}{\\rowcolor{gpgray}}",
      "\\begin{longtable}{lllrrrrrr}",
      paste0(
        "\\caption{\\large MCMC diagnostics by parameter for ",
        distribution_name,
        " simulation scenarios by fitted model.}",
        "\\label{tab:mcmc-diagnostic-summary-", label_suffix, "}\\\\"
      ),
      "\\toprule",
      paste0(
        "Scenario & Model & Parameter & Mean R-hat & Max R-hat & ",
        "Min bulk ESS & Med. bulk ESS & Min tail ESS & Med. tail ESS \\\\"
      ),
      "\\midrule",
      "\\endfirsthead",
      paste0(
        "\\multicolumn{9}{l}{\\textit{Continued from previous page}} \\\\"
      ),
      "\\toprule",
      paste0(
        "Scenario & Model & Parameter & Mean R-hat & Max R-hat & ",
        "Min bulk ESS & Med. bulk ESS & Min tail ESS & Med. tail ESS \\\\"
      ),
      "\\midrule",
      "\\endhead",
      "\\bottomrule",
      "\\endfoot",
      row_lines,
      "\\end{longtable}",
      "\\endgroup"
    )
  }

  tex_lines <- c(
    make_distribution_table("Log-logistic", "loglogistic"),
    "",
    make_distribution_table("Weibull", "weibull")
  )

  writeLines(tex_lines, output_file)
  invisible(output_file)
}

##########


# Per-scenario analysis (accumulates cross-scenario comparison data)
##########

all_metrics_combined <- list()
all_runtime_combined <- list()
all_runtime_plot_data_combined <- list()
all_plot_data_combined <- list()
all_sd_plot_data_combined <- list()
all_ic_combined <- list()
all_mcmc_diagnostic_combined <- list()

if (run_standard_analysis) {
for (sc in scenario_list) {
  scenario_folder <- sc$folder
  scenario_label <- sc$label

  message("\n========================================")
  message("Processing: ", scenario_label)
  message("========================================")

  sim_out_list <- lapply(sc$files, readRDS)
  results <- unlist(lapply(sim_out_list, function(x) x$results), recursive = FALSE)
  true_params <- sim_out_list[[1]]$true_params

  lmm_metrics <- tibble(parameter = character(), true = numeric(), bias = numeric(), rmse = numeric())

  # Metrics
  bp1_metrics <- make_bp_metrics("bp1", results, true_params)
  bp2_metrics <- make_bp_metrics("bp2", results, true_params)
  gp_metrics <- make_bp_metrics("gp",  results, true_params)

  model_metrics <- list(
    LMM = lmm_metrics %>% mutate(model = "LMM", coverage = NA_real_),
    bp1 = bp1_metrics,
    bp2 = bp2_metrics,
    gp = gp_metrics
  )

  message("--- Combined metrics | ", scenario_label, " ---")
  combined_metrics <- bind_rows(model_metrics[include_models]) %>%
    dplyr::select(model, parameter, true, sample_mean, bias, rmse, coverage)
  print(combined_metrics %>% mutate(across(where(is.numeric), ~ round(.x, 4))), n = Inf)

  all_metrics_combined[[scenario_folder]] <- combined_metrics %>%
    mutate(scenario = scenario_label)

  # Plot data
  bp_plot_data <- bind_rows(lapply(plot_models, extract_bp_plot_data,
                                   results = results, true_params = true_params))
  bp_plot_data$model <- factor(bp_plot_data$model, levels = plot_models)

  bp_sd_plot_data <- bind_rows(lapply(c("bp2", "gp"), extract_bp_sd_plot_data,
                                      results = results, true_params = true_params))
  if (nrow(bp_sd_plot_data) > 0) {
    bp_sd_plot_data$model <- factor(bp_sd_plot_data$model, levels = plot_models)
  }

  aggregate_estimates_table <- bp_plot_data %>%
    group_by(model, parameter) %>%
    summarise(
      n_sim = sum(!is.na(estimate)),
      estimate_mean = mean(estimate, na.rm = TRUE),
      sample_mean = ifelse(any(!is.na(sample_mean)), mean(sample_mean, na.rm = TRUE), NA_real_),
      true = mean(true, na.rm = TRUE),
      std_error = ifelse(n_sim > 1, sd(estimate, na.rm = TRUE) / sqrt(n_sim), NA_real_),
      lower_95_ci = estimate_mean - 1.96 * std_error,
      upper_95_ci = estimate_mean + 1.96 * std_error,
      .groups = "drop"
    ) %>%
    rename(estimate = estimate_mean) %>%
    arrange(model, parameter)

  # Clean results table
  clean_results_table <- aggregate_estimates_table %>%
    dplyr::select(model, parameter, n_sim, estimate, true, std_error, lower_95_ci, upper_95_ci) %>%
    mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
    arrange(model, parameter)

  message("--- Clean results table | ", scenario_label, " ---")
  print(clean_results_table, n = Inf)

  if (nrow(bp_sd_plot_data) > 0) {
    aggregate_sd_estimates_table <- bp_sd_plot_data %>%
      group_by(model, parameter) %>%
      summarise(
        n_sim = sum(!is.na(estimate)),
        estimate_mean = mean(estimate, na.rm = TRUE),
        sample_mean = NA_real_,
        true = if (all(is.na(true))) NA_real_ else mean(true, na.rm = TRUE),
        std_error = ifelse(n_sim > 1, sd(estimate, na.rm = TRUE) / sqrt(n_sim), NA_real_),
        lower_95_ci = estimate_mean - 1.96 * std_error,
        upper_95_ci = estimate_mean + 1.96 * std_error,
        .groups = "drop"
      ) %>%
      rename(estimate = estimate_mean) %>%
      arrange(model, parameter)

    full_results_with_sd_table <- bind_rows(
      aggregate_estimates_table,
      aggregate_sd_estimates_table
    ) %>%
      dplyr::select(model, parameter, n_sim, estimate, true, std_error,
                    lower_95_ci, upper_95_ci) %>%
      mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
      arrange(model, factor(
        parameter,
        levels = c("beta_0", "beta_1", "beta_2", "alpha_AFT", "log_AF",
                   "sigma_long", "sd_intercept", "sd_slope")
      ))

    message("--- Full results table with BP/GP standard deviation terms | ",
            scenario_label, " ---")
    print(full_results_with_sd_table, n = Inf)
  } else {
    full_results_with_sd_table <- tibble()
  }

  bias_rmse_coverage_table <- aggregate_estimates_table %>%
    dplyr::select(model, parameter, estimate, lower_95_ci, upper_95_ci) %>%
    left_join(
      combined_metrics %>% dplyr::select(model, parameter, bias, rmse, coverage),
      by = c("model", "parameter")
    ) %>%
    mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
    arrange(model, parameter)

  message("--- Bias/RMSE/Coverage table | ", scenario_label, " ---")
  print(bias_rmse_coverage_table, n = Inf)

  out_prefix <- if (scenario_folder == "default") "" else paste0(scenario_folder, "_")

  # Information criteria
  ic_summary <- make_ic_summary(results, include_models) %>%
    mutate(across(where(is.numeric), ~ round(.x, 4)))

  if (nrow(ic_summary) > 0) {
    message("--- Information criteria | ", scenario_label, " ---")
    print(ic_summary, n = Inf)

    paste0("<!-- Scenario: ", scenario_label, " -->\n") %>%
      cat(file = here("Simulation-Workflow", "Results",
                      paste0(out_prefix, "information_criteria_table.md")))
    kable(ic_summary, format = "markdown") %>%
      cat(file = here("Simulation-Workflow", "Results",
                      paste0(out_prefix, "information_criteria_table.md")), append = TRUE)

    all_ic_combined[[scenario_folder]] <- ic_summary %>%
      mutate(scenario = scenario_label)
  }

  mcmc_diagnostic_table <- make_mcmc_diagnostic_table(results, c("bp2", "gp")) %>%
    mutate(across(where(is.numeric), ~ round(.x, 4)))
  mcmc_diagnostic_long <- bind_rows(lapply(c("bp2", "gp"),
                                           extract_mcmc_diagnostics,
                                           results = results))

  if (nrow(mcmc_diagnostic_table) > 0) {
    message("--- MCMC diagnostics | ", scenario_label, " ---")
    print(mcmc_diagnostic_table, n = Inf)

    paste0("<!-- Scenario: ", scenario_label, " -->\n") %>%
      cat(file = here("Simulation-Workflow", "Results",
                      paste0(out_prefix, "mcmc_diagnostic_table.md")))
    kable(mcmc_diagnostic_table, format = "markdown") %>%
      cat(file = here("Simulation-Workflow", "Results",
                      paste0(out_prefix, "mcmc_diagnostic_table.md")), append = TRUE)

    all_mcmc_diagnostic_combined[[scenario_folder]] <- mcmc_diagnostic_table %>%
      mutate(scenario = scenario_label)
  }

  if (nrow(mcmc_diagnostic_long) > 0) {
    mcmc_bulk_ess_histogram <- make_mcmc_ess_histogram(
      mcmc_diagnostic_long, scenario_label,
      ess_col = "ess_bulk",
      ess_label = "Bulk ESS"
    )
    mcmc_tail_ess_histogram <- make_mcmc_ess_histogram(
      mcmc_diagnostic_long, scenario_label,
      ess_col = "ess_tail",
      ess_label = "Tail ESS"
    )
    mcmc_rhat_histogram <- make_mcmc_ess_histogram(
      mcmc_diagnostic_long, scenario_label,
      ess_col = "rhat",
      ess_label = "R-hat"
    )

    if (!is.null(mcmc_bulk_ess_histogram)) {
      ggsave(
        filename = here("Writeup", "figures",
                        paste0(out_prefix, "mcmc_bulk_ess_histograms.pdf")),
        plot = mcmc_bulk_ess_histogram,
        width = 15,
        height = 11
      )
    }
    if (!is.null(mcmc_tail_ess_histogram)) {
      ggsave(
        filename = here("Writeup", "figures",
                        paste0(out_prefix, "mcmc_tail_ess_histograms.pdf")),
        plot = mcmc_tail_ess_histogram,
        width = 15,
        height = 11
      )
    }
    if (!is.null(mcmc_rhat_histogram)) {
      ggsave(
        filename = here("Writeup", "figures",
                        paste0(out_prefix, "mcmc_rhat_histograms.pdf")),
        plot = mcmc_rhat_histogram,
        width = 15,
        height = 11
      )
    }
  }

  paste0("<!-- Scenario: ", scenario_label, " -->\n") %>%
    cat(file = here("Simulation-Workflow", "Results",
                    paste0(out_prefix, "clean_results_table.md")))
  kable(clean_results_table, format = "markdown") %>%
    cat(file = here("Simulation-Workflow", "Results",
                    paste0(out_prefix, "clean_results_table.md")), append = TRUE)

  paste0("<!-- Scenario: ", scenario_label, " -->\n") %>%
    cat(file = here("Simulation-Workflow", "Results",
                    paste0(out_prefix, "bias_rmse_coverage_table.md")))
  kable(bias_rmse_coverage_table, format = "markdown") %>%
    cat(file = here("Simulation-Workflow", "Results",
                    paste0(out_prefix, "bias_rmse_coverage_table.md")), append = TRUE)

  if (nrow(full_results_with_sd_table) > 0) {
    paste0("<!-- Scenario: ", scenario_label, " -->\n") %>%
      cat(file = here("Simulation-Workflow", "Results",
                      paste0(out_prefix, "full_results_with_sd_table.md")))
    kable(full_results_with_sd_table, format = "markdown") %>%
      cat(file = here("Simulation-Workflow", "Results",
                      paste0(out_prefix, "full_results_with_sd_table.md")),
          append = TRUE)
  }

  all_plot_data_combined[[scenario_folder]] <- bp_plot_data %>%
    mutate(scenario = scenario_label)

  if (nrow(bp_sd_plot_data) > 0) {
    all_sd_plot_data_combined[[scenario_folder]] <- bp_sd_plot_data %>%
      mutate(scenario = scenario_label)
  }

  # Runtime
  runtime_models <- intersect(include_models, c("LMM", "bp1", "bp2", "gp", "jm2"))
  runtime_long <- map_dfr(runtime_models, function(m) {
    rts <- vapply(results, function(r) {
      rt <- r[[m]]$runtime
      if (is.null(rt)) NA_real_ else as.numeric(rt)
    }, numeric(1))
    tibble(model = m, runtime_s = rts)
  }) %>% filter(!is.na(runtime_s))

  runtime_summary <- runtime_long %>%
    group_by(model) %>%
    summarise(n = sum(!is.na(runtime_s)),
              mean_s = mean(runtime_s, na.rm = TRUE),
              median_s = median(runtime_s, na.rm = TRUE),
              sd_s = sd(runtime_s, na.rm = TRUE),
              min_s = min(runtime_s, na.rm = TRUE),
              max_s = max(runtime_s, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
    arrange(mean_s)

  message("--- Runtime summary | ", scenario_label, " ---")
  print(runtime_summary)

  paste0("<!-- Scenario: ", scenario_label, " -->\n") %>%
    cat(file = here("Simulation-Workflow", "Results",
                    paste0(out_prefix, "runtime_summary.md")))
  kable(runtime_summary, format = "markdown") %>%
    cat(file = here("Simulation-Workflow", "Results",
                    paste0(out_prefix, "runtime_summary.md")), append = TRUE)

  all_runtime_combined[[scenario_folder]] <- runtime_summary %>%
    mutate(scenario = scenario_label)

  all_runtime_plot_data_combined[[scenario_folder]] <- runtime_long %>%
    mutate(scenario = scenario_label)

}
}

##########


# PH/LMM comparison outputs
##########

if (length(ph_scenario_list) > 0) {
  message("\n========================================")
  message("PH/LMM cross-scenario comparisons")
  message("========================================")

  ph_lmm_metrics_all <- list()
  ph_lmm_plot_all <- list()
  ph_runtime_all <- list()
  ph_availability_all <- list()

  for (base_sf in names(ph_scenario_list)) {
    sc <- ph_scenario_list[[base_sf]]
    sim_out_list <- lapply(sc$files, readRDS)
    results <- unlist(lapply(sim_out_list, function(x) x$results), recursive = FALSE)
    true_params <- sim_out_list[[1]]$true_params

    scenario_number <- as.integer(sub("^s([0-9]+).*$", "\\1", base_sf))
    distribution <- if (grepl("_ll(?:_|$)", base_sf, perl = TRUE)) {
      "Log-Logistic"
    } else {
      "Weibull"
    }

    lmm_metrics <- make_bp_metrics("LMM", results, true_params)
    lmm_estimates <- extract_bp_plot_data("LMM", results, true_params)

    ph_lmm_metrics_all[[base_sf]] <- lmm_metrics %>%
      mutate(
        scenario = sc$label,
        scenario_short = paste0("S", scenario_number, " ",
                                if_else(distribution == "Log-Logistic", "LL", "WB")),
        scenario_number = scenario_number,
        distribution = distribution
      )

    ph_lmm_plot_all[[base_sf]] <- lmm_estimates %>%
      mutate(
        scenario = sc$label,
        scenario_short = paste0("S", scenario_number, " ",
                                if_else(distribution == "Log-Logistic", "LL", "WB")),
        scenario_number = scenario_number,
        distribution = distribution
      )

    ph_runtime_all[[base_sf]] <- map_dfr(c("LMM", "jm2"), function(model_name) {
      tibble(
        model = model_name,
        runtime_s = vapply(results, function(r) {
          model_result <- r[[model_name]]
          if (is.null(model_result) || is.null(model_result$runtime)) {
            return(NA_real_)
          }
          as.numeric(model_result$runtime)
        }, numeric(1))
      )
    }) %>%
      filter(is.finite(runtime_s)) %>%
      mutate(
        scenario = sc$label,
        scenario_short = paste0("S", scenario_number, " ",
                                if_else(distribution == "Log-Logistic", "LL", "WB")),
        scenario_number = scenario_number,
        distribution = distribution
      )

    ph_availability_all[[base_sf]] <- tibble(
      scenario = sc$label,
      scenario_number = scenario_number,
      distribution = distribution,
      n_replicates = length(results),
      n_lmm_summaries = sum(vapply(
        results,
        function(r) !is.null(r$LMM$summary),
        logical(1)
      )),
      n_jm2_summaries = sum(vapply(
        results,
        function(r) !is.null(r$jm2$summary),
        logical(1)
      )),
      n_jm2_runtimes = sum(vapply(
        results,
        function(r) !is.null(r$jm2$runtime) && is.finite(r$jm2$runtime),
        logical(1)
      ))
    )
  }

  ph_lmm_metrics <- bind_rows(ph_lmm_metrics_all)
  ph_lmm_plot_data <- bind_rows(ph_lmm_plot_all)
  ph_runtime_data <- bind_rows(ph_runtime_all)
  ph_summary_availability <- bind_rows(ph_availability_all) %>%
    arrange(scenario_number, distribution)

  ph_lmm_performance_summary <- ph_lmm_plot_data %>%
    group_by(scenario_number, distribution, scenario, parameter) %>%
    summarise(
      n = sum(is.finite(estimate)),
      true = first(true),
      mean_estimate = mean(estimate, na.rm = TRUE),
      empirical_sd = sd(estimate, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(
      ph_lmm_metrics %>%
        dplyr::select(scenario_number, distribution, parameter,
                      bias, rmse, coverage),
      by = c("scenario_number", "distribution", "parameter")
    ) %>%
    arrange(scenario_number, distribution, parameter) %>%
    mutate(across(where(is.numeric), ~ round(.x, 4)))

  ph_runtime_summary <- ph_runtime_data %>%
    mutate(model = recode(model, LMM = "LMM", jm2 = "PH joint model")) %>%
    group_by(scenario_number, distribution, scenario, model) %>%
    summarise(
      n = sum(is.finite(runtime_s)),
      mean_seconds = mean(runtime_s, na.rm = TRUE),
      median_seconds = median(runtime_s, na.rm = TRUE),
      q25_seconds = quantile(runtime_s, 0.25, na.rm = TRUE),
      q75_seconds = quantile(runtime_s, 0.75, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(scenario_number, distribution, model) %>%
    mutate(across(where(is.numeric), ~ round(.x, 2)))

  message("--- PH-file LMM performance summary ---")
  print(ph_lmm_performance_summary, n = Inf)
  message("--- LMM versus PH runtime summary ---")
  print(ph_runtime_summary, n = Inf)
  message("--- PH summary availability ---")
  print(ph_summary_availability, n = Inf)

  writeLines(
    kable(ph_lmm_performance_summary, format = "markdown"),
    here("Simulation-Workflow", "Results",
         "cross_scenario_ph_lmm_performance_summary.md")
  )
  writeLines(
    kable(ph_runtime_summary, format = "markdown"),
    here("Simulation-Workflow", "Results",
         "cross_scenario_ph_lmm_runtime_summary.md")
  )
  writeLines(
    kable(ph_summary_availability, format = "markdown"),
    here("Simulation-Workflow", "Results",
         "cross_scenario_ph_summary_availability.md")
  )

  if (all(ph_summary_availability$n_jm2_summaries == 0)) {
    warning(
      "The _ph result files contain jm2 runtimes but no stored jm2 parameter summaries; ",
      "PH coefficient-performance comparisons cannot be produced from these RDS files."
    )
  }
}

##########


# Cross-scenario comparison
##########

if (length(all_metrics_combined) > 1) {
  message("\n========================================")
  message("Cross-scenario comparison")
  message("========================================")

  all_scenarios_metrics <- bind_rows(all_metrics_combined) %>%
    dplyr::select(scenario, model, parameter, true, bias, rmse, coverage) %>%
    mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
    arrange(scenario, model, parameter)

  print(all_scenarios_metrics, n = Inf)

  kable(all_scenarios_metrics, format = "markdown") %>%
    cat(file = here("Simulation-Workflow", "Results", "all_scenarios_metrics.md"))

  cross_scenario_runtime_summary <- bind_rows(all_runtime_combined) %>%
    filter(model %in% c("bp2", "gp")) %>%
    transmute(
      scenario,
      model = recode(model, bp2 = "BP", gp = "GP"),
      n,
      mean_seconds = round(mean_s, 2),
      median_seconds = round(median_s, 2),
      sd_seconds = round(sd_s, 2),
      min_seconds = round(min_s, 2),
      max_seconds = round(max_s, 2)
    )

  if (exists("ph_runtime_data") && nrow(ph_runtime_data) > 0) {
    ph_lmm_runtime_summary <- ph_runtime_data %>%
      filter(model %in% c("LMM", "jm2")) %>%
      mutate(model = recode(model, LMM = "LMM", jm2 = "PH")) %>%
      group_by(scenario_number, distribution, scenario, model) %>%
      summarise(
        n = sum(is.finite(runtime_s)),
        mean_seconds = round(mean(runtime_s, na.rm = TRUE), 2),
        median_seconds = round(median(runtime_s, na.rm = TRUE), 2),
        sd_seconds = round(sd(runtime_s, na.rm = TRUE), 2),
        min_seconds = round(min(runtime_s, na.rm = TRUE), 2),
        max_seconds = round(max(runtime_s, na.rm = TRUE), 2),
        .groups = "drop"
      ) %>%
      dplyr::select(scenario, model, n, mean_seconds, median_seconds,
                    sd_seconds, min_seconds, max_seconds)

    cross_scenario_runtime_summary <- bind_rows(
      ph_lmm_runtime_summary,
      cross_scenario_runtime_summary
    )
  }

  cross_scenario_runtime_summary <- cross_scenario_runtime_summary %>%
    mutate(
      scenario_number = as.integer(regmatches(scenario, regexpr("[0-9]+", scenario))),
      distribution = ifelse(
        grepl("log.logistic|Log-Logistic", scenario, ignore.case = TRUE),
        "Log-logistic", "Weibull"
      ),
      model = factor(model, levels = c("LMM", "PH", "BP", "GP"))
    ) %>%
    arrange(scenario_number, distribution, model) %>%
    dplyr::select(-scenario_number, -distribution)

  message("--- Cross-scenario runtime summary (seconds) ---")
  print(cross_scenario_runtime_summary, n = Inf)

  kable(cross_scenario_runtime_summary, format = "markdown") %>%
    cat(file = here(
      "Simulation-Workflow", "Results", "cross_scenario_runtime_summary.md"
    ))

  all_runtime_plot_data <- bind_rows(all_runtime_plot_data_combined)
  all_runtime_plot_data$model <- factor(all_runtime_plot_data$model, levels = plot_models)
  all_runtime_plot_data$scenario <- factor(
    all_runtime_plot_data$scenario,
    levels = sapply(scenario_list, `[[`, "label")
  )
  all_runtime_plot_data <- all_runtime_plot_data %>%
    mutate(
      model_type = ifelse(grepl("^gp", model, ignore.case = TRUE), "GP", "BP"),
      dgm_type = ifelse(
        grepl("log.logistic|\\bll\\b", scenario, ignore.case = TRUE),
        "loglogistic", "weibull"
      ),
      scenario_num = as.integer(regmatches(scenario, regexpr("[0-9]+", scenario))),
      fill_group = paste(model_type, dgm_type, paste0("S", scenario_num), sep = "_")
    )

  short_scenario_labels <- function(x) {
    scenario_num <- regmatches(x, regexpr("[0-9]+", x))
    dgm <- ifelse(
      grepl("log.logistic|Log-Logistic|\\bll\\b", x, ignore.case = TRUE),
      "LL", "WB"
    )
    paste0("S", scenario_num, ": ", dgm)
  }

  cross_scenario_runtime_boxplot <- ggplot(
    all_runtime_plot_data,
    aes(x = scenario, y = runtime_s / 60, fill = fill_group)
  ) +
    geom_boxplot(
      alpha = 0.72, width = 0.65, outlier.size = 0.6,
      position = position_dodge(width = 0.75)
    ) +
    scale_fill_manual(values = boxplot_cols, name = "Model / DGM / Scenario") +
    scale_x_discrete(labels = short_scenario_labels) +
    labs(
      title = "Runtime Across Scenarios",
      x = "Scenario", y = "Runtime (minutes)"
    ) +
    theme_bw() +
    theme(
      legend.position = "right",
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14, colour = "black"),
      axis.text.x = element_text(size = 14, colour = "black"),
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 14),
      strip.text = element_text(size = 16),
      panel.grid.major.x = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 18),
      plot.subtitle = element_text(hjust = 0.5, size = 15)
    )

  ggsave(
    filename = here("Writeup", "figures", "cross_scenario_runtime_boxplot.pdf"),
    plot = cross_scenario_runtime_boxplot, width = 14, height = 8
  )

  if (length(all_ic_combined) > 0) {
    all_scenarios_ic <- bind_rows(all_ic_combined) %>%
      dplyr::select(scenario, model, n_ic, mean_looic, sd_looic,
                    mean_waic, sd_waic, mean_elpd_loo,
                    mean_elpd_waic, mean_p_loo, mean_p_waic) %>%
      arrange(scenario, mean_looic)

    print(all_scenarios_ic, n = Inf)

    kable(all_scenarios_ic, format = "markdown") %>%
      cat(file = here("Simulation-Workflow", "Results", "all_scenarios_information_criteria.md"))

    cross_scenario_loo_waic_summary <- all_scenarios_ic %>%
      transmute(
        Scenario = as.integer(regmatches(scenario, regexpr("[0-9]+", scenario))),
        DGM = ifelse(
          grepl("log.logistic", scenario, ignore.case = TRUE),
          "Log-logistic", "Weibull"
        ),
        Model = recode(model, bp2 = "BP", gp = "GP"),
        n = n_ic,
        `Mean LOOIC` = round(mean_looic, 2),
        `SD LOOIC` = round(sd_looic, 2),
        `Mean WAIC` = round(mean_waic, 2),
        `SD WAIC` = round(sd_waic, 2)
      ) %>%
      arrange(Scenario, DGM, Model)

    message("--- Cross-scenario LOO-CV and WAIC summary ---")
    print(cross_scenario_loo_waic_summary, n = Inf)

    kable(cross_scenario_loo_waic_summary, format = "markdown") %>%
      cat(file = here(
        "Simulation-Workflow", "Results", "cross_scenario_loo_waic_summary.md"
      ))

    cross_scenario_ic_plot_data <- all_scenarios_ic %>%
      transmute(
        scenario = factor(
          scenario,
          levels = sapply(scenario_list, `[[`, "label")
        ),
        model = factor(model, levels = c("bp2", "gp")),
        LOOIC = mean_looic,
        LOOIC_sd = sd_looic,
        WAIC = mean_waic,
        WAIC_sd = sd_waic
      ) %>%
      pivot_longer(
        cols = c(LOOIC, WAIC),
        names_to = "criterion",
        values_to = "estimate"
      ) %>%
      mutate(
        sd = ifelse(criterion == "LOOIC", LOOIC_sd, WAIC_sd),
        model_label = recode(as.character(model), bp2 = "BP", gp = "GP"),
        criterion = factor(
          criterion,
          levels = c("LOOIC", "WAIC"),
          labels = c("LOO-CV (LOOIC)", "WAIC")
        )
      )

    cross_scenario_loo_waic_plot <- ggplot(
      cross_scenario_ic_plot_data,
      aes(x = scenario, y = estimate, colour = model_label, group = model_label)
    ) +
      geom_errorbar(
        aes(ymin = estimate - sd, ymax = estimate + sd),
        width = 0.15,
        linewidth = 0.45,
        position = position_dodge(width = 0.45)
      ) +
      geom_point(
        size = 2.2,
        position = position_dodge(width = 0.45)
      ) +
      facet_wrap(~ criterion, scales = "free_y", ncol = 1) +
      scale_colour_manual(
        values = c(BP = "#D73027", GP = "#2166AC"),
        name = "Model"
      ) +
      scale_x_discrete(
        labels = short_scenario_labels
      ) +
      labs(
        title = "LOO-CV and WAIC Across Scenarios",
        subtitle = "Points show means; error bars show +/- 1 SD across simulation replicates",
        x = "Scenario",
        y = "Information criterion (lower is better)"
      ) +
      theme_bw() +
      theme(
        legend.position = "right",
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14, colour = "black"),
        axis.text.x = element_text(size = 14, colour = "black"),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 14),
        strip.background = element_rect(fill = "grey90", colour = "grey40"),
        strip.text = element_text(size = 16),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 15)
      )

    ggsave(
      filename = here("Writeup", "figures", "cross_scenario_loo_waic_plot.pdf"),
      plot = cross_scenario_loo_waic_plot,
      width = 14,
      height = 10
    )
  }

  if (length(all_mcmc_diagnostic_combined) > 0) {
    all_scenarios_mcmc_diagnostics <- bind_rows(all_mcmc_diagnostic_combined) %>%
      dplyr::select(scenario, model, variable, n, mean_rhat, max_rhat,
                    min_bulk_ess, median_bulk_ess,
                    min_tail_ess, median_tail_ess) %>%
      arrange(scenario, model, variable)

    message("--- Cross-scenario MCMC diagnostics by parameter ---")
    print(all_scenarios_mcmc_diagnostics, n = Inf)

    kable(all_scenarios_mcmc_diagnostics, format = "markdown") %>%
      cat(file = here("Simulation-Workflow", "Results",
                      "all_scenarios_mcmc_diagnostics.md"))

    write_mcmc_diagnostics_parameter_tex(
      all_scenarios_mcmc_diagnostics,
      here("Writeup", "tables",
           "cross_scenario_mcmc_diagnostics_by_parameter.tex")
    )

    cross_scenario_mcmc_diagnostic_summary <- all_scenarios_mcmc_diagnostics %>%
      group_by(scenario, model) %>%
      summarise(
        n_parameters = n_distinct(variable),
        max_rhat = max(max_rhat, na.rm = TRUE),
        min_bulk_ess = min(min_bulk_ess, na.rm = TRUE),
        min_tail_ess = min(min_tail_ess, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        scenario_number = as.integer(regmatches(scenario, regexpr("[0-9]+", scenario))),
        distribution = ifelse(
          grepl("log.logistic|Log-Logistic", scenario, ignore.case = TRUE),
          "Log-logistic", "Weibull"
        ),
        model = recode(model, bp2 = "BP", gp = "GP"),
        across(c(max_rhat, min_bulk_ess, min_tail_ess), ~ round(.x, 3))
      ) %>%
      arrange(scenario_number, distribution, model) %>%
      dplyr::select(
        scenario, model, n_parameters, max_rhat, min_bulk_ess, min_tail_ess
      )

    message("--- Cross-scenario MCMC diagnostic summary ---")
    print(cross_scenario_mcmc_diagnostic_summary, n = Inf)

    kable(cross_scenario_mcmc_diagnostic_summary, format = "markdown") %>%
      cat(file = here("Simulation-Workflow", "Results",
                      "cross_scenario_mcmc_diagnostic_summary.md"))

    write_mcmc_diagnostic_summary_tex(
      all_scenarios_mcmc_diagnostics,
      here("Writeup", "tables", "mcmc_diagnostic_summary.tex")
    )
  }

  # Cross-scenario plots
  all_plot_data <- bind_rows(all_plot_data_combined)
  all_plot_data$model    <- factor(all_plot_data$model,    levels = plot_models)
  all_plot_data$scenario <- factor(all_plot_data$scenario, levels = sapply(scenario_list, `[[`, "label"))
  all_plot_data$parameter <- factor(all_plot_data$parameter, levels = primary_parameter_order)

  main_table_plot_data <- all_plot_data
  if (exists("ph_lmm_plot_data") && nrow(ph_lmm_plot_data) > 0) {
    main_table_plot_data <- bind_rows(
      all_plot_data,
      ph_lmm_plot_data %>% dplyr::select(any_of(names(all_plot_data)))
    )
  }

  table2_1_expanded_saft_lmm_distribution_estimates <- main_table_plot_data %>%
    mutate(
      scenario_num = as.integer(regmatches(scenario, regexpr("[0-9]+", scenario))),
      distribution = ifelse(
        grepl("log.logistic", scenario, ignore.case = TRUE),
        "Log-logistic", "Weibull"
      ),
      model_label = recode(as.character(model), LMM = "LMM", bp2 = "BP", gp = "GP"),
      parameter = factor(
        parameter,
        levels = primary_parameter_order
      )
    ) %>%
    filter(
      model_label %in% c("LMM", "BP", "GP"),
      !is.na(parameter)
    ) %>%
    group_by(scenario_num, parameter, distribution, model_label) %>%
    summarise(
      true = first(true[!is.na(true)]),
      estimate = mean(estimate, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(distribution_model = paste(distribution, model_label)) %>%
    dplyr::select(scenario_num, parameter, true, distribution_model, estimate) %>%
    pivot_wider(names_from = distribution_model, values_from = estimate) %>%
    arrange(scenario_num, parameter) %>%
    transmute(
      Scenario = scenario_num,
      Parameter = recode(
        as.character(parameter),
        beta_0 = "beta_0",
        beta_1 = "beta_1",
        beta_2 = "beta_2",
        alpha_AFT = "alpha",
        log_AF = "gamma"
      ),
      True = round(true, 3),
      `Log-logistic LMM` = round(`Log-logistic LMM`, 3),
      `Log-logistic BP` = round(`Log-logistic BP`, 3),
      `Log-logistic GP` = round(`Log-logistic GP`, 3),
      `Weibull LMM` = round(`Weibull LMM`, 3),
      `Weibull BP` = round(`Weibull BP`, 3),
      `Weibull GP` = round(`Weibull GP`, 3)
    )

  message("--- Table 2.1: estimates by data-generation distribution ---")
  print(table2_1_expanded_saft_lmm_distribution_estimates, n = Inf)

  format_main_result <- function(x) {
    ifelse(is.na(x), "--", formatC(x, format = "f", digits = 3))
  }

  main_results_rows <- table2_1_expanded_saft_lmm_distribution_estimates %>%
    mutate(
      Scenario = if_else(duplicated(Scenario), "", as.character(Scenario)),
      Parameter = recode(
        Parameter,
        beta_0 = "$\\beta_0$",
        beta_1 = "$\\beta_1$",
        beta_2 = "$\\beta_2$",
        alpha = "$\\alpha$",
        gamma = "$\\gamma$"
      ),
      across(-c(Scenario, Parameter), format_main_result)
    )

  main_results_body <- pmap_chr(main_results_rows, function(...) {
    values <- unname(unlist(list(...)))
    paste0(paste(values, collapse = " & "), " \\\\")
  })
  scenario_end_rows <- which(table2_1_expanded_saft_lmm_distribution_estimates$Parameter == "gamma")
  main_results_body[scenario_end_rows[-length(scenario_end_rows)]] <- paste0(
    main_results_body[scenario_end_rows[-length(scenario_end_rows)]],
    "\n\\addlinespace"
  )

  main_results_tex <- c(
    "\\begin{table}[!htbp]",
    "\\centering",
    "\\scriptsize",
    "\\caption{LMM, BP, and GP estimates by data-generation distribution.}",
    "\\label{tab:saft-lmm-distribution-estimates}",
    "\\begin{tabular}{llrrrrrrr}",
    "\\toprule",
    "Scen. & Par. & True & \\multicolumn{3}{c}{Log-logistic} & \\multicolumn{3}{c}{Weibull} \\\\",
    "\\cmidrule(lr){4-6} \\cmidrule(lr){7-9}",
    " & & & LMM & BP & GP & LMM & BP & GP \\\\",
    "\\midrule",
    main_results_body,
    "\\bottomrule",
    "\\end{tabular}",
    "\\vspace{0.5ex}",
    "\\begin{minipage}{0.97\\linewidth}",
    "\\scriptsize Notes: Log-logistic = log-logistic(1.20, 23.0); Weibull = Weibull(0.90, 38.0). LMM denotes the standalone longitudinal mixed model; BP and GP denote the Bernstein-polynomial and Gaussian-process sAFT joint models. Dashes indicate parameters not estimated by the LMM. Entries are means across 1,000 simulation replicates per scenario and fitted model.",
    "\\end{minipage}",
    "\\end{table}"
  )
  writeLines(
    main_results_tex,
    here("Writeup", "tables", "main_results_table.tex")
  )

  if (length(all_sd_plot_data_combined) > 0) {
    all_sd_plot_data <- bind_rows(all_sd_plot_data_combined)
    all_sd_plot_data$model <- factor(all_sd_plot_data$model, levels = plot_models)
    all_sd_plot_data$scenario <- factor(
      all_sd_plot_data$scenario,
      levels = sapply(scenario_list, `[[`, "label")
    )

    full_main_table_plot_data <- bind_rows(main_table_plot_data, all_sd_plot_data)

    main_results_with_sd_table <- full_main_table_plot_data %>%
      mutate(
        scenario_num = as.integer(regmatches(scenario, regexpr("[0-9]+", scenario))),
        distribution = ifelse(
          grepl("log.logistic", scenario, ignore.case = TRUE),
          "Log-logistic", "Weibull"
        ),
        model_label = recode(as.character(model), LMM = "LMM", bp2 = "BP", gp = "GP"),
        parameter = factor(
          parameter,
          levels = c("beta_0", "beta_1", "beta_2", "alpha_AFT", "log_AF",
                     "sigma_long", "sd_intercept", "sd_slope")
        )
      ) %>%
      filter(
        model_label %in% c("LMM", "BP", "GP"),
        !is.na(parameter)
      ) %>%
      group_by(scenario_num, parameter, distribution, model_label) %>%
      summarise(
        true = if (all(is.na(true))) NA_real_ else first(true[!is.na(true)]),
        estimate = mean(estimate, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(distribution_model = paste(distribution, model_label)) %>%
      dplyr::select(scenario_num, parameter, true, distribution_model, estimate) %>%
      pivot_wider(names_from = distribution_model, values_from = estimate) %>%
      arrange(scenario_num, parameter) %>%
      transmute(
        Scenario = scenario_num,
        Parameter = recode(
          as.character(parameter),
          beta_0 = "beta_0",
          beta_1 = "beta_1",
          beta_2 = "beta_2",
          alpha_AFT = "alpha",
          log_AF = "gamma",
          sigma_long = "sigma[e]",
          sd_intercept = "sigma_b0",
          sd_slope = "sigma_b1"
        ),
        True = round(true, 3),
        `Log-logistic LMM` = round(`Log-logistic LMM`, 3),
        `Log-logistic BP` = round(`Log-logistic BP`, 3),
        `Log-logistic GP` = round(`Log-logistic GP`, 3),
        `Weibull LMM` = round(`Weibull LMM`, 3),
        `Weibull BP` = round(`Weibull BP`, 3),
        `Weibull GP` = round(`Weibull GP`, 3)
      )

    message("--- Full main results table with BP/GP standard deviation terms ---")
    print(main_results_with_sd_table, n = Inf)

    kable(main_results_with_sd_table, format = "markdown") %>%
      cat(file = here("Simulation-Workflow", "Results",
                      "main_results_with_sd_table.md"))

    main_results_with_sd_rows <- main_results_with_sd_table %>%
      mutate(
        Scenario = if_else(duplicated(Scenario), "", as.character(Scenario)),
        Parameter = recode(
          Parameter,
          beta_0 = "$\\beta_0$",
          beta_1 = "$\\beta_1$",
          beta_2 = "$\\beta_2$",
          alpha = "$\\alpha$",
          gamma = "$\\gamma$",
          sigma = "$\\sigma$",
          sigma_b0 = "$\\sigma_{b0}$",
          sigma_b1 = "$\\sigma_{b1}$"
        ),
        across(-c(Scenario, Parameter), format_main_result)
      )

    main_results_with_sd_body <- pmap_chr(main_results_with_sd_rows, function(...) {
      values <- unname(unlist(list(...)))
      paste0(paste(values, collapse = " & "), " \\\\")
    })
    full_scenario_end_rows <- which(!duplicated(
      main_results_with_sd_table$Scenario,
      fromLast = TRUE
    ))
    if (length(full_scenario_end_rows) > 1) {
      main_results_with_sd_body[full_scenario_end_rows[-length(full_scenario_end_rows)]] <- paste0(
        main_results_with_sd_body[full_scenario_end_rows[-length(full_scenario_end_rows)]],
        "\n\\addlinespace"
      )
    }

    main_results_with_sd_tex <- c(
      "\\begin{table}[!htbp]",
      "\\centering",
      "\\scriptsize",
      "\\caption{LMM, BP, and GP estimates by data-generation distribution, including BP and GP standard deviation terms.}",
      "\\label{tab:saft-lmm-distribution-estimates-with-sd}",
      "\\begin{tabular}{llrrrrrrr}",
      "\\toprule",
      "Scen. & Par. & True & \\multicolumn{3}{c}{Log-logistic} & \\multicolumn{3}{c}{Weibull} \\\\",
      "\\cmidrule(lr){4-6} \\cmidrule(lr){7-9}",
      " & & & LMM & BP & GP & LMM & BP & GP \\\\",
      "\\midrule",
      main_results_with_sd_body,
      "\\bottomrule",
      "\\end{tabular}",
      "\\vspace{0.5ex}",
      "\\begin{minipage}{0.97\\linewidth}",
      "\\scriptsize Notes: Log-logistic = log-logistic(1.20, 23.0); Weibull = Weibull(0.90, 38.0). LMM denotes the standalone longitudinal mixed model; BP and GP denote the Bernstein-polynomial and Gaussian-process sAFT joint models. Dashes indicate parameters not estimated by the model or true values not stored in the simulation output. Entries are means across simulation replicates per scenario and fitted model.",
      "\\end{minipage}",
      "\\end{table}"
    )
    writeLines(
      main_results_with_sd_tex,
      here("Writeup", "tables", "main_results_with_sd_table.tex")
    )
  }

  table2_2_saft_bp_gp_distribution_details <- all_plot_data %>%
    mutate(
      Scenario = as.integer(regmatches(scenario, regexpr("[0-9]+", scenario))),
      DGM = ifelse(
        grepl("log.logistic", scenario, ignore.case = TRUE),
        "Log-logistic", "Weibull"
      ),
      Model = recode(as.character(model), bp2 = "BP", gp = "GP"),
      Parameter = factor(
        parameter,
        levels = c("beta_0", "beta_1", "beta_2", "alpha_AFT", "log_AF")
      )
    ) %>%
    filter(Model %in% c("BP", "GP"), !is.na(Parameter)) %>%
    group_by(Scenario, Model, Parameter, DGM) %>%
    summarise(
      Bias = mean(estimate - true, na.rm = TRUE),
      ESD = sd(estimate, na.rm = TRUE),
      PSD = mean(posterior_sd, na.rm = TRUE),
      RMSE = sqrt(mean((estimate - true)^2, na.rm = TRUE)),
      CP = mean(ci_lower <= true & ci_upper >= true, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = DGM,
      values_from = c(Bias, ESD, PSD, RMSE, CP),
      names_glue = "{DGM} { .value }"
    ) %>%
    arrange(Scenario, Model, Parameter) %>%
    transmute(
      Scenario,
      Model,
      Parameter = recode(
        as.character(Parameter),
        beta_0 = "beta_0",
        beta_1 = "beta_1",
        beta_2 = "beta_2",
        alpha_AFT = "alpha",
        log_AF = "gamma"
      ),
      `Log-logistic Bias` = round(`Log-logistic Bias`, 4),
      `Log-logistic ESD` = round(`Log-logistic ESD`, 4),
      `Log-logistic PSD` = round(`Log-logistic PSD`, 4),
      `Log-logistic RMSE` = round(`Log-logistic RMSE`, 4),
      `Log-logistic CP` = round(`Log-logistic CP`, 3),
      `Weibull Bias` = round(`Weibull Bias`, 4),
      `Weibull ESD` = round(`Weibull ESD`, 4),
      `Weibull PSD` = round(`Weibull PSD`, 4),
      `Weibull RMSE` = round(`Weibull RMSE`, 4),
      `Weibull CP` = round(`Weibull CP`, 3)
    )

  message("--- Table 2.2: detailed BP and GP performance by distribution ---")
  print(table2_2_saft_bp_gp_distribution_details, n = Inf, width = Inf)

  # Ordered colour vector for x-axis text: ll scenarios in colour_set1, wb in colour_set1_dark
  scenario_x_colours <- scenario_colour_map[levels(all_plot_data$scenario)]

  all_truth <- all_plot_data %>%
    filter(!is.na(true)) %>%
    distinct(parameter, true)

  parameter_labeller <- labeller(
    parameter = as_labeller(
      c(
        beta_0 = "beta[0]",
        beta_1 = "beta[1]",
        beta_2 = "beta[2]",
        alpha_AFT = "alpha",
        log_AF = "gamma"
      ),
      label_parsed
    )
  )

  cross_scenario_boxplot <- ggplot(
    all_plot_data,
    aes(x = scenario, y = estimate, fill = scenario)
  ) +
    geom_boxplot(alpha = 0.72, width = 0.6, outlier.size = 0.7) +
    geom_hline(
      data = all_truth, aes(yintercept = true),
      linetype = "dashed", colour = "red", linewidth = 0.55,
      inherit.aes = FALSE
    ) +
    facet_wrap(~ parameter, scales = "free_y", ncol = 3,
               labeller = parameter_labeller) +
    scale_fill_manual(values = scenario_colour_map) +
    labs(
      title = "Mean estimates across scenarios",
      x = "Scenario", y = "Mean estimates", fill = "Scenario"
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 30, hjust = 1, colour = "black"),
      strip.background = element_rect(fill = "grey90", colour = "grey40"),
      strip.text = element_text(size = 13),
      panel.grid.major.x = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )

  ggsave(filename = here("Writeup", "figures", "cross_scenario_boxplot.pdf"),
         plot = cross_scenario_boxplot, width = 12, height = 8)

  # Combined cross-scenario boxplot showing both BP and GP models side-by-side.
  # fill_group encodes model-type × DGM × scenario number so each box gets a
  # unique colour drawn from the four palette blocks in boxplot_cols.
  combined_plot_data <- all_plot_data %>%
    mutate(
      model_type   = ifelse(grepl("^gp", model, ignore.case = TRUE), "GP", "BP"),
      dgm_type     = ifelse(grepl("log.logistic|\\bll\\b", scenario, ignore.case = TRUE),
                            "loglogistic", "weibull"),
      scenario_num = as.integer(regmatches(scenario,
                                           regexpr("[0-9]+", scenario))),
      fill_group   = paste(model_type, dgm_type,
                           paste0("S", scenario_num), sep = "_")
    )

  cross_scenario_boxplot_combined <- ggplot(
    combined_plot_data,
    aes(x = scenario, y = estimate, fill = fill_group)
  ) +
    geom_boxplot(alpha = 0.85, width = 0.6, outlier.size = 0.7, outlier.alpha = 0.25,
                 position = position_dodge(width = 0.75)) +
    geom_hline(
      data = all_truth, aes(yintercept = true),
      linetype = "dashed", colour = "red", linewidth = 0.55,
      inherit.aes = FALSE
    ) +
    facet_wrap(~ parameter, scales = "free_y", ncol = 3,
               labeller = parameter_labeller) +
    scale_fill_manual(values = boxplot_cols, name = "Model / DGM / Scenario") +
    labs(
      title = "Mean estimates across scenarios (BP+GP sAFT-JM and LMM)",
      x = "Scenario", y = "Mean estimates"
    ) +
    theme_bw() +
    theme(
      legend.position  = "right",
      axis.text.x      = element_text(angle = 30, hjust = 1, colour = "black"),
      strip.background = element_rect(fill = "grey90", colour = "grey40"),
      strip.text = element_text(size = 13),
      panel.grid.major.x = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )

  ggsave(filename = here("Writeup", "figures", "cross_scenario_boxplot_combined.pdf"),
         plot = cross_scenario_boxplot_combined, width = 14, height = 9)

  scenario_specific_plot_data <- combined_plot_data
  if (exists("ph_lmm_plot_data") && nrow(ph_lmm_plot_data) > 0) {
    scenario_specific_plot_data <- bind_rows(
      scenario_specific_plot_data,
      ph_lmm_plot_data %>%
        filter(parameter %in% c("beta_0", "beta_1", "beta_2")) %>%
        mutate(
          model_type = "LMM",
          dgm_type = ifelse(distribution == "Log-Logistic", "loglogistic", "weibull"),
          scenario_num = scenario_number,
          fill_group = paste(model_type, dgm_type, paste0("S", scenario_num), sep = "_")
        ) %>%
        dplyr::select(any_of(names(scenario_specific_plot_data)))
    )
  }

  for (scenario_i in 1:5) {
    scenario_plot_data <- scenario_specific_plot_data %>%
      filter(scenario_num == scenario_i)

    if (nrow(scenario_plot_data) == 0) next

    scenario_truth <- scenario_plot_data %>%
      filter(!is.na(true)) %>%
      distinct(parameter, true)

    scenario_boxplot_combined <- ggplot(
      scenario_plot_data,
      aes(x = scenario, y = estimate, fill = fill_group)
    ) +
      geom_boxplot(alpha = 0.85, width = 0.6, outlier.size = 0.7, outlier.alpha = 0.25,
                   position = position_dodge(width = 0.75)) +
      geom_hline(
        data = scenario_truth, aes(yintercept = true),
        linetype = "dashed", colour = "red", linewidth = 0.55,
        inherit.aes = FALSE
      ) +
      facet_wrap(~ parameter, scales = "free_y", ncol = 3,
                 labeller = parameter_labeller) +
      scale_fill_manual(values = boxplot_cols, name = "Model / DGM / Scenario") +
      scale_x_discrete(
        labels = function(x) {
          ifelse(grepl("Log-Logistic", x), "Loglogistic", "Weibull")
        }
      ) +
      labs(
        title = paste0("Mean estimates in Scenario ", scenario_i,
                       " (BP+GP sAFT-JM and LMM)"),
        x = "Scenario", y = "Mean estimates"
      ) +
      theme_bw() +
      theme(
        legend.position  = "right",
        axis.title       = element_text(size = 13),
        axis.text        = element_text(size = 11, colour = "black"),
        axis.text.x      = element_text(angle = 30, hjust = 1, colour = "black", size = 11),
        legend.title     = element_text(size = 12),
        legend.text      = element_text(size = 11),
        strip.background = element_rect(fill = "grey90", colour = "grey40"),
        strip.text       = element_text(size = 15),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 16)
      )

    ggsave(
      filename = here("Writeup", "figures",
                      paste0("scenario_", scenario_i, "_boxplot_combined.pdf")),
      plot = scenario_boxplot_combined, width = 14, height = 9
    )
    ggsave(
      filename = here("Writeup", "pictures",
                      paste0("scenario_", scenario_i, "_boxplot_combined.pdf")),
      plot = scenario_boxplot_combined, width = 14, height = 9
    )
  }

  if (exists("all_sd_plot_data") && nrow(all_sd_plot_data) > 0) {
    sd_combined_plot_data <- all_sd_plot_data %>%
      mutate(
        model_type = ifelse(grepl("^gp", model, ignore.case = TRUE), "GP", "BP"),
        dgm_type = ifelse(grepl("log.logistic|\\bll\\b", scenario, ignore.case = TRUE),
                          "loglogistic", "weibull"),
        scenario_num = as.integer(regmatches(scenario, regexpr("[0-9]+", scenario))),
        fill_group = paste(model_type, dgm_type, paste0("S", scenario_num), sep = "_"),
        parameter = factor(
          parameter,
          levels = c("sigma_long", "sd_intercept", "sd_slope")
        )
      ) %>%
      filter(!is.na(parameter))

    sd_parameter_labeller <- labeller(
      parameter = as_labeller(
        c(
          sigma_long = "sigma[e]",
          sd_intercept = "sigma[b0]",
          sd_slope = "sigma[b1]"
        ),
        label_parsed
      )
    )

    for (scenario_i in 1:5) {
      scenario_sd_plot_data <- sd_combined_plot_data %>%
        filter(scenario_num == scenario_i)

      if (nrow(scenario_sd_plot_data) == 0) next

      scenario_sd_truth <- scenario_sd_plot_data %>%
        filter(!is.na(true)) %>%
        distinct(parameter, true)

      scenario_sd_boxplot_combined <- ggplot(
        scenario_sd_plot_data,
        aes(x = scenario, y = estimate, fill = fill_group)
      ) +
        geom_boxplot(
          alpha = 0.85,
          width = 0.6,
          outlier.size = 0.7,
          outlier.alpha = 0.25,
          position = position_dodge(width = 0.75)
        ) +
        geom_hline(
          data = scenario_sd_truth,
          aes(yintercept = true),
          linetype = "dashed",
          colour = "red",
          linewidth = 0.55
        ) +
        facet_wrap(~ parameter, scales = "free_y", ncol = 3,
                   labeller = sd_parameter_labeller) +
        scale_fill_manual(values = boxplot_cols, name = "Model / DGM / Scenario") +
        scale_x_discrete(
          labels = function(x) {
            ifelse(grepl("Log-Logistic", x), "Loglogistic", "Weibull")
          }
        ) +
        labs(
          title = paste0("Standard Deviation Parameters in Scenario ", scenario_i),
          x = "Scenario",
          y = "Posterior Mean"
        ) +
        theme_bw() +
        theme(
          legend.position = "right",
          axis.title = element_text(size = 13),
          axis.text = element_text(size = 11, colour = "black"),
          axis.text.x = element_text(angle = 30, hjust = 1, colour = "black", size = 11),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 11),
          strip.background = element_rect(fill = "grey90", colour = "grey40"),
          strip.text = element_text(size = 15),
          panel.grid.major.x = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 16)
        )

      sd_boxplot_file <- here(
        "Writeup", "figures",
        paste0("scenario_", scenario_i, "_sd_boxplot_combined.pdf")
      )
      tryCatch(
        ggsave(
          filename = sd_boxplot_file,
          plot = scenario_sd_boxplot_combined,
          width = 14,
          height = 9
        ),
        error = function(e) {
          fallback_file <- sub("\\.pdf$", "_with_truth.pdf", sd_boxplot_file)
          message("Could not write ", sd_boxplot_file, "; writing ", fallback_file)
          ggsave(
            filename = fallback_file,
            plot = scenario_sd_boxplot_combined,
            width = 14,
            height = 9
          )
        }
      )
    }
  }

}

##########


