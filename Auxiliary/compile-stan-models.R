library(tidyverse)
library(survival)
library(MASS)
library(future)
library(furrr)
library(simsurv)
library(rstan)
library(brms)
library(cmdstanr)

# Set CmdStan path globally for main session
set_cmdstan_path("~/cmdstan/cmdstan-2.37.0")

# Bernstein Polynomials (original)
cmdstan_model("~/WPP/Stan/JM/Bernstein-Polynomials-JM-Hist.stan",
              cpp_options = list(stan_threads = TRUE),
              force_recompile = TRUE)

# Bernstein Polynomials (new)
cmdstan_model("~/WPP/Stan/JM/bernstein-polynomials.stan",
              cpp_options = list(stan_threads = TRUE),
              force_recompile = TRUE)

# Gaussian Processes 
cmdstan_model("~/WPP/Stan/JM/gaussian-process.stan",
              cpp_options = list(stan_threads = TRUE),
              force_recompile = TRUE)






