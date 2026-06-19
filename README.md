# sAFTjm

Joint modelling of longitudinal and time-to-event outcomes with semiparametric accelerated failure time joint models.

This repository contains the simulation and analysis code used to compare:

- `LMM`: a standalone longitudinal mixed model.
- `bp1`: the original Bernstein-polynomial sAFT joint model.
- `bp2`: the current Bernstein-polynomial sAFT joint model.
- `gp`: the Gaussian-process sAFT joint model.
- `jm2`: an optional proportional-hazards joint model fit with `JMbayes2`.

## Directory Layout

- `Simulation-Workflow/`: data generation, model fitting, HPC batch execution, and result formatting.
- `Auxiliary/`: supporting scripts for Stan compilation and baseline-hazard plots.
- `Stan/`: Stan model files for the BP and GP sAFT joint models.
- `local-fit.R`: single-replicate local runner for checking that the LMM, BP, and GP workflow runs outside the HPC array.

Run scripts from this `sAFTjm-main/` project root. Generated outputs stay inside this folder, with simulation results and logs under `Simulation-Workflow/`, plots under `Figures/`, and non-figure images under `Pictures/`.

## Typical Workflow

1. Compile the Stan models if needed using `Auxiliary/compile-stan-models.R` or the matching Slurm script.
2. Run a local smoke test with `local-fit.R`.
3. Run the simulation grid through `Simulation-Workflow/run_hpc_simulate.slurm`.
4. Summarise completed `.rds` outputs with `Simulation-Workflow/result-formatting.R`.

See the READMEs in `Simulation-Workflow/` and `Auxiliary/` for file-level details.
