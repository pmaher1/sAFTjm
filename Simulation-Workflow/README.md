# Simulation Workflow

This directory contains the main simulation pipeline for the sAFT joint-model comparison. Note that this was constructed for alignment with high-performance computing (HPC), in particular UQ's HPC system [Bunya](https://rcc.uq.edu.au/systems/high-performance-computing/bunya).

## Files

- `functions.R`: shared implementation for data generation and model fitting. It defines `simulate_joint_dataset()` and `fit_one_rep()`, builds Stan data for the BP and GP models, fits selected models, extracts posterior summaries, and computes LOO/WAIC for Stan fits.
- `datagen.R`: scenario registry and data-generating parameter setup. It maps scenario IDs such as `s1_ll_n300` to treatment effects, baseline-hazard families, sample sizes, censoring settings, and shared simulation constants.
- `hpc-simulate.R`: batch runner used by Slurm. It reads command-line arguments, sources `datagen.R`, calls `fit_one_rep()` across replicate indices, and saves batch `.rds` files under `Simulation-Workflow/Results/<scenario_id>/`.
- `run_hpc_simulate.slurm`: Slurm array script. It expands the treatment, baseline, and sample-size grid into scenario IDs, maps array task IDs to scenario/batch pairs, and launches `hpc-simulate.R`.
- `result-formatting.R`: post-processing script for completed simulation outputs. It reads saved `.rds` batches, builds cross-scenario summaries, diagnostics, information-criteria summaries, runtime summaries, figures, and supporting images.

## Outputs

Generated simulation batches are written under:

```text
Simulation-Workflow/Results/
```

Slurm logs are written under:

```text
Simulation-Workflow/Logs/
```

Figures are written under:

```text
Figures/
```

Supporting images are written under:

```text
Pictures/
```

## Notes

- `hpc-simulate.R` accepts model names after the worker argument. Supported models are defined in `functions.R`.
- `run_hpc_simulate.slurm` currently defines the production scenario grid. Keep its array size aligned with the number of scenarios times `N_BATCHES`.
- Batch files are named `sAFTjm-<scenario_id>-<aft_mode>-<batch_id>.rds`.
- Run scripts from the `sAFTjm-main/` project root; Stan files are read from `Stan/`.
