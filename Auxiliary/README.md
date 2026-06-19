# Auxiliary Scripts

This directory contains support scripts that are useful for setup and reporting, but are not the main simulation runner.

## Files

- `compile-stan-models.R`: compiles the Stan models used by the workflow through `cmdstanr`. It targets the original BP model, the current BP model, and the GP model.
- `compile-stan-models.slurm`: Slurm wrapper for compiling Stan models on the cluster. It loads R, sets single-threaded BLAS/OpenMP environment variables, and runs the compile script.
- `baseline_haz_comparison.R`: creates baseline-hazard comparison plots for the Weibull and log-logistic data-generating distributions. It also contains helper code for reconstructing fitted BP/GP baseline hazards when the required posterior summaries are available in saved `.rds` outputs.

## Outputs

The plotting script writes figures to:

```text
Figures/
```

The compile scripts create or refresh CmdStan executables next to the Stan model files.

## Notes

- Run these scripts from the `sAFTjm-main/` project root; Stan files are read from `Stan/`.
- The baseline reconstruction section depends on what variables were saved in the simulation summaries; if the fitted baseline quantities are absent, only the data-generating baseline plots can be produced.
