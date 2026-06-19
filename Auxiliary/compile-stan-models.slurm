#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=0:10:00
#SBATCH --partition=general
#SBATCH --account=a_uqccr
#SBATCH --job-name=compile-stan

module purge
module load r/4.4.0-gfbf-2023a

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo "Compiling Stan models..."
Rscript Auxiliary/compile-stan-models.R
echo "Done."

