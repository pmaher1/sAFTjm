#!/bin/bash
#SBATCH --job-name=WPP-sim
#SBATCH --array=1-500
#SBATCH --exclude=bun091
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=6GB
#SBATCH --time=06:00:00
#SBATCH --output=WPP/Simulation-Workflow/Logs/task-%04a.out
#SBATCH --error=WPP/Simulation-Workflow/Logs/task-%04a.err
#SBATCH --mail-type=FAIL

N_PER_BATCH=20
N_WORKERS=${SLURM_CPUS_PER_TASK}
N_BATCHES=50

# Scenario grid: treatment effect x baseline hazard x sample size.
# Keep #SBATCH --array in sync with ${#SCENARIOS[@]} * N_BATCHES.
TREATMENT_SCENARIOS=(s1 s2 s3 s4 s5)
BASELINE_SCENARIOS=(ll wb)
SAMPLE_SIZE_SCENARIOS=(n300)

SCENARIOS=()
for SAMPLE_SIZE_ID in "${SAMPLE_SIZE_SCENARIOS[@]}"; do
    for BASELINE_ID in "${BASELINE_SCENARIOS[@]}"; do
        for TREATMENT_ID in "${TREATMENT_SCENARIOS[@]}"; do
            SCENARIOS+=("${TREATMENT_ID}_${BASELINE_ID}_${SAMPLE_SIZE_ID}")
        done
    done
done

# Derive 0-based indices from the 1-based array task ID.
TASK_IDX=$(( SLURM_ARRAY_TASK_ID - 1 ))
SCENARIO_IDX=$(( TASK_IDX / N_BATCHES ))
BATCH_ID=$(( TASK_IDX % N_BATCHES + 1 ))
if [ "${SCENARIO_IDX}" -ge "${#SCENARIOS[@]}" ]; then
    echo "Array task ${SLURM_ARRAY_TASK_ID} exceeds configured scenario grid" >&2
    exit 1
fi
SCENARIO_ID=${SCENARIOS[$SCENARIO_IDX]}

module load r/4.4.0-gfbf-2023a

# Each jm2 chain is already an independent process. Prevent BLAS/OpenMP from
# creating extra threads inside every chain and oversubscribing the 16 CPUs.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

mkdir -p "${SLURM_SUBMIT_DIR}/WPP/Simulation-Workflow/Logs"

cd "${SLURM_SUBMIT_DIR}"

Rscript WPP/Simulation-Workflow/hpc-simulate.R \
    ${BATCH_ID} \
    ${N_PER_BATCH} \
    ${N_WORKERS} \
    bp2 \
    ${SCENARIO_ID}
