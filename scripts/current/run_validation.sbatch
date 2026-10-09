#!/usr/bin/env bash
#SBATCH --job-name=spike_validation
#SBATCH --partition=batch
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=results/logs/%x_%j.out
#SBATCH --error=results/logs/%x_%j.err

set -euo pipefail
trap 'echo "ERROR @ line $LINENO: $BASH_COMMAND" >&2' ERR
REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${REPO_ROOT}"
module load "${MATLAB_MODULE:-trytonp/matlab/R2025b}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
echo "Job: ${SLURM_JOB_ID:-NA}; host: $(hostname); date: $(date -Is)"
git rev-parse HEAD
git status --short
matlab -nodisplay -nosplash -batch "run('scripts/current/run_validation.m')"
