#!/usr/bin/env bash
#SBATCH --job-name=phase_b2_4
#SBATCH --partition=batch
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=results/logs/%x_%j.out
#SBATCH --error=results/logs/%x_%j.err

set -euo pipefail
trap 'echo "ERROR @ line $LINENO: $BASH_COMMAND" >&2' ERR

mkdir -p results/logs

echo "=== START ==="
echo "JobID: ${SLURM_JOB_ID:-NA}"
echo "Node:  ${SLURMD_NODENAME:-NA}"
echo "Start: $(date)"
echo "CWD:   $(pwd)"
echo

module load trytonp/matlab/R2025b

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

SOFT_TIMEOUT_S=$((72*3600 - 10*60))

if command -v timeout >/dev/null 2>&1; then
  timeout --preserve-status --signal=TERM --kill-after=30s "${SOFT_TIMEOUT_S}" \
    matlab -nodisplay -nosplash -batch "run('scripts/current/run_phase_beta_2p2_4.m')"
else
  matlab -nodisplay -nosplash -batch "run('scripts/current/run_phase_beta_2p2_4.m')"
fi

echo
echo "End: $(date)"
echo "=== END ==="
