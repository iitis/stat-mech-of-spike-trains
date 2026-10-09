#!/usr/bin/env bash
#SBATCH --partition=batch
#SBATCH --time=3-00:00:00
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

RUN_SCRIPT_IN="${1:?Usage: sbatch --job-name=phase_b0_2 scripts/current/submit_convergence_1M.sh scripts/current/run_convergence_beta_0_2.m}"

if [[ -f "${RUN_SCRIPT_IN}" ]]; then
  RUN_SCRIPT="${RUN_SCRIPT_IN}"
elif [[ -f "scripts/current/${RUN_SCRIPT_IN}" ]]; then
  RUN_SCRIPT="scripts/current/${RUN_SCRIPT_IN}"
else
  echo "Cannot find run script: ${RUN_SCRIPT_IN}" >&2
  echo "Tried:" >&2
  echo "  ${REPO_ROOT}/${RUN_SCRIPT_IN}" >&2
  echo "  ${REPO_ROOT}/scripts/current/${RUN_SCRIPT_IN}" >&2
  exit 2
fi

mkdir -p results/logs

echo "=== START ==="
echo "JobID: ${SLURM_JOB_ID:-NA}"
echo "Node:  ${SLURMD_NODENAME:-NA}"
echo "Start: $(date)"
echo "CWD:   $(pwd)"
echo "Script: ${RUN_SCRIPT}"
echo

module load trytonp/matlab/R2025b

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

SOFT_TIMEOUT_S=$((72*3600 - 10*60))

if command -v timeout >/dev/null 2>&1; then
  timeout --preserve-status --signal=TERM --kill-after=30s "${SOFT_TIMEOUT_S}" \
    matlab -nodisplay -nosplash -batch "run('${RUN_SCRIPT}')"
else
  matlab -nodisplay -nosplash -batch "run('${RUN_SCRIPT}')"
fi

echo
echo "End: $(date)"
echo "=== END ==="
