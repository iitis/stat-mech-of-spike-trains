#!/usr/bin/env bash
#SBATCH --job-name=phase
#SBATCH --partition=batch
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

# Głośne błędy: linia + komenda
trap 'echo "ERROR @ line $LINENO: $BASH_COMMAND" >&2' ERR

mkdir -p logs

echo "=== START ==="
echo "JobID: ${SLURM_JOB_ID:-NA}"
echo "Node:  ${SLURMD_NODENAME:-NA}"
echo "Start: $(date)"
echo "CWD:   $(pwd)"
echo

RUNNER="./run_run_phase_diagram_v0220_bin.sh"
MCRROOT="/apl/trytonp/matlab/R2025b"

ARGS=(
  Gphase30x13_betalin6to9_h3
  100
  0.3
  5
  0.0 4.5 0.15
  500
  500
  3
  6.0 9.0 0.25
  100000
  10000
  heaviside3
  1
)


SOFT_TIMEOUT_S=$((72*3600 - 10*60))   # 261600

echo "Runner: ${RUNNER}"
echo "MCR:    ${MCRROOT}"
echo "Args:   ${ARGS[*]}"
echo "Soft timeout: ${SOFT_TIMEOUT_S}s"
echo

test -x "${RUNNER}" || { echo "ERROR: brak wykonywalnego ${RUNNER}"; exit 2; }
test -d "${MCRROOT}" || { echo "ERROR: MCRROOT nie istnieje: ${MCRROOT}"; exit 2; }

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

if command -v timeout >/dev/null 2>&1; then
  timeout --preserve-status --signal=TERM --kill-after=30s "${SOFT_TIMEOUT_S}" \
    "${RUNNER}" "${MCRROOT}" "${ARGS[@]}"
else
  echo "WARN: brak 'timeout' w systemie, lecę bez miękkiego limitu."
  "${RUNNER}" "${MCRROOT}" "${ARGS[@]}"
fi

echo
echo "End: $(date)"
echo "=== END ==="
