#!/usr/bin/env bash
#SBATCH --job-name=missing_exprnd
#SBATCH --partition=batch
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=results/logs/%x_%j.out
#SBATCH --error=results/logs/%x_%j.err

set -euo pipefail
trap 'echo "ERROR @ line $LINENO: $BASH_COMMAND" >&2' ERR

mkdir -p results/logs

echo "=== START missing exprnd rerun ==="
echo "JobID: ${SLURM_JOB_ID:-NA}"
echo "Node:  ${SLURMD_NODENAME:-NA}"
echo "Start: $(date)"
echo "CWD:   $(pwd)"
echo

module load trytonp/matlab/R2025b

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

matlab -nodisplay -nosplash -batch "run('scripts/current/run_missing_exprnd_points.m')"

echo
echo "End: $(date)"
echo "=== END missing exprnd rerun ==="
