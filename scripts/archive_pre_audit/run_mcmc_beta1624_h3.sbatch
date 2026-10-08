#!/usr/bin/env bash
#SBATCH --job-name=phase
#SBATCH --partition=batch
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --array=1-30%10
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -euo pipefail

trap 'echo "ERROR @ line $LINENO: $BASH_COMMAND" >&2' ERR

mkdir -p logs results

echo "=== START ==="
echo "JobID:        ${SLURM_JOB_ID:-NA}"
echo "ArrayJobID:   ${SLURM_ARRAY_JOB_ID:-NA}"
echo "ArrayTaskID:  ${SLURM_ARRAY_TASK_ID:-NA}"
echo "Node:         ${SLURMD_NODENAME:-NA}"
echo "Start:        $(date)"
echo "CWD:          $(pwd)"
echo

RUNNER="./run_run_phase_diagram_v0208_bin.sh"
MCRROOT="/apl/trytonp/matlab/R2025b"

test -x "${RUNNER}" || { echo "ERROR: brak wykonywalnego ${RUNNER}"; exit 2; }
test -d "${MCRROOT}" || { echo "ERROR: MCRROOT nie istnieje: ${MCRROOT}"; exit 2; }

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

# Oryginalna siatka:
# sdel: 2.0, 2.5, 3.0, 3.5, 4.0, 4.5
# beta exponent/log10 beta: 1.6, 1.8, 2.0, 2.2, 2.4
#
# beta = 10 ^ beta_exp
SDELS=(2.0 2.5 3.0 3.5 4.0 4.5)
BETA_EXPS=(1.6 1.8 2.0 2.2 2.4)

NS=${#SDELS[@]}
NB=${#BETA_EXPS[@]}
NTOTAL=$((NS * NB))

TASK_ID="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID not set}"

if (( TASK_ID < 1 || TASK_ID > NTOTAL )); then
  echo "ERROR: SLURM_ARRAY_TASK_ID=${TASK_ID}, ale dozwolone 1..${NTOTAL}"
  exit 2
fi

# Mapowanie task_id -> ib,is
# task_id=1  -> ib=1,is=1
# task_id=2  -> ib=1,is=2
# ...
# task_id=6  -> ib=1,is=6
# task_id=7  -> ib=2,is=1
ZERO_BASED=$((TASK_ID - 1))
IB0=$((ZERO_BASED / NS))
IS0=$((ZERO_BASED % NS))

IB=$((IB0 + 1))
IS=$((IS0 + 1))

SDEL="${SDELS[$IS0]}"
BETA_EXP="${BETA_EXPS[$IB0]}"

BETA_VALUE="$(awk -v x="${BETA_EXP}" 'BEGIN { printf "%.10g", 10^x }')"

RUN_NAME="Gphase6x5_beta16to24_h3_ib${IB}_is${IS}_bexp${BETA_EXP}_sdel${SDEL}"
OUT_FILE="results/${RUN_NAME}_phase.mat"

echo "Grid:"
echo "  NB=${NB}"
echo "  NS=${NS}"
echo "  NTOTAL=${NTOTAL}"
echo
echo "This task:"
echo "  TASK_ID=${TASK_ID}"
echo "  ib=${IB}/${NB}"
echo "  is=${IS}/${NS}"
echo "  beta_exp=${BETA_EXP}"
echo "  beta≈${BETA_VALUE}"
echo "  sdel=${SDEL}"
echo "  out_file=${OUT_FILE}"
echo

# Resume / skip: jeśli wynik już istnieje, nie licz ponownie.
if [[ -s "${OUT_FILE}" ]]; then
  echo "SKIP: wynik już istnieje: ${OUT_FILE}"
  echo "=== END ==="
  exit 0
fi

# Tutaj zwężamy zakresy do pojedynczego punktu:
#
# Oryginał:
#   sdel_min sdel_max sdel_step = 2.0 4.5 0.5
#   beta_min beta_max beta_step = 1.6 2.4 0.2
#
# Teraz:
#   sdel_min=sdel_max=SDEL
#   beta_min=beta_max=BETA_EXP
#
ARGS=(
  "results/${RUN_NAME}"
  500
  0.3
  20
  "${SDEL}" "${SDEL}" 0.5
  1000
  2000
  3
  "${BETA_EXP}" "${BETA_EXP}" 0.2
  200000
  20000
  heaviside3
  1
)

# Dla joba 24h zostawiamy 10 minut bufora na miękkie zamykanie.
SOFT_TIMEOUT_S=$((24*3600 - 10*60))

echo "Runner:       ${RUNNER}"
echo "MCR:          ${MCRROOT}"
echo "Args:         ${ARGS[*]}"
echo "Soft timeout: ${SOFT_TIMEOUT_S}s"
echo

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