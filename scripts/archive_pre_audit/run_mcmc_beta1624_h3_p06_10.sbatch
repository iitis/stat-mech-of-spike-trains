#!/usr/bin/env bash
#SBATCH --job-name=phase_p06_10
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

mkdir -p logs results

echo "=== START ==="
echo "JobID: ${SLURM_JOB_ID:-NA}"
echo "Node:  ${SLURMD_NODENAME:-NA}"
echo "Start: $(date)"
echo "CWD:   $(pwd)"
echo

RUNNER="./run_run_phase_diagram_v0208_bin.sh"
MCRROOT="/apl/trytonp/matlab/R2025b"

test -x "${RUNNER}" || { echo "ERROR: brak wykonywalnego ${RUNNER}"; exit 2; }
test -d "${MCRROOT}" || { echo "ERROR: MCRROOT nie istnieje: ${MCRROOT}"; exit 2; }

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

# Oryginalna siatka:
# sdel:      2.0, 2.5, 3.0, 3.5, 4.0, 4.5
# beta_exp:  1.6, 1.8, 2.0, 2.2, 2.4
#
# Kolejność punktów:
# 1=(beta_exp=1.6,sdel=2.0), ..., 6=(1.6,4.5),
# 7=(1.8,2.0), ..., 30=(2.4,4.5)

SDELS=(2.0 2.5 3.0 3.5 4.0 4.5)
BETA_EXPS=(1.6 1.8 2.0 2.2 2.4)

NS=${#SDELS[@]}
NB=${#BETA_EXPS[@]}
NTOTAL=$((NS * NB))

START_POINT=6
END_POINT=10

# Miękki timeout na pojedynczy punkt.
# Przy jobie 3-dniowym i 5 punktach daje to margines na nierówne czasy.
POINT_TIMEOUT_S=$((14*3600))

echo "Runner: ${RUNNER}"
echo "MCR:    ${MCRROOT}"
echo "Grid:   NB=${NB}, NS=${NS}, NTOTAL=${NTOTAL}"
echo "Chunk:  points ${START_POINT}..${END_POINT}"
echo "Point soft timeout: ${POINT_TIMEOUT_S}s"
echo

for POINT_ID in $(seq "${START_POINT}" "${END_POINT}"); do
  echo
  echo "============================================================"
  echo "POINT ${POINT_ID}/${NTOTAL}"
  echo "Start point: $(date)"
  echo "============================================================"

  ZERO_BASED=$((POINT_ID - 1))
  IB0=$((ZERO_BASED / NS))
  IS0=$((ZERO_BASED % NS))

  IB=$((IB0 + 1))
  IS=$((IS0 + 1))

  SDEL="${SDELS[$IS0]}"
  BETA_EXP="${BETA_EXPS[$IB0]}"
  BETA_VALUE="$(awk -v x="${BETA_EXP}" 'BEGIN { printf "%.10g", 10^x }')"

  RUN_NAME="Gphase6x5_beta16to24_h3_p${POINT_ID}_ib${IB}_is${IS}_bexp${BETA_EXP}_sdel${SDEL}"
  OUT_PREFIX="results/${RUN_NAME}"
  OUT_FILE="${OUT_PREFIX}_phase.mat"

  echo "ib=${IB}/${NB}"
  echo "is=${IS}/${NS}"
  echo "beta_exp=${BETA_EXP}"
  echo "beta≈${BETA_VALUE}"
  echo "sdel=${SDEL}"
  echo "out_file=${OUT_FILE}"

  if [[ -s "${OUT_FILE}" ]]; then
    echo "SKIP: wynik już istnieje: ${OUT_FILE}"
    continue
  fi

  ARGS=(
    "${OUT_PREFIX}"
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

  echo "Args: ${ARGS[*]}"
  echo

  if command -v timeout >/dev/null 2>&1; then
    timeout --preserve-status --signal=TERM --kill-after=30s "${POINT_TIMEOUT_S}" \
      "${RUNNER}" "${MCRROOT}" "${ARGS[@]}"
  else
    echo "WARN: brak 'timeout' w systemie, lecę bez miękkiego limitu."
    "${RUNNER}" "${MCRROOT}" "${ARGS[@]}"
  fi

  echo "End point: $(date)"
done

echo
echo "End: $(date)"
echo "=== END ==="
