#!/usr/bin/env bash
set -euo pipefail

mkdir -p logs results

sbatch run_mcmc_beta1624_h3_p01_05.sbatch
sbatch run_mcmc_beta1624_h3_p06_10.sbatch
sbatch run_mcmc_beta1624_h3_p11_15.sbatch
sbatch run_mcmc_beta1624_h3_p16_20.sbatch
sbatch run_mcmc_beta1624_h3_p21_25.sbatch
sbatch run_mcmc_beta1624_h3_p26_30.sbatch
