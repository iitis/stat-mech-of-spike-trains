# Spike-train Monte Carlo simulations

This directory contains MATLAB code for Monte Carlo sampling of a fixed-spike-count spike-train ensemble with refractory constraints.

The code was copied from the older `matlab_compilation` directory into a cleaner post-audit structure. The goal is to make the simulations easier to test, rerun, document, and share with collaborators or reviewers.

## Directory structure

- `src/` contains MATLAB source functions.
- `scripts/` contains executable MATLAB and SLURM scripts.
- `tests/` contains sanity checks and audit tests.
- `configs/` contains parameter definitions for production runs.
- `results/` contains new post-audit simulation outputs.
- `figures/` contains generated figures.
- `docs/` contains notes on the method, audit, and reproduction.
- `archive_pre_audit/` contains selected metadata from the older pre-audit setup.

## Basic MATLAB setup

From the repository root, start MATLAB or run MATLAB in batch mode and add the source directory:

```matlab
addpath(genpath('src'));
```

If you also want to run scripts directly from the `scripts/` directory, add:

```matlab
addpath(genpath('scripts'));
```

## Quick audit tests

Before running production simulations, run the audit tests:

```matlab
run('tests/run_all_tests.m');
```

At the moment this file is a placeholder. During the audit it should be extended to run checks such as:

- `test_initialize_mcmc`
- `test_flip_states`
- `test_walk_states`
- `test_transfer_delay_sign`
- `test_apply_operator_fft_ifft`

## Production phase-diagram runs

The production phase diagram is intended to be generated in three linear beta batches:

- beta range `0.0` to `2.0`
- beta range `2.2` to `4.0`
- beta range `4.2` to `6.0`

The beta grid is linear, not logarithmic.

Typical batch scripts should live in `scripts/` and write output into:

```text
results/raw/beta_0_2/
results/raw/beta_2_4/
results/raw/beta_4_6/
```

## Important diagnostics

For each grid point, the post-audit code should store diagnostic information, including:

- final spike count,
- whether the fixed spike-count constraint was satisfied,
- whether the refractory constraint was satisfied,
- final energy,
- acceptance rate,
- random seed,
- code version string.

A simulation point should not be marked as valid unless both the fixed spike-count and refractory checks pass.

## Audit status

This directory was created as a clean post-audit copy of the older `matlab_compilation` folder. Old results are not used for submitted figures unless explicitly regenerated with the audited code.

Planned critical fixes include:

- fixing `initialize_mcmc` output handling from `ref_period`,
- adding `validate_state` checks,
- replacing the relocation move by a symmetric proposal or adding a Hastings correction,
- adding diagnostics to the MCMC core,
- testing the transfer-operator delay sign with an impulse test,
- rerunning phase diagrams after the audit.

## Notes on compiled binaries

Compiled MATLAB binaries are optional. They may be useful for deployment on a cluster, but they are not required for reproducing the source-code runs. The source MATLAB code in `src/` and scripts in `scripts/` should be treated as the primary reproducible version.

## Pre-audit files

The older `matlab_compilation` directory is retained separately for provenance. Selected metadata may be copied into `archive_pre_audit/`, but old result files should not be mixed with post-audit results.
