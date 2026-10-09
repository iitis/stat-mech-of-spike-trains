# Reproducibility and validation

See `README.md` for experiment parameters and basic commands.

## Data provenance

The main map uses 100,000 iterations per point. Seven points failed during the original calculation and were rerun separately:

| Beta | Delay standard-deviation parameter | Replacement run seed |
| --- | --- | --- |
| 0.2 | 4.5 | 100201 |
| 1.8 | 4.5 | 101801 |
| 2.0 | 4.5 | 102001 |
| 2.4 | 4.5 | 102401 |
| 4.0 | 4.5 | 104001 |
| 4.4 | 4.5 | 104401 |
| 6.0 | 4.5 | 106001 |

Replacement calculations are in `results/raw/missing_exprnd_fix/`. This historical name is not a diagnosis of every failure. Raw TXT files retain the failed entries; processed TXT exports represent the reconstructed map.

The merge script rebuilds three batches under `results/processed/main_map/`, preserving raw inputs. It validates replacement coordinates, parameters, and diagnostics, and copies activity, energy, validity flags, snapshots, and graph/delay fields. `source_file_grid`, `seed_grid`, and `reconstruction_history` identify their origins. New provenance paths are relative to the repository; legacy metadata in raw inputs are retained.

A base seed initializes a whole batch, not each individual point. Replacement points use separate graph/delay realizations and random streams. Reordering or splitting a run changes random-number consumption. A fresh simulation is distinct from reconstruction of historical results.

The coarser 1M grid was added to assess convergence of the 100k results. It is not the source of the main entropy figures. A quantitative convergence statement requires specified observables and tolerances, with explicit treatment of graph/delay variability. Neither experiment is disorder averaged, and sampler invariant tests do not establish equilibration or mixing of production chains.

## Static checks

From the repository root:

```bash
python3 tests/check_repository_static.py
```

The checker uses Python's standard library and Git. It checks UTF-8, Python syntax, README commands, all six experiment configurations, SLURM targets, sampler test wiring, reconstruction inputs/output paths, and `git diff --check`. If Bash is available it also runs `bash -n` on every shell/SLURM script. It does not validate MATLAB semantics.

## Validation on cluster

Copy the complete working directory, including new files. Run static checks first, then submit from the repository root:

```bash
mkdir -p results/logs
sbatch scripts/current/run_validation.sbatch
```

The wrapper defaults to partition `batch` and module `trytonp/matlab/R2025b`, as in the original scripts. Confirm their availability on the target cluster. To override them:

```bash
sbatch --partition=YOUR_PARTITION --export=ALL,MATLAB_MODULE=YOUR_MATLAB_MODULE scripts/current/run_validation.sbatch
```

On an allocated compute node with MATLAB loaded, the equivalent command is:

```bash
matlab -nodisplay -nosplash -batch "run('scripts/current/run_validation.m')" > results/logs/validation.log 2>&1
```

The driver runs eight sampler tests, a small smoke simulation, reconstruction, and `test_main_map_reconstruction`. It repeats reconstruction and verification, then regenerates entropy/IPR figures. No production-map runs are launched. It overwrites smoke outputs, reconstructed outputs, and generated figures, while preserving original main-map and replacement inputs.

The reconstruction test compares non-snapshot numerical/graph arrays against raw data with exactly seven replacements, checks replacement snapshots, per-point provenance, history lengths, and every TXT row. It does not scan every untouched snapshot. A MATLAB exception terminates batch execution with failure; success prints `VALIDATION PASSED`.

Keep full stdout and stderr, including MATLAB version, platform, Git commit, working-tree status, stage markers, and stack traces. Logs are `results/logs/spike_validation_JOBID.out` and `.err`. If MATLAB never starts, inspect scheduler state and exit code as well.