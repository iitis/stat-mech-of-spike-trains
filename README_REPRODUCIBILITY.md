# Reproducibility notes

This repository contains the Monte Carlo and mean-field numerical code used for the manuscript

**Statistical Mechanics of Spike Trains: Self-Consistent Activity Ensembles and Phase Transitions**.

## Repository structure

- `src/` — core MATLAB routines for the Monte Carlo sampler, transfer operator, delay sampling, and state validation.
- `scripts/current/` — scripts used for the current numerical phase-map calculations and plotting.
- `tests/` — MATLAB tests checking move invariants, refractory constraints, transfer-operator delay conventions, and sampler behavior.
- `figures/` — processed figure outputs.
- `results/raw/` — raw or semi-processed numerical outputs used for the current figures.

## Monte Carlo phase map

Each point of the Monte Carlo phase map corresponds to one quenched graph-and-delay realization and one finite Markov chain. The map should therefore be interpreted as a single-realization finite-size phase portrait, not as a disorder-averaged phase diagram.

## Notes on stochasticity

The simulations use random graph and delay realizations. Reproducing bitwise-identical outputs may require matching MATLAB version, random seed handling, and job-splitting conventions. The repository is intended to reproduce the qualitative phase-map pipeline and the processed numerical figures accompanying the manuscript.

## Large data

The repository contains the data currently needed to reproduce the manuscript figures. Larger raw simulation outputs, if needed, should be archived separately in a versioned data repository such as Zenodo.
