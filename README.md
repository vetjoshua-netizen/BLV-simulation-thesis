# BLV Individual-Based Simulation

This repository contains an individual-based stochastic simulation model for bovine leukaemia virus (BLV) transmission in a dairy herd. It was developed for a thesis project to compare a baseline scenario with three intervention mechanisms:

- corrective management, which reduces the horizontal force of infection;
- test-and-remove;
- preventive vaccination.

The model combines scheduled herd events with stochastic events simulated using the Gillespie stochastic simulation algorithm (SSA).

## Requirements

- R (tested with R 4.6.1)
- R package: `R6`
- R package: `ggplot2` for `simulation.run.R` and its plots

Install the required packages in R if necessary:

```r
install.packages(c("R6", "ggplot2"))
```

## Files

| File | Purpose |
| --- | --- |
| `simulation.R` | Core model, parameter definitions, intervention logic, simulation engine, and replicate runner. |
| `analysis.R` | Functions for summarising simulation outputs. |
| `simulation.run.R` | Main scenario script. Edit parameters here, run simulations, and create summary CSV files and plots. |
| `test.baseline.R` | Baseline-model checks, including reproducibility and consistency checks. |
| `test.interventions.R` | Focused checks for corrective management, vaccination, and test-and-remove. |
| `simulation_annotated.R` | Annotated version of the model code for explanation and reference. |

## Running a scenario

Open `simulation.run.R` in RStudio and edit the scenario name, model parameters, and number of replicates as needed.

```r
n_runs <- 1000L
master_seed <- 12345L

results <- run_strategy(
  p,
  n_runs = n_runs,
  master_seed = master_seed
)
```

The script creates intervention summaries and annual records in CSV format. It also generates population plots when run interactively.

## Reproducibility

`master_seed` controls the full set of replicate seeds. `run_strategy()` uses the master seed to generate a unique seed for each replicate, stores those seeds in `results$replicate_seeds`, and runs each replicate with its assigned seed.

Use the same `master_seed` and `n_runs` across intervention scenarios when making comparisons. This ensures that every scenario uses the same set of replicate seeds.

## Important interpretation note

The model is stochastic, so individual replicates can have different final prevalence values even when they use identical parameters. Compare scenarios using summaries across replicates, such as mean or median prevalence and the probability of eradication.

In this code, `eradicated` is an operational surveillance outcome after the required number of clean surveillance rounds. It does not necessarily mean that the herd is infection-free.

## Tests

Run the checks from the `simulation` directory:

```r
source("test.baseline.R")
source("test.interventions.R")
```

## Notes

Several epidemiological parameter defaults in `simulation.R` are placeholders. Confirm all parameter values and assumptions against the literature and the final thesis methodology before reporting results.
