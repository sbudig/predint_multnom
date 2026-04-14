# Supplementary Information / Reproducible Research Files

Manuscript title: "Prediction intervals for overdispersed multinomial data with application to historical controls"

Authors: Sören Budig, Frank Schaarschmidt, Max Menssen

Contact: budig@cell.uni-hannover.de

The folder `Code_and_Data` contains the code, example data, figures, and simulation results required to reproduce the analyses and graphics for the manuscript.

The project was developed in R on Windows. The main scripts rely on the following R packages: `here`, `dplyr`/`tidyverse`, `VGAM`, `MCMCpack`, `future.apply`, `mvtnorm`, `cmdstanr`, `patchwork`, `viridis`, `RColorBrewer`, `ggsci`, `wesanderson`, `forcats`, and `xtable`.

To run the Bayesian methods from source, a working CmdStan installation is required in addition to the `cmdstanr` package.

The folder contains the following subfolders and files:

## `./Code/`

This folder contains the code for the simulation study, the example analyses, and the shared function library used throughout the project.

- `predint_mult_source.R`: main source file with functions for data generation, model fitting, dispersion estimation, interval construction, and simulation helpers.
- `dirichlet_multinomial_gamma.stan`, `dirichlet_multinomial_cauchy.stan`, `dirichlet_multinomial_beta_rho.stan`: Stan model definitions for the Bayesian approaches.
- `dirichlet_multinomial_gamma`, `dirichlet_multinomial_cauchy`, `dirichlet_multinomial_beta_rho` and the corresponding `.exe` files: compiled CmdStan model artifacts that can be regenerated from the `.stan` files.

### `./Code/Example_Analysis/`

This folder contains the scripts and datasets for the applied examples discussed in the manuscript. Note that the dart example is additional and was not part of the manuscript.

- `Example_Analysis_patho.R`: reproduces the pathology example analysis.
- `Example_Analysis_dart.R`: reproduces the dartboard example analysis.
- `df_pat.csv`: pathology example dataset.
- `df_dart.csv`: dartboard example dataset.
- `example_patho_all_ints.rds`: saved output object from the pathology example.

### `./Code/Simulation_Study/`

This folder contains the scripts needed to generate and process the simulation study results.

- `predint_mult_sim.R`: main simulation script. It sources the shared functions, compiles the Stan models, defines the simulation settings, and writes intermediate `.rds` results to the `Results/intermediate_results` folders.
- `predint_mult_processing.R`: reads the intermediate `.rds` files, summarizes coverage and interval-width measures, and produces the figures used in the manuscript.

## Reproduction Notes

- Run the scripts from the `Code_and_Data` project root so that `here::here()` resolves paths correctly.
- Use `Code/Simulation_Study/predint_mult_sim.R` to generate simulation output.
- Use `Code/Simulation_Study/predint_mult_processing.R` after the simulations have finished to summarize the results and recreate the figures.
- Use the scripts in `Code/Example_Analysis/` to reproduce the applied examples.
- Bayesian methods are computationally intensive and require Stan model compilation.






