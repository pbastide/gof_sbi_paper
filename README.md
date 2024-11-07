# Repository to reproduce the analyses in Le Mailloux, Bastide, Marin, Estoup 2024+

## Content

Please refer to the `README` files in each sub-directory for more information on the analyses.
All the `R` scripts assume that the working directory is the root directory
(location of the `.Rproj` and `.here` files).

* `data`: 
  * `simulation_toy_gaussian_laplace.R`: R script to generate the toy "Laplace - Gaussian" dataset
  * `gaussianlaplace` repository containing the toy "Laplace - Gaussian" dataset

* `diyabc_simulation_scripts`: Simulation scripts for diyabc
  * `sim_depindep.R`: simulation script for the "Dep-Indep" example
  * `deptindep_resim_nocond`: header files for "Dep-Indep" simulation
  * `sim_human.R`: simulation script for the "Human-like" example
  * `human_resim`: header files for "Human-like" simulation

* `gof_sbi_paper.Rproj`: Root `Rproj` file.

* `human_analysis`: files for the empirical human data analysis.

* `simulation_study`: script to reproduce the simulation study.
  * `gof_simulaiton_study.R`: Script for simulation studies.
  * `gof_simulaiton_study_plots.R`: Script to plot the results.
  * `gof_lof_knn_example_plot.R`: Script to plot an example analysis.

## Requirements

* `DIYABC`: https://github.com/diyabc/diyab.
  Version 1.1.54 or above.
  For simulation of summary statistics from population genetics models.
  Pre-build releases are available from https://github.com/diyabc/diyabc/releases

* `R`: https://cran.r-project.org/index.html.
  Version 4.4. or above.

* `abcgof`: https://github.com/pbastide/abcgof
  Version v0.0.1.
  Can be installed from `R` with: `devtools::install_github("pbastide/abcgof", ref = "v0.0.1")`
