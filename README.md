# Repository to reproduce the analyses in Le Mailloux, Bastide, Marin, Estoup 2024+

## Content

* `diyabc_reftables`: simulated reference tables. Please see dedicated README file.
  * `gaussianlaplace`: toy "Laplace - Gaussian" dataset
    * `simulation_toy_gaussian_laplace.R`: R script to generate the toy "Laplace - Gaussian" dataset
  * `depindep`: Dep-Indep dataset
    * `01_Generate_reftableRF.bin_for_depindep_fast_NoConditionMediumMoins.sh`: Shell script to generate bin reftable
    * `02_Rscript_to_generate_parameter_files_for_each_scenario_separetly.R`: R script to translate bin to txt reftables
    * `03_depindep_fast_NoConditionSmallT_reftableRF_110000_param.sh`: shell script to generate replicated reftables
    * `DSIM_indep_dep_4pop_SNPind_10indsPerPop.snp` SNP data file
    * `headerRF_for_depindep_fast_NoConditionMediumMoins.txt` header file for DIYABC
  * `humanfast`: Human like dataset
    * `01_Generate_reftableRF.bin_for_diyabc_generated_reftableRF_human_fast.sh`: Shell script to generate bin reftable
    * `02_Rscript_to_generate_parameter_files_for_each_scenario_separetly.R`: R script to translate bin to txt reftables
    * `03_human_fast_reftableRF_110000_param.sh`: shell script to generate replicated reftables
    * `human_snp_all22chr_hudson_10IndPerPop.snp` SNP data file
    * `headerRF_for_diyabc_generated_reftableRF_human_fast.txt` header file for DIYABC
  * `humanreal`: Human like dataset
    * `01_Generate_reftableRF.bin_for_diyabc_generated_reftableRF_human_12000snp_maf_hudson_reftableRF_allS_11000PerScen.sh`: Shell script to generate bin reftable (pre-gof)
    * `02_Rscript_to_generate_parameter_files_for_each_scenario_separetly.R`: R script to translate bin to txt reftables (pre-gof)
    * `03_human_real_data_12000snp_maf_hudson_reftableRF_allS_11000PerScen_param.sh`: shell script to generate replicated reftables (pre-gof)
    * `04_Generate_reftableRF.bin_for_diyabc_generated_reftableRF_human_12000snp_maf_hudson_reftableRF_S2S3_110000PerScen`: Shell script to generate bin reftable (post-gof)
    * `05_Rscript_to_generate_parameter_files_for_each_scenario_separetly.R`: R script to translate bin to txt reftables (post-gof)
    * `06_human_real_data_12000snp_maf_hudson_reftableRF_S2S3only_110000PerScen_param.sh`: shell script to generate replicated reftables (post-gof)
    * `human_snp_all22chr_maf_hudson.snp` SNP data file
    * `headerRF_for_diyabc_generated_reftableRF_human_12000snp_maf_hudson_reftableRF_allS_11000PerScen.txt` header file for DIYABC (pre-gof)
    * `headerRF_for_diyabc_generated_reftableRF_human_12000snp_maf_hudson_reftableRF_S2S3_110000PerScen.txt` header file for DIYABC (pre-gof)

* `diyabc_resimulation_scripts`: Simulation scripts for diyabc
  * `sim_depindep.R`: simulation script for the "Dep-Indep" example
  * `depindep_resim`: header files for "Dep-Indep" simulation
  * `sim_humanfast.R`: simulation script for the "Human-like" example
  * `humanfast_resim`: header files for "Human-like" simulation
  * `sim_humanreal.R`: simulation script for the "Human" real example
  * `humanrealt_resim`: header files for "Human" real simulation

* `gof_sbi_paper.Rproj`: Root `Rproj` file.

* `human_analysis`: files for the empirical human data analysis.
  * `GOF_PRE_POST_HUMAN.R` R script to perform GOF analyses on the real human dataset
  
* `human_data`: data files for Human analysis
    * `statobsRF_maf_hudson_SNP_1_12000.txt` Observed dataset SNP positions 1 to 120000
    * `statobsRF_maf_hudson_SNP_12001_24000.txt` Observed dataset for last 12000 positions

* `simulation_study`: script to reproduce the simulation study.
  * `01_gof_simulaiton_study.R`: Script for simulation studies.
  * `02_gof_simulaiton_study_plots.R`: Script to plot the results.
  * `03_gof_lof_knn_example_plot.R`: Script to plot an example analysis.

## Requirements

* `DIYABC`: https://github.com/diyabc/diyab.
  Version 1.1.51 or above.
  For simulation of summary statistics from population genetics models.
  Pre-build releases are available from https://github.com/diyabc/diyabc/releases/tag/v1.1.51

* `R`: https://cran.r-project.org/index.html.
  Version 4.4. or above.
  All the `R` scripts assume that the working directory is the root directory
  (location of the `.Rproj` and `.here` files).

* `abcgof`: https://github.com/pbastide/abcgof
  Version v0.0.1.
  Can be installed from `R` with: `devtools::install_github("pbastide/abcgof", ref = "v0.0.1")`
