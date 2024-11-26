# Reference Table Simulation Scripts

Files to simulate reference tables from three examples.
Each folder contains all the files needed for the simulations.

## Toy Gaussian Laplace `gaussianlaplace`

R file `simulation_toy_gaussian_laplace.R`
simulates reference tables from the toy Laplace-Gaussian example.

## Dep-Indep Example `depindep`

Run the files in the following order:
  * Run shell script `01_Generate_reftableRF.bin_for_depindep_fast_NoConditionMediumMoins.sh` 
  to generate the reference table.
  * Run R script `02_Rscript_to_generate_parameter_files_for_each_scenario_separetly.R`
  to translate the `bin` file into parameter and sumstats text files.
  * Run `03_depindep_fast_NoConditionSmallT_reftableRF_110000_param.sh` 
  to generate replicates of the reference table.
  
## Human Like "Fast" Example `humanfast`

Run the files in the following order:
  * Run shell script `01_Generate_reftableRF.bin_for_diyabc_generated_reftableRF_human_fast.sh` 
  to generate the reference table.
  * Run R script `02_Rscript_to_generate_parameter_files_for_each_scenario_separetly.R`
  to translate the `bin` file into parameter and sumstats text files.
  * Run `03_human_fast_reftableRF_110000_param.sh` 
  to generate replicates of the reference table.
  
## Human Like Real Example `humanreal`

For the PRE-GOF analysis, run the scripts in the following order:
  * Run shell script `01_Generate_reftableRF.bin_for_diyabc_generated_reftableRF_human_12000snp_maf_hudson_reftableRF_allS_11000PerScen.sh` 
  to generate the reference table.
  * Run R script `02_Rscript_to_generate_parameter_files_for_each_scenario_separetly.R`
  to translate the `bin` file into parameter and sumstats text files.
  * Run `03_human_real_data_12000snp_maf_hudson_reftableRF_allS_11000PerScen_param.sh` 
  to generate replicates of the reference table.
  
For the POST-GOF analysis, run the scripts in the following order:
  * Run shell script `04_Generate_reftableRF.bin_for_diyabc_generated_reftableRF_human_12000snp_maf_hudson_reftableRF_S2S3_110000PerScen.sh` 
  to generate the reference table.
  * Run R script `05_Rscript_to_generate_parameter_files_for_each_scenario_separetly.R`
  to translate the `bin` file into parameter and sumstats text files.
  * Run `06_human_real_data_12000snp_maf_hudson_reftableRF_S2S3only_110000PerScen_param.sh` 
  to generate replicates of the reference table.

## Requirements

Please see general README file for software requirement.
