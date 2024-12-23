# Reference Table Simulation Scripts

Files to simulate reference tables from three example studies based on the analysis of simulated pseudo-observed (pods) datasets.
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

Files to simulate reference tables for the Human Like Real Example `humanreal` 
corresponding to pre-inference and post-inference GOF analyses of a real human SNP observed dataset.
Each folder contains all the files needed for the simulations.
  
For the pre-inference GOF analysis, run the scripts in the following order:
  * Run shell script `01_Generate_reftableRF.bin_for_diyabc_generated_reftableRF_human_12000snp_maf_hudson_reftableRF_allS_11000PerScen.sh` 
  to generate the reference table.
  * Run R script `02_Rscript_to_generate_parameter_files_for_each_scenario_separetly.R`
  to translate the `bin` file into parameter and sumstats text files.
  * Run `03_human_real_data_12000snp_maf_hudson_reftableRF_allS_11000PerScen_param.sh` 
  to generate replicates of the reference table.
  
For the post-inference GOF analysis, run the scripts in the following order:
  * Run shell script `04_Generate_reftableRF.bin_for_diyabc_generated_reftableRF_human_12000snp_maf_hudson_reftableRF_S2S3_110000PerScen.sh` 
  to generate the reference table.
  * Run R script `05_Rscript_to_generate_parameter_files_for_each_scenario_separetly.R`
  to translate the `bin` file into parameter and sumstats text files.
  * Run `06_human_real_data_12000snp_maf_hudson_reftableRF_S2S3only_110000PerScen_param.sh` 
  to generate replicates of the reference table.

> **Note** *for advanced users*: running the scripts
> `03_human_real_data_12000snp_maf_hudson_reftableRF_allS_11000PerScen_param.sh` and
> `06_human_real_data_12000snp_maf_hudson_reftableRF_S2S3only_110000PerScen_param.sh`
> can be avoided, because GOF pre-inference analyses do not use replicates of simulated points,
> and because the `abcgof` R package is able to produce "on the fly" replicate points from the selected points
> for the post-inference GOF.
> The input files can therfore more directly be the sumstats text files produced in the previous steps 
> based on the scripts
> `02_Rscript_to_generate_parameter_files_for_each_scenario_separetly.R` (pre-inference GOF) or
> `05_Rscript_to_generate_parameter_files_for_each_scenario_separetly.R` (post-inference GOF).
> In that case, the subsequent analysis using script
> `human_analysis/GOF_PRE_POST_HUMAN.R` should load the datafiles
> `human_12000snp_maf_hudson_reftableRF_allS_11000PerScen_param_S*.txt` instead of 
> `yellow_points_from_human_12000snp_maf_hudson_reftableRF_allS_11000PerScen_param_S*.txt`
> (pre-inference GOF),
> and the datafiles
> `human_12000snp_maf_hudson_reftableRF_S2S3_110000PerScen_param_S[2-3].txt` instead of
> `yellow_points_from_human_12000snp_maf_hudson_reftableRF_S2S3_110000PerScen_param_S[2-3].txt`.
> For practical reasons, we however used the scripts
> `03_human_real_data_12000snp_maf_hudson_reftableRF_allS_11000PerScen_param.sh` and
> `06_human_real_data_12000snp_maf_hudson_reftableRF_S2S3only_110000PerScen_param.sh`
> as described above in the computation processed in the associated manuscript.

## Requirements

Please see general README file for software requirement.
