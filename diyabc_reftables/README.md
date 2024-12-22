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
  
Files to simulate reference tables for the Human Like Real Example `humanreal` corresponding to pre-inference and post-inference GOF analyses of a real human SNP observed dataset.
Each folder contains all the files needed for the simulations.
  
For the pre-inference GOF analysis, run the scripts in the following order:
  * Run shell script `01_Generate_reftableRF.bin_for_diyabc_generated_reftableRF_human_12000snp_maf_hudson_reftableRF_allS_11000PerScen.sh` 
  to generate the reference table.
  * Run R script `02_Rscript_to_generate_parameter_files_for_each_scenario_separetly.R`
  to translate the `bin` file into parameter and sumstats text files.
  * Run `03_human_real_data_12000snp_maf_hudson_reftableRF_allS_11000PerScen_param.sh` 
  to generate replicates of the reference table.
  
  ## REMARK: For advanced users (only), use of the script above `03_human_real_data_12000snp_maf_hudson_reftableRF_allS_11000PerScen_param.sh` can be avoided
  ## because GOF pre-inference analyses do not use replicates of simulated points. Therefore, for a pre-inference GOF analysis using the abcgof.R package, the input file 
  ## can (more directly) be the sumstats text files produced in the previous step based on the script `02_Rscript_to_generate_parameter_files_for_each_scenario_separetly.R`
  ## i.e. the human_12000snp_maf_hudson_reftableRF_allS_11000PerScen_param_S*.txt files, 
  ## instead of the files yellow_points_from_human_12000snp_maf_hudson_reftableRF_allS_11000PerScen_param_S*.txt. For practical reasons, we however used the script 
  ## `03_human_real_data_12000snp_maf_hudson_reftableRF_allS_11000PerScen_param.sh` in the computation processed in the associated manuscript.
  
For the post-inference GOF analysis, run the scripts in the following order:
  * Run shell script `04_Generate_reftableRF.bin_for_diyabc_generated_reftableRF_human_12000snp_maf_hudson_reftableRF_S2S3_110000PerScen.sh` 
  to generate the reference table.
  * Run R script `05_Rscript_to_generate_parameter_files_for_each_scenario_separetly.R`
  to translate the `bin` file into parameter and sumstats text files.
  * Run `06_human_real_data_12000snp_maf_hudson_reftableRF_S2S3only_110000PerScen_param.sh` 
  to generate replicates of the reference table.
  
  ## REMARK: For advanced users (only), use of the script above `06_human_real_data_12000snp_maf_hudson_reftableRF_S2S3only_110000PerScen_param.sh` can be avoided
   ## because for GOF post-inference analyses using the abcgof.R package, the package is able to produce ‘on the fly’ replicate (green) points from the selected (yellow) points. 
   ## The input file hence can (more directly be) human_12000snp_maf_hudson_reftableRF_S2S3_110000PerScen_param_S2.txt (or human_12000snp_maf_hudson_reftableRF_S2S3_110000PerScen_param_S3.txt)
   ## obtained by running the previous script `05_Rscript_to_generate_parameter_files_for_each_scenario_separetly.R`,
   ## instead of the file yellow_points_from_human_12000snp_maf_hudson_reftableRF_S2S3_110000PerScen_param_S2.txt (or yellow_points_from_human_12000snp_maf_hudson_reftableRF_S2S3_110000PerScen_param_S3).txt.
   ## For practical reasons, we however used the script `06_human_real_data_12000snp_maf_hudson_reftableRF_S2S3only_110000PerScen_param.sh` in the computation processed in the associated manuscript.

## Requirements

Please see general README file for software requirement.
