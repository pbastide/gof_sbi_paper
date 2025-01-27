# Copyright (C) {2024} {GLM, PB, AE, JMM}
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

# Generate simulate new sumstats from a matrix of parameters
# This script needs to have diyabc installed on the machine
# Please refer to https://github.com/diyabc/diyabc for installation instructions
# Pre-build releases are available from https://github.com/diyabc/diyabc/releases
#
library(here)
#######################################################################################################################
## Functions
#######################################################################################################################

DIYABC_PATH <- "/path/to/diyabc/executable"
DIYABC_CMD <- paste0(DIYABC_PATH, " -p ./")

# function to match parameter orders
get_parameter_order <- function(scenario = paste0(1:6)) {
  #scenario <- paste0(scen.ref)
  param_order <- switch(scenario,
                        "1" = c("scenario", "N1", "N2", "N3", "N4", "t1", "t2", "d3", "Nbn3", "d4", "Nbn4", "N34", "t3", "d34", "Nbn34", "t4", "Na"),
                        "2" = c("scenario", "N1", "N2", "N3", "N4", "t1", "ra", "t2", "d3", "Nbn3", "d4", "Nbn4", "N34", "t3", "d34", "Nbn34", "t4", "Na"),
                        "3" = c("scenario", "N1", "N2", "N3", "N4", "t1", "ra", "t2", "d3", "Nbn3", "d4", "Nbn4", "N34", "t3", "d34", "Nbn34", "t4", "Na"),
                        "4" = c("scenario", "N1", "N2", "N3", "N4", "t11", "t22", "d3", "Nbn3", "t33", "d4", "Nbn4", "t44", "Na"),
                        "5" = c("scenario", "N1", "N2", "N3", "N4", "t11", "ra", "t22", "d3", "Nbn3", "t33", "d4", "Nbn4", "t44", "Na"),
                        "6" = c("scenario", "N1", "N2", "N3", "N4", "t11", "ra", "t22", "d3", "Nbn3", "t33", "d4", "Nbn4", "t44", "Na"))
  return(param_order)
}

# function to reorder parameters
reorder_param_table <- function(params_table, scenario) {
  new_params <- params_table[, match(get_parameter_order(paste0(scenario)), colnames(params_table))]
  return(new_params)
}

setup_sim <- function(path_to_headers, datestamp, seed, ncores) {
  cwd <- getwd()
  on.exit(setwd(cwd))
  path_to_sim <- file.path(cwd,"sim")
  if (!dir.exists(path_to_sim)) {
    dir.create(path_to_sim)
    setwd(path_to_sim)
    # copy headers and maf
    files_to_copy <- list.files(path_to_headers, full.names = TRUE)
    file.copy(files_to_copy, path_to_sim)
    # Generate RNG
    system(paste0(DIYABC_CMD, " -n \"t:", ncores, ";c:1;s:", seed, ";f:f\""))

  }
  return(path_to_sim)
}

simulate_human <- function(params_table, scenario, path_to_sim, ncores) {
  # re-order table
  new_params_table <- reorder_param_table(params_table, scenario)
  # tmpext <- sub("/", "", tempfile("", "", ""))
  # newparamname <- paste0("reordered_param_human_", tmpext)
  newparamname <- "reordered_param_human.txt"
  write.table(file = here(path_to_sim, newparamname), new_params_table, quote = F, col.names = T, row.names = F)
  #write.table(file = paste0(path_to_sim,"/", newparamname), new_params_table, quote = F, col.names = T, row.names = F)

  # do sim
  cwd <- getwd()
  on.exit(setwd(cwd))
  setwd(path_to_sim)
  # outputsumstats <- paste0("generated_sumstats_human.", tmpext)
  outputsumstats <- "generated_sumstats_human.txt"
  system(paste0(DIYABC_CMD, " -o ", newparamname, " -i ", outputsumstats, " -g 100 -m -t ", ncores))
  # read result
  new_sumstats <- read.table(outputsumstats, header = FALSE)
  # delete temp files
  unlink(c(newparamname, outputsumstats))
  return(new_sumstats)
}







# #######################################################################################################################
# ## Test
# #######################################################################################################################
# library(here)
#
# # Files names
# inputfile <- here("data", "dummy", "test_human_params_scenario_2.txt")
# path_to_headers <- here("data", "2024-03-18_human_resim")
# seed <- 1289
# ncores <- 4
#
# # setup simulations
# path_to_sim <- setup_sim(path_to_headers, seed, ncores)
#
# # parameters
# scenario <- 2
# params_table <- read.table(inputfile, header = TRUE)
#
# # sim
# new_sumstats <- simulate_human(params_table, scenario, path_to_sim, ncores)
#
# new_sumstats
