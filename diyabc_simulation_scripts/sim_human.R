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


#######################################################################################################################
## Functions
#######################################################################################################################
## Change here to the local path of your install of diyabc
DIYABC_PATH <- "/path/to/diyabc/executable"
DIYABC_PATH <- "/home/pbastide/warm/diyabc/build/src-JMC-C++/general"
DIYABC_CMD <- paste0(DIYABC_PATH, " -p ./")

#' @title Get parameter order
#'
#' @description
#' This function gives the correct parameter order for the Human-like example.
#'
#' @return The vector of ordered parameters
#'
get_parameter_order <- function(scenario = paste0(1:6)) {
  scenario <- match.arg(scenario)
  param_order <- switch(scenario,
                        "1" = c("scenario", "N1", "N2", "N3", "N4", "t1", "t2", "d3", "Nbn3", "d4", "Nbn4", "N34", "t3", "d34", "Nbn34", "t4", "Na"),
                        "2" = c("scenario", "N1", "N2", "N3", "N4", "t1", "ra", "t2", "d3", "Nbn3", "d4", "Nbn4", "N34", "t3", "d34", "Nbn34", "t4", "Na"),
                        "3" = c("scenario", "N1", "N2", "N3", "N4", "t1", "ra", "t2", "d3", "Nbn3", "d4", "Nbn4", "N34", "t3", "d34", "Nbn34", "t4", "Na"),
                        "4" = c("scenario", "N1", "N2", "N3", "N4", "t11", "t22", "d3", "Nbn3", "t33", "d4", "Nbn4", "t44", "Na"),
                        "5" = c("scenario", "N1", "N2", "N3", "N4", "t11", "ra", "t22", "d3", "Nbn3", "t33", "d4", "Nbn4", "t44", "Na"),
                        "6" = c("scenario", "N1", "N2", "N3", "N4", "t11", "ra", "t22", "d3", "Nbn3", "t33", "d4", "Nbn4", "t44", "Na"))
  return(param_order)
}

#' @title Re-order parameters
#'
#' @description
#' This function re-orders parameters to the correct order.
#'
#' @param params_table a matrix of parameters, with named columns.
#'
#' @return The re-ordered `params_table` matrix
#'
reorder_param_table <- function(params_table, scenario) {
  new_params <- params_table[, match(get_parameter_order(paste0(scenario)), colnames(params_table))]
  return(new_params)
}

#' @title Set-up the simulation environment
#'
#' @description
#' This function sets up the simulation environment to prepare for
#' re-simulation according to specified parameter values.
#'
#' @param path_to_headers local path to headers used in the simulation
#' @param datestamp string for date stamp
#' @param seed integer for setting the seed
#' @param ncores the number of cores to used in simulations
#'
#' @return `path_to_sim` the path to the set-up simulation environment.
#'
setup_sim <- function(path_to_headers, datestamp, seed, ncores) {
  cwd <- getwd()
  on.exit(setwd(cwd))
  path_to_sim <- file.path(cwd, "results", paste0(datestamp, "_sim_human"))
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

#' @title Simulate the Human-like model
#'
#' @description
#' This function simulates new summary statistics from the Human like model
#' for each line of the parameter matrix.
#'
#' @param params_table a matrix of parameters, with named columns.
#' @param path_to_sim the path to the set-up simulation environment.
#' @param ncores the number of cores to used in simulations.
#'
#' @return `new_sumstats` the simulated summary statistics.
#'
simulate_human <- function(params_table, scenario, path_to_sim, ncores) {
  # re-order table
  new_params_table <- reorder_param_table(params_table, scenario)
  tmpext <- sub("/", "", tempfile("", "", ""))
  newparamname <- paste0("reordered_param_human_", tmpext)
  write.table(file = here(path_to_sim, newparamname), new_params_table, quote = F, col.names = T, row.names = F)
  # do sim
  cwd <- getwd()
  on.exit(setwd(cwd))
  setwd(path_to_sim)
  outputsumstats <- paste0("generated_sumstats_human_", tmpext)
  system(paste0(DIYABC_CMD, " -o ", newparamname, " -i ", outputsumstats, " -g 100 -m -t ", ncores))
  # read result
  new_sumstats <- read.table(outputsumstats, header = FALSE)
  # delete temp files
  unlink(c(newparamname, outputsumstats))
  return(new_sumstats)
}
