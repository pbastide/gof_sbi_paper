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

# Script to run simulations to asses the power of the GoF tests.
#
# This script runs the pre-inference and post-inference GoF tests on 3 datasets:
# - The toy "gaussianlaplace" example;
# - The simple Population Genetics "depindep" example;
# - The more complex "humanfast" example.
# Computations can be intensive, especially for the post-inference GoF test in the Human-like example.
#
# REQUIREMENTS:
#   - DIYABC: version 1.1.54 or above. https://github.com/diyabc/diyabc/releases.
#   - R: version 4.4. or above.
#   - abcgof: version v0.0.1.
#
# REFERENCE TABLES:
#   Must be simulated using DIYABC.
#   Please see file `diyabc_reftables/README.md` for more details on the simulation procedure.
#
# NEW POINT SIMULATIONS:
# The post-inference tests require the simulation of new particles from the reference model.
#   - "gaussianlaplace": packages used in "diyabc_reftables/simulation_toy_gaussian_laplace.R" must be installed on the machine
#   - "depindep":
#     - Simulation script: "diyabc_resimulation_scripts/sim_depindep.R"
#     - Reference tables: "diyabc_resimulation_scripts/depindep_resim"
#   - "human":
#     - Simulation script: "diyabc_resimulation_scripts/humanfast.R"
#     - Reference tables: "diyabc_resimulation_scripts/humanfast_resim"
#
# OUTPUTS:
#   - Power and calibration results
#   - To be plotted by script `02_gof_simulaiton_study_plots.R`


################################################################################
## Load abcgof package
## Version v0.0.1 of abcgof must be used to run this script
devtools::install_github("pbastide/abcgof", ref = "v0.0.1")
library(abcgof)
################################################################################
## Seed for replication
myseed <- 1974
set.seed(myseed)
################################################################################
## Parallel computations
ncores <- 8
################################################################################
## Paths are relative to the root of the repository, containing the ".Rproj" file.
library(here)
PATH <- here()
datestamp_day <- format(Sys.time(), "%Y-%m-%d")

################################################################################
## Parameters
# number of pseudo-observations (pods) for power computation
n.test <- 1000
# test level for power computation
level <- 0.05
# total number of points in the reference dataset for the pre-inference GoFs
n.ref.prior.all <- c(500, 1000, 2000, 3000, 4000, 5000)
# total number of points in the reference dataset for the post-inference GoFs
n.ref.hp.all <- c(50*1000, 100*1000)
# number of points kept in the posterior for the post-inference GoFs
n.post.all <- c(500, 1000, 2000)
# proportion of points to use for calibration in the post-inference GoFs
split <- 0.5
# values of k on which to compute the kNN and LOF scores
k_all <- 1:20
# range of k values for the "max-LOF" score computation
k_range <- c(5, 20)

################################################################################
## Parameters
# gaussianlaplace example: must match the date of the "diyabc_reftables/simulation_toy_gaussian_laplace.R" script.
datestampe_gaussianlaplace <- "2024-03-19"

################################################################################
## utility function to add power and parameters to a result
add_power_parameters <- function(res, level, dataname, test_in_alternative, method, nref, npost) {
  ## Power
  res <- abcgof::compute_power(res, level = level)
  res <- abcgof::compute_power_sd(res)
  ## parameters
  res$dataset <- dataname
  res$test_in_alternative <- test_in_alternative
  res$method <- method
  res$nref <- nref
  res$npost <- npost
  return(res)
}

################################################################################
## Main Computation Loop
for (dataname in c("gaussianlaplace", "depindep", "human")) {
  ## In the human case, reduce number of pods to reduce the computation burden
  if (dataname == "human") n.test <- 500
  ## choose whether we compute the power of the test (test_in_alternative = TRUE),
  ## or if we check the calibration (test_in_alternative = FALSE)
  for (test_in_alternative in c(TRUE, FALSE)) {
    ################################################################################
    ################################################################################
    ## DATASETS
    ################################################################################
    ################################################################################

    ################################################################################
    ## Load toy gaussian laplace dataset
    if (dataname == "gaussianlaplace") {
      ## reference table from the Laplace model
      refname <- "laplace"
      ## PODs from the Gaussian model (power computation) or the same Laplace model (test calibration)
      testname <- ifelse(test_in_alternative, "gaussian", "laplace")

      ## Choose dataset
      lmomname <- "salmu"
      data_name <- paste0("2024-03-19", "_toy_gof_freq_lmom_", lmomname, "_ref_", refname, "_test_", testname)
      results_name <- paste0(datestamp_day, "_gof_resim_", dataname, "_lmom_", lmomname, "_ref_", refname, "_test_", testname)
      ## Load dataset
      all_data <- readRDS(here("diyabc_reftables", "gaussianlaplace", paste0(data_name, "_data.rds")))

      ## Format data
      # reference table
      data.ref <- all_data$data.ref
      # replicated reference table (for post-inference GoF tests)
      data.ref.replica <- all_data$data.ref.replica
      # PODs
      data.test <- all_data$data.test.obs
      # Replicated PODs (for post-inference GoF tests)
      data.test.replica <- all_data$data.test.new
      # Parameters associated with the reference table
      param.ref.all <- all_data$param.ref
      # Simulation function used for the reference table (for post-inference GoF tests)
      sim.fun.ref <- all_data$sim.fun.ref
      lmonfun <- all_data$lmonfun
      path_to_sim <- ""
      rm(all_data)

      ## transformation of parameters (for loclin and ridge regression)
      param_transform <- rep("logit", 2)
      param_lower_bound <- all_data$param_lower_bound
      param_upper_bound <- all_data$param_upper_bound
      if (refname == "laplace") {
        param_lower_bound[2] <- param_lower_bound[2] / sqrt(2)
        param_upper_bound[2] <- param_upper_bound[2] / sqrt(2)
      }
      names(param_transform) <- names(param_lower_bound) <- names(param_upper_bound) <- colnames(param.ref.all)
      trans_and_back <- abcgof:::check_param_transform(param.ref.all, param_transform, param_lower_bound, param_upper_bound)
    }

    ################################################################################
    ## Load Dep-Indep dataset
    if (dataname == "depindep") {
      ## reference table from scenario 2
      refname <- "scenario2"
      scen.ref <- 2
      ## PODs from the scenario 1 (power computation) or the same scenario 2 (test calibration)
      testname <- ifelse(test_in_alternative, "scenario1", "scenario2")
      scen.test <- ifelse(test_in_alternative, 1, 2)
      ## name of the result files
      results_name <- paste0(datestamp_day, "_gof_resim_", dataname, "_ref_", refname, "_test_", testname)

      ## Load dataset
      # names
      folder <- here("diyabc_reftables", "depindep")
      PATH.yellow.S1 <- file.path(folder, "yellow_points_from_depindep_fast_NoConditionMediumMoinsT_reftableRF_110000_param_S1.txt")
      PATH.yellow.S2 <- file.path(folder, "yellow_points_from_depindep_fast_NoConditionMediumMoinsT_reftableRF_110000_param_S2.txt")
      PATH.green.S1 <- file.path(folder, "green_points_from_depindep_fast_NoConditionMediumMoinsT_reftableRF_110000_param_S1.txt")
      PATH.green.S2 <- file.path(folder, "green_points_from_depindep_fast_NoConditionMediumMoinsT_reftableRF_110000_param_S2.txt")
      PATH.param.S1 <- file.path(folder, "depindep_fast_NoConditionMediumMoinsT_reftableRF_110000_param_S1.txt")
      PATH.param.S2 <- file.path(folder, "depindep_fast_NoConditionMediumMoinsT_reftableRF_110000_param_S2.txt")
      # load data
      data1.S1 <- as.matrix(read.table(PATH.yellow.S1, header = FALSE))
      data1.S2 <- as.matrix(read.table(PATH.yellow.S2, header = FALSE))
      data2.S1 <- as.matrix(read.table(PATH.green.S1, header = FALSE)) #replica
      data2.S2 <- as.matrix(read.table(PATH.green.S2, header = FALSE))
      params.S1 <- as.matrix(read.table(PATH.param.S1, header = TRUE))
      params.S2 <- as.matrix(read.table(PATH.param.S2, header = TRUE))
      # sanity checks
      stopifnot(all(dim(data1.S1) == dim(data2.S1)))
      stopifnot(all(dim(data1.S2) == dim(data2.S2)))
      stopifnot(all(nrow(data1.S1) == nrow(params.S1)))
      stopifnot(all(nrow(data1.S2) == nrow(params.S2)))

      ## re-draw lines at random
      sel_ind <- sample(1:101000, 101000, replace = FALSE)

      ## Choose reference table
      if (scen.ref == 2) {
        # reference table
        data.ref <- data1.S2[sel_ind[1:100000], 2:131]
        # replicated reference table (for post-inference GoF tests)
        data.ref.replica <- data2.S2[sel_ind[1:100000], 2:131]
        # Parameters associated with the reference table
        param.ref.all <- params.S2[sel_ind[1:100000], -1] # remove scenario
      } else if (scen.ref == 1) {
        # reference table
        data.ref <- data1.S1[sel_ind[1:100000], 2:131]
        # replicated reference table (for post-inference GoF tests)
        data.ref.replica <- data2.S1[sel_ind[1:100000], 2:131]
        # Parameters associated with the reference table
        param.ref.all <- params.S1[sel_ind[1:100000], -1] # remove scenario
      }
      ## Choose PODs
      if (scen.test == 2) {
        # PODs
        data.test <- data1.S2[sel_ind[100001:101000], 2:131]
        # Replicated PODs (for post-inference GoF tests)
        data.test.replica <- data2.S2[sel_ind[100001:101000], 2:131]
      } else if (scen.test == 1) {
        # PODs
        data.test <- data1.S1[sel_ind[100001:101000], 2:131]
        # Replicated PODs (for post-inference GoF tests)
        data.test.replica <- data2.S1[sel_ind[100001:101000], 2:131]
      }

      ## Simulation function used for the reference table (for post-inference GoF tests)
      # Source simulation function
      source("diyabc_resimulation_scripts/sim_depindep.R")
      path_to_headers <- here("diyabc_resimulation_scripts", "depindep_resim")
      seed <- 1289
      ncores_sim <- ncores
      datestamp_simu <- datestamp_day
      if (!test_in_alternative) datestamp_simu <- paste0(datestamp_simu, "_H0")
      path_to_sim <- setup_sim(path_to_headers, datestamp_simu, seed, ncores_sim)
      # Define function
      sim.fun.ref <- function(params, ncores_sim, path_to_sim) {
        # tmpseed <- round(runif(1, 1, 5000))
        # path_to_sim <- setup_sim(path_to_headers, datestamp_day, tmpseed, ncores_sim)
        params <- round(params)
        params <- cbind(scen.ref, params)
        colnames(params)[1] <- "scenario"
        new_sumstats <- simulate_depindep(params, path_to_sim, ncores_sim)
        # unlink(path_to_sim, recursive = TRUE)
        return(as.matrix(new_sumstats)[, -1])
      }

      ## transformation of parameters (for loclin and ridge regression)
      param_transform <- rep("logit", 7)
      param_lower_bound <- c(1000, 1000, 1000, 1000, 1, 61, 121) - 0.49
      param_upper_bound <- c(10000, 10000, 10000, 10000, 60, 120, 180) + 0.49
      names(param_transform) <- names(param_lower_bound) <- names(param_upper_bound) <- colnames(param.ref.all)
      trans_and_back <- abcgof:::check_param_transform(param.ref.all, param_transform, param_lower_bound, param_upper_bound)
    }

    ################################################################################
    ## Load Human-like dataset
    if (dataname == "human") {
      ## reference table from scenario 2
      refname <- "scenario2"
      scen.ref <- 2
      ## PODs from the scenario 3 (power computation) or the same scenario 2 (test calibration)
      testname <- ifelse(test_in_alternative, "scenario3", "scenario2")
      scen.test <- ifelse(test_in_alternative, 3, 2)
      ## name of the result files
      results_name <- paste0(datestamp_day, "_gof_resim_", dataname, "_ref_", refname, "_test_", testname)

      ## Load dataset
      # names
      folder <- here("diyabc_reftables", "humanfast")
      file.scen3 = "yellow_points_from_human_fast_reftableRF_110000_param_S3.txt"
      file.scen2 = "yellow_points_from_human_fast_reftableRF_110000_param_S2.txt"
      file.replica.scen3 = "green_points_from_human_fast_reftableRF_110000_param_S3.txt"
      file.replica.scen2 = "green_points_from_human_fast_reftableRF_110000_param_S2.txt"
      file.param.scen3 <- file.path(folder, "human_fast_reftableRF_110000_param_S3.txt")
      file.param.scen2 <- file.path(folder, "human_fast_reftableRF_110000_param_S2.txt")

      ## re-draw lines at random
      sel_ind <- sample(1:101000, 101000, replace = FALSE)

      ## Load data
      # reference table
      file.ref <- file.scen2
      data.ref <- as.matrix(read.table(file.path(folder, file.ref), header = FALSE))[sel_ind[1:100000], 2:131]
      # replicated reference table (for post-inference GoF tests)
      file.ref.rep <- file.replica.scen2
      data.ref.replica <- as.matrix(read.table(file.path(folder, file.ref.rep), header = FALSE))[sel_ind[1:100000], 2:131]
      # Parameters associated with the reference table
      file.param <- file.param.scen2
      param.ref.all <- as.matrix(read.table(file.param, header = TRUE))[sel_ind[1:100000], -1]
      # PODS
      file.test <- ifelse(test_in_alternative, file.scen3, file.scen2)
      data.test <- as.matrix(read.table(file.path(folder, file.test), header = FALSE))[sel_ind[100001:101000], 2:131]
      # Replicated PODs
      file.test.rep <- ifelse(test_in_alternative, file.replica.scen3, file.replica.scen2)
      data.test.replica <- as.matrix(read.table(file.path(folder, file.test.rep), header = FALSE))[sel_ind[100001:101000], 2:131]

      ## Simulation function used for the reference table (for post-inference GoF tests)
      # Source simulation function
      source("diyabc_resimulation_scripts/sim_humanfast")
      path_to_headers <- here("diyabc_resimulation_scripts", "humanfast_resim")
      seed <- 1289
      ncores_sim <- ncores
      datestamp_simu <- datestamp_day
      if (length(n.ref.hp.all) == 1) datestamp_simu <- paste0(datestamp_simu, "_", n.ref.hp.all)
      if (!test_in_alternative) datestamp_simu <- paste0(datestamp_simu, "_H0")
      path_to_sim <- setup_sim(path_to_headers, datestamp_simu, seed, ncores_sim)
      # Define function
      sim.fun.ref <- function(params, ncores_sim, path_to_sim) {
        params[, !grepl("ra", colnames(params))] <- round(params)[, !grepl("ra", colnames(params))]
        params <- cbind(scen.ref, params)
        colnames(params)[1] <- "scenario"
        new_sumstats <- simulate_human(params, scen.ref, path_to_sim, ncores_sim)
        return(as.matrix(new_sumstats)[, -1])
      }

      ## transformation of parameters (for loclin and ridge regression)
      # Parameter bounds
      N1_bounds <- c(1000.0,100000.0)
      N2_bounds <- c(1000.0,100000.0)
      N3_bounds <- c(1000.0,100000.0)
      N4_bounds <- c(1000.0,100000.0)
      t1_bounds <- c(1.0,30.0)
      t2_bounds <- c(100.0,10000.0)
      d3_bounds <- c(0.0,50.0)
      Nbn3_bounds <- c(5,500)
      d4_bounds <- c(0.0,50.0)
      Nbn4_bounds <- c(5,500)
      N34_bounds <- c(1000.0,100000.0)
      t3_bounds <- c(100.0,10000.0)
      d34_bounds <- c(0.0,50.0)
      Nbn34_bounds <- c(5,500)
      t4_bounds <- c(100.0,10000.0)
      Na_bounds <- c(100.0,10000.0)
      ra_bounds <- c(0.05,0.95)
      t11_bounds <- c(1.0,30.0)
      t22_bounds <- c(100.0,10000.0)
      t33_bounds <- c(100.0,10000.0)
      t44_bounds <- c(100.0,10000.0)
      param_bounds <- cbind(N1_bounds, N2_bounds, N3_bounds, N4_bounds,
                            t1_bounds, ra_bounds, t2_bounds,
                            d3_bounds, Nbn3_bounds,
                            d4_bounds, Nbn4_bounds,
                            N34_bounds, t3_bounds,
                            d34_bounds, Nbn34_bounds,
                            t4_bounds, Na_bounds)
      param_lower_bound <- param_bounds[1, ]
      param_upper_bound <- param_bounds[2, ]
      # Constraint that t4 > t3 > t2
      ind_times_order <- c(7, 13, 16)
      # Logit transform except for t4, t3, t2 (handled separately)
      param_transform <- rep("logit", 17)
      param_transform[ind_times_order] <- "none"
      names(param_transform) <- names(param_lower_bound) <- names(param_upper_bound) <- colnames(param.ref.all)
      # Add offset for integer parameters
      param_lower_bound[!grepl("ra", names(param_lower_bound))] <- param_lower_bound[!grepl("ra", names(param_lower_bound))] - 0.49
      param_upper_bound[!grepl("ra", names(param_upper_bound))] <- param_upper_bound[!grepl("ra", names(param_upper_bound))] + 0.49
      # parameter by parameter transform
      trans_and_back_no_order <- abcgof:::check_param_transform(param.ref.all, param_transform, param_lower_bound, param_upper_bound)
      # Handle constraint t4 > t3 > t2
      trans_order <- function(x, ind, alpha = 0.25) {
        y <- x
        for (ii in 2:length(ind)) {
          y[, ind[ii]] <- x[, ind[ii]] - 1 - x[, ind[ii-1]]
          y[y[, ind[ii]] == 0, ind[ii]] <- y[y[, ind[ii]] == 0, ind[ii]] + alpha
        }
        z <- y
        z[, ind[1]] <- abcgof:::logit(y[, ind[1]],
                                      param_lower_bound[ind[1]],
                                      param_upper_bound[ind[1]] - 2)
        for (ii in 2:length(ind)) {
          z[, ind[ii]] <- abcgof:::logitVec(y[, ind[ii]],
                                            0,
                                            param_upper_bound[ind[ii]] - (length(ind) - ii + 1) - x[, ind[ii-1]])
        }
        return(z)
      }
      trans_back <- function(z, ind, alpha = 0.25) {
        x <- y <- z
        y[, ind[1]] <- abcgof:::logistic(z[, ind[1]],
                                         param_lower_bound[ind[1]],
                                         param_upper_bound[ind[1]] - 2)
        x[, ind[1]] <- y[, ind[1]]
        for (ii in 2:length(ind)) {
          y[, ind[ii]] <- abcgof:::logisticVec(z[, ind[ii]],
                                               0,
                                               param_upper_bound[ind[ii]] - (length(ind) - ii + 1) - x[, ind[ii-1]])
          x[, ind[ii]] <- y[, ind[ii]] + x[, ind[ii-1]] + 1
        }
        # rounding
        x[, ind] <- round(x[, ind])
        return(x)
      }
      # Put together all transformations
      trans_par_order <- function(x) {
        y <- trans_and_back_no_order$transform(x)
        y <- trans_order(y, ind_times_order)
        colnames(y) <- colnames(x)
        return(y)
      }
      backtrans_par_order <- function(y) {
        x <- trans_back(y, ind_times_order)
        x <- trans_and_back_no_order$back_transform(x)
        colnames(x) <- colnames(y)
        return(x)
      }
      trans_and_back <- list(transform = trans_par_order,
                             back_transform = backtrans_par_order)
    }

    ################################################################################
    ################################################################################
    ## Normalisation
    ################################################################################
    ################################################################################

    # normalize data by standard deviation
    norm <- sd
    data_norm <- norm_stats(rbind(data.test, data.test.replica, data.ref.replica),
                            data.ref,
                            norm = norm)

    # save normalized data
    data.test.normalized.all <- as.matrix(data_norm$target[1:nrow(data.test), ])
    data.test.replica.normalized.all <- as.matrix(data_norm$target[(nrow(data.test)+1):(nrow(data.test)+nrow(data.test.replica)), ])
    data.ref.replica.normalized.all <- as.matrix(data_norm$target[(nrow(data.test)+nrow(data.test.replica)+1):(nrow(data.test)+nrow(data.test.replica)+nrow(data.ref.replica)), ])
    data.ref.normalized.all <- as.matrix(data_norm$sumstat)
    data.ref.normvec <- data_norm$normvec
    # delete unormalized data
    rm(data.test)
    rm(data_norm)
    rm(data.ref.replica)
    rm(data.ref)

    ################################################################################
    ################################################################################
    # Prior GOF
    ################################################################################
    ################################################################################

    for (n.ref.prior in n.ref.prior.all) {
      set.seed(1289)
      ## file management
      results_name_nref <- paste0(results_name, "_nrep_", n.ref.prior)
      res_file <-  here("results", paste0(results_name_nref, "_prior.rds"))
      if (!file.exists(res_file)) { # only run if result file does not exist
        ## Subset of the data for reference
        ind_ref <- sample(1:nrow(data.ref.normalized.all), n.ref.prior, replace = FALSE)
        data.ref.normalized <- data.ref.normalized.all[ind_ref, ]
        data.ref.replica.normalized <- data.ref.replica.normalized.all[ind_ref, ]
        param.ref <- param.ref.all[ind_ref, ]
        ## Subset of the data for test (pods)
        ind_test <- sample(1:nrow(data.test.normalized.all), n.test, replace = FALSE)
        data.test.normalized <- data.test.normalized.all[ind_test, ]
        ## Split reference table between calibration and reference
        n.calib.prior <- floor(split * n.ref.prior)
        calibid <- sample(seq_len(nrow(data.ref.normalized)), size = n.calib.prior)
        data.norm.calib.prior <- data.ref.normalized[calibid, ]
        data.norm.ref.prior <- data.ref.normalized[-calibid, ]
        data.norm.ref.rep.prior <- data.ref.replica.normalized[-calibid, ] # replicate for post gof
        ## Precomputation of distances for improved speed
        Tab.dbscan <- dbscan::kNN(data.norm.ref.prior, max(k_all)) #Pre-computation
        dist <- Tab.dbscan$dist
        id <- Tab.dbscan$id
        ## Prior GoF
        gof.prior <- abcgof::gof.fit(data.test = data.test.normalized,
                                     data.ref = data.norm.ref.prior,
                                     data.calib = data.norm.calib.prior,
                                     k = k_all,
                                     k_range = k_range,
                                     score = c("lof", "kNN"),
                                     ncores = ncores,
                                     kdist = dist, kid = id)
        ## Power computation
        gof.prior <- add_power_parameters(gof.prior, level, dataname, test_in_alternative, "prior", n.ref.prior, n.ref.prior)
        ## save the results
        saveRDS(gof.prior, file = res_file)
      }
    }

    ################################################################################
    ################################################################################
    # Localized and Post-Inference GOF
    ################################################################################
    ################################################################################

    ################################################################################
    # Loop over n.ref
    for (n.ref in n.ref.hp.all) {
      set.seed(1289)
      results_name_nref <- paste0(results_name, "_nrep_", n.ref)

      ################################################################################
      # Subset of the data to keep only n.ref lines
      ind_ref <- sample(1:nrow(data.ref.normalized.all), n.ref, replace = FALSE)
      data.ref.normalized <- data.ref.normalized.all[ind_ref, ]
      data.ref.replica.normalized <- data.ref.replica.normalized.all[ind_ref, ]
      param.ref <- param.ref.all[ind_ref, ]

      ind_test <- sample(1:nrow(data.test.normalized.all), n.test, replace = FALSE)
      data.test.normalized <- data.test.normalized.all[ind_test, ]
      data.test.replica.normalized <- data.test.replica.normalized.all[ind_test, ]

      gof.post <- NULL
      gof.freq <- NULL
      gof.hp.rej <- NULL
      gof.hp.loclin <- NULL
      gof.hp.ridge <- NULL

      ####################################################################
      ## loop over n.post the number of particles kept after localization

      for (n.post in n.post.all) {
        eps <- n.post / n.ref
        results_name_nref_eps <- paste0(results_name_nref, "_eps_", eps*100)

        ################################################################################
        ### Localized prior GOF
        ################################################################################

        res_file <-  here("results", paste0(results_name_nref_eps, "_post.rds"))

        if (!file.exists(res_file)) { # only run if result file does not exist
          # Select calibration points
          n.calib.even <- n.ref * eps * split
          calib.ref.id <- sample(seq_len(nrow(data.ref.normalized)), size = n.calib.even, replace = FALSE)
          data.norm.calib.prior.even <- data.ref.normalized[calib.ref.id, ]
          data.norm.ref.post.even <- data.ref.normalized[-calib.ref.id, ]
          data.norm.ref.rep.post.even <- data.ref.replica.normalized[-calib.ref.id, ]
          # find epsilon such that the number of points kept in the posterior is n.ref.even
          # ie  n.ref.even == eps_post * (n.ref - n.calib.even)
          eps_post <- eps * split / (1 - eps * split)

          gof.post <- abcgof::gof.fit.post(data.test = data.test.normalized,
                                           data.calib = data.norm.calib.prior.even,
                                           data.ref = data.norm.ref.post.even,
                                           data.ref.replica = data.norm.ref.rep.post.even,
                                           k = k_all,
                                           k_range = k_range,
                                           eps = eps_post,
                                           score = c("lof", "kNN"),
                                           ncores = ncores)
          # Power
          gof.post <- add_power_parameters(gof.post, level, dataname, test_in_alternative, "post", n.ref, n.post)
          # save
          saveRDS(gof.post, file = res_file)
        }

        ################################################################################
        ### Pre-computed post inference holdout GOF
        ################################################################################

        res_file <-  here("results", paste0(results_name_nref_eps, "_freq.rds"))

        if (!file.exists(res_file)) { # only run if result file does not exist
          gof.freq <- abcgof::gof.fit.freq(data.test.obs = data.test.normalized,
                                           data.test.new = data.test.replica.normalized,
                                           data.ref = data.ref.normalized,
                                           data.ref.replica = data.ref.replica.normalized,
                                           k = k_all,
                                           k_range = k_range,
                                           eps = eps,
                                           split = split,
                                           score = c("lof", "kNN"),
                                           ncores = ncores)
          # Power
          gof.freq <- add_power_parameters(gof.freq, level, dataname, test_in_alternative, "freq", n.ref, n.post)
          # save
          saveRDS(gof.freq, file = res_file)
        }

        ################################################################################
        ### Holdout Posterior GOF with re-simulations
        ################################################################################

        res_file <-  here("results", paste0(results_name_nref_eps, "_hp_rej.rds"))

        if (!file.exists(res_file)) { # only run if result file does not exist
          gof.hp.rej <- abcgof::gof.fit.holdout(data.test.obs = data.test.normalized,
                                                data.test.new = data.test.replica.normalized,
                                                data.param = param.ref,
                                                data.ref = data.ref.normalized,
                                                sim.fun = sim.fun.ref,
                                                method = "rejection",
                                                k = k_all,
                                                k_range = k_range,
                                                eps = eps,
                                                split = split,
                                                normvec = data.ref.normvec,
                                                score = c("lof", "kNN"),
                                                ncores = ifelse(dataname == "gaussianlaplace", ncores, 1),
                                                ncores_sim = ncores,
                                                path_to_sim = path_to_sim)
          # Power
          gof.hp.rej <- add_power_parameters(gof.hp.rej, level, dataname, test_in_alternative, "hp.rej", n.ref, n.post)
          # save
          saveRDS(gof.hp.rej, file = res_file)
        }

        ################################################################################
        ### Holdout Posterior GOF - with re-simulations and local linear adjustment
        ################################################################################

        if (n.ref * eps > ncol(data.test.normalized) + 1) { # only if enough points for the regression
          res_file <-  here("results", paste0(results_name_nref_eps, "_hp_loclin.rds"))

          if (!file.exists(res_file)) { # only run if result file does not exist
            gof.hp.loclin <- abcgof::gof.fit.holdout(data.test.obs = data.test.normalized,
                                                     data.test.new = data.test.replica.normalized,
                                                     data.param = param.ref,
                                                     data.ref = data.ref.normalized,
                                                     sim.fun = sim.fun.ref,
                                                     method = "loclinear",
                                                     kernel = "epanechnikov",
                                                     trans.fun = trans_and_back$transform,
                                                     back.trans.fun = trans_and_back$back_transform,
                                                     k = k_all,
                                                     k_range = k_range,
                                                     eps = eps,
                                                     split = split,
                                                     normvec = data.ref.normvec,
                                                     score = c("lof", "kNN"),
                                                     ncores = ifelse(dataname == "gaussianlaplace", ncores, 1),
                                                     ncores_sim = ncores,
                                                     path_to_sim = path_to_sim)
            # Power
            gof.hp.loclin <- add_power_parameters(gof.hp.loclin, level, dataname, test_in_alternative, "hp.loclin", n.ref, n.post)
            # save
            saveRDS(gof.hp.loclin, file = res_file)
          }
        }

        ################################################################################
        ### Holdout Posterior GOF with re-simulations and ridge linear adjustment
        ################################################################################

        res_file <-  here("results", paste0(results_name_nref_eps, "_hp_ridge.rds"))

        if (!file.exists(res_file)) { # only run if result file does not exist
          gof.hp.ridge <- abcgof::gof.fit.holdout(data.test.obs = data.test.normalized,
                                                  data.test.new = data.test.replica.normalized,
                                                  data.param = param.ref,
                                                  data.ref = data.ref.normalized,
                                                  sim.fun = sim.fun.ref,
                                                  method = "ridge",
                                                  kernel = "epanechnikov",
                                                  trans.fun = trans_and_back$transform,
                                                  back.trans.fun = trans_and_back$back_transform,
                                                  k = k_all,
                                                  k_range = k_range,
                                                  eps = eps,
                                                  split = split,
                                                  normvec = data.ref.normvec,
                                                  score = c("lof", "kNN"),
                                                  ncores = ifelse(dataname == "gaussianlaplace", ncores, 1),
                                                  ncores_sim = ncores,
                                                  path_to_sim = path_to_sim)
          # Power
          gof.hp.ridge <- add_power_parameters(gof.hp.ridge, level, dataname, test_in_alternative, "hp.ridge", n.ref, n.post)
          # save
          saveRDS(gof.hp.ridge, file = res_file)
        }

        ################################################################################
        ### Holdout Posterior GOF with re-simulations - not using the replicate of the data
        ################################################################################

        res_file <-  here("results", paste0(results_name_nref_eps, "_hp_rej_norep.rds"))

        if (!file.exists(res_file)) { # only run if result file does not exist
          gof.hp.rej.norep <- abcgof::gof.fit.holdout(data.test.obs = data.test.normalized,
                                                      data.test.new = data.test.normalized,
                                                      data.param = param.ref,
                                                      data.ref = data.ref.normalized,
                                                      sim.fun = sim.fun.ref,
                                                      method = "rejection",
                                                      k = k_all,
                                                      k_range = k_range,
                                                      eps = eps,
                                                      split = split,
                                                      normvec = data.ref.normvec,
                                                      score = c("lof", "kNN"),
                                                      ncores = ifelse(dataname == "gaussianlaplace", ncores, 1),
                                                      ncores_sim = ncores,
                                                      path_to_sim = path_to_sim)
          # Power
          gof.hp.rej.norep <- add_power_parameters(gof.hp.rej.norep, level, dataname, test_in_alternative, "hp.rej.norep", n.ref, n.post)
          # save
          saveRDS(gof.hp.rej.norep, file = res_file)
        }

        ################################################################################
        ### Holdout Posterior GOF - with resimulations and local linear adj - norep
        ################################################################################

        if (n.ref * eps > ncol(data.test.normalized) + 1) { # only if enough points for the regression
          res_file <-  here("results", paste0(results_name_nref_eps, "_hp_loclin_norep.rds"))

          if (!file.exists(res_file)) { # only run if result file does not exist
            gof.hp.loclin.norep <- abcgof::gof.fit.holdout(data.test.obs = data.test.normalized,
                                                           data.test.new = data.test.normalized,
                                                           data.param = param.ref,
                                                           data.ref = data.ref.normalized,
                                                           sim.fun = sim.fun.ref,
                                                           method = "loclinear",
                                                           kernel = "epanechnikov",
                                                           trans.fun = trans_and_back$transform,
                                                           back.trans.fun = trans_and_back$back_transform,
                                                           k = k_all,
                                                           k_range = k_range,
                                                           eps = eps,
                                                           split = split,
                                                           normvec = data.ref.normvec,
                                                           score = c("lof", "kNN"),
                                                           ncores = ifelse(dataname == "gaussianlaplace", ncores, 1),
                                                           ncores_sim = ncores,
                                                           path_to_sim = path_to_sim)
            # Power
            gof.hp.loclin.norep <- add_power_parameters(gof.hp.loclin.norep, level, dataname, test_in_alternative, "hp.loclin.norep", n.ref, n.post)
            # save
            saveRDS(gof.hp.loclin.norep, file = res_file)
          }
        }

        ################################################################################
        ### Holdout Posterior GOF - with resimulations and local linear adj ridge - norep
        ################################################################################

        res_file <-  here("results", paste0(results_name_nref_eps, "_hp_ridge_norep.rds"))

        if (!file.exists(res_file)) { # only run if result file does not exist
          gof.hp.ridge.norep <- abcgof::gof.fit.holdout(data.test.obs = data.test.normalized,
                                                        data.test.new = data.test.normalized,
                                                        data.param = param.ref,
                                                        data.ref = data.ref.normalized,
                                                        sim.fun = sim.fun.ref,
                                                        method = "ridge",
                                                        kernel = "epanechnikov",
                                                        trans.fun = trans_and_back$transform,
                                                        back.trans.fun = trans_and_back$back_transform,
                                                        k = k_all,
                                                        k_range = k_range,
                                                        eps = eps,
                                                        split = split,
                                                        normvec = data.ref.normvec,
                                                        score = c("lof", "kNN"),
                                                        ncores = ifelse(dataname == "gaussianlaplace", ncores, 1),
                                                        ncores_sim = ncores,
                                                        path_to_sim = path_to_sim)
          # Power
          gof.hp.ridge.norep <- add_power_parameters(gof.hp.ridge.norep, level, dataname, test_in_alternative, "hp.ridge.norep", n.ref, n.post)
          # save
          saveRDS(gof.hp.ridge.norep, file = res_file)
        }
      }
    }
  }
}
