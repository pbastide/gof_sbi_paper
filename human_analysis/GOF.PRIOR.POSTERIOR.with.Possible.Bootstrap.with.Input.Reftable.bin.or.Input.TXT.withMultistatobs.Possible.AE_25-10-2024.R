# Copyright (C) {2024} {GLM, PB, JMM, AE}
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
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
# This script runs the pre-inference and post-inference GoF tests on the Human empirical example
#
# REQUIREMENTS:
#   - DIYABC: version 1.1.54 or above. https://github.com/diyabc/diyabc/releases.
#   - R: version 4.4. or above.
#   - abcgof: version v0.0.1. (see below)
#
# DATA:
#   - TODO
#
# NEW POINT SIMULATIONS:
# The post-inference tests require the simulation of new particles from the reference model.
#   - "human": TODO
#
# OUTPUTS:
# Pre and post inference GoF results

################################################################################
library(abcgof)
library(abcrf)
################################################################################
# Fonction to convert times in seconds to time in HH:MM:SS
convert_time <- function(time_in_seconds) {
  hours <- floor(time_in_seconds / 3600)
  minutes <- floor((time_in_seconds %% 3600) / 60)
  seconds <- (time_in_seconds %% 3600) %% 60
  return(sprintf("%02d:%02d:%05.2f", hours, minutes, seconds))
}
################################################################################
## location of data, simulation and results
swd = getwd() # swd = source work directory
folder <- swd
################################################################################
# For controling starting point and bootsrapping
myseed <- 1162
set.seed(myseed)
################################################################################
n.boot=100 #
################################################################################
# Parallel computations
# ncores for Prior GoF
ncores.prior = 32
# ncores for post-inference GoF
ncores.posterior = 1
# ncores for simulations with DIYABC
ncores_sim = 32
################################################################################
# Which type of analysis ?
GOF.PRIOR = TRUE
GOF.POSTERIOR = TRUE
post.rej = TRUE
post.rej.loclin = FALSE
post.rej.ridge = FALSE
################################################################################
# Parameters
n.ref.prior.all <- c(500, 1000, 2000, 3000, 4000, 5000) # total number of points in the reference dataset for the post and freq gof
n.ref.hp.all <- c(50*1000, 100*1000) # total number of points in the reference dataset for the post and freq gof
n.post.all <- c(500, 1000, 2000) # number of points kept in the posterior # On ne parle plus de epsilon mais equivalent
split <- 0.5 # split of the posterior for the frequentist test
k_all <- 1:20 # k vector
k_range <- c(5, 20) #### range effectif continue de k values

################################################################################
# Output
datestamp_day <- "2024-08-03"
res.file.name.prior <- paste0("empirical_human_results_", datestamp_day, "_PRIOR.txt")
res.file.name.posterior = paste0("empirical_human_results_", datestamp_day, "_POSTERIOR_rej.txt")
################################################################################

#####################################################################################
########## DATA #####################################################################
#####################################################################################
## Human real data files

file.dataobs.replica = "statobsRF_maf_hudson_SNP_1_12000.txt"
file.dataobs = "statobsRF_maf_hudson_SNP_12001_24000.txt" # For "GOF posterior" only
data.obs <- as.matrix(read.table(file.dataobs,header=TRUE), drop = FALSE)
data.obs.replica <- as.matrix(read.table(file.dataobs.replica,header=TRUE), drop = FALSE)

# For GOF PRIOR analysis
n.scenario.prior.gof = 6
prior.gof.file.ref = "yellow_points_from_human_12000snp_maf_hudson_reftableRF_allS_11000PerScen_param_S"
#prior.gof.file.ref = "green_points_from_human_12000snp_maf_hudson_reftableRF_allS_11000PerScen_param_S"

# For any GOF POSTERIOR analysis
scen.ref <- 2 # ID number of the analysed scenario
scenario <- paste0(scen.ref)
posterior.gof.file.ref = "yellow_points_from_human_12000snp_maf_hudson_reftableRF_S2S3_110000PerScen_param_S2.txt"
#posterior.gof.file.ref = "green_points_from_human_12000snp_maf_hudson_reftableRF_S2S3_110000PerScen_param_S2.txt"
posterior.gof.file.param <- "human_12000snp_maf_hudson_reftableRF_S2S3_110000PerScen_param_S2.txt" # Pas optimal
param.ref.all <- as.matrix(read.table(file.path(folder, posterior.gof.file.param), header = TRUE))[1:100000, -1]
#############################################################################################

################################################################################
################################################################################
# START OF MAIN COMPUTATION SECTION
################################################################################
################################################################################

if (GOF.PRIOR==TRUE) {
  ################################################################################
  ################################################################################
  # Prior GOF
  ################################################################################
  ################################################################################
  sink(file = res.file.name.prior)
  cat("#################### PRIOR GOF #################","\n")
  cat("myseed =",myseed,"\n")
  cat("k_all =",k_all,"\n")
  cat("k_range =",k_range,"\n")
  cat("\n")
  cat("Name of the reference table:", name.refTable, "\n")
  cat("Name of the observed dataset:", name.observed.dataset, "\n")
  cat("Number of simulations loaded from the reference table:", refTable$nrec, "\n")
  cat("Number of scenarios (i.e. models) in the reference table:", refTable$nscen, "\n")
  cat("Number of simulations available for each scenario from the loaded reference table:", refTable$nrecscen, "\n")
  cat("Number of parameters recovered from the reference table:", refTable$nparam, "\n")
  cat("Number of summary statistics in the reference table:", ncol(refTable$stats), "\n")
  cat("\n")
  sink()
  n.scenario.prior.gof=N.scen

  for (i.scen in 1:n.scenario.prior.gof) {

    ##### ALA HUMAN reftable txt JAS #######
    prior.gof.file.ref.i.scen = paste0(prior.gof.file.ref,i.scen,".txt")
    #prior.gof.file.ref.i.scen.replica = paste0(prior.gof.file.ref.replica,i.scen,".txt")
    #prior.gof.file.param.i.scen = paste0(prior.gof.file.param,i.scen,".txt")
    data.ref <- as.matrix(read.table(file.path(folder, prior.gof.file.ref.i.scen), header = FALSE))[1:10000, 2:131]
    #data.ref.replica <- as.matrix(read.table(file.path(folder, prior.gof.file.ref.i.scen.replica), header = FALSE))[1:10000, 2:131]

    cat("\n")
    cat("##########################################################","\n")
    cat("PRIOR GOF analyzed scenario =",i.scen,"\n")
    cat("\n")

    sink(file = res.file.name.prior, append = TRUE)
    cat("\n")
    cat("##########################################################","\n")
    cat("PRIOR GOF analyzed scenario =",i.scen,"\n")
    sink()

    for (n.ref.prior in n.ref.prior.all) {
      #res_file <-  here("results", paste0(results_name_nref, "_prior.rds"))
      data.ref.n.ref.prior <- data.ref[1:n.ref.prior,]
      #param.ref = param.ref.all[1:n.ref.prior,]

      ## GOF computation
      # Capturer le temps de début
      start_time <- proc.time()

      gof.prior <- gfit(target=data.obs,
                        sumstat=data.ref.n.ref.prior,
                        nb.replicate = n.ref.prior*split,
                        score = c("lof", "kNN"),
                        k = k_all,
                        k_range = k_range,
                        norm = sd,
                        ncores = ncores.prior,
                        nboot = n.boot)

      # Capturer le temps de fin
      end_time <- proc.time()
      # Calculer la durée d'exécution
      execution_time <- end_time - start_time

      # Convertir le temps écoulé en heures, minutes et secondes
      execution_time_converted <- convert_time(execution_time["elapsed"])

      summary.lof = summary(gof.prior, score = "lof", k = "max", level = 0.95)
      summary.kNN = summary(gof.prior, score = "kNN", k = 1, level = 0.95)

      # Capture des résumés sous forme de vecteurs de chaînes de caractères
      summary.lof_text <- capture.output(print(summary.lof))
      summary.kNN_text <- capture.output(print(summary.kNN))

      cat("--------------> n.sim =", n.ref.prior,"- n.calib =", n.ref.prior*split,"\n")
      print(execution_time_converted)
      cat("\n")
      print(summary.lof)
      cat("\n")
      print(summary.kNN)
      cat("\n")

      # Écrire le contenu combiné dans un fichier txt
      sink(file = res.file.name.prior, append = TRUE )
      cat("--------------> n.sim =", n.ref.prior,"- n.calib =", n.ref.prior*split,"\n")
      print(execution_time_converted)
      cat("\n")
      print(summary.lof)
      cat("\n")
      print(summary.kNN)
      cat("\n")
      sink()
    }
  }
}

if (GOF.POSTERIOR ==TRUE) {
  ################################################################################
  ################################################################################
  # Posterior GOF (rej, loclin et ridge)
  ################################################################################
  ################################################################################

  sink(file = res.file.name.posterior)

  cat("#################### POSTERIOR GOF - analysed scenario ",scen.ref,"#################","\n")
  cat("myseed =",myseed,"\n")
  cat("k_all =",k_all,"\n")
  cat("k_range =",k_range,"\n")
  cat("\n")
  sink()

  data.ref <- as.matrix(read.table(file.path(folder, posterior.gof.file.ref), header = FALSE))[1:100000, 2:131]
  #data.ref.replica <- as.matrix(read.table(file.path(folder, posterior.gof.file.ref.i.scen.replica), header = FALSE))[1:10000, 2:131]
  colnames(data.ref) = colnames(data.obs)
  #colnames(data.ref.replica) = colnames(data.obs)
  param.ref.all <- as.matrix(read.table(file.path(folder, posterior.gof.file.param.i.scen), header = TRUE))[1:100000, -1]


  ################################################################################
  ################################################################################
  # DEBUT PREPARATION OF SIMULATION
  ################################################################################
  ################################################################################
  # Prepare simulation function
  source("sim_human_AE-01-05-2024.R")
  path_to_headers <- swd

  seed <- myseed
  ncores_sim <- ncores_sim
  datestamp_simu <- datestamp_day
  path_to_sim <- setup_sim(path_to_headers, datestamp_simu, myseed, ncores_sim)

  #  Ajouter scen.ref et rep.run comme paramètres pour la fonction, afin que leur utilisation soit explicite et contrôlée de manière externe.
  sim.fun.ref <- function(params, ncores_sim, path_to_sim, scen.ref) {
    # Liste des motifs à exclure pour arrondir (colnames qui contiennent "raa", "µmic_1", "pmic_1", ou "snimic_1")
    cols_to_exclude <- c("raa", "µmic_1", "pmic_1", "snimic_1")
    # Appliquer !grepl pour exclure les colonnes correspondant aux motifs ci-dessus
    cols_to_round <- !grepl(paste(cols_to_exclude, collapse = "|"), colnames(params))
    params[, cols_to_round] <- round(params[, cols_to_round])
    # Ajouter une colonne de scénario référence au début des paramètres
    params <- cbind(scenario = scen.ref, params)
    # Simuler les statistiques sommaires
    new_sumstats <- simulate_human(params, scen.ref, path_to_sim, ncores_sim)
    # Retourner les statistiques sommaires en excluant la première colonne
    return(as.matrix(new_sumstats)[, -1])
  }

  ########### Transformation of parameters (si regression et si human!!!) ###############
  # All param: any order
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

  # Param dans "le bon ordre"
  param_bounds <- cbind(N1_bounds, N2_bounds, N3_bounds, N4_bounds,
                        t1_bounds, ra_bounds, t2_bounds,
                        d3_bounds, Nbn3_bounds,
                        d4_bounds, Nbn4_bounds,
                        N34_bounds, t3_bounds,
                        d34_bounds, Nbn34_bounds,
                        t4_bounds, Na_bounds)

  param_lower_bound <- param_bounds[1, ]
  param_upper_bound <- param_bounds[2, ]

  ind_times_order <- c(7, 13, 16) # necessaire car t4<t3<t2
  param_transform <- rep("logit", 17)
  param_transform[ind_times_order] <- "none"

  names(param_transform) <- names(param_lower_bound) <- names(param_upper_bound) <- colnames(param.ref.all) # Tous le meme nom

  # Pour les parametres entier (pas du type ra) = un peu avant ou apres pour pouvoir avoir inf et sup apres logit
  param_lower_bound[!grepl("ra", names(param_lower_bound))] <- param_lower_bound[!grepl("ra", names(param_lower_bound))] - 0.49
  param_upper_bound[!grepl("ra", names(param_upper_bound))] <- param_upper_bound[!grepl("ra", names(param_upper_bound))] + 0.49

  trans_and_back_no_order <- gofabcpkg:::check_param_transform(param.ref.all, param_transform, param_lower_bound, param_upper_bound)

  trans_order <- function(x, ind, alpha = 0.25) {
    y <- x
    for (ii in 2:length(ind)) {
      y[, ind[ii]] <- x[, ind[ii]] - 1 - x[, ind[ii-1]]
      y[y[, ind[ii]] == 0, ind[ii]] <- y[y[, ind[ii]] == 0, ind[ii]] + alpha
    }
    z <- y
    z[, ind[1]] <- gofabcpkg:::logit(y[, ind[1]],
                                     param_lower_bound[ind[1]],
                                     param_upper_bound[ind[1]] - 2)
    for (ii in 2:length(ind)) {
      z[, ind[ii]] <- gofabcpkg:::logitVec(y[, ind[ii]],
                                           0,
                                           param_upper_bound[ind[ii]] - (length(ind) - ii + 1) - x[, ind[ii-1]])
    }
    return(z)
  }
  trans_back <- function(z, ind, alpha = 0.25) {
    x <- y <- z
    y[, ind[1]] <- gofabcpkg:::logistic(z[, ind[1]],
                                        param_lower_bound[ind[1]],
                                        param_upper_bound[ind[1]] - 2)
    x[, ind[1]] <- y[, ind[1]]
    for (ii in 2:length(ind)) {
      y[, ind[ii]] <- gofabcpkg:::logisticVec(z[, ind[ii]],
                                              0,
                                              param_upper_bound[ind[ii]] - (length(ind) - ii + 1) - x[, ind[ii-1]])
      x[, ind[ii]] <- y[, ind[ii]] + x[, ind[ii-1]] + 1
    }
    # rounding
    x[, ind] <- round(x[, ind])
    return(x)
  }
  colnames(param.ref.all)[ind_times_order]
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


  ##################################


  ################################################################################
  ################################################################################
  # END PREPARATION OF SIMULATION
  ################################################################################
  ################################################################################

  ################################################################################
  ################################################################################
  # OUTPUT
  ################################################################################
  ################################################################################

  sink(file = res.file.name.posterior, append = TRUE)
  cat("#####################################","\n")
  cat("Analysed scenario = S",scen.ref,"\n")
  cat("Reftable file name =",posterior.gof.file.ref,"\n")
  #cat("refscen.replica =",posterior.gof.file.ref.i.scen.replica,"\n")
  cat("dataobs =",file.dataobs,"\n")
  cat("dataobs.replica =",file.dataobs.replica,"\n")
  cat("myseed =",myseed,"\n")
  cat("k_all =",k_all,"\n")
  cat("k_range =",k_range,"\n")
  cat("\n")
  sink()

  cat("\n")
  cat("#################### posterior GOF #################","\n")
  cat("refscen =",posterior.gof.file.ref,"\n")
  #cat("refscen.replica =",posterior.gof.file.ref.i.scen.replica,"\n")
  cat("dataobs =",file.dataobs,"\n")
  cat("dataobs.replica =",file.dataobs.replica,"\n")
  cat("\n")

  ################################################################################
  # Loop over n.ref
  for (n.ref in n.ref.hp.all) {

    data.ref.n.ref <- data.ref[1:n.ref,]
    param.ref.n.ref = param.ref.all[1:n.ref,]

    ####################################################################
    ## loop over n.post and hence eps

    for (n.post in n.post.all) {
      eps <- n.post / n.ref

      if (post.rej == TRUE) {
        ########## Rejection ##############################################
        # Capturer le temps de début
        start_time <- proc.time()

        gof.hp.rej <- hpgfit(target = data.obs,
                             target.replica = data.obs.replica ,
                             param = param.ref.n.ref,
                             sumstat = data.ref.n.ref,
                             sim.fun= sim.fun.ref,
                             method = c("rejection"),
                             kernel = c("epanechnikov"),
                             lambda = c(0.0001, 0.001, 0.01),
                             #param_transform = "none",
                             # trans.fun = trans_and_back$transform,
                             # back.trans.fun = trans_and_back$back_transform,
                             # param_lower_bound = param_lower_bound,
                             # param_upper_bound = param_upper_bound,
                             score = c("lof", "kNN"),
                             k = k_all,
                             k_range = k_range,
                             eps = eps,
                             split = 0.5,
                             norm = sd,
                             ncores = ncores.posterior,
                             nboot = n.boot,
                             ncores_sim = ncores_sim,
                             path_to_sim = path_to_sim,
                             scen.ref = scen.ref)

        # Capturer le temps de fin
        end_time <- proc.time()
        # Calculer la durée d'exécution
        execution_time <- end_time - start_time

        # Convertir le temps écoulé en heures, minutes et secondes
        execution_time_converted <- convert_time(execution_time["elapsed"])

        summary.lof.95 = summary(gof.hp.rej, score = "lof", k = "max", level = 0.95)
        summary.kNN.95 = summary(gof.hp.rej, score = "kNN", k = 1, level = 0.95)
        summary.lof.90 = summary(gof.hp.rej, score = "lof", k = "max", level = 0.90)
        summary.kNN.90 = summary(gof.hp.rej, score = "kNN", k = 1, level = 0.90)

        cat("--------------> REJ: n.sim =", n.ref,"  n.post = ",n.post, "  eps = ",eps, "  n.calib =", n.post*split,"\n")
        print(execution_time_converted)
        cat("\n")
        print(summary.lof.95)
        cat("\n")
        print(summary.kNN.95)
        cat("\n")
        # cat("\n")
        # print(summary.lof.90)
        # cat("\n")
        # print(summary.kNN.90)
        # cat("\n")

        # Écrire le contenu combiné dans un fichier txt
        sink(file = res.file.name.posterior, append = TRUE )
        cat("\n")
        cat("--------------> REJ: n.sim =", n.ref,"  n.post = ",n.post, "  eps = ",eps, "  n.calib =", n.post*split,"\n")
        print(execution_time_converted)
        cat("\n")
        print(summary.lof.95)
        cat("\n")
        print(summary.kNN.95)
        cat("\n")
        # cat("\n")
        # print(summary.lof.90)
        # cat("\n")
        # print(summary.kNN.90)
        # cat("\n")
        sink()

        if (post.rej.loclin == TRUE) {
          ########## rej.loclin ##############################################
          start_time <- proc.time()

          gof.hp.loclin  <- hpgfit(target = data.obs,
                                   target.replica = data.obs.replica ,
                                   param = param.ref.n.ref,
                                   sumstat = data.ref.n.ref,
                                   sim.fun= sim.fun.ref,
                                   method = "loclinear",
                                   kernel = c("epanechnikov"),
                                   lambda = c(0.0001, 0.001, 0.01),
                                   #param_transform = "none",
                                   trans.fun = trans_and_back$transform,
                                   back.trans.fun = trans_and_back$back_transform,
                                   # param_lower_bound = param_lower_bound,
                                   # param_upper_bound = param_upper_bound,
                                   score = c("lof", "kNN"),
                                   k = k_all,
                                   k_range = k_range,
                                   eps = eps,
                                   split = 0.5,
                                   norm = sd,
                                   ncores = ncores.posterior,
                                   nboot = n.boot,
                                   ncores_sim = ncores_sim,
                                   path_to_sim = path_to_sim,
                                   scen.ref = scen.ref)
          #							 ,
          #                             n.post = n.post) ### AE?

          # Capturer le temps de fin
          end_time <- proc.time()
          # Calculer la durée d'exécution
          execution_time <- end_time - start_time

          # Convertir le temps écoulé en heures, minutes et secondes
          execution_time_converted <- convert_time(execution_time["elapsed"])

          summary.lof.95 = summary(gof.hp.loclin, score = "lof", k = "max", level = 0.95)
          summary.kNN.95 = summary(gof.hp.loclin, score = "kNN", k = 1, level = 0.95)
          summary.lof.90 = summary(gof.hp.loclin, score = "lof", k = "max", level = 0.90)
          summary.kNN.90 = summary(gof.hp.loclin, score = "kNN", k = 1, level = 0.90)

          cat("--------------> LOCLIN: n.sim =", n.ref,"  n.post = ",n.post, "  eps = ",eps, "  n.calib =", n.post*split,"\n")
          print(execution_time_converted)
          cat("\n")
          print(summary.lof.95)
          cat("\n")
          print(summary.kNN.95)
          cat("\n")
          cat("\n")
          print(summary.lof.90)
          cat("\n")
          print(summary.kNN.90)
          cat("\n")

          # Écrire le contenu combiné dans un fichier txt
          sink(file = res.file.name.posterior, append = TRUE )
          cat("\n")
          cat("--------------> LOCLIN: n.sim =", n.ref,"  n.post = ",n.post, "  eps = ",eps, "  n.calib =", n.post*split,"\n")
          print(execution_time_converted)
          cat("\n")
          print(summary.lof.95)
          cat("\n")
          print(summary.kNN.95)
          cat("\n")
          cat("\n")
          print(summary.lof.90)
          cat("\n")
          print(summary.kNN.90)
          cat("\n")
          sink()
        }

        if (post.rej.ridge == TRUE) {
          ########## rej.ridge ##############################################
          start_time <- proc.time()

          gof.hp.ridge  <- hpgfit(target = data.obs,
                                  target.replica = data.obs.replica ,
                                  param = param.ref.n.ref,
                                  sumstat = data.ref.n.ref,
                                  sim.fun= sim.fun.ref,
                                  method = "ridge",
                                  kernel = c("epanechnikov"),
                                  lambda = c(0.0001, 0.001, 0.01),
                                  #param_transform = "none",
                                  trans.fun = trans_and_back$transform,
                                  back.trans.fun = trans_and_back$back_transform,
                                  # param_lower_bound = param_lower_bound,
                                  # param_upper_bound = param_upper_bound,
                                  score = c("lof", "kNN"),
                                  k = k_all,
                                  k_range = k_range,
                                  eps = eps,
                                  split = 0.5,
                                  norm = sd,
                                  ncores = ncores.posterior,
                                  nboot = n.boot,
                                  ncores_sim = ncores_sim,
                                  path_to_sim = path_to_sim,
                                  scen.ref = scen.ref)

          # Capturer le temps de fin
          end_time <- proc.time()
          # Calculer la durée d'exécution
          execution_time <- end_time - start_time

          # Convertir le temps écoulé en heures, minutes et secondes
          execution_time_converted <- convert_time(execution_time["elapsed"])

          summary.lof.95 = summary(gof.hp.ridge, score = "lof", k = "max", level = 0.95)
          summary.kNN.95 = summary(gof.hp.ridge, score = "kNN", k = 1, level = 0.95)
          summary.lof.90 = summary(gof.hp.ridge, score = "lof", k = "max", level = 0.90)
          summary.kNN.90 = summary(gof.hp.ridge, score = "kNN", k = 1, level = 0.90)

          cat("--------------> RIDGE: n.sim =", n.ref,"  n.post = ",n.post, "  eps = ",eps, "  n.calib =", n.post*split,"\n")
          print(execution_time_converted)
          cat("\n")
          print(summary.lof.95)
          cat("\n")
          print(summary.kNN.95)
          cat("\n")
          cat("\n")
          print(summary.lof.90)
          cat("\n")
          print(summary.kNN.90)
          cat("\n")

          # Écrire le contenu combiné dans un fichier txt
          sink(file = res.file.name.posterior, append = TRUE )
          cat("\n")
          cat("--------------> RIDGE: n.sim =", n.ref,"  n.post = ",n.post, "  eps = ",eps, "  n.calib =", n.post*split,"\n")
          print(execution_time_converted)
          cat("\n")
          print(summary.lof.95)
          cat("\n")
          print(summary.kNN.95)
          cat("\n")
          cat("\n")
          print(summary.lof.90)
          cat("\n")
          print(summary.kNN.90)
          cat("\n")
          sink()
        }
      }
    }
  }
}

