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

# Script simulate the toy "gaussianlaplace" dataset
#
# This script simulates reference tables using the toy Gaussian - Laplace models.
#
# It saves all the simulated data in the "reftables/gaussianlaplace" repository.

################################################################################
## Load packages
library(Lmoments)
library(VGAM)
library(lmom)
library(gofabcpkg)
library(doParallel)
library(foreach)

################################################################################
## File management
library(here)
datestamp_day <- format(Sys.time(), "%Y-%m-%d")
dir.create(here("reftables", "gaussianlaplace"))

################################################################################
## Parameters
n.ref = 100*1000
n.test = 1000
n.calib = 1000
d = 350
rmax = 20 #nbr de L-moments
lmomnameall <- c("salmu", "Lmoments") ### samlmu ou Lmoments ???
refnameall <- c("gaussian", "laplace")

################################################################################
## Simulation
for (lmomname in lmomnameall) {
  for (refname in refnameall) {
    for (testname in refnameall) {
      set.seed(1289)
      # name of result
      results_name <- paste0(datestamp_day, "_toy_gof_freq_lmom_", lmomname, "_ref_", refname, "_test_", testname)
      # l moment function
      if (rmax != 20) stop("The value of rmax changed.")
      lmonfun <- switch(lmomname,
                        salmu = function(vec) {
                          lmom::samlmu(vec, nmom = 20)
                        },
                        Lmoments = function(vec) {
                          Lmoments::Lmoments(vec, rmax = 20)
                        })
      # sim functions
      simgauss <- function(n) {
        mu <- runif(n, -10, 10)
        sigma <- runif(n, 1, 4)
        param <- cbind(mu, sigma)
        colnames(param) <- c("mu", "sigma")
        if (d != 350) stop("The value of d changed.")
        sim.fun <- switch(lmomname,
                          salmu = function(params, ...) {
                            t(apply(params, 1, function(oneparam) lmom::samlmu(stats::rnorm(350, mean = oneparam[1], sd = oneparam[2]), nmom = 20)))
                          },
                          Lmoments = function(params, ...) {
                            t(apply(params, 1, function(oneparam) Lmoments::Lmoments(stats::rnorm(350, mean = oneparam[1], sd = oneparam[2]), rmax = 20)))
                          })
        sim <- sim.fun(param)
        simrep <- sim.fun(param)
        return(list(param = param, sim = sim, simrep = simrep, sim.fun = sim.fun))
      }
      simlaplace <- function(n) {
        mu <- runif(n, -10, 10)
        sigma <- runif(n, 1, 4)
        param <- cbind(mu, sigma/sqrt(2))
        colnames(param) <- c("mu", "scale")
        if (d != 350) stop("The value of d changed.")
        sim.fun <- switch(lmomname,
                          salmu = function(params, ...) {
                            t(apply(params, 1, function(oneparam) lmom::samlmu(VGAM::rlaplace(350, location = oneparam[1], scale = oneparam[2]), nmom = 20)))
                          },
                          Lmoments = function(params, ...) {
                            t(apply(params, 1, function(oneparam) Lmoments::Lmoments(VGAM::rlaplace(350, location = oneparam[1], scale = oneparam[2]), rmax = 20)))
                          })
        sim <- sim.fun(param)
        simrep <- sim.fun(param)
        return(list(param = param, sim = sim, simrep = simrep, sim.fun = sim.fun))
      }
      simref <- switch(refname,
                       gaussian = simgauss,
                       laplace = simlaplace)
      simtest <- switch(testname,
                        gaussian = simgauss,
                        laplace = simlaplace)

      # Ref
      refdata <- simref(n.ref)
      param.ref <- refdata$param
      sim.fun.ref <- refdata$sim.fun
      data.ref <- refdata$sim
      data.ref.replica <- refdata$simrep
      data.ref.calib <- simref(n.calib)$sim

      # Test
      testdata <- simtest(n.test)
      data.test.obs <- testdata$sim
      data.test.new <- testdata$simrep
      param.test <- testdata$param

      ## Plot PCA
      trainall <- rbind(data.ref, data.test.obs)
      res.pca <- FactoMineR::PCA(trainall, graph = FALSE)
      colo <-c(rep(1,n.ref),rep(2,n.test))
      mix <- c(sample(1:n.ref, 2000), (n.ref+1):(n.ref+n.test))
      pdf(file = here("reftables", "gaussianlaplace", paste0(results_name, "_pca.pdf")))
      plot(res.pca$ind$coord[mix,c(1,2)],col=colo[mix],pch="*")
      dev.off()

      ## Save dataset
      saveRDS(list(data.ref = data.ref,
                   param.ref = param.ref,
                   sim.fun.ref = sim.fun.ref,
                   data.ref.replica = data.ref.replica,
                   data.ref.calib = data.ref.calib,
                   data.test.obs = data.test.obs,
                   data.test.new = data.test.new,
                   param.test = param.test,
                   lmonfun = lmonfun,
                   datestamp_day = datestamp_day,
                   param_upper_bound = c(10, 4),
                   param_lower_bound = c(-10, 1)),
              file =  here("reftables", "gaussianlaplace", paste0(results_name, "_data.rds")))
    }
  }
}
