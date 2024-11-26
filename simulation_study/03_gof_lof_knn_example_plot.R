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

# Script to analyse and plot a simple example
#
# REQUIREMENTS:
#   - R: version 4.4. or above.
#   - abcgof: version v0.0.1. (see below)
#
# DATA:
#   - "depindep": TODO
#
# OUTPUTS:
#   - PDF plots of the results

library(abcgof)
library(here)
library(ggplot2)
library(tidyr)

################################################################################
# Load Data
set.seed(1974)
datestamp_day <- "2024-04-24"
PATH <- here()

## Select data
dataname <- "depindep"
refname <- "scenario2"
testname <- "scenario1"
results_name <- paste0(datestamp_day, "_gof_resim_", dataname, "_ref_", refname, "_test_", testname)
scen.ref <- 2
scen.test <- 1

## Read data
folder <- file.path(PATH, "..", "..", "gof-abc-data",
                    "20240401_data_fast",
                    "indep_dep_4pop_SNP_prior_Uniform_WINDOWS_NC_MediumMoinsT_GOOD")
PATH.yellow.S1 <- file.path(folder, "yellow_points_from_depindep_fast_NoConditionMediumMoinsT_reftableRF_110000_param_S1.txt")
PATH.yellow.S2 <- file.path(folder, "yellow_points_from_depindep_fast_NoConditionMediumMoinsT_reftableRF_110000_param_S2.txt")
data1.S1 <- as.matrix(read.table(PATH.yellow.S1, header = FALSE))
data1.S2 <- as.matrix(read.table(PATH.yellow.S2, header = FALSE))

## subset data
n.ref.prior <- 1000
n.test <- 1000
results_name_nref <- paste0(results_name, "_nrep_", n.ref.prior)
res_file <-  here("results", paste0(results_name_nref, "_prior_lof_kNN_comparison"))

sel_ind <- sample(1:(n.ref.prior+n.test), n.ref.prior+n.test, replace = FALSE)
# ref
data.ref <- data1.S2[sel_ind[1:n.ref.prior], 2:131]
# test
data.test <- data1.S1[sel_ind[(n.ref.prior+1):(n.ref.prior+n.test)], 2:131]

################################################################################
# Normalization and split
data_norm <- norm_stats(data.test, data.ref, norm = sd)
data.norm.test <- as.matrix(data_norm$target)

n.calib.prior <- floor(0.5 * n.ref.prior)
calibid <- sample(seq_len(nrow(data_norm$sumstat)), size = n.calib.prior)
data.norm.calib <- data_norm$sumstat[calibid, ]
data.norm.ref <- data_norm$sumstat[-calibid, ]

################################################################################
# Prior GOF
gof.prior <- abcgof::gof.fit(data.test = data.norm.test,
                             data.ref = data.norm.ref,
                             data.calib = data.norm.calib,
                             k = 1:20,
                             k_range = c(5, 20),
                             score = c("lof", "kNN"),
                             ncores = 6)

## points that are rejected by lof, but not by kNN
sel1 <- gof.prior$lof$pval[, "max"] <= 0.005 & gof.prior$kNN$pval[, "1"] > 0.2
sum(sel1)
gof.prior$lof$score.test[sel1, "max"]
gof.prior$kNN$score.test[sel1, "1"]

## points that are rejected by neither
sel2 <- gof.prior$lof$pval[, "max"] > 0.28 & gof.prior$kNN$pval[, "1"] > 0.59
sel2 <- sel2 & (gof.prior$lof$pval[, "max"] <= 0.3 | gof.prior$kNN$pval[, "1"] <= 0.6)
sum(sel2)
gof.prior$lof$pval[sel2, "max"]
gof.prior$kNN$pval[sel2, "1"]

sel <- sel1 | sel2

################################################################################
# plot data

## PCA
trainall <- rbind(data.norm.test, data.norm.calib, data.norm.ref)
res.pca <- FactoMineR::PCA(trainall, graph = FALSE)
precent_explain <- round(res.pca$svd$vs^2 / sum(res.pca$svd$vs^2) * 100, 2)

## format for plotting
pca_plot <- as.data.frame(res.pca$ind$coord[c(sel, rep(TRUE, nrow(data.norm.calib)), rep(TRUE, nrow(data.norm.ref))), ])
pca_plot$dataset <- c(rep("test points", sum(sel)), rep("calibration points", nrow(data.norm.calib)), rep("reference points", nrow(data.norm.ref)))
pca_plot$status <- c("C", "B", rep("A", nrow(data.norm.ref)+nrow(data.norm.calib)))

pca_test <- pca_plot[1:2, ]
pca_test$plof <- paste0("p[lof]==", gof.prior$lof$pval[sel, "max"])
pca_test$pknn <- paste0("p[kNN]==", gof.prior$kNN$pval[sel, "1"])

## plot
p <- ggplot(pca_plot, aes(x = Dim.1, y = Dim.3, color = dataset, shape = status, size = dataset)) +
  geom_point(alpha = 0.8) +
  scale_colour_manual(values = c("calibration points" = "#E69F00",   # palette.colors(3)
                                 "reference points" = "#56B4E9",
                                 "test points" = "#000000"),
                      guide = guide_legend(override.aes = list(size = 1)),
                      name = element_blank(),
                      breaks = c("reference points", "calibration points")) +
  scale_shape(guide = "none") +
  scale_size_manual(values = c("calibration points" = 0.9,
                               "reference points" = 0.9,
                               "test points" = 1.2),
                    guide = "none") +
  geom_text(data = pca_test[1, ],
            aes(hjust = "right", vjust = "bottom", label = plof),
            nudge_x = -0.6, nudge_y = -0.2, parse = TRUE, size = 10, size.unit = "pt") +
  geom_text(data = pca_test[1, ],
            aes(hjust = "right", vjust = "top", label = pknn),
            nudge_x = -0.6, nudge_y = -0.2, parse = TRUE, size = 10, size.unit = "pt") +
  geom_text(data = pca_test[2, ],
            aes(hjust = "right", vjust = "bottom", label = plof),
            nudge_x = 0.6, nudge_y = 1.5, parse = TRUE, size = 10, size.unit = "pt") +
  geom_text(data = pca_test[2, ],
            aes(hjust = "right", vjust = "top", label = pknn),
            nudge_x = 0.6, nudge_y = 1.5, parse = TRUE, size = 10, size.unit = "pt") +
  xlab(paste0("PCA dimension 1 (", precent_explain[1], "%)")) +
  ylab(paste0("PCA dimension 3 (", precent_explain[3], "%)")) +
  theme_bw() +
  theme(text = element_text(size = 10),
        # title = element_text(size = 7),
        # panel.grid.minor = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.99, 0.01),
        legend.justification = c("right", "bottom"),
        legend.key.size = unit(10, 'pt'),
        # plot.margin = unit(c(0.1,0.1,-0.31,-0.5), "cm"),
        # axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
  )
p
colorblindr::cvd_grid(p)

columnwidth <- 3.25 # * 72 # 1 point / 72 = 1 inch
ggsave(filename = here("results", paste0(results_name_nref, "_lof_kNN_comparison.pdf")),
       plot = p,
       width = columnwidth,
       height = columnwidth,
       unit = "in")


## plot simplified
pca_plot$dataset[1:2] <- c("test point 2", "test point 1")
pca_test$dataset <- c("test point 2", "test point 1")
pca_plot$dataset <- factor(pca_plot$dataset)
levels(pca_plot$dataset) <- c("simulated points", "simulated points", "test point 1", "test point 2")
cols_pvalues <- c("#D55E00", "#009E73")
p <- ggplot(pca_plot, aes(x = Dim.1, y = Dim.3, color = dataset, shape = dataset, size = dataset)) +
  geom_point(alpha = 0.8) +
  scale_colour_manual(values = c("#56B4E9", "#000000", "#000000"),
                      guide = guide_legend(override.aes = list(size = 2)),
                      name = element_blank(),
                      #breaks = c("simulated points", "test points", "test points"),
                      labels = c("simulated points", "test point 1", "test point 2")) +
  # scale_shape(guide = "none") +
  scale_shape_manual(values = c(16, 17, 15),
                     #guide = guide_legend(override.aes = list(size = 1)),
                     name = element_blank(),
                     #breaks = c("A", "B", "C"),
                     labels = c("simulated points", "test point 1", "test point 2")) +
  scale_size_manual(values = c("simulated points" = 0.9,
                               "test point 1" = 1.2,
                               "test point 2" = 1.2),
                    guide = "none") +
  geom_text(data = pca_test[1, ],
            aes(hjust = "right", vjust = "bottom", label = plof), color = cols_pvalues[2],
            nudge_x = -0.6, nudge_y = -0.2, parse = TRUE, size = 10, size.unit = "pt") +
  geom_text(data = pca_test[1, ],
            aes(hjust = "right", vjust = "top", label = pknn), color = cols_pvalues[1],
            nudge_x = -0.6, nudge_y = -0.2, parse = TRUE, size = 10, size.unit = "pt") +
  geom_text(data = pca_test[2, ],
            aes(hjust = "right", vjust = "bottom", label = plof), color = cols_pvalues[1],
            nudge_x = 0.6, nudge_y = 1.5, parse = TRUE, size = 10, size.unit = "pt") +
  geom_text(data = pca_test[2, ],
            aes(hjust = "right", vjust = "top", label = pknn), color = cols_pvalues[1],
            nudge_x = 0.6, nudge_y = 1.5, parse = TRUE, size = 10, size.unit = "pt") +
  xlab(paste0("PCA dimension 1 (", precent_explain[1], "%)")) +
  ylab(paste0("PCA dimension 3 (", precent_explain[3], "%)")) +
  theme_bw() +
  theme(text = element_text(size = 10),
        # title = element_text(size = 7),
        # panel.grid.minor = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.99, 0.01),
        legend.justification = c("right", "bottom"),
        legend.key.size = unit(10, 'pt'),
        # plot.margin = unit(c(0.1,0.1,-0.31,-0.5), "cm"),
        # axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
  )
p
colorblindr::cvd_grid(p)

columnwidth <- 3.25 # * 72 # 1 point / 72 = 1 inch
ggsave(filename = here("result_figures", paste0(results_name_nref, "_lof_kNN_comparison_blue.pdf")),
       plot = p,
       width = columnwidth,
       height = columnwidth,
       unit = "in")


