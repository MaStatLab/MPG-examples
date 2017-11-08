options(stringsAsFactors = FALSE)

library(cremid)
library(MASS)
library(dplyr)
library(mclust)
library(DPpackage)
library(ROCR)
library(ggplot2)

setwd("~/Dropbox/Duke/Thesis/locally_tied_stick_breaking/examples/MPG-examples/comparisons/")
source("sample_data.R")
source("helpers.R")
source("plots.R")

# Number of mixture components for MPG.
K <- 100 # 100
# Number of observation per group.
n <- 100 # 100
# Number of simulations.
n.test <- 100 # 100

# Shift 1
shift.1 <- dplyr::bind_rows(lapply(1:n.test, 
                                   RunAnalysis, 
                                   SampleShift1Data, 
                                   "Local Shift", 
                                   n, 
                                   K))

# Shift All
shift.all <- dplyr::bind_rows(lapply(1:n.test, 
                                     RunAnalysis, 
                                     SampleShiftAllData, 
                                     "Global Shift", 
                                     n, 
                                     K))

# Weight 2
weight.2 <- dplyr::bind_rows(lapply(1:n.test, 
                                    RunAnalysis, 
                                    SampleWeight2Data, 
                                    "Local Weight Change", 
                                    n, 
                                    K))

# Weight All
# Note: error with seed = 79. This is why 79 is skipped.
tests <- c(1:78, 80:(n.test + 1))
weight.all <- dplyr::bind_rows(lapply(tests, 
                                      RunAnalysis, 
                                      SampleWeightAllData, 
                                      "Global Weight Change", 
                                      n, 
                                      K))

# L1 distance plots
all.output <- dplyr::bind_rows(shift.1, shift.all, weight.2, weight.all)

distance.name <- "distance"
pdf(paste0(distance.name, ".pdf"), height = 5, width = 12)
theme_set(theme_gray(base_size = 18))
DistancePlot(all.output) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
dev.off()

setEPS()
postscript(paste0(distance.name, ".eps"), height = 5, width = 12)
theme_set(theme_gray(base_size = 18))
DistancePlot(all.output) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
dev.off()

# ROC curves
roc.name <- "roc"
pdf(paste0(roc.name, ".pdf"), height = 12, width = 12)
par(mfrow = c(2, 2))
RocPlot(shift.1, legend = TRUE)
RocPlot(shift.all, legend = FALSE)
RocPlot(weight.2, legend = FALSE)
RocPlot(weight.all, legend = FALSE)
dev.off()

setEPS()
postscript(paste0(roc.name, ".eps"), height = 12, width = 12)
par(mfrow = c(2, 2))
RocPlot(shift.1, legend = TRUE)
RocPlot(shift.all, legend = FALSE)
RocPlot(weight.2, legend = FALSE)
RocPlot(weight.all, legend = FALSE)
dev.off()


