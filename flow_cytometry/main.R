options(stringsAsFactors = FALSE)

library(MPG)
library(dplyr)
library(mvtnorm)
library(mclust)

source("helpers.R")

# Number of observations per group for inference.
n.obs <- 5000
# Number of obserations per group for out of sample comparison.
n.testing <- 1000
# Dimensions of the original data to use.
dims <- 5:10
# Number of mixture components for MPG.
K <- 100
path <- "~/Dropbox/Duke/Thesis/Flowcytometry/Equapol"

# Comparison between 3 technical replicates within a lab
# to see the amount of experimental variation.           
name.null.1 <- "EQAPOL_101_EP4_K69038GP_Unstim_02064.csv"
name.null.2 <- "EQAPOL_101_EP4_K69038GP_Unstim_02065.csv"
name.null.3 <- "EQAPOL_101_EP4_K69038GP_Unstim_02066.csv"

null.data <- GetData(path, name.null.1, name.null.2, name.null.3)
null.training <- SampleData(null.data, n.obs, dims, seed = 1)

mpg.null.adaptive <- RunMPG(null.training, K, "adaptive")
save(mpg.null.adaptive, file = "mpg_null_adaptive.rda")
mpg.null.fixed <- RunMPG(null.training, K, "fixed")
save(mpg.null.fixed, file = "mpg_null_fixed.rda")
mc.null <- RunMClust(null.training)
save(mc.null, file = "mc_null.rda")

null.testing <- SampleData(null.data, n.testing, dims, seed = 1492)
(mpg.null.score.fixed <- GetMPGScore(null.testing, mpg.null.fixed))
save(mpg.null.score.fixed, file = "mpg_null_score_fixed.rda")
(mpg.null.score.adaptive <- GetMPGScore(null.testing, mpg.null.adaptive))
save(mpg.null.score.adaptive, file = "mpg_null_score_adaptive.rda")
(mc.null.score <- GetMCScore(null.testing, mc.null))
save(mc.null.score, file = "mc_null_score.rda")

flow.null.name.adaptive <- "flow_null_adaptive"
pdf(paste0(flow.null.name.adaptive, ".pdf"), 
    pointsize = 18, width = 15, height = 5)
par(mfrow = c(1, 3))
PlotHistograms(mpg.null.adaptive)
dev.off()

setEPS()
postscript(paste0(flow.null.name.adaptive, ".eps"), pointsize = 18, 
           height = 5, width = 15)
par(mfrow = c(1, 3))
PlotHistograms(mpg.null.adaptive)
dev.off()

flow.null.name.fixed <- "flow_null_fixed"
pdf(paste0(flow.null.name.fixed, ".pdf"), 
    pointsize = 18, width = 15, height = 5)
par(mfrow = c(1, 3))
PlotHistograms(mpg.null.fixed)
dev.off()

setEPS()
postscript(paste0(flow.null.name.fixed, ".eps"), pointsize = 18, 
           height = 5, width = 15)
par(mfrow = c(1, 3))
PlotHistograms(mpg.null.fixed)
dev.off()

# Same subject different stimulation conditions
# to see changes with stimulation. 
name.alt.1 <- "EQAPOL_101_EP4_K69038GP_Unstim_02064.csv"
name.alt.2 <- "EQAPOL_101_EP4_K69038GP_CEF_02072.csv"
name.alt.3 <- "EQAPOL_101_EP4_K69038GP_CMVpp65_02069.csv"

alt.data <- GetData(path, name.alt.1, name.alt.2, name.alt.3)
alt.training <- SampleData(alt.data, n.obs, dims, seed = 2)

mpg.alt.adaptive <- RunMPG(alt.training, K, "adaptive")
save(mpg.alt.adaptive, file = "mpg_alt_adaptive.rda")
mpg.alt.fixed <- RunMPG(alt.training, K, "fixed")
save(mpg.alt.fixed, file = "mpg_alt_fixed.rda")
mc.alt <- RunMClust(alt.training)
save(mc.alt, file = "mc_alt.rda")

alt.testing <- SampleData(alt.data, n.testing, dims, seed = 1493)
(mpg.alt.score.fixed <- GetMPGScore(alt.testing, mpg.alt.fixed))
save(mpg.alt.score.fixed, file = "mpg_alt_score_fixed.rda")
(mpg.alt.score.adaptive <- GetMPGScore(alt.testing, mpg.alt.adaptive))
save(mpg.alt.score.adaptive, file = "mpg_alt_score_adaptive.rda")
(mc.alt.score <- GetMCScore(alt.testing, mc.alt))
save(mc.alt.score, file = "mc_alt_score.rda")

flow.alt.name.adaptive <- "flow_alt_adaptive"
pdf(paste0(flow.alt.name.adaptive, ".pdf"), 
    pointsize = 18, width = 15, height = 5)
par(mfrow = c(1, 3))
PlotHistograms(mpg.alt.adaptive)
dev.off()

setEPS()
postscript(paste0(flow.alt.name.adaptive, ".eps"), pointsize = 18, 
           height = 5, width = 15)
par(mfrow = c(1, 3))
PlotHistograms(mpg.alt.adaptive)
dev.off()

flow.alt.name.fixed <- "flow_alt_fixed"
pdf(paste0(flow.alt.name.fixed, ".pdf"), 
    pointsize = 18, width = 15, height = 5)
par(mfrow = c(1, 3))
PlotHistograms(mpg.alt.fixed)
dev.off()

setEPS()
postscript(paste0(flow.alt.name.fixed, ".eps"), pointsize = 18, 
           height = 5, width = 15)
par(mfrow = c(1, 3))
PlotHistograms(mpg.alt.fixed)
dev.off()