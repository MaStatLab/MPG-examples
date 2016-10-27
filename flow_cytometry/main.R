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

mpg.null <- RunMPG(null.training, K)
save(mpg.null, file = "mpg_null.rda")
mc.null <- RunMClust(null.training)
save(mc.null, file = "mc_null.rda")

null.testing <- SampleData(null.data, n.testing, dims, seed = 1492)
(mpg.null.score <- GetMPGScore(null.testing, mpg.null))
save(mpg.null.score, file = "mpg_null_score.rda")
(mc.null.score <- GetMCScore(null.testing, mc.null))
save(mc.null.score, file = "mc_null_score.rda")

flow.null.name <- "flow_null"
pdf(paste0(flow.null.name, ".pdf"), pointsize = 18, width = 15, height = 5)
par(mfrow = c(1, 3))
PlotHistograms(mpg.null)
dev.off()

setEPS()
postscript(paste0(flow.null.name, ".eps"), pointsize = 18, 
           height = 5, width = 15)
par(mfrow = c(1, 3))
PlotHistograms(mpg.null)
dev.off()

# Same subject different stimulation conditions
# to see changes with stimulation. 
name.alt.1 <- "EQAPOL_101_EP4_K69038GP_Unstim_02064.csv"
name.alt.2 <- "EQAPOL_101_EP4_K69038GP_CEF_02072.csv"
name.alt.3 <- "EQAPOL_101_EP4_K69038GP_CMVpp65_02069.csv"

alt.data <- GetData(path, name.alt.1, name.alt.2, name.alt.3)
alt.training <- SampleData(alt.data, n.obs, dims, seed = 2)

mpg.alt <- RunMPG(alt.training, K)
save(mpg.alt, file = "mpg_alt.rda")
mc.alt <- RunMClust(alt.training)
save(mc.alt, file = "mc_alt.rda")

alt.testing <- SampleData(alt.data, n.testing, dims, seed = 1493)
(mpg.alt.score <- GetMPGScore(alt.testing, mpg.alt))
save(mpg.alt.score, file = "mpg_alt_score.rda")
(mc.alt.score <- GetMCScore(alt.testing, mc.alt))
save(mc.alt.score, file = "mc_alt_score.rda")

flow.alt.name <- "flow_alt"
pdf(paste0(flow.alt.name, ".pdf"), pointsize = 18, width = 15, height = 5)
par(mfrow = c(1, 3))
PlotHistograms(mpg.alt)
dev.off()

setEPS()
postscript(paste0(flow.alt.name, ".eps"), pointsize = 18, 
           height = 5, width = 15)
par(mfrow = c(1, 3))
PlotHistograms(mpg.alt)
dev.off()