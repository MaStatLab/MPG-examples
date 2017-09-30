options(stringsAsFactors = FALSE)

library(cremid)
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
# Number of mixture components for Cremid.
K <- 100
path <- "~/Dropbox/Duke/Thesis/Flowcytometry/Equapol"

# Comparison between 3 technical replicates within a lab
# to see the amount of experimental variation.           
name.null.1 <- "EQAPOL_101_EP4_K69038GP_Unstim_02064.csv"
name.null.2 <- "EQAPOL_101_EP4_K69038GP_Unstim_02065.csv"
name.null.3 <- "EQAPOL_101_EP4_K69038GP_Unstim_02066.csv"

null.data <- GetData(path, name.null.1, name.null.2, name.null.3)
null.training <- SampleData(null.data, n.obs, dims, seed = 1)

cremid.null.fixed <- RunCremid(null.training, K, "fixed")
save(cremid.null.fixed, file = "cremid_null_fixed.rda")
mc.null <- RunMClust(null.training)
save(mc.null, file = "mc_null.rda")


n.sim <- 100
cremid.null.score.fixed <- rep(NA, n.sim)
mc.null.score <- rep(NA, n.sim)

for (i in 1:n.sim) {
  print(i)
  null.testing <- SampleData(null.data, n.testing, dims, seed = i + 1)
  cremid.null.score.fixed[i] <- GetCremidScore(null.testing, cremid.null.fixed)
  mc.null.score[i] <- GetMCScore(null.testing, mc.null)
}

mean(cremid.null.score.fixed)
sd(cremid.null.score.fixed)

mean(mc.null.score)
sd(mc.null.score)


save(cremid.null.score.fixed, file = "cremid_null_score_fixed.rda")
save(mc.null.score, file = "mc_null_score.rda")

flow.null.name.fixed.hist <- "flow_null_fixed_hist"
pdf(paste0(flow.null.name.fixed.hist, ".pdf"), 
    pointsize = 18, width = 15, height = 5)
par(mfrow = c(1, 3))
PlotHistograms(cremid.null.fixed, xlim = c(0, 1))
dev.off()

setEPS()
postscript(paste0(flow.null.name.fixed.hist, ".eps"), pointsize = 18, 
           height = 5, width = 15)
par(mfrow = c(1, 3))
PlotHistograms(cremid.null.fixed, xlim = c(0, 1))
dev.off()

flow.null.name.fixed.free.hist <- "flow_null_fixed_free_hist"
pdf(paste0(flow.null.name.fixed.free.hist, ".pdf"), 
    pointsize = 18, width = 15, height = 5)
par(mfrow = c(1, 3))
PlotHistograms(cremid.null.fixed)
dev.off()

setEPS()
postscript(paste0(flow.null.name.fixed.free.hist, ".eps"), pointsize = 18, 
           height = 5, width = 15)
par(mfrow = c(1, 3))
PlotHistograms(cremid.null.fixed)
dev.off()

flow.null.name.fixed.scatter <- "flow_null_fixed_scatter"
pdf(paste0(flow.null.name.fixed.scatter, ".pdf"), 
    pointsize = 18, width = 15, height = 5)
par(mfrow = c(1, 3))
ScatterPlotHyper(cremid.null.fixed)
dev.off()

setEPS()
postscript(paste0(flow.null.name.fixed.scatter, ".eps"), pointsize = 18, 
           height = 5, width = 15)
par(mfrow = c(1, 3))
ScatterPlotHyper(cremid.null.fixed)
dev.off()


# Same subject different stimulation conditions
# to see changes with stimulation. 
name.alt.1 <- "EQAPOL_101_EP4_K69038GP_Unstim_02064.csv"
name.alt.2 <- "EQAPOL_101_EP4_K69038GP_CEF_02072.csv"
name.alt.3 <- "EQAPOL_101_EP4_K69038GP_CMVpp65_02069.csv"

alt.data <- GetData(path, name.alt.1, name.alt.2, name.alt.3)
alt.training <- SampleData(alt.data, n.obs, dims, seed = 2)

cremid.alt.fixed <- RunCremid(alt.training, K, "fixed")
save(cremid.alt.fixed, file = "cremid_alt_fixed.rda")
mc.alt <- RunMClust(alt.training)
save(mc.alt, file = "mc_alt.rda")


n.sim <- 100
cremid.alt.score.fixed <- rep(NA, n.sim)
mc.alt.score <- rep(NA, n.sim)

for (i in 1:n.sim) {
  print(i)
  alt.testing <- SampleData(alt.data, n.testing, dims, seed = i + 1)
  cremid.alt.score.fixed[i] <- GetCremidScore(alt.testing, cremid.alt.fixed)
  mc.alt.score[i] <- GetMCScore(alt.testing, mc.alt)
}

mean(cremid.alt.score.fixed)
sd(cremid.alt.score.fixed)

mean(mc.alt.score)
sd(mc.alt.score)

save(cremid.alt.score.fixed, file = "cremid_alt_score_fixed.rda")
save(mc.alt.score, file = "mc_alt_score.rda")

flow.alt.name.fixed.hist <- "flow_alt_fixed_hist"
pdf(paste0(flow.alt.name.fixed.hist, ".pdf"), 
    pointsize = 18, width = 15, height = 5)
par(mfrow = c(1, 3))
PlotHistograms(cremid.alt.fixed, xlim = c(0, 1))
dev.off()

setEPS()
postscript(paste0(flow.alt.name.fixed.hist, ".eps"), pointsize = 18, 
           height = 5, width = 15)
par(mfrow = c(1, 3))
PlotHistograms(cremid.alt.fixed, xlim = c(0, 1))
dev.off()

flow.alt.name.fixed.free.hist <- "flow_alt_fixed_free_hist"
pdf(paste0(flow.alt.name.fixed.free.hist, ".pdf"), 
    pointsize = 18, width = 15, height = 5)
par(mfrow = c(1, 3))
PlotHistograms(cremid.alt.fixed)
dev.off()

setEPS()
postscript(paste0(flow.alt.name.fixed.free.hist, ".eps"), pointsize = 18, 
           height = 5, width = 15)
par(mfrow = c(1, 3))
PlotHistograms(cremid.alt.fixed)
dev.off()

flow.alt.name.fixed.scatter <- "flow_alt_fixed_scatter"
pdf(paste0(flow.alt.name.fixed.scatter, ".pdf"), 
    pointsize = 18, width = 15, height = 5)
par(mfrow = c(1, 3))
ScatterPlotHyper(cremid.alt.fixed)
dev.off()

setEPS()
postscript(paste0(flow.alt.name.fixed.scatter, ".eps"), pointsize = 18, 
           height = 5, width = 15)
par(mfrow = c(1, 3))
ScatterPlotHyper(cremid.alt.fixed)
dev.off()