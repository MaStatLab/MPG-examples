options(stringsAsFactors = FALSE)

library(cremid)
library(dplyr)
library(mvtnorm)
library(mclust)
library(ggplot2)

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
prior <- list(K = K, truncation_type = "fixed", shared_alpha = FALSE)
cremid.null.fixed <- RunCremid(null.training, prior = prior)
save(cremid.null.fixed, file = "cremid_null_fixed.rda")
# load("cremid_null_fixed.rda")
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

cremid.alt.fixed <- RunCremid(alt.training, prior = prior)
save(cremid.alt.fixed, file = "cremid_alt_fixed.rda")
# load("cremid_alt_fixed.rda")
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

# Sensitivity analysis
# prior.1 <- list(K = K, 
#                 truncation_type = "fixed", 
#                 tau_varphi = c(1, 1),
#                 shared_alpha = FALSE)
# cremid.null.fixed.1 <- RunCremid(null.training, prior = prior.1)
# save(cremid.null.fixed.1, file = "cremid_null_fixed_1.rda")
# cremid.alt.fixed.1 <- RunCremid(alt.training, prior = prior.1)
# save(cremid.alt.fixed.1, file = "cremid_alt_fixed_1.rda")
# 
# prior.2 <- list(K = K, 
#                 truncation_type = "fixed", 
#                 tau_varphi = c(2, 2),
#                 shared_alpha = FALSE)
# cremid.null.fixed.2 <- RunCremid(null.training, prior = prior.2)
# save(cremid.null.fixed.2, file = "cremid_null_fixed_2.rda")
# cremid.alt.fixed.2 <- RunCremid(alt.training, prior = prior.2)
# save(cremid.alt.fixed.2, file = "cremid_alt_fixed_2.rda")
# 
# prior.3 <- list(K = K, 
#                 truncation_type = "fixed", 
#                 tau_rho = c(1, 1),
#                 shared_alpha = FALSE)
# cremid.null.fixed.3 <- RunCremid(null.training, prior = prior.3)
# save(cremid.null.fixed.3, file = "cremid_null_fixed_3.rda")
# cremid.alt.fixed.3 <- RunCremid(alt.training, prior = prior.3)
# save(cremid.alt.fixed.3, file = "cremid_alt_fixed_3.rda")
# 
# prior.4 <- list(K = K, 
#                 truncation_type = "fixed", 
#                 tau_rho = c(2, 2),
#                 shared_alpha = FALSE)
# cremid.null.fixed.4 <- RunCremid(null.training, prior = prior.4)
# save(cremid.null.fixed.4, file = "cremid_null_fixed_4.rda")
# cremid.alt.fixed.4 <- RunCremid(alt.training, prior = prior.4)
# save(cremid.alt.fixed.4, file = "cremid_alt_fixed_4.rda")
# 
# prior.5 <- list(K = K, 
#                 truncation_type = "fixed", 
#                 epsilon_range = c(1E-10, 0.5),
#                 shared_alpha = FALSE)
# cremid.null.fixed.5 <- RunCremid(null.training, prior = prior.5)
# save(cremid.null.fixed.5, file = "cremid_null_fixed_5.rda")
# cremid.alt.fixed.5 <- RunCremid(alt.training, prior = prior.5)
# save(cremid.alt.fixed.5, file = "cremid_alt_fixed_5.rda")
# 
# prior.6 <- list(K = K, 
#                 truncation_type = "fixed", 
#                 epsilon_range = c(1E-10, 0.25),
#                 shared_alpha = FALSE)
# cremid.null.fixed.6 <- RunCremid(null.training, prior = prior.6)
# save(cremid.null.fixed.6, file = "cremid_null_fixed_6.rda")
# cremid.alt.fixed.6 <- RunCremid(alt.training, prior = prior.6)
# save(cremid.alt.fixed.6, file = "cremid_alt_fixed_6.rda")
# 
# prior.7 <- list(K = K, 
#                 truncation_type = "fixed", 
#                 shared_alpha = TRUE)
# cremid.null.fixed.7 <- RunCremid(null.training, prior = prior.7)
# save(cremid.null.fixed.7, file = "cremid_null_fixed_7.rda")
# cremid.alt.fixed.7 <- RunCremid(alt.training, prior = prior.7)
# save(cremid.alt.fixed.7, file = "cremid_alt_fixed_7.rda")
# 
load("cremid_alt_fixed_4.rda")
load("cremid_alt_fixed_3.rda")
load("cremid_alt_fixed_2.rda")
load("cremid_alt_fixed_1.rda")
load("cremid_alt_fixed.rda")
load("cremid_alt_fixed_5.rda")
load("cremid_alt_fixed_6.rda")
load("cremid_alt_fixed_7.rda")

load("cremid_null_fixed_4.rda")
load("cremid_null_fixed_3.rda")
load("cremid_null_fixed_2.rda")
load("cremid_null_fixed_1.rda")
load("cremid_null_fixed.rda")
load("cremid_null_fixed_5.rda")
load("cremid_null_fixed_6.rda")
load("cremid_null_fixed_7.rda")


alt.string <- "Different stimulation conditions"
null.string <- "Control study"



# Varphi

varphi.alt <- bind_rows(
  ExtractPosterior(cremid.alt.fixed$chain$varphi,
                   "varphi",
                   cremid.alt.fixed$prior$tau_varphi[1]),
  ExtractPosterior(cremid.alt.fixed.1$chain$varphi,
                   "varphi",
                   cremid.alt.fixed.1$prior$tau_varphi[1]),
  ExtractPosterior(cremid.alt.fixed.2$chain$varphi,
                   "varphi",
                   cremid.alt.fixed.2$prior$tau_varphi[1])) %>%
  dplyr::mutate(version = alt.string)

varphi.null <- bind_rows(
  ExtractPosterior(cremid.null.fixed$chain$varphi,
                   "varphi",
                   cremid.null.fixed$prior$tau_varphi[1]),
  ExtractPosterior(cremid.null.fixed.1$chain$varphi,
                   "varphi",
                   cremid.null.fixed.1$prior$tau_varphi[1]),
  ExtractPosterior(cremid.null.fixed.2$chain$varphi,
                   "varphi",
                   cremid.null.fixed.2$prior$tau_varphi[1])) %>%
  dplyr::mutate(version = null.string)

pdf("varphi_sensitivity.pdf",
    pointsize = 18, width = 12, height = 5)
dplyr::bind_rows(varphi.null, varphi.alt) %>%
DensityPlot(., expression(paste(a[varphi], " = ",
                                          b[varphi])), expression(varphi)) +
  xlim(0, 0.75) +
  facet_wrap(~ version)
dev.off()

# Rho

rho.alt <- bind_rows(
  ExtractPosterior(cremid.alt.fixed$chain$rho,
                   "rho",
                   cremid.alt.fixed$prior$tau_rho[1]),
  ExtractPosterior(cremid.alt.fixed.3$chain$rho,
                   "rho",
                   cremid.alt.fixed.3$prior$tau_rho[1]),
  ExtractPosterior(cremid.alt.fixed.4$chain$rho,
                   "rho",
                   cremid.alt.fixed.4$prior$tau_rho[1])) %>%
  dplyr::mutate(version = alt.string)

rho.null <- bind_rows(
  ExtractPosterior(cremid.null.fixed$chain$rho,
                   "rho",
                   cremid.null.fixed$prior$tau_rho[1]),
  ExtractPosterior(cremid.null.fixed.3$chain$rho,
                   "rho",
                   cremid.null.fixed.3$prior$tau_rho[1]),
  ExtractPosterior(cremid.null.fixed.4$chain$rho,
                   "rho",
                   cremid.null.fixed.4$prior$tau_rho[1])) %>%
  dplyr::mutate(version = null.string)

pdf("rho_sensitivity.pdf",
    pointsize = 18, width = 12, height = 5)
dplyr::bind_rows(rho.null, rho.alt) %>%
  DensityPlot(., expression(paste(a[rho], " = ",
                                        b[rho], " = ")), expression(rho)) +
  facet_wrap(~ version, scales = "free") +
  geom_blank(aes(x = 0.0)) +
  geom_blank(aes(x = 1.0))
dev.off()

# Epsilon

epsilon.alt <- bind_rows(
  ExtractPosterior(cremid.alt.fixed$chain$epsilon,
                   "epsilon",
                   cremid.alt.fixed$prior$epsilon_range[2]),
  ExtractPosterior(cremid.alt.fixed.5$chain$epsilon,
                   "epsilon",
                   cremid.alt.fixed.5$prior$epsilon_range[2]),
  ExtractPosterior(cremid.alt.fixed.6$chain$epsilon,
                   "epsilon",
                   cremid.alt.fixed.6$prior$epsilon_range[2])) %>%
  dplyr::mutate(version = alt.string)

epsilon.null <- bind_rows(
  ExtractPosterior(cremid.null.fixed$chain$epsilon,
                   "epsilon",
                   cremid.null.fixed$prior$epsilon_range[2]),
  ExtractPosterior(cremid.null.fixed.5$chain$epsilon,
                   "epsilon",
                   cremid.null.fixed.5$prior$epsilon_range[2]),
  ExtractPosterior(cremid.null.fixed.6$chain$epsilon,
                   "epsilon",
                   cremid.null.fixed.6$prior$epsilon_range[2])) %>%
  dplyr::mutate(version = null.string)

pdf("epsilon_sensitivity.pdf", 
    pointsize = 18, width = 12, height = 5)
dplyr::bind_rows(epsilon.alt, epsilon.null) %>%
  DensityPlot(., expression(b[epsilon]), expression(epsilon)) +
  facet_wrap( ~ version) + 
  xlim(0, 0.5)
dev.off()  

# Alpha 

alpha.alt <- bind_rows(
  ExtractPosterior(cremid.alt.fixed$chain$alpha[, 1],
                   "alpha",
                   "0"),
  ExtractPosterior(cremid.alt.fixed$chain$alpha[, 2],
                   "alpha",
                   "1"),
  ExtractPosterior(cremid.alt.fixed.7$chain$alpha[, 1],
                   "alpha",
                   "shared")) %>%
  dplyr::mutate(version = alt.string)

alpha.null <- bind_rows(
  ExtractPosterior(cremid.null.fixed$chain$alpha[, 1],
                   "alpha",
                   "0"),
  ExtractPosterior(cremid.null.fixed$chain$alpha[, 2],
                   "alpha",
                   "1"),
  ExtractPosterior(cremid.null.fixed.7$chain$alpha[, 1],
                   "alpha",
                   "shared")) %>%
  dplyr::mutate(version = null.string)

pdf("alpha_sensitivity.pdf", 
    pointsize = 18, width = 12, height = 5)
dplyr::bind_rows(alpha.alt, alpha.null) %>%
  DensityPlot(., expression(alpha), expression(alpha)) +
  facet_wrap( ~ version) +
  xlim(0, 40)
dev.off()

# Number of mixture components

c.alt <- Clusters(cremid.alt.fixed) %>%
  dplyr::mutate(version = alt.string,
                alpha = "not shared")
c.null <- Clusters(cremid.null.fixed) %>%
  dplyr::mutate(version = null.string,
                alpha = "not shared")
c.alt.5 <- Clusters(cremid.alt.fixed.7) %>%
  dplyr::mutate(version = alt.string,
                alpha = "shared")
c.null.5 <- Clusters(cremid.null.fixed.7) %>%
  dplyr::mutate(version = null.string,
                alpha = "shared")


cluster.data <- bind_rows(c.alt, c.null, c.alt.5, c.null.5)

pdf("z_sensitivity.pdf", 
    pointsize = 18, width = 12, height = 5)
cluster.data %>%
  ggplot(aes(x = shared)) +
  geom_histogram(aes(fill = alpha)) +
  facet_grid(~ version)
dev.off()


