options(stringsAsFactors = FALSE)

library(MASS)
library(cremid)

setwd(file.path("~/Dropbox/Duke/Thesis/locally_tied_stick_breaking",
                "examples/MPG-examples/calibration"))
source("../comparisons/sample_data.R")

# Number of mixture components.
K <- 4
# Number of dimensions.
p <- 4
# Number of samples per group.
n <- 1000
# Number of groups.
G <- 3

# Define mixture probabilities.
probs <- matrix(NA, ncol = K, nrow = G)
probs[1, ] <- c(0.16, 0.8, 0.02, 0.02) 
probs[2, ] <- c(0.09, 0.8, 0.09, 0.02) 
probs[3, ] <- c(0.02, 0.8, 0.16, 0.02) 

# Define mixture means.
mu <- array(NA, dim = c(K, G, p))
mu[1, , ] <- matrix(rep(c(1, 9), G * p / 2), ncol = p, byrow = TRUE)
mu[2, , ] <- matrix(rep(c(8, 8), G * p / 2), ncol = p, byrow = TRUE)
mu[3, , ] <- matrix(rep(c(1, 1), G * p / 2), ncol = p, byrow = TRUE)
mu[4, , ] <- matrix(rep(c(7, 1), G * p / 2), ncol = p, byrow = TRUE)
# Introduce shifts. 
mu[1, 2, 2] <- mu[1, 1, 2] - 1
mu[1, 3, 2] <- mu[1, 1, 2] - 2
mu[4, 2, 1] <- mu[4, 1, 1] + 1
mu[4, 3, 1] <- mu[4, 1, 1] + 2

# Define Covariances.
Sigma <- array(NA, dim = c(K, p, p))
Sigma[1, , ] <- diag(rep(1, p))
Sigma[2, , ] <- diag(rep(2, p))
Sigma[3, , ] <- diag(rep(0.2, p))
Sigma[4, , ] <- diag(rep(0.1, p))

# Sample data.
set.seed(1)
data <- SampleData(n, probs, mu, Sigma)

# Reshuffle data.
n.tot <- n * G
idx <- sample(n.tot, n.tot)
Y <- data$Y[idx, ]
C <- data$C[idx]

# Run model.
prior <- list(K = 100, shared_alpha = FALSE)
ans <- Fit(Y, C, prior = prior)

# Calibrate.
cal <- Calibrate(ans)

ScatterPlot <- function(Y, C, c, dim.1, dim.2, title) {
  plot(Y[C == c, ], 
       pch = 1, 
       main = paste(title, c, sep = " "), 
       xlim = range(Y[, dim.1]), ylim = range(Y[, dim.2]), 
       col = "darkgrey",
       xaxt = "n", yaxt = "n",
       ylab = "", xlab = "")
  axis(side = 1, at = c(1, 8))
  axis(side = 2, at = c(1, 8))
  abline(h = 8, lty = 2)
  abline(v = 8, lty = 2)
  abline(h = 1, lty = 2)
  abline(v = 1, lty = 2)
}

FinalPlots <- function() {
  sample.name <- "sample"
  cal.sample.name <- "calibrated sample"
  ScatterPlot(Y, C, 1, 1, 2, sample.name)
  ScatterPlot(Y, C, 2, 1, 2, sample.name)
  ScatterPlot(Y, C, 3, 1, 2, sample.name)
  ScatterPlot(cal$Y_cal, C, 1, 1, 2, cal.sample.name)
  ScatterPlot(cal$Y_cal, C, 2, 1, 2, cal.sample.name)
  ScatterPlot(cal$Y_cal, C, 3, 1, 2, cal.sample.name)
}

file.name <- "calibration"
pdf(paste0(file.name, ".pdf"), height = 8, width = 12, pointsize = 18)
par(mfrow=c(2, 3))
par(mar = c(4, 4, 4, 4))
FinalPlots()
dev.off()

setEPS()
postscript(paste0(file.name, ".eps"), height = 8, width = 12, pointsize = 18)
par(mfrow=c(2, 3))
par(mar = c(4, 4, 4, 4))
FinalPlots()
dev.off()
