MakePriorMuller <- function(y, J) {
  return(list(pe1 = 0.1,
              pe0 = 0.1,
              ae = 1,
              be = 1,
              a0 = rep(1, J + 1),
              b0 = rep(1, J + 1),
              nu = 9,
              tinv = 0.25 * var(y),
              m0 = apply(y, 2, mean),
              S0 = var(y),
              nub = 9,
              tbinv = var(y)))
}

CremidTesting <- function(ans) {
  return(mean(ans$chain$rho) * mean(ans$chain$varphi) )
}

EvalTrueDensity <- function(x, mu, sigma, probs) {
  
  K <- dim(mu)[1]
  G <- dim(mu)[2]
  p <- dim(mu)[3]
  
  output <- array(0, dim = c(nrow(x), G, p))   
  
  for (k in 1:K) {
    for (g in 1:G) {
      for (dim in 1:p) {
        output[, g, dim] <- output[, g, dim] + 
          probs[g, k] * dnorm(x[, dim], mu[k, g, dim], sqrt(sigma[dim, k]))
      }
    }
  }
  return(output)
}

EvalCremidDensity <- function(x, ans) {
  G <- dim(ans$chain$mu)[1]
  K <- ans$prior$K
  p <- ncol(ans$data$Y)
  niter <- ans$mcmc$nsave
  
  output <- array(0, dim = c(nrow(x), G, p))   
  Omega <- ans$chain$Omega
  mu <- ans$chain$mu  
  probs <- ans$chain$w
  
  for (it in 1:niter) {
    for (k in 1:K) {
      if (sum(ans$chain$w[, k, it]) > G / 1E4) {
        ix <- seq(p * (k - 1) + 1, p * k)
        sigma <- sqrt(diag(solve(Omega[, ix, it])))
        for (g in 1:G) {
          m <- mu[g, ix, it]
          for (dim in 1:p) {
            output[, g, dim] <- output[, g, dim] +
              probs[g, k, it] * dnorm(x[, dim], m[dim], sigma[dim])
          }
        }        
      }
    }
  }
  return(output / niter)
}

EvalMullerDensity <- function(fit) {
  output <- array(0, dim = c(fit$ngrid, fit$nstudy, fit$nvar))    
  for (i in 1:fit$ngrid) {
    output[i, , ] <- 
      fit$densms[1:fit$nstudy, i + (1:fit$nvar - 1) * fit$ngrid]
  }
  return(output)
}

EvalMCDensity <- function(x, mc) {
  
  G <- length(mc)
  p <- nrow(mc[[1]]$mu)
  
  output <- array(0, dim = c(nrow(x), G, p))   
  
  for (g in 1:G) {
    probs <- mc[[g]]$probs
    mu <- mc[[g]]$mu
    sigma <- mc[[g]]$sigma
    K <- ncol(mu)
    for (k in 1:K) {
      for (dim in 1:p) {
        output[, g, dim] <- output[, g, dim] + 
          probs[k] * dnorm(x[, dim], mu[dim, k], sqrt(sigma[dim, k]))
      }
    }
  }
  return(output)
}

MeasureL1Distance <- function(data.true, data.estimated) {
  return(mean(abs(data.true - data.estimated)))
}

MClust <- function(data) {
  mc <- lapply(1:data$G, function(c) {
    temp <- densityMclust(data$Y[data$C == c, ])
    return(list(probs = temp$parameters$pro,
                mu = temp$parameters$mean,
                sigma = apply(temp$parameters$variance$sigma, 3, diag)))
  })
}

RunAnalysis <- function(seed, f, name, n, K) {
  tryCatch(RunSingleDataSet(seed, f, name, n, K), 
           error = function(e) cat("ERROR\n"))
}

RunSingleDataSet <- function(seed, f, name, n, K) {
  
  cat(sprintf("Seed %d\n", seed))
  
  mcmc <- list(nburn = 5000,
               nsave = 1000,
               nskip = 1,
               ndisplay = 10000)

  # Generate data
  data <- match.fun(f)(n, seed)
  cat("\tMuller\n")
  
  # Run Muller
  muller.1 <- HDPMdensity(y = data$Y,
                          study = data$C,
                          prior = MakePriorMuller(data$Y, data$G),
                          mcmc = mcmc,
                          status = TRUE)  
  muller.0 <- HDPMdensity(y = data$Y,
                          study = data$C0,
                          prior = MakePriorMuller(data$Y, data$G),
                          mcmc = mcmc,
                          status = TRUE)
  
  # Run CREMID
  cat("\tCremid\n")
  fixed.prior <- list(K = K, truncation_type = "fixed", shared_alpha = FALSE)
  cremid.fixed.1 <- Fit(data$Y, data$C, prior = fixed.prior, mcmc = mcmc)
  cremid.fixed.0 <- Fit(data$Y, data$C0, prior = fixed.prior, mcmc = mcmc)
  
  cat("\tMclust\n")
  # Run Mclust
  mc <- MClust(data)
  
  # Density estimation
  x <- muller.1$grid
  y.true <- EvalTrueDensity(x, data$mu, data$sigma, data$probs)
  y.muller <- EvalMullerDensity(muller.1)
  y.cremid.fixed <- EvalCremidDensity(x, cremid.fixed.1)
  y.mc <- EvalMCDensity(x, mc)
  
  muller.df <- data.frame(seed = seed,
                          name = name,
                          method = "HDPM",
                          null = mean(muller.0$coef["eps"]),
                          alternative = mean(muller.1$coef["eps"]),
                          distance = MeasureL1Distance(y.true, y.muller))
  
  cremid.fixed.df <- data.frame(seed = seed,
                             name = name,
                             method = "cremid",
                             null = CremidTesting(cremid.fixed.0),
                             alternative = CremidTesting(cremid.fixed.1),
                             distance = MeasureL1Distance(y.true, y.cremid.fixed))

  mc.df <- data.frame(seed = seed,
                      name = name,
                      method = "Mclust",
                      null = NA,
                      alternative = NA,
                      distance = MeasureL1Distance(y.true, y.mc))
  
  return(dplyr::bind_rows(muller.df, cremid.fixed.df, mc.df))
}


