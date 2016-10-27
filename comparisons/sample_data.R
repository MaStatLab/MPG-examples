SampleShift1Data <- function(n, seed = 1) {
  
  # Number of mixture components
  K <- 4
  # Number of dimensions
  p <- 4
  # Number of groups
  G <- 3
  mu <- array(NA, dim = c(K, G, p))
  Sigma <- array(NA, dim = c(K, p, p))
  probs <- matrix(rep(c(0.3, 0.3, 0.2, 0.2), G), ncol = K, byrow = TRUE) 
    
  set.seed(seed)
  
  for (k in 1:K) {
    mu[k, , ] <- matrix(rep(runif(p, min = 0, max = 10), G), ncol = p, 
                        byrow = TRUE)
  }
  
  mu[1, 2, 2] <- mu[1, 1, 2] + 0.5
  mu[1, 3, 2] <- mu[1, 1, 2] + 1
  
  Sigma[1, , ] <- diag(rep(0.7, p)) + 0.2
  Sigma[2, , ] <- diag(rep(1, p)) + 1
  Sigma[3, , ] <- diag(rep(.5, p)) - 0.1
  Sigma[4, , ] <- diag(rep(.1, p))
  
  return(SampleData(n, probs, mu, Sigma))
}

SampleShiftAllData <- function(n, seed = 1) {
  
  K <- 4
  p <- 4
  G <- 3
  mu <- array(NA, dim = c(K, G, p))
  Sigma <- array(NA, dim = c(K, p, p) )
  
  probs <- matrix(rep(c( 0.3, 0.3, 0.2, 0.2), G), ncol = K, byrow = TRUE) 
    
  set.seed(seed)
  
  for (k in 1:K) {
    mu[k, , ] <- matrix(rep(runif(p, min = 0, max = 10), G), ncol = p, 
                        byrow = TRUE)
  }
    
  mu[, 2, ] <- mu[, 1, ] + 0.1
  mu[, 3, ] <- mu[, 1, ] - 0.1
  
  Sigma[1, , ] <- diag(rep(0.7, p)) + 0.2
  Sigma[2, , ] <- diag(rep(1, p)) + 1
  Sigma[3, , ] <- diag(rep(0.5, p)) - 0.1
  Sigma[4, , ] <- diag(rep(0.1, p))
  
  return(SampleData(n, probs, mu, Sigma))
}

SampleWeight2Data <- function(n, seed = 1) {
  
  K <- 4
  p <- 4
  G <- 3
  mu <- array(NA, dim = c(K, G, p))
  Sigma <- array(NA, dim = c(K, p, p) )
  
  probs <- matrix(NA, ncol = K, nrow = G, byrow = TRUE) 
  for (j in 1:G) {
    probs[j, ] <- c(0.19 - (j - 1) * 0.08, 0.7, 0.01 + (j - 1) * 0.08, 0.1) 
  }
  
  set.seed(seed)
  
  for (k in 1:K) {
    mu[k, , ] <- matrix(rep(runif(p, min = 0, max = 10), G), ncol = p, 
                        byrow = TRUE)
  }
  
  Sigma[1, , ] <- diag(rep(1, p))
  Sigma[2, , ] <- diag(rep(2, p))
  Sigma[3, , ] <- diag(rep(0.2, p))
  Sigma[4, , ] <- diag(rep(0.1, p))  
  
  return(SampleData(n, probs, mu, Sigma))
}

SampleWeightAllData <- function(n, seed = 1) {
  
  set.seed(seed)
  
  K <- 8
  mu1 <- rep(0, K)
  G <- 3
  Sigma1 <- diag(0.5, K) + 0.1
  probs <- matrix(NA, nrow = G, ncol = K)
  for (j in 1:G) {
    temp <- mvrnorm(1, mu1, Sigma1)
    probs[j, ] <- exp(temp) / sum(exp(temp))
  }
  
  p <- 4
  mu <- array(NA, dim = c(K, G, p))
  Sigma <- array(NA, dim = c(K, p, p) )
    
  for (k in 1:K) {
    mu[k, , ] <- matrix(rep(runif(p, min = 0, max = 10), G), ncol = p, 
                        byrow = TRUE)
  }
  
  Sigma[1, , ] <- diag(rep(1, p))
  Sigma[2, , ] <- diag(rep(2, p))
  Sigma[3, , ] <- diag(rep(0.2, p))
  for (k in 4:K) {
    Sigma[k, , ] <- diag(rep(0.1, p))  
  }
    
  return(SampleData(n, probs, mu, Sigma))
}

SampleData <- function(n, probs, mu, Sigma) {
  
  K <- dim(mu)[1]
  G <- dim(mu)[2]
  p <- dim(mu)[3]
  
  sigma <- ExtractMarginalVariances(Sigma)
  
  Z <- ComputeZ(G, K, n, probs)
  mix.com.sum <- ComputeMixComSum(Z, G, K)
  
  data <- InternalDraws(n, mix.com.sum, Z, mu, Sigma)
  
  Y <- ComputeY(data)
  C <- ComputeC(data)  
  C0 <- ComputeC0(C, n, G)
  
  return(list(Y = Y,
              C = C,
              C0 = C0,
              G = G,
              mu = mu,
              sigma = sigma,
              probs = probs))  
}

ExtractMarginalVariances <- function(Sigma) {
  p <- dim(Sigma)[2]
  K <- dim(Sigma)[1]
  sigma <- matrix(unlist(lapply(1:K, function(x) {
    diag(Sigma[x, , ])
  })), ncol = K)
  return(sigma)
}

ComputeY <- function(data) {
  G <- dim(data)[1] 
  n <- dim(data)[2]
  p <- dim(data)[3]
  Y <- matrix(NA, nrow = G * n, ncol = p)
  for (j in 1:G) {
    Y[((j - 1) * n + 1):(j * n), ] <- data[j, , ]
  }
  return(Y)
}

ComputeC <- function(data) {
  G <- dim(data)[1] 
  n <- dim(data)[2]
  C <- rep(NA, G * n)
  for (j in 1:G) {
    C[((j - 1) * n + 1):(j * n)] <- j 
  }
  return(C)
}

ComputeC0 <- function(C, n, G) {
  return(C[sample(1:(n * G), n * G)])
}

ComputeZ <- function(G, K, n, probs) {
  
  Z <- matrix(NA, nrow = G, ncol = K)
  for (j in 1:G) {
    Z[j, ] <- table(factor(sample(K, 
                                  n,
                                  replace = TRUE,
                                  prob = probs[j, ]),
                           levels = 1:K))
  }
  return(Z)
  
}

ComputeMixComSum <- function(Z, G, K) {
  mix.com.sum <- matrix(0, nrow = G, ncol = K + 1)
  for (j in 1:G) {
    mix.com.sum[j, 2:(K + 1)] <- cumsum(Z[j, ])
  } 
  return(mix.com.sum)
}

InternalDraws <- function(n, mix.com.sum, Z, mu, Sigma) {
  
  p <- dim(mu)[3]
  K <- ncol(Z)
  G <- nrow(mix.com.sum)
  data <- array(NA, dim = c(G, n, p))
  
  for (k in 1:K) {
    for (j in 1:G) {
      if (mix.com.sum[j, k] < mix.com.sum[j, k + 1]) {
        data[j, (mix.com.sum[j, k] + 1):mix.com.sum[j, k + 1], ] <- 
          mvrnorm(Z[j, k], mu[k, j, ], Sigma[k, , ])
      }
    }
  }
  return(data)
}