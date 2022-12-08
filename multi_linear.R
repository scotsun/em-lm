library(tidyverse)


reverse_conditional <- function(params, y, X) {
  mu <- params$mu
  sigma_x_sq <- params$sigma_x_sq
  sigma_sq <- params$sigma_sq
  b <- params$b
  
  dist_params <- list()
  for (j in 2:ncol(X)) {
    temp_param <- list(
      ksi = (b[j] * (y - X[, -j] %*% b[-j]) * sigma_x_sq[j-1] + mu[j-1] * sigma_sq) /
        (b[j] ^ 2 * sigma_x_sq[j-1] + sigma_sq),
      phi_sq = sigma_x_sq[j-1] * sigma_sq / (b[j] ^ 2 * sigma_x_sq[j-1] + sigma_sq)
    )
    dist_params[[j-1]] <- temp_param
  }
  return(dist_params)
}

.pseudo_x <- function(params, y, X) {
  mu <- params$mu
  sigma_x_sq <- params$sigma_x_sq
  sigma_sq <- params$sigma_sq
  b <- params$b
  
  M <- (!is.na(X)) + 0
  X[!M] <- 1000 # arbitrary value to mask NAs
  
  conditional_params <- reverse_conditional(params, y, X)
  # column-bind the intercept data
  ksis <- cbind(1, sapply(conditional_params, function(x) x$ksi))
  out <- M * X + (1 - M) * ksis
  return(out)
}

pseduo_sum_x <- function(params, y, X) {
  pseduo_x_complete <- .pseudo_x(params, y, X)
  return(t(apply(pseduo_x_complete, 2, sum)))
}

pseudo_sum_x_sq <- function(params, y, X) {
  mu <- params$mu
  sigma_x_sq <- params$sigma_x_sq
  sigma_sq <- params$sigma_sq
  b <- params$b
  
  n <- length(y)
  p <- ncol(X)
  out <- array(dim = c(p, p, n))
  
  M <- (!is.na(X)) + 0
  
  conditional_params <- reverse_conditional(params, y, X)
  # include intercept data (constant 1)
  phi_sqs <- c(0, sapply(conditional_params, function(x) x$phi_sq))
  E_X <- .pseudo_x(params, y, X)
  E_XX <- t(E_X) %*% E_X + diag(apply(1 - M, 2, sum) * phi_sqs)
  return(E_XX)
}

m_step <- function(params, y, X) {
  mu <- params$mu
  sigma_x_sq <- params$sigma_x_sq
  sigma_sq <- params$sigma_sq
  b <- params$b
  
  n <- length(y)
  p <- ncol(X) - 1
  
  conditional_params <- reverse_conditional(params, y, X)
  phi_sqs <- sapply(conditional_params, function(x) x$phi_sq)
  
  pseduo_x <- .pseudo_x(params, y, X)
  S <- pseudo_sum_x_sq(params, y, X)
  
  
  new_mu <- rep(NA, p)
  new_sigma_x_sq <- rep(NA, p)
  for (j in 1:p) {
    new_mu[j] <- sum(pseduo_x[, j+1]) / n
    new_sigma_x_sq[j] <- (S[j+1,j+1] - mu[j] * sum(pseduo_x[, j+1]))/n + mu[j]^2
  }
  new_sigma_sq <- (t(y) %*% y - 2 * t(b) %*% t(pseduo_x) %*% y + t(b) %*% S %*% b) / n
  new_b <- solve(S) %*% t(pseduo_x) %*% y
  return(list(
    mu = new_mu,
    sigma_x_sq = new_sigma_x_sq,
    sigma_sq = drop(new_sigma_sq),
    b = drop(new_b)
  ))
}

set.seed(823)
x1 <- rnorm(200, mean = 10, sd = 5)
x2 <- rnorm(200, mean = 0, sd = 1)
y <- rnorm(200, mean = 5 + 1 * x1 + 3 * x2, sd = 2)
p1 <- 1 / (1 + exp(-(2 + 1 * x2)))
m1 <- rbinom(200, 1, p1)
x1 <- ifelse(m1 == 1, x1, NA)
X <- cbind(1, x1, x2)

# lm(y ~ 0 + X) %>% summary()

params <- list(
  mu = c(10, 10),
  sigma_x_sq = c(1, 1),
  sigma_sq = 1,
  b = c(0, 0, 0)
)

for (i in seq_len(40)) {
  params <- m_step(params, y, X)
  # print(round(unlist(params), 3))
}
mean(!is.na(x1))
params$b




