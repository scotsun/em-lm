library(tidyverse)
library(progress)

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

reverse_cond <- function(params, y, X) {
  # decompose params
  mu <- params$mu
  sigma_x_sq <- params$sigma_x_sq
  sigma_sq <- params$sigma_sq
  b <- params$b
  # cond_var & cond_mean
  cond_var <- solve(b[-1] %*% t(b[-1]) / sigma_sq + solve(diag(sigma_x_sq)))
  cond_mean <- (((y - b[1])/sigma_sq) %o% b[-1] + (rep(1, 200) %o% (solve(diag(sigma_x_sq)) %*% mu))[,,1]) %*% cond_var
}

pseudo_x <- function(params, y, X) {
  M <- (!is.na(X)) + 0
  # arbitrary value to mask NAs
  X[!M] <- -999
  # decompose params
  mu <- params$mu
  sigma_x_sq <- params$sigma_x_sq
  sigma_sq <- params$sigma_sq
  b <- params$b
  # reverse conditional distribution f(x|y) mean
  conditional_params <- reverse_conditional(params, y, X)
  # column-bind the intercept data
  ksis <- cbind(1, sapply(conditional_params, function(x) x$ksi))
  out <- M * X + (1 - M) * ksis
  return(out)
}

.pseduo_sum_x <- function(params, y, X) {
  pseduo_x_complete <- pseudo_x(params, y, X)
  return(t(apply(pseduo_x_complete, 2, sum)))
}

pseudo_x_sq <- function(param, y, X) {
  n <- length(y)
  M <- (!is.na(X)) + 0
  # arbitrary value to mask NAs
  X[!M] <- -1
  # decompose params
  mu <- params$mu
  sigma_x_sq <- params$sigma_x_sq
  sigma_sq <- params$sigma_sq
  b <- params$b
  # reverse conditional distribution f(x|y) mean & var
  conditional_params <- reverse_conditional(params, y, X)
  ksis <- cbind(1, sapply(conditional_params, function(x) x$ksi))
  phi_sqs <- c(0, sapply(conditional_params, function(x) x$phi_sq))
  phi_sqs_ <- matrix(rep(phi_sqs, n), nrow = n, byrow = TRUE)
  out <- M * X^2 + (1 - M) * (ksis^2 + phi_sqs_)
  return(out)
}

.pseudo_sum_x_sq <- function(params, y, X) {
  M <- (!is.na(X)) + 0
  mu <- params$mu
  sigma_x_sq <- params$sigma_x_sq
  sigma_sq <- params$sigma_sq
  b <- params$b
  
  n <- length(y)
  p <- ncol(X)
  out <- array(dim = c(p, p, n))
  
  conditional_params <- reverse_conditional(params, y, X)
  # include intercept data (constant 1)
  phi_sqs <- c(0, sapply(conditional_params, function(x) x$phi_sq))
  E_X <- pseudo_x(params, y, X)
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
  
  pseduo_x_ <- pseudo_x(params, y, X)
  S <- .pseudo_sum_x_sq(params, y, X)
  
  new_mu <- rep(NA, p)
  new_sigma_x_sq <- rep(NA, p)
  for (j in 1:p) {
    new_mu[j] <- sum(pseduo_x_[, j+1]) / n
    new_sigma_x_sq[j] <- (S[j+1,j+1] - 2 * mu[j] * sum(pseduo_x_[, j+1]))/n + mu[j]^2
  }
  new_sigma_sq <- (t(y) %*% y - 2 * t(b) %*% t(pseduo_x_) %*% y + t(b) %*% S %*% b) / n
  new_b <- solve(S) %*% t(pseduo_x_) %*% y
  return(list(
    mu = new_mu,
    sigma_x_sq = new_sigma_x_sq,
    sigma_sq = drop(new_sigma_sq),
    b = drop(new_b)
  ))
}

em <- function(params, y, X, maxit = 50, tol = 1e-5) {
  # iterative estimation
  for (i in seq_len(maxit)) {
    prev_params <- params
    params <- m_step(params, y, X)
    err <- max(abs(unlist(prev_params) - unlist(params)))
    if (err < tol) {
      return(params)
    }
  }
  cat("error > tolerance, error:", err, "\n")
  return(params)
}

em_vcov <- function(params, y, X) {
  num_params <- length(unlist(params))
  p <- ncol(X) - 1 # remove intercept column
  n <- length(y)
  M <- (!is.na(X)) + 0
  # decompose params
  mu <- params$mu
  sigma_x_sq <- params$sigma_x_sq
  sigma_sq <- params$sigma_sq
  b <- params$b
  # augmented data
  pseudo_x_ <- pseudo_x(params, y, X)
  pseudo_x_sq_ <- pseudo_x_sq(params, y, X)
  # inference: 
  conditional_params <- reverse_conditional(params, y, X)
  phi_sqs <- c(0, sapply(conditional_params, function(x) x$phi_sq))
  # score
  score <- matrix(NA, nrow = n, ncol = num_params)
  colnames(score) <- c(paste("mu", seq_len(p), sep = "_"),
                       paste("sigma", seq_len(p), "sq", sep = "_"),
                       "sigma_sq",
                       paste("b", seq(0,p), sep = "_"))
  for (i in 1:n) {
    # mu
    score[i, 1:p] <- (pseudo_x_[i, 2:(1 + p)] - mu) / sigma_x_sq
    # sigma_x_sq
    score[i, (1+p):(2*p)] <- -1 / (2 * sigma_x_sq) + 
      1 / (2 * sigma_x_sq ^ 2) * (pseudo_x_sq_[i, 2:(1 + p)] - 2 * mu * pseudo_x_[i, 2:(1 + p)] + mu)
    # sigma
    e_xxt <- pseudo_x_[i,] %*% t(pseudo_x_[i,]) + diag((1 - M[i,]) * phi_sqs)
    score[i, (2*p + 1)] <-  -1 / (2 * sigma_sq) + 
      1 / (2 * sigma_sq ^ 2) * (y[i]^2 - 2 * t(b) %*% pseudo_x_[i,] * y[i] + b %*% e_xxt %*% b)
    # b
    score[i, (2*p + 2):num_params] <- -1 / sigma_sq * (-pseudo_x_[i,] * y[i] + e_xxt %*% b)
  }
  # score -to-> fisher information
  fisher_info <- solve(n * var(score))
  return(fisher_info)
}

em_theoretical_ci <- function(est_params, y, X) {
  p <- ncol(X) - 1
  out <- data.frame(
    key = c(
      paste("mu", seq_len(p), sep = "_"),
      paste("sigma", seq_len(p), "sq", sep = "_"),
      "sigma_sq",
      paste("b", seq(0, p), sep = "_")
    ),
    value = unname(unlist(est_params))
  )
  se <- sqrt(diag(em_vcov(est_params, y, X)))
  out %>% mutate(lower_bd = value + qnorm(0.025) * se,
                 upper_db = value + qnorm(0.975) * se,
                 value = NULL) %>% 
    arrange(key) %>% 
    tibble()
}

em_bootstrapping <- function(B, params, y, X, ...) {
  pb <- progress_bar$new(
    format = "  downloading [:bar] :percent eta: :eta",
    total = B, clear = FALSE, width= 60)
  n <- length(y)
  p <- ncol(X) - 1
  boot_ests <- matrix(NA, nrow = B, ncol = length(unlist(params)))
  colnames(boot_ests) <- c(
    paste("mu", seq_len(p), sep = "_"),
    paste("sigma", seq_len(p), "sq", sep = "_"),
    "sigma_sq",
    paste("b", seq(0, p), sep = "_")
  )
  for (b in seq_len(B)) {
    pb$tick()
    boot_idx <- sample(n, n, replace = TRUE)
    boot_y <- y[boot_idx]
    boot_X <- X[boot_idx,]
    boot_ests[b,] <- unlist(em(params, boot_y, boot_X, ...))
  }
  return(boot_ests)
}

boot_confint <- function(boot_ests) {
  boot_ests %>% as.data.frame() %>% gather() %>% group_by(key) %>%
    summarise(lower_bd = quantile(value, 0.025),
              upper_bd = quantile(value, 0.975))
}

set.seed(824)
x1 <- rnorm(200, mean = 10, sd = 5)
x2 <- rnorm(200, mean = 0, sd = 1)
x3 <- rnorm(200, mean = 10, sd = 5)
y <- rnorm(200, mean = 5 + 1 * x1 + 3 * x2 + 0 * x3, sd = 2)
p1 <- 1 / (1 + exp(-(2 + 1 * x2)))
p3 <- 0.8
m1 <- rbinom(200, 1, p1)
m3 <- rbinom(200, 1, p3)
x1 <- ifelse(m1 == 1, x1, NA)
x3 <- ifelse(m3 == 1, x3, NA)
X <- cbind(1, x1, x2, x3)

params <- list(
  mu = c(10, 0, 10),
  sigma_x_sq = c(25, 1, 20),
  sigma_sq = 1,
  b = c(0, 0, 0, 0)
)

for (i in 1:20) {
  params <- m_step(params, y, X)
  print(round(unlist(params)))
}

est_param <- em(params, y, X)
em_bootstrapping(5000, params, y, X) %>% boot_confint()
em_theoretical_ci(est_param, y, X)
