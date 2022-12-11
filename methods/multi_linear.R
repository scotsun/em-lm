library(tidyverse)
library(progress)

reverse_cond <- function(params, y, X) {
  n <- length(y)
  # decompose params
  mu <- params$mu
  Sigma_x <- params$Sigma_x
  sigma_sq <- params$sigma_sq
  b <- params$b
  # cond_var & cond_mean
  cond_var <-
    solve(b[-1] %*% t(b[-1]) / sigma_sq + solve(Sigma_x))
  cond_mean <-
    (((y - b[1]) / sigma_sq) %o% b[-1] + (rep(1, n) %o% (solve(Sigma_x) %*% mu))[, , 1]) %*% cond_var
  return(list(ksis = cond_mean,
              phi = cond_var))
}

pseudo_x <- function(params, y, X) {
  M <- (!is.na(X)) + 0
  # arbitrary value to mask NAs
  X[!M] <- -999
  # reverse conditional distribution f(x|y) mean
  conditional_params <- reverse_cond(params, y, X)
  # column-bind the intercept data
  ksis_aug <- cbind(1, conditional_params$ksis)
  out <- M * X + (1 - M) * ksis_aug
  return(out)
}

.pseduo_sum_x <- function(params, y, X) {
  pseduo_x_ <- pseudo_x(params, y, X)
  return(apply(pseduo_x_, 2, sum))
}

pseudo_x_sq <- function(params, y, X) {
  n <- length(y)
  p <- ncol(X) - 1
  M <- (!is.na(X)) + 0
  # reverse conditional distribution f(x|y) var-cov matrix
  conditional_params <- reverse_cond(params, y, X)
  phi_aug <- cbind(0, rbind(0, conditional_params$phi))
  pseudo_x_ <- pseudo_x(params, y, X)
  # for each observation, E~[x_i x_i'] is a matrix
  # E~[x_i x_i'] = E~[x_i] E~[x_i]' + M_i M_i' phi_aug
  out <- array(NA, dim = c(p + 1, p + 1, n))
  for (i in 1:n) {
    out[,,i] <- pseudo_x_[i,] %*% t(pseudo_x_[i,]) + 
      (1-M[i,]) %*% t((1-M[i,])) * phi_aug
  }
  return(out)
}

.pseudo_sum_x_sq <- function(params, y, X) {
  pseudo_x_sq_ <- pseudo_x_sq(params, y, X)
  return(rowSums(pseudo_x_sq_, dims = 2))
}

m_step <- function(params, y, X) {
  # decompose params
  mu <- params$mu
  sigma_x_sq <- params$sigma_x_sq
  sigma_sq <- params$sigma_sq
  b <- params$b
  # get dims
  n <- length(y)
  p <- ncol(X) - 1
  
  pseudo_x_ <- pseudo_x(params, y, X)
  S <- .pseudo_sum_x_sq(params, y, X)
  
  new_mu <- apply(pseudo_x(params, y, X)[,-1], 2, mean)
  new_Sigma_x <- S[2:(p+1), 2:(p+1)] / n - mu %*% t(mu)
  new_sigma_sq <- (t(y) %*% y - 2 * t(b) %*% t(pseudo_x_) %*% y + t(b) %*% S %*% b) / n
  new_b <- solve(S) %*% t(pseudo_x_) %*% y
  return(list(
    mu = new_mu,
    Sigma_x = new_Sigma_x,
    sigma_sq = drop(new_sigma_sq),
    b = drop(new_b)
  ))
}

em <- function(params, y, X, maxit = 200, tol = 1e-5) {
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
  p <- ncol(X) - 1
  # p+1 beta, p mu, p sigma_x_sq (only consider var), 1 sigma_sq from error
  num_params <- 3 * p + 2 
  n <- length(y)
  M <- (!is.na(X)) + 0
  # decompose params
  mu <- params$mu
  Sigma_x <- params$Sigma_x
  sigma_sq <- params$sigma_sq
  b <- params$b
  # augmented data
  pseudo_x_ <- pseudo_x(params, y, X)
  pseudo_x_sq_ <- pseudo_x_sq(params, y, X)
  # inference: score
  score <- matrix(NA, nrow = n, ncol = num_params)
  colnames(score) <- c(paste("mu", seq_len(p), sep = "_"),
                       paste("sigma", seq_len(p), "sq", sep = "_"),
                       "sigma_sq",
                       paste("b", seq(0,p), sep = "_"))
  for (i in 1:n) {
    # mu
    score[i, 1:p] <- solve(Sigma_x) %*% (pseudo_x_[i, -1] - mu)
    # Sigma_x: diagonal elem
    score[i, (1 + p):(2 * p)] <-
      diag(-1 / 2 * (solve(Sigma_x) - solve(Sigma_x) %*% (
        pseudo_x_sq_[-1,-1, i] - params$mu %o% pseudo_x_[i,-1] - t(params$mu %o% pseudo_x_[i,-1]) + mu %o% mu
      ) %*% solve(Sigma_x)))
    # sigma
    score[i, (2*p + 1)] <-  -1 / (2 * sigma_sq) + 
      1 / (2 * sigma_sq ^ 2) * (y[i]^2 - 2 * t(b) %*% pseudo_x_[i,] * y[i] + b %*% pseudo_x_sq_[,,i] %*% b)
    # b
    score[i, (2*p + 2):num_params] <- -1 / sigma_sq * (-pseudo_x_[i,] * y[i] + pseudo_x_sq_[,,i] %*% b)
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
    value = c(
      est_params$mu,
      diag(est_params$Sigma_x),
      est_params$sigma_sq,
      est_params$b
    )
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
  boot_ests <- matrix(NA, nrow = B, ncol = 3 * p + 2)
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
    est_params <- em(params, boot_y, boot_X, ...)
    boot_ests[b, ] <- c(est_params$mu,
                        diag(est_params$Sigma_x),
                        est_params$sigma_sq,
                        est_params$b)
  }
  return(boot_ests)
}

boot_confint <- function(boot_ests) {
  boot_ests %>% as.data.frame() %>% gather() %>% group_by(key) %>%
    summarise(lower_bd = quantile(value, 0.025),
              upper_bd = quantile(value, 0.975))
}
