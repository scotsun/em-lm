library(tidyverse)
library(progress)

pseduo_x <- function(params, x, y, m) {
  # mask NA with an arbitrary value
  x <- replace_na(x, -999)
  # decompose params
  mu_1 <- params$mu_1
  sigma_1_sq <- params$sigma_1_sq
  sigma_sq <- params$sigma_sq
  b1 <- params$b1
  b2 <- params$b2
  # reverse conditional distribution f(x|y) mean
  ksi <- (b2 * (y - b1) * sigma_1_sq + mu_1 * sigma_sq) / (b2^2 * sigma_1_sq + sigma_sq)
  return(m * x + (1 - m) * ksi)
}

.pseudo_sum_x <- function(params, x, y, m, likelihood_mode) {
  pseduo_x_ <- pseduo_x(params, x, y, m) 
  switch (likelihood_mode,
    "full" = return(sum(pseduo_x_)),
    "observed" = return(sum(m * pseduo_x_))
  )
}

pseduo_x_sq <- function(params, x, y, m) {
  # mask NA with an arbitrary value
  x <- replace_na(x, -999)
  # decompose params
  mu_1 <- params$mu_1
  sigma_1_sq <- params$sigma_1_sq
  sigma_sq <- params$sigma_sq
  b1 <- params$b1
  b2 <- params$b2
  # reverse conditional distribution f(x|y) mean & var
  ksi <- (b2 * (y - b1) * sigma_1_sq + mu_1 * sigma_sq) / (b2^2 * sigma_1_sq + sigma_sq)
  phi_sq <- sigma_1_sq * sigma_sq / (b2^2 * sigma_1_sq + sigma_sq)
  return(m * x^2 + (1 - m) * (ksi^2 + phi_sq))
}

.pseduo_sum_x_sq <- function(params, x, y, m, likelihood_mode) {
  pseduo_x_sq_ <- pseduo_x_sq(params, x, y, m) 
  switch (likelihood_mode,
          "full" = return(sum(pseduo_x_sq_)),
          "observed" = return(sum(m * pseduo_x_sq_))
  )
}


m_step <- function(params, x, y, m, likelihood_mode) {
  mu_1 <- params$mu_1
  sigma_1_sq <- params$sigma_1_sq
  sigma_sq <- params$sigma_sq
  b1 <- params$b1
  b2 <- params$b2
  
  A <- .pseudo_sum_x(params, x, y, m, "full")
  B <- .pseduo_sum_x_sq(params, x, y, m, "full")
  A.observed <- .pseudo_sum_x(params, x, y, m, "observed")
  B.observed <- .pseduo_sum_x_sq(params, x, y, m, "observed")
  C <- sum(pseduo_x(params, x, y, m) * y)
  
  n <- length(m)
  if (likelihood_mode == "full") {
    new_mu_1 <- A / n
    new_sigma_1_sq <- 2 / n * (B / 2 - mu_1 * A + n / 2 * mu_1 ^ 2)
  } else if (likelihood_mode == "observed") {
    new_mu_1 <- A.observed / sum(m)
    new_sigma_1_sq <- 2 / sum(m) * (B.observed / 2 - mu_1 * A.observed + n / 2 * mu_1 ^ 2)
  }
  new_sigma_sq <- 2 / n * (sum(y ^ 2) / 2 -
                             b1 * sum(y) -
                             b2 * C +
                             n / 2 * b1 ^ 2 +
                             b1 * b2 * A +
                             1 / 2 * b2 ^ 2 * B)
  new_bs <- solve(matrix(c(n, A, A, B), nrow = 2, ncol = 2, byrow = TRUE)) %*% c(sum(y), C)
  return(list(
    mu_1 = new_mu_1,
    sigma_1_sq = new_sigma_1_sq,
    sigma_sq = new_sigma_sq,
    b1 = new_bs[1],
    b2 = new_bs[2]
  ))
}

em <- function(params, x, y, m, likelihood_mode, maxit = 50, tol = 1e-5) {
  # iterative estimation
  for (i in seq_len(maxit)) {
    prev_params <- params
    params <- m_step(params, x, y, m, likelihood_mode)
    err <- max(abs(unlist(prev_params) - unlist(params)))
    if (err < tol) {
      return(params)
    }
  }
  cat("error > tolerance, error:", err, "\n")
  return(params)
}

em_vcov <- function(params, x, y, m) {
  # decompose params
  mu_1 <- params$mu_1
  sigma_1_sq <- params$sigma_1_sq
  sigma_sq <- params$sigma_sq
  b1 <- params$b1
  b2 <- params$b2
  # augmented data
  pseduo_x_ <- pseduo_x(params, x, y, m)
  pseduo_x_sq_ <- pseduo_x_sq(params, x, y, m)
  # inference: score
  fs_b1 <- 1 / sigma_sq * (y - b1 - b2 * pseduo_x_)
  fs_b2 <-
    1 / sigma_sq * (y * pseduo_x_ - b1 * pseduo_x_ - b2 * pseduo_x_sq_)
  fs_sigma_sq <-
    -1 / (2 * sigma_sq) + 1 / (2 * sigma_sq ^ 2) * (y ^ 2 - 2 * b1 * y - 2 * b2 * y * pseduo_x_ +
                                                      b1 ^ 2 + 2 * b1 * b2 * pseduo_x_ + b2 ^ 2 * pseduo_x_sq_)
  fs_mu_1 <- 1 / sigma_1_sq * (pseduo_x_ - mu_1)
  fs_sigma_1_sq <-
    -1 / (2 * sigma_1_sq) + 1 / (2 * sigma_1_sq ^ 2) * (pseduo_x_sq_ - 2 * mu_1 * pseduo_x_ + mu_1 ^ 2)
  score <- cbind(fs_mu_1, fs_sigma_1_sq, fs_sigma_sq, fs_b1, fs_b2)
  colnames(score) <- names(params)
  # score -to-> fisher information
  n <- length(y)
  fisher_info <- solve(n * var(score))
  return(fisher_info)
}

em_theoretical_ci <- function(est_params, x, y, m) {
  out <- est_params %>% as.data.frame() %>% gather()
  se <- sqrt(diag(em_vcov(est_params, x, y, m)))
  out %>% mutate(lower_bd = value + qnorm(0.025) * se,
                 upper_db = value + qnorm(0.975) * se,
                 value = NULL) %>% 
    arrange(key) %>% 
    tibble()
}

em_bootstrapping <- function(B, params, x, y, m, ...) {
  pb <- progress_bar$new(
    format = "  downloading [:bar] :percent eta: :eta",
    total = B, clear = FALSE, width= 60)
  n <- length(y)
  boot_ests <- matrix(NA, nrow = B, ncol = length(params))
  colnames(boot_ests) <- names(params)
  for (b in seq_len(B)) {
    pb$tick()
    boot_idx <- sample(n, n, replace = TRUE)
    boot_x <- x[boot_idx]
    boot_y <- y[boot_idx]
    boot_m <- m[boot_idx]
    boot_ests[b,] <- unlist(em(params, boot_x, boot_y, boot_m, ...))
  }
  return(boot_ests)
}

boot_confint <- function(boot_ests) {
  boot_ests %>% as.data.frame() %>% gather() %>% group_by(key) %>%
    summarise(lower_bd = quantile(value, 0.025),
              upper_bd = quantile(value, 0.975))
}
