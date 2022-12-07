library(tidyverse)

sigma_to_theory_rSq <- function(beta, phi, theta, sigma) {
  var_matrix <- matrix(
    c(phi ^ 2, 0,
      0, theta * (1 - theta)),
    byrow = TRUE,
    nrow = 2,
    ncol = 2
  )
  var_explained <- t(beta) %*% var_matrix %*% beta
  return(var_explained / (var_explained + sigma ^ 2))
}

theory_rSq_to_sigma <- function(beta, phi, theta, rSq) {
  var_matrix <- matrix(
    c(phi ^ 2, 0,
      0, theta * (1 - theta)),
    byrow = TRUE,
    nrow = 2,
    ncol = 2
  )
  var_explained <- t(beta) %*% var_matrix %*% beta
  return(drop(sqrt(var_explained * (1 - rSq) / rSq)))
}

generate_full_data <-
  function(beta0,
           beta,
           phi,
           theta,
           sigma,
           n,
           misspecified = FALSE,
           curvature_scale = NULL) {
    x <- rnorm(n, 0, phi)
    z <- rbinom(n, 1, theta)
    epsilon <- rnorm(n, 0, sigma)
    if (misspecified) {
      x_star <- curvature_scale * (x - phi) * (x + phi) * x + x
      y <- beta0 + beta[1] * x_star + beta[2] * z + epsilon
    } else {
      y <- beta0 + beta[1] * x + beta[2] * z + epsilon
    }
    df <- data.frame(x = x, z = as.factor(z), y = y)
    return(df)
  }

logistic_scale_to_latent_rSq <- function(gamma, sigma, s, df) {
  var_matrix_lgr <- df %>%
    select(x, y) %>%
    cov()
  var_explained_lgr <- t(gamma) %*% var_matrix_lgr %*% gamma
  return(var_explained_lgr / (var_explained_lgr + 1 / 3 * s ^ 2 * pi ^ 2))
}

latent_rSq_to_logistic_scale <- function(gamma, rSq, df) {
  var_matrix_lgr <- df %>%
    select(x, y) %>% cov() %>% unname()
  var_explained_lgr <- t(gamma) %*% var_matrix_lgr %*% gamma
  return(sqrt(3 * var_explained_lgr * (1 - rSq) / rSq) / pi)
}

generate_completeness <- function(gamma0, gamma, df, s) {
  n <- dim(df)[1]
  # latent variable (sigmoid + gamma0) in LGR
  sigmoid <-
    gamma[1] * df$x + gamma[2] * df$y + rlogis(n, location = 0, scale = s)
  completeness <- (sigmoid > -gamma0 + 0)
  return(completeness)
}

generate_incomplete_data <- function(df, completeness) {
  df <- df %>% mutate(z = ifelse(completeness == 1, z, NA))
  return(df)
}

