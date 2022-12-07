library(tidyverse)

.pseudo_sum_x <- function(params, x, y, m, likelihood_mode) {
  x <- replace_na(x, -999)
  
  mu_1 <- params$mu_1
  sigma_1_sq <- params$sigma_1_sq
  sigma_sq <- params$sigma_sq
  b1 <- params$b1
  b2 <- params$b2
  
  ksi <- (b2 * (y - b1) * sigma_1_sq + mu_1 * sigma_sq) / (b2^2 * sigma_1_sq + sigma_sq)
  
  switch (likelihood_mode,
    "full" = return(sum(m * x) + sum((1 - m) * ksi)),
    "observed" = return(sum(m * x))
  )
}


.pseduo_sum_x_sq <- function(params, x, y, m, likelihood_mode) {
  x <- replace_na(x, -999)
  
  mu_1 <- params$mu_1
  sigma_1_sq <- params$sigma_1_sq
  sigma_sq <- params$sigma_sq
  b1 <- params$b1
  b2 <- params$b2
  
  n <- length(x)
  ksi <- (b2 * (y - b1) * sigma_1_sq + mu_1 * sigma_sq) / (b2^2 * sigma_1_sq + sigma_sq)
  phi_sq <- rep(sigma_1_sq * sigma_sq / (b2^2 * sigma_1_sq + sigma_sq), n)
  
  switch (likelihood_mode,
          "full" = return(sum(m * x^2) + sum((1 - m) * (ksi^2 + phi_sq))),
          "observed" = return(sum(m * x^2))
  )
}

.pseduo_sum_xy <- function(params, x, y, m) {
  x <- replace_na(x, -999)
  
  mu_1 <- params$mu_1
  sigma_1_sq <- params$sigma_1_sq
  sigma_sq <- params$sigma_sq
  b1 <- params$b1
  b2 <- params$b2
  
  ksi <- (b2 * (y - b1) * sigma_1_sq + mu_1 * sigma_sq) / (b2^2 * sigma_1_sq + sigma_sq)
  out <- sum(m * x * y) + sum((1 - m) * ksi * y)
  return(out)
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
  C <- .pseduo_sum_xy(params, x, y, m)
  
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


set.seed(823)
x <- rnorm(200, mean = 10, sd = 5)
y <- rnorm(200, mean = 10 + 2 * x, sd = 1)
p <- 0.8
m <- rbinom(200, 1, p)
x <- ifelse(m == 1, x, NA)


params <- list(
  mu_1 = 1,
  sigma_1_sq = 1,
  sigma_sq = 1,
  b1 = 0,
  b2 = 0
)


print(mean(m))
for (i in seq_len(20)) {
  params <- m_step(params, x, y, m, "full")
}
print(round(unlist(params), 3))
for (i in seq_len(20)) {
  params <- m_step(params, x, y, m, "observed")
}
print(round(unlist(params), 3))
print(c(sum((residuals(lm(y ~ x)))^2)/(length(y) - 2), coef(lm(y ~ x))))


