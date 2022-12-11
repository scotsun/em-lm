source("./methods/multi_linear.R")
source("./methods/complete_case.R")

n <- 200
p <- 3

gamma0 <- -12
obs_rate <- 0.80
noise_scale <- 25
filename <- paste0("gamma", gamma0, "p", obs_rate * 100, "sigmaSq", noise_scale, ".rds")

message("simulation for ", filename)

# simulation
N <- 1000
param_names <- c(paste("mu", seq_len(p), sep = "_"),
                 paste("sigma", seq_len(p), "sq", sep = "_"),
                 "sigma_sq",
                 paste("b", seq(0,p), sep = "_"))
estimates <- matrix(NA, nrow = N * 2, ncol = 3 * p + 2)
colnames(estimates) <- param_names

ses <- matrix(NA, nrow = N * 2, ncol = 3 * p + 2)
colnames(ses) <- param_names

coverage <- matrix(NA, nrow = N * 2, ncol = 3 * p + 2)
colnames(coverage) <- param_names

pb <- progress_bar$new(
  format = "  downloading [:bar] :percent eta: :eta",
  total = N,
  clear = FALSE,
  width = 60
)
for (i in seq_len(N)) {
  pb$tick()
  # simulated data generation
  x1 <- rnorm(n, mean = 10, sd = 5)
  x2 <- rnorm(n, mean = 0, sd = 1)
  x3 <- rnorm(n, mean = 10, sd = 5)
  y <- rnorm(n, mean = 10 + 1 * x1 + 3 * x2 + 0 * x3, sd = sqrt(noise_scale))
  p1 <- 1 / (1 + exp(-(gamma0 + 1 * y)))
  p3 <- obs_rate
  m1 <- rbinom(n, 1, p1)
  m3 <- rbinom(n, 1, p3)
  x1 <- ifelse(m1 == 1, x1, NA)
  x3 <- ifelse(m3 == 1, x3, NA)
  X <- cbind(1, x1, x2, x3)
  rm(x1, x2, x3, p1, p3, m1, m3)
  # truth
  truth <- c(10, 0, 10, 25, 1, 25, noise_scale, 10, 1, 3, 0)
  # initial params
  params <- list(
    mu = c(10, 10, 10),
    Sigma_x = diag(c(1, 1, 1)),
    sigma_sq = 1,
    b = c(0, 0, 0, 0)
  )
  # em and log estimates & ses
  ## full
  params <- em(params, y, X)
  estimates[i, ] <- c(params$mu, diag(params$Sigma_x), params$sigma_sq, params$b)
  ses[i, ] <- em_vcov(params, y, X) %>% diag() %>% sqrt()
  coverage[i, ] <-
    (truth >= estimates[i, ] + qnorm(0.025) * ses[i, ]) &
    (truth <= estimates[i, ] + qnorm(0.975) * ses[i, ])
  ## complete-case
  cc <- complete_case_lm(y, X)
  estimates[(N + i), ] <- c(cc$point_est$mu, diag(cc$point_est$Sigma_x), cc$point_est$sigma_sq, cc$point_est$b)
  ses[(N + i), ] <- cc$se
  coverage[(N + i), ] <-
    (truth >= estimates[(N + i), ] + qnorm(0.025) * ses[(N + i), ]) &
    (truth <= estimates[(N + i), ] + qnorm(0.975) * ses[(N + i), ])
}
bias <- estimates - outer(rep(1, 2 * N), truth)

estimates <- estimates %>% as.data.frame() %>%
  mutate(estimator_type = c(rep("full", N), rep("cc", N)))
ses <- ses %>% as.data.frame() %>%
  mutate(estimator_type = c(rep("full", N), rep("cc", N)))
bias <- bias %>% as.data.frame() %>%
  mutate(estimator_type = c(rep("full", N), rep("cc", N)))
coverage <- coverage %>% as.data.frame() %>%
  mutate(estimator_type = c(rep("full", N), rep("cc", N)))

bias %>%
  group_by(estimator_type) %>%
  summarise_all(.funs = mean)

estimates %>%
  group_by(estimator_type) %>%
  summarise_all(.funs = mean)

estimates %>%
  group_by(estimator_type) %>%
  summarise_all(.funs = var)

ses %>%
  group_by(estimator_type) %>%
  summarise_all(.funs = mean)

coverage %>%
  group_by(estimator_type) %>%
  summarise_all(.funs = mean)

saveRDS(list(estimates = estimates,
             ses = ses,
             bias = bias,
             coverage = coverage),
        paste0("./simulation_rlt/mult_linear/", filename))
