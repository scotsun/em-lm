source("./methods/simple_linear.R")
source("./methods/complete_case.R")

# tuning parameters
n <- 200

mu_1 <- 10
sigma_1_sq <- 25
sigma_sq <- 25
b1 <- 10
b2 <- 2
truth <- c(mu_1, sigma_1_sq, sigma_sq, b1, b2)

obs_rate <- 0.80

filename <- paste0("p", 100 * obs_rate, "sigmaSq", sigma_sq, ".rds")

message("simulation for ", filename)

# simulation
N <- 5000
# c(rep("full", N), rep("observed", N), rep("cc", N))
estimates <- matrix(NA, nrow = N * 3, ncol = 5)
colnames(estimates) <-
  c("mu_1", "sigma_1_sq", "sigma_sq", "b1", "b2")

ses <- matrix(NA, nrow = N * 3, ncol = 5)
colnames(ses) <- c("mu_1", "sigma_1_sq", "sigma_sq", "b1", "b2")

coverage <- matrix(NA, nrow = N * 3, ncol = 5)
colnames(coverage) <-
  c("mu_1", "sigma_1_sq", "sigma_sq", "b1", "b2")

pb <- progress_bar$new(
  format = "  downloading [:bar] :percent eta: :eta",
  total = N,
  clear = FALSE,
  width = 60
)
for (i in seq_len(N)) {
  pb$tick()
  # simulated data generation
  x <- rnorm(n, mean = mu_1, sd = sqrt(sigma_1_sq))
  y <- rnorm(n, mean = b1 + b2 * x, sd = sqrt(sigma_sq))
  p <- obs_rate
  m <- rbinom(n, 1, p)
  x <- ifelse(m == 1, x, NA)
  # initial params
  params <- list(
    mu_1 = 1,
    sigma_1_sq = 1,
    sigma_sq = 1,
    b1 = 0,
    b2 = 0
  )
  # em and log estimates & ses
  ## full
  params_full <- em(params, x, y, m, likelihood_mode = "full")
  estimates[i, ] <- unlist(params_full)
  ses[i, ] <- em_vcov(params_full, x, y, m) %>% diag() %>% sqrt()
  coverage[i, ] <-
    (truth >= estimates[i, ] + qnorm(0.025) * ses[i, ]) &
    (truth <= estimates[i, ] + qnorm(0.975) * ses[i, ])
  ## observed
  params_observed <-
    em(params, x, y, m, likelihood_mode = "observed")
  estimates[(N + i), ] <- unlist(params_observed)
  ses[(N + i), ] <-
    em_vcov(params_observed, x, y, m) %>% diag() %>% sqrt()
  coverage[(N + i), ] <-
    (truth >= estimates[(N + i), ] + qnorm(0.025) * ses[(N + i), ]) &
    (truth <= estimates[(N + i), ] + qnorm(0.975) * ses[(N + i), ])
  ## complete-case
  cc <- complete_case_lm(y, cbind(1, x))
  estimates[(2 * N + i), ] <- unlist(cc$point_est)
  ses[(2 * N + i), ] <- cc$se
  coverage[(2 * N + i), ] <-
    (truth >= estimates[(2 * N + i), ] + qnorm(0.025) * ses[(2 * N + i), ]) &
    (truth <= estimates[(2 * N + i), ] + qnorm(0.975) * ses[(2 * N + i), ])
}
bias <- estimates - outer(rep(1, 3 * N), truth)

estimates <- estimates %>% as.data.frame() %>%
  mutate(estimator_type = c(rep("full", N), rep("observed", N), rep("cc", N)))
ses <- ses %>% as.data.frame() %>%
  mutate(estimator_type = c(rep("full", N), rep("observed", N), rep("cc", N)))
bias <- bias %>% as.data.frame() %>%
  mutate(estimator_type = c(rep("full", N), rep("observed", N), rep("cc", N)))
coverage <- coverage %>% as.data.frame() %>%
  mutate(estimator_type = c(rep("full", N), rep("observed", N), rep("cc", N)))

saveRDS(list(estimates = estimates,
             ses = ses,
             bias = bias,
             coverage = coverage),
        paste0("./simulation_rlt/simple_linear/", filename))


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
